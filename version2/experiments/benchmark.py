"""Bounded, reproducible full-kernel and separately labeled query benchmarks.

Each construction sample, inference configuration and memory pass runs in a fresh
subprocess. Timed sections exclude interpreter startup/imports. No partial
distribution is accepted when a resource guard refuses a workload.
"""
from __future__ import annotations

import argparse
import csv
import ctypes
import gc
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import time
import tracemalloc

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "results"
CAP = 256 * 1024**2
SEED = 20260905
METHODS = ("dense", "csr", "hash", "streamed_kernel")
# Must precede NumPy/SciPy imports in every child.
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS"):
    os.environ[_name] = "1"
sys.path.insert(0, str(BASE))
if (BASE / ".deps").exists():
    sys.path.insert(0, str(BASE / ".deps"))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def process_memory():
    """Windows process absolute counters; not an operation-only peak estimate."""
    if os.name != "nt":
        return None
    from ctypes import wintypes
    class Counters(ctypes.Structure):
        _fields_ = [("cb",wintypes.DWORD),("PageFaultCount",wintypes.DWORD),
                    ("PeakWorkingSetSize",ctypes.c_size_t),("WorkingSetSize",ctypes.c_size_t),
                    ("QuotaPeakPagedPoolUsage",ctypes.c_size_t),("QuotaPagedPoolUsage",ctypes.c_size_t),
                    ("QuotaPeakNonPagedPoolUsage",ctypes.c_size_t),("QuotaNonPagedPoolUsage",ctypes.c_size_t),
                    ("PagefileUsage",ctypes.c_size_t),("PeakPagefileUsage",ctypes.c_size_t)]
    result = Counters(); result.cb = ctypes.sizeof(result)
    current_process = ctypes.windll.kernel32.GetCurrentProcess
    current_process.restype = wintypes.HANDLE
    getter = ctypes.windll.psapi.GetProcessMemoryInfo
    getter.argtypes = [wintypes.HANDLE,ctypes.POINTER(Counters),wintypes.DWORD]
    if not getter(current_process(),ctypes.byref(result),result.cb):
        return None
    return {"working_set_bytes":int(result.WorkingSetSize),
            "absolute_peak_working_set_bytes":int(result.PeakWorkingSetSize),
            "pagefile_usage_bytes":int(result.PagefileUsage),
            "absolute_peak_pagefile_bytes":int(result.PeakPagefileUsage)}


def physical_memory():
    if os.name != "nt":
        return None
    class MemoryStatus(ctypes.Structure):
        _fields_ = [("length",ctypes.c_ulong),("load",ctypes.c_ulong),
                    ("total_physical",ctypes.c_ulonglong),("available_physical",ctypes.c_ulonglong),
                    ("total_pagefile",ctypes.c_ulonglong),("available_pagefile",ctypes.c_ulonglong),
                    ("total_virtual",ctypes.c_ulonglong),("available_virtual",ctypes.c_ulonglong),
                    ("available_extended_virtual",ctypes.c_ulonglong)]
    value = MemoryStatus(); value.length = ctypes.sizeof(value)
    if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(value)):
        return {"total_physical_bytes":int(value.total_physical),
                "available_physical_bytes":int(value.available_physical),
                "memory_load_percent":int(value.load)}
    return None


def recursive_object_bytes(obj, seen=None):
    """Unique retained built-in object graph; a Python footprint, not numeric payload."""
    seen = set() if seen is None else seen
    identity = id(obj)
    if identity in seen:
        return 0
    seen.add(identity)
    size = sys.getsizeof(obj)
    if isinstance(obj,dict):
        size += sum(recursive_object_bytes(k,seen)+recursive_object_bytes(v,seen) for k,v in obj.items())
    elif isinstance(obj,(list,tuple,set,frozenset)):
        size += sum(recursive_object_bytes(v,seen) for v in obj)
    return size


def estimates(model):
    return {"dense":8*model.G*model.U,
            "csr":2*(16*model.nnz+8*(model.U+model.G+2)),
            "hash":220*model.nnz+100*model.U,
            "streamed_kernel":8*model.G}


def footprint(kernel, method):
    if method == "dense":
        return {"numeric_array_payload_bytes":int(kernel.nbytes),
                "array_dtypes":{"data":str(kernel.dtype)},
                "hash_recursive_python_object_bytes":None,
                "stored_entries":int(kernel.size)}
    if method == "csr":
        return {"numeric_array_payload_bytes":int(kernel.data.nbytes+kernel.indices.nbytes+kernel.indptr.nbytes),
                "array_payload_detail_bytes":{k:int(getattr(kernel,k).nbytes) for k in ("data","indices","indptr")},
                "array_dtypes":{k:str(getattr(kernel,k).dtype) for k in ("data","indices","indptr")},
                "hash_recursive_python_object_bytes":None,"stored_entries":int(kernel.nnz)}
    if method == "hash":
        return {"numeric_array_payload_bytes":None,"array_dtypes":None,
                "hash_recursive_python_object_bytes":recursive_object_bytes(kernel),
                "stored_entries":sum(len(v) for v in kernel.values())}
    return {"numeric_array_payload_bytes":0,"hash_recursive_python_object_bytes":None,
            "stored_entries":0,"note":"No full kernel materialized; model metadata and output vector excluded."}


def metadata(np, scipy, threadpoolctl):
    return {"utc":time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime()),
            "python":sys.version,"python_executable":sys.executable,"platform":platform.platform(),
            "machine":platform.machine(),"processor":platform.processor(),
            "processor_identifier":os.environ.get("PROCESSOR_IDENTIFIER"),
            "logical_cpu_count":os.cpu_count(),"physical_memory":physical_memory(),
            "numpy":np.__version__,"numpy_path":np.__file__,
            "scipy":scipy.__version__,"scipy_path":scipy.__file__,
            "threadpoolctl":threadpoolctl.__version__,"threadpools":threadpoolctl.threadpool_info(),
            "core_sha256":digest(BASE/"genetics"/"core.py"),"harness_sha256":digest(__file__),
            "pid":os.getpid(),"gc_enabled":gc.isenabled()}


def worker(job):
    import numpy as np
    import scipy
    import scipy.sparse  # Eager import: excluded from construction/memory measurements.
    import threadpoolctl
    from genetics.core import InheritanceModel, ResourceLimitError, local_cross

    phase, n = job["phase"], job["loci"]
    method = job.get("method","factored_log_query")
    with threadpoolctl.threadpool_limits(limits=1):
        meta = metadata(np,scipy,threadpoolctl)
        if any(p.get("num_threads",1)!=1 for p in meta["threadpools"]):
            raise RuntimeError("Could not enforce one thread in a detected numerical threadpool")
        model = InheritanceModel([2]*n,max_bytes=CAP)
        info = {"phase":phase,"loci":n,"method":method,"metadata":meta,
                "genotype_states":model.G,"unordered_pair_columns":model.U,
                "theoretical_nnz":model.nnz,"budget_bytes":CAP}
        if phase == "preflight":
            info["status"] = "preflight_only"
            info["estimates_bytes"] = estimates(model)
            outcomes = []
            for kind in METHODS:
                estimate = info["estimates_bytes"][kind]
                try:
                    # Exactly the existing model guard, with no kernel allocation.
                    model._guard(estimate,kind+" full kernel" if kind!="streamed_kernel" else "population result")
                except ResourceLimitError as exc:
                    outcomes.append({"method":kind,"status":"resource_limit_refused",
                                     "estimate_bytes":estimate,"reason":str(exc),"measured":False})
                else:
                    outcomes.append({"method":kind,"status":"not_measured_scope_limit",
                                     "estimate_bytes":estimate,
                                     "reason":"Beyond preregistered complete-workload range n=1..5; estimate fits guard.",
                                     "measured":False})
            info["outcomes"] = outcomes
            return info
        if phase.startswith("query"):
            hetero = model.encode((1,)*n)
            target = 0
            expected = -n*math.log(4)
            def query():
                return model.log_probability(hetero,hetero,target)
            info.update({"query":"log probability of all-00 offspring from two all-01 parents",
                         "expected_log_probability":expected,"full_output_enumerated":False,
                         "positive_full_cross_support":3**n})
            if phase == "query_memory":
                # Fresh process; no query warmup or allocation trace during timing runs.
                local_cross.cache_clear(); gc.collect()
                before = process_memory()
                tracemalloc.start()
                result = query()
                retained,peak = tracemalloc.get_traced_memory(); tracemalloc.stop()
                info.update({"status":"measured","log_probability":result,
                             "tracemalloc_retained_bytes":retained,"tracemalloc_peak_bytes":peak,
                             "process_memory_before":before,"process_memory_after":process_memory(),
                             "tracemalloc_scope":"one query; existing model and integer inputs excluded"})
            else:
                repeats = job["repeats"]
                batch = {10:100,50:50,100:20,1000:5}.get(n,10)
                query()
                times=[]
                for _ in range(repeats):
                    start=time.perf_counter()
                    for _ in range(batch): result=query()
                    times.append(time.perf_counter()-start)
                model_times=[]
                for _ in range(job["build_repeats"]):
                    gc.collect(); start=time.perf_counter()
                    temporary=InheritanceModel([2]*n,max_bytes=CAP)
                    model_times.append(time.perf_counter()-start)
                    del temporary
                info.update({"status":"measured","warmup_calls":1,"calls_per_timed_batch":batch,
                             "raw_batch_seconds":times,"raw_seconds_per_query":[t/batch for t in times],
                             "raw_model_construction_seconds":model_times,"log_probability":result})
            info["absolute_log_error"] = abs(info["log_probability"]-expected)
            if info["absolute_log_error"]>1e-9: raise AssertionError(info)
            return info

        rng=np.random.default_rng(SEED+n)
        frequencies=rng.random(model.G)+np.finfo(float).eps
        frequencies/=frequencies.sum()
        info.update({"input_seed":SEED+n,"input_sha256":hashlib.sha256(frequencies.tobytes()).hexdigest(),
                     "input_support":int(np.count_nonzero(frequencies)),
                     "guard_estimate_bytes":estimates(model)[method],
                     "full_dense_entries":model.G*model.U,
                     "density":model.nnz/(model.G*model.U)})
        def build():
            return None if method=="streamed_kernel" else model.kernel(method)
        if phase == "build":
            if method=="streamed_kernel":
                info.update({"status":"not_applicable","seconds":None,
                             "reason":"No full representation construction; common model setup excluded."})
                return info
            # Warm interpreter/build path once; clear local-cross cache per cold build.
            warm=build(); del warm; local_cross.cache_clear(); gc.collect()
            start=time.perf_counter(); kernel=build(); duration=time.perf_counter()-start
            info.update({"status":"measured","seconds":duration,"warmup_builds":1,
                         "local_cross_cache":"cleared before timed construction",
                         "footprint":footprint(kernel,method)})
            return info
        if phase == "memory":
            local_cross.cache_clear(); gc.collect()
            before=process_memory(); tracemalloc.start()
            kernel=build()
            result=model.next_generation(frequencies,kernel)
            retained,peak=tracemalloc.get_traced_memory(); tracemalloc.stop()
            info.update({"status":"measured","tracemalloc_retained_bytes":retained,
                         "tracemalloc_peak_bytes":peak,"process_memory_before":before,
                         "process_memory_after":process_memory(),"footprint":footprint(kernel,method),
                         "tracemalloc_scope":"kernel construction plus one full update; imports, model and input vector excluded",
                         "output_mass":float(result.sum())})
            return info
        kernel=build()
        model.next_generation(frequencies,kernel)  # one untimed warmup
        times=[]
        for _ in range(job["repeats"]):
            start=time.perf_counter()
            result=model.next_generation(frequencies,kernel)
            times.append(time.perf_counter()-start)
        if not np.isfinite(result).all() or (result<0).any() or abs(result.sum()-1)>1e-12:
            raise AssertionError("Output violated probability invariants")
        info.update({"status":"measured","raw_seconds":times,"warmup_updates":1,
                     "timed_operation":"one full next_generation including pair-weight formation and input validation",
                     "local_cross_cache":"warm from construction and untimed update",
                     "output":result.tolist(),"output_mass":float(result.sum()),
                     "footprint":footprint(kernel,method)})
        return info


def launch(job, directory, resume=False, timeout=240):
    key=f"{job['phase']}_n{job['loci']}_{job.get('method','query')}_{job.get('sample',0)}"
    cache=directory/(key+".json")
    core_hash=digest(BASE/"genetics"/"core.py")
    script_hash=digest(__file__)
    if resume and cache.exists():
        previous=json.loads(cache.read_text(encoding="utf-8"))
        if previous.get("job")==job and previous.get("core_sha256")==core_hash and previous.get("harness_sha256")==script_hash:
            print("Reused "+key,flush=True)
            return previous["result"]
    start=time.perf_counter()
    try:
        proc=subprocess.run([sys.executable,"-X","utf8",str(Path(__file__).resolve()),
                             "--worker",json.dumps(job)],capture_output=True,text=True,encoding="utf-8",
                            timeout=timeout,creationflags=subprocess.CREATE_NO_WINDOW if os.name=="nt" else 0)
        if proc.returncode:
            result={"phase":job["phase"],"loci":job["loci"],"method":job.get("method","factored_log_query"),
                    "status":"failed","returncode":proc.returncode,"stderr":proc.stderr[-6000:],"stdout":proc.stdout[-2000:]}
        else:
            result=json.loads(proc.stdout)
    except subprocess.TimeoutExpired:
        result={"phase":job["phase"],"loci":job["loci"],"method":job.get("method","factored_log_query"),
                "status":"timeout_not_measured","timeout_seconds":timeout}
    result["worker_wall_seconds_including_startup"]=time.perf_counter()-start
    wrapper={"job":job,"core_sha256":core_hash,"harness_sha256":script_hash,"result":result}
    cache.write_text(json.dumps(wrapper,indent=2)+"\n",encoding="utf-8")
    print(f"{key}: {result['status']} ({result['worker_wall_seconds_including_startup']:.2f}s worker wall)",flush=True)
    return result


def stats(values):
    return {"median":statistics.median(values),"minimum":min(values),"maximum":max(values),
            "mean":statistics.mean(values),"stdev":statistics.stdev(values) if len(values)>1 else 0}


def main(args):
    OUT.mkdir(parents=True,exist_ok=True)
    directory=OUT/"benchmark_workers"; directory.mkdir(exist_ok=True)
    start_hash=digest(BASE/"genetics"/"core.py")
    started=time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime())
    records=[]
    for n in args.loci:
        for method in METHODS:
            if method!="streamed_kernel":
                for sample in range(args.build_repeats):
                    records.append(launch({"phase":"build","loci":n,"method":method,"sample":sample},directory,args.resume,args.timeout))
            records.append(launch({"phase":"inference","loci":n,"method":method,"repeats":args.repeats},directory,args.resume,args.timeout))
            records.append(launch({"phase":"memory","loci":n,"method":method},directory,args.resume,args.timeout))
    preflight=launch({"phase":"preflight","loci":6},directory,args.resume,args.timeout)
    records.append(preflight)
    for n in args.query_loci:
        records.append(launch({"phase":"query_timing","loci":n,"repeats":args.repeats,
                               "build_repeats":args.build_repeats},directory,args.resume,args.timeout))
        records.append(launch({"phase":"query_memory","loci":n},directory,args.resume,args.timeout))

    summaries=[]; accuracy=[]
    for n in args.loci:
        measurements=[r for r in records if r["loci"]==n and r["phase"]=="inference" and r["status"]=="measured"]
        if not measurements: continue
        reference=next((r for r in measurements if r["method"]=="dense"),measurements[0])
        for result in measurements:
            error=max(abs(a-b) for a,b in zip(result["output"],reference["output"]))
            accuracy.append({"loci":n,"method":result["method"],"reference_method":reference["method"],
                             "max_absolute_difference":error,"same_input":result["input_sha256"]==reference["input_sha256"],
                             "passed":error<=1e-12 and result["input_sha256"]==reference["input_sha256"]})
            construction=[r["seconds"] for r in records if r["loci"]==n and r["method"]==result["method"] and r["phase"]=="build" and r["status"]=="measured"]
            memory=next((r for r in records if r["loci"]==n and r["method"]==result["method"] and r["phase"]=="memory"),{})
            summaries.append({"workload":"complete_random_mating_update","loci":n,"method":result["method"],
                              "G":result["genotype_states"],"U":result["unordered_pair_columns"],
                              "nnz":result["theoretical_nnz"],"density":result["density"],
                              "construction_status":("measured" if construction else
                                                     "not_applicable" if result["method"]=="streamed_kernel" else "failed_or_missing"),
                              "construction_seconds":stats(construction) if construction else None,
                              "inference_seconds":stats(result["raw_seconds"]),
                              "construction_repetitions":len(construction),"inference_repetitions":len(result["raw_seconds"]),
                              "numeric_array_payload_bytes":result["footprint"]["numeric_array_payload_bytes"],
                              "hash_recursive_python_object_bytes":result["footprint"].get("hash_recursive_python_object_bytes"),
                              "tracemalloc_peak_bytes":memory.get("tracemalloc_peak_bytes"),
                              "tracemalloc_status":memory.get("status","missing"),
                              "output_mass":result["output_mass"],"max_absolute_difference_from_dense":error})
    queries=[]
    for r in records:
        if r["phase"]=="query_timing" and r["status"]=="measured":
            memory=next((m for m in records if m["phase"]=="query_memory" and m["loci"]==r["loci"]),{})
            queries.append({"workload":"one_conditional_log_probability_query","loci":r["loci"],
                            "query_seconds":stats(r["raw_seconds_per_query"]),
                            "model_construction_seconds":stats(r["raw_model_construction_seconds"]),
                            "calls_per_timed_batch":r["calls_per_timed_batch"],
                            "repetitions":len(r["raw_batch_seconds"]),"log_probability":r["log_probability"],
                            "absolute_log_error":r["absolute_log_error"],"full_output_enumerated":False,
                            "tracemalloc_peak_bytes":memory.get("tracemalloc_peak_bytes")})
    failures=[r for r in records if r["status"] in ("failed","timeout_not_measured")]
    stable=start_hash==digest(BASE/"genetics"/"core.py")
    report={"benchmark_started_utc":started,"benchmark_finished_utc":time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime()),
            "status":"passed" if not failures and stable and all(a["passed"] for a in accuracy) else "incomplete_or_failed",
            "core_unchanged_during_run":stable,"core_sha256":start_hash,"harness_sha256":digest(__file__),
            "budget_bytes":CAP,"input_seed_base":SEED,"full_loci":args.loci,"query_loci":args.query_loci,
            "methodology":{
                "complete_workload":"One exact no-pruning random-mating update from the identical full-support seeded genotype vector for each n. Biallelic independent-segregation loci. All G(G+1)/2 parental columns included.",
                "construction":"Three fresh-worker samples per materialized representation, one warmup construction each, cold local-cross cache before each timed build. Imports, model setup and input creation excluded.",
                "inference":"Seven perf_counter samples after one warmup; validated end-to-end model.next_generation includes input validation, pair-weight formation and underflow/reachability checks. These are not pure matrix-vector timings. No allocation tracing during timing.",
                "memory":"Separate fresh untimed worker traces representation construction plus one update. tracemalloc is tracked allocator peak, not total OS memory. Windows absolute process peaks include imports. Retained numeric ndarray payload and recursive hash Python object footprint are different measures and must not be presented as interchangeable.",
                "threads":"Environment variables and threadpoolctl limit detected numerical libraries to one thread. Python pair loops are serial. Host background activity is uncontrolled; no CPU affinity pinning.",
                "streamed_kernel":"Python streaming full pair-enumeration baseline with the same mathematics as the original ABO/ABO+Rh random-mating listings. Uses unordered pairs with off-diagonal weight 2; not a MATLAB timing and not an execution of the paper's truncated square matrices.",
                "queries":"Separate workload: one selected offspring log probability from two all-heterozygous parents. This does not construct or return the full population distribution and is not a full-kernel speedup comparison.",
                "preflight":"n=6 allocation guards checked without allocation. Fitting representations beyond n=5 are scope skips, not fabricated resource refusals.",
                "timing_limits":"Microbenchmarks on one available Windows CPU. Results are implementation- and workload-specific; no GPU or MATLAB timings."},
            "full_workload_summary":summaries,"query_workload_summary":queries,
            "accuracy_checks":accuracy,"preflight":preflight,"failed_workers":failures,"raw_workers":records}
    (OUT/"benchmark.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    flat=[]
    for row in summaries:
        flat.append({k:v for k,v in row.items() if k not in ("construction_seconds","inference_seconds")}
                    | {"construction_median_seconds":row["construction_seconds"]["median"] if row["construction_seconds"] else None,
                       "inference_median_seconds":row["inference_seconds"]["median"]})
    if flat:
        with (OUT/"benchmark.csv").open("w",newline="",encoding="utf-8") as fh:
            writer=csv.DictWriter(fh,fieldnames=list(flat[0]));writer.writeheader();writer.writerows(flat)
    (OUT/"benchmark_queries.json").write_text(json.dumps(queries,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":report["status"],"complete_rows":len(summaries),"query_rows":len(queries),
                      "failed_workers":len(failures),"all_accuracy_checks_pass":all(a["passed"] for a in accuracy)},indent=2),flush=True)
    if report["status"]!="passed": raise SystemExit(1)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker",help=argparse.SUPPRESS)
    parser.add_argument("--loci",type=int,nargs="+",default=[1,2,3,4,5])
    parser.add_argument("--query-loci",type=int,nargs="+",default=[10,50,100,1000])
    parser.add_argument("--repeats",type=int,default=7)
    parser.add_argument("--build-repeats",type=int,default=3)
    parser.add_argument("--timeout",type=int,default=240)
    parser.add_argument("--resume",action="store_true")
    options=parser.parse_args()
    if options.worker:
        print(json.dumps(worker(json.loads(options.worker))))
    else:
        if any(n<1 or n>5 for n in options.loci): parser.error("Complete workloads are bounded to n=1..5")
        if options.repeats<1 or options.build_repeats<1: parser.error("Repetition counts must be positive")
        main(options)
