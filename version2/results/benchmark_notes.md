# Benchmark interpretation and evidence

The final run is `benchmark.json`, `benchmark.csv` and `benchmark_queries.json`, completed on 2026-09-05 from 12:55:02 to 12:56:22 UTC. It passed all 20 complete-workload comparisons and four query workloads, with no failed workers. The maximum absolute difference between complete-update implementations was 6.938893903907228e-17. The core source hash did not change during the run. The raw JSON stores every timing sample, worker metadata, input fingerprint, accuracy comparison and separate memory measurement. `benchmark_workers/` contains the individual worker checkpoints.

The earlier files named `benchmark_initial_memory_import_included.*` are **superseded**: their CSR memory pass included the first lazy import of `scipy.sparse`. The final harness imports this module before all timing and tracing. The initial files are retained as provenance and must not supply final memory claims.

## Workload and timing limits

Every complete workload uses n biallelic loci, all 3^n genotype states and all G(G+1)/2 unordered parental pairs. Within each n, every method receives exactly the same seeded full-support genotype-frequency vector. The timed operation is a **validated complete population update**, including pair-weight construction, input validation and the engine's underflow/reachability checks. It is not a pure matrix-vector multiplication. The streamed baseline enumerates all unordered pair kernels and retains the mathematics of the original full random-mating MATLAB listings; it is a Python implementation and provides no MATLAB performance measurement.

Construction uses three samples in fresh subprocesses. Each has an untimed warmup build, then clears the local-cross cache before the timed construction. Inference uses seven samples after one warmup update in a separate fresh worker for each method and n. Tracing runs in another fresh worker and is never enabled during timing. Detected numerical threadpools were limited to one thread. The Python loops are serial. CPU affinity was not pinned, and host background activity was not controlled. These measurements describe this implementation, host and input distribution; they are not universal algorithmic speed guarantees.

The range n=1..5 was selected in the task/harness before measurement. There was no external or public study preregistration; the word 'preregistered' in the raw n=6 scope-skip message refers only to that internally declared bounded range.

## Representative complete workload: five loci

G=243, U=29,646, nonzero inheritance entries=381,250, dense entries=7,203,978, density=5.2922149401%.

| Method | Median construction (ms), 3 repeats | Median validated update (ms), 7 repeats | Tracked allocation peak, separate untimed pass (bytes) |
|---|---:|---:|---:|
| Dense | 744.5774 | 12.5597 | 58,648,729 |
| CSR | 772.9773 | 9.5350 | 13,079,664 |
| Hash | 739.6091 | 119.8532 | 29,929,827 |
| Streamed kernel | No full kernel constructed | 831.9384 | 511,475 |

Construction must be accounted for when interpreting the advantage of retained kernels. Update timings use already-built dense/CSR/hash representations; streaming recomputes each cross during the update. Adding medians from separate phases would be an estimate of combined cost, not an independently measured cold end-to-end run.

The dense array payload is **57,631,824 bytes**. CSR's three numeric array payloads total **4,575,976 bytes**. The hash representation's recursively counted Python object graph is **29,427,332 bytes**. These are different retained-memory measures: numeric payload excludes object headers, while the hash graph includes unique reachable built-in objects (including any shared/interned objects reached by the graph). Neither is the same as peak live allocations or process resident memory. The streamed full-kernel payload is zero because no full kernel is retained; that does not imply zero model, input, output or temporary memory.

`tracemalloc` peaks cover kernel construction plus one update, excluding imports and the already-created model/input vector. They describe traced allocations, not total OS working set. The raw Windows process counters are absolute fresh-process peaks and include imports; do not subtract or relabel them as isolated algorithmic memory without an additional measurement design.

## Six-locus preflight

The 256 MiB per-operation representation estimate cap is **268,435,456 bytes**. No n=6 full update was timed.

| Method | Guard estimate (bytes) | Outcome |
|---|---:|---|
| Dense | 1,551,807,720 | Refused by existing allocation guard before allocation |
| CSR | 186,769,056 | Fits guard; not measured beyond bounded experiment scope |
| Hash | 1,281,296,000 | Refused by existing allocation guard before allocation |
| Streamed kernel | 5,832 for output vector | Fits output guard; not measured beyond bounded experiment scope |

The cap is an operation estimate, not an operating-system memory limit. In particular, the output-vector estimate for streaming does not claim to bound all Python metadata or elapsed work. A refusal is an observed software preflight result; the refused memory requirement is an estimate, not a measured allocation or out-of-memory event.

## Separate selected-outcome query workload

Both parents are heterozygous at every locus; the requested child is homozygous for allele 0 at every locus. Its exact probability is 4^(-n). The benchmark requests the log probability to retain very small positive outcomes without floating-point underflow. It does not enumerate or return the 3^n-state offspring distribution.

| Loci | Median time per query (ms) | Timed calls per batch | Separate traced peak (bytes) |
|---|---:|---:|---:|
| 10 | 0.011803 | 100 | 12,425 |
| 50 | 0.049378 | 50 | 24,609 |
| 100 | 0.105590 | 20 | 38,801 |
| 1,000 | 1.246460 | 5 | 266,069 |

Each row uses seven batches after one warmup query; raw batch and per-query times are preserved. At 1,000 loci the computed log probability is -1386.2943611198643, with absolute error 2.6147972675971687e-11 relative to `-n*log(4)`. The larger elapsed growth includes mixed-radix Python integer decoding and cannot, by these four measurements alone, establish a formal asymptotic complexity result. Query results must not be plotted as full-population inference speedups against the complete-kernel methods.

## Hardware and libraries

A read-only `Get-CimInstance Win32_Processor` query confirmed **12th Gen Intel Core i7-12850HX**, 16 physical cores and 24 logical processors. The queried maximum-clock field was 2,100 MHz; it is not a measured operating or turbo frequency. Windows reported 137,108,377,600 usable physical-memory bytes. The benchmark used Windows 11 build 26200, NumPy 2.5.2 and SciPy 1.18.1 from the project dependency directory. NumPy's detected OpenBLAS 0.3.34.0.0 pool was configured to one thread. Complete Python, platform, library paths, source hashes and per-worker process/memory metadata appear in the raw report.
