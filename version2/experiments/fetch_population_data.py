"""Verify or retrieve the public Ensembl/1000 Genomes rs334 aggregate snapshot.

No output files are replaced. By default the command validates a retained
provenance file offline and compares its rows with the frozen CSV. Pass --live
to query the public API; --output-dir must name a NEW directory to save a new
snapshot. Genotype calls are not clinically confirmed phenotypes.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import math
from pathlib import Path
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
URL = "https://rest.ensembl.org/variation/human/rs334?population_genotypes=1;pops=1;content-type=application/json"
RELEASE_URL = "https://rest.ensembl.org/info/data?content-type=application/json"
PREFIX = "1000GENOMES:phase_3:"
GROUPS = {
    "AFR": ["ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI"],
    "AMR": ["CLM", "MXL", "PEL", "PUR"],
    "EAS": ["CDX", "CHB", "CHS", "JPT", "KHV"],
    "EUR": ["CEU", "FIN", "GBR", "IBS", "TSI"],
    "SAS": ["BEB", "GIH", "ITU", "PJL", "STU"],
}
FIELDS = ["dataset_id", "population", "superpopulation", "locus", "allele_1",
          "allele_2", "hom_1", "het", "hom_2", "n", "observation_status", "source_id"]


def nonnegative_int(value, label):
    if type(value) is not int or value < 0:
        raise ValueError(f"Invalid nonnegative integer {label}: {value!r}")
    return value


def extract_rows(payload, source_id):
    """Cross-check genotype frequencies against independently supplied allele counts.

    Zero completion is permitted only when both arrays account for the same 2N
    alleles. This denotes an unobserved genotype in this call set, never an
    impossible biological state. A new allele or changed denominator fails closed.
    """
    calls = [x for x in payload["population_genotypes"] if x["population"].startswith(PREFIX)]
    alleles = [x for x in payload.get("populations", payload.get("allele_counts", []))
               if x["population"].startswith(PREFIX)]
    mapping = [x for x in payload["mappings"] if x["assembly_name"] == "GRCh38"
               and x["seq_region_name"] == "11" and x["start"] == 5227002
               and x["end"] == 5227002 and x["strand"] == 1]
    if not mapping:
        raise ValueError("Expected GRCh38 positive-strand rs334 mapping is absent")
    codes = sorted(code for codes in GROUPS.values() for code in codes)
    allowed = set(codes) | set(GROUPS) | {"ALL"}
    if {x["population"][len(PREFIX):] for x in calls} != allowed:
        raise ValueError("Unexpected or incomplete phase 3 population inventory")

    def counts_for(code):
        selected = [x for x in calls if x["population"] == PREFIX + code]
        result = {"TT": 0, "AT": 0, "AA": 0}
        n = sum(nonnegative_int(x["count"], "genotype count") for x in selected)
        if n == 0:
            raise ValueError(f"No genotype calls for {code}")
        for call in selected:
            key = "".join(sorted(call["genotype"].replace("|", "/").split("/")))
            if key not in result:
                raise ValueError(f"Unexpected genotype {key} in {code}")
            result[key] += call["count"]
            if not math.isclose(call["frequency"], call["count"] / n, rel_tol=0, abs_tol=1e-10):
                raise ValueError(f"Genotype frequency mismatch in {code}")
        observed_alleles = {"T": 0, "A": 0}
        for allele in [x for x in alleles if x["population"] == PREFIX + code]:
            if allele["allele"] not in observed_alleles:
                raise ValueError(f"Unexpected allele in {code}")
            observed_alleles[allele["allele"]] += nonnegative_int(allele["allele_count"], "allele count")
            if not math.isclose(allele["frequency"], allele["allele_count"] / (2*n), rel_tol=0, abs_tol=1e-10):
                raise ValueError(f"Allele frequency mismatch in {code}")
        expected_alleles = {"T": 2*result["TT"] + result["AT"],
                           "A": 2*result["AA"] + result["AT"]}
        if observed_alleles != expected_alleles or sum(observed_alleles.values()) != 2*n:
            raise ValueError(f"Allele/genotype denominator mismatch in {code}")
        return result

    counts = {code: counts_for(code) for code in allowed}
    for group, members in {**GROUPS, "ALL": codes}.items():
        if any(counts[group][key] != sum(counts[code][key] for code in members)
               for key in ("TT", "AT", "AA")):
            raise ValueError(f"Component/aggregate disagreement in {group}")
    if sum(counts["ALL"].values()) != 2504:
        raise ValueError("Phase 3 total differs from the selected 2,504-person panel")
    rows = []
    for code in codes:
        c = counts[code]
        rows.append(dict(zip(FIELDS, ["1000genomes_phase3_rs334", code,
                    next(g for g, members in GROUPS.items() if code in members),
                    "rs334", "T", "A", c["TT"], c["AT"], c["AA"], sum(c.values()),
                    "published_genotype_call_counts_zero_completed", source_id])))
    return rows


def csv_bytes(rows):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDS, quoting=csv.QUOTE_ALL, lineterminator="\r\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def read_json_url(url):
    request = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "Genetics-Version-II-aggregate-validation/1.0"})
    with urllib.request.urlopen(request, timeout=30) as response:
        raw = response.read(5_000_001)
    if len(raw) > 5_000_000:
        raise ValueError("Unexpectedly large API response")
    return json.loads(raw), raw


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--live", action="store_true", help="Query current Ensembl aggregates instead of reading the frozen provenance")
    mode.add_argument("--from-provenance", "--snapshot", dest="from_provenance", type=Path, default=ROOT / "data/observed_genotypes.provenance.json")
    parser.add_argument("--output-dir", type=Path, help="Save into a NEW directory; existing directories are refused")
    parser.add_argument("--compare", type=Path, default=ROOT / "data/observed_genotypes.csv")
    args = parser.parse_args(argv)
    if args.output_dir and args.output_dir.exists():
        parser.error("--output-dir must not already exist; frozen snapshots are never overwritten")
    if args.live:
        payload, raw = read_json_url(URL)
        release, _ = read_json_url(RELEASE_URL)
        now = datetime.now(timezone.utc)
        source_id = "ensembl_rs334_phase3_" + now.strftime("%Y%m%d")
        provenance = {
            "schema_version": 1, "dataset_id": "1000genomes_phase3_rs334",
            "accessed_date": now.date().isoformat(), "retrieved_at_utc": now.isoformat(),
            "source_id": source_id, "source_url": URL, "ensembl_release": release["releases"],
            "primary_publication_doi": "10.1038/nature15393",
            "http_response_sha256": hashlib.sha256(raw).hexdigest(), "raw_http_response_saved": False,
            "mappings": payload["mappings"],
            "population_genotypes": [x for x in payload["population_genotypes"] if x["population"].startswith(PREFIX)],
            "allele_counts": [x for x in payload["populations"] if x["population"].startswith(PREFIX)],
            "component_superpopulation_map": GROUPS,
            "scope": "Aggregate reference-panel genotype calls; no clinical prediction validation. Zero-completed absent entries are not biologically impossible.",
            "reuse": {"policy_url": "https://www.internationalgenome.org/IGSR_disclaimer/", "individual_level_data_used": False},
        }
    else:
        provenance = json.loads(args.from_provenance.read_text(encoding="utf-8-sig"))
        payload, source_id = provenance, provenance["source_id"]
    rows = extract_rows(payload, source_id)
    generated = csv_bytes(rows)
    provenance["csv_sha256"] = hashlib.sha256(generated).hexdigest()
    old_rows = list(csv.DictReader(args.compare.open(encoding="utf-8-sig", newline="")))
    # Retrieval date/source ID can differ without the biological counts changing.
    value_fields = [x for x in FIELDS if x != "source_id"]
    values = lambda rs: [[str(r[x]) for x in value_fields] for r in sorted(rs, key=lambda r: r["population"])]
    matched = values(old_rows) == values(rows)
    if args.output_dir:
        args.output_dir.mkdir(parents=True, exist_ok=False)
        (args.output_dir / "observed_genotypes.csv").write_bytes(generated)
        (args.output_dir / "observed_genotypes.provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"rows": len(rows), "n": sum(r["n"] for r in rows),
                      "genotypes_TT_AT_AA": [sum(r[k] for r in rows) for k in ("hom_1", "het", "hom_2")],
                      "matches_frozen_counts": matched, "csv_sha256": provenance["csv_sha256"],
                      "output_dir": str(args.output_dir) if args.output_dir else None}, indent=2))
    return 0 if matched else 2


if __name__ == "__main__":
    raise SystemExit(main())
