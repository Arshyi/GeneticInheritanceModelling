"""Verify or retrieve the public Ensembl/1000 Genomes rs12913832 aggregate snapshot.

Same contract as fetch_population_data.py: no output file is replaced. By default
the command validates the retained provenance offline against the frozen CSV.
Pass --live to query the public API; --output-dir must name a NEW directory.

rs12913832 lies in intron 86 of HERC2 and is the dominant common determinant of
blue versus brown iris colour. These are genotype calls from a reference panel.
The panel is not phenotyped: no eye colour is observed anywhere in this dataset,
and nothing here validates an eye-colour prediction for any individual.

Strand note. On the GRCh38 forward strand the alleles are A and G, and the
blue-associated allele is G. Much of the literature reports this variant on the
opposite strand as C/T, where C is the blue-associated allele. The orientation
used here is checked against the recorded mapping AND against the known
population pattern, rather than assumed from an allele label.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
URL = ("https://rest.ensembl.org/variation/human/rs12913832"
       "?population_genotypes=1;pops=1;content-type=application/json")
RELEASE_URL = "https://rest.ensembl.org/info/data?content-type=application/json"
PREFIX = "1000GENOMES:phase_3:"
GROUPS = {
    "AFR": ["ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI"],
    "AMR": ["CLM", "MXL", "PEL", "PUR"],
    "EAS": ["CDX", "CHB", "CHS", "JPT", "KHV"],
    "EUR": ["CEU", "FIN", "GBR", "IBS", "TSI"],
    "SAS": ["BEB", "GIH", "ITU", "PJL", "STU"],
}
# GRCh38 primary assembly placement of rs12913832.
CHROM, POSITION = "15", 28120472
FIELDS = ["dataset_id", "population", "superpopulation", "locus", "allele_1",
          "allele_2", "hom_1", "het", "hom_2", "n", "observation_status", "source_id"]


def nonnegative_int(value, label):
    if type(value) is not int or value < 0:
        raise ValueError(f"Invalid nonnegative integer {label}: {value!r}")
    return value


def extract_rows(payload, source_id):
    """Genotype counts cross-checked against independently supplied allele counts.

    Zero completion is permitted only when both arrays account for the same 2N
    alleles. A zero denotes a genotype unobserved in this call set, never an
    impossible state. A new allele or a changed denominator fails closed.
    """
    calls = [x for x in payload["population_genotypes"] if x["population"].startswith(PREFIX)]
    alleles = [x for x in payload.get("populations", []) if x["population"].startswith(PREFIX)]
    mapping = [x for x in payload["mappings"]
               if x["assembly_name"] == "GRCh38" and x["seq_region_name"] == CHROM
               and x["start"] == POSITION and x["end"] == POSITION and x["strand"] == 1]
    if not mapping:
        raise ValueError("Expected GRCh38 positive-strand rs12913832 mapping is absent")
    codes = sorted(code for group in GROUPS.values() for code in group)
    allowed = set(codes) | set(GROUPS) | {"ALL"}
    if {x["population"][len(PREFIX):] for x in calls} != allowed:
        raise ValueError("Unexpected or incomplete phase 3 population inventory")

    def counts_for(code):
        selected = [x for x in calls if x["population"] == PREFIX + code]
        result = {"AA": 0, "AG": 0, "GG": 0}
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
        observed = {"A": 0, "G": 0}
        for allele in [x for x in alleles if x["population"] == PREFIX + code]:
            if allele["allele"] not in observed:
                raise ValueError(f"Unexpected allele {allele['allele']} in {code}")
            observed[allele["allele"]] += nonnegative_int(allele["allele_count"], "allele count")
            if not math.isclose(allele["frequency"], allele["allele_count"] / (2 * n),
                                rel_tol=0, abs_tol=1e-10):
                raise ValueError(f"Allele frequency mismatch in {code}")
        if observed["A"] + observed["G"] != 2 * n:
            raise ValueError(f"Allele total does not match 2N in {code}")
        if 2 * result["AA"] + result["AG"] != observed["A"]:
            raise ValueError(f"Genotype and allele counts disagree in {code}")
        status = "complete" if all(result.values()) else "zero_completed_unobserved"
        return result, n, status

    rows = []
    for code in codes:
        result, n, status = counts_for(code)
        superpop = next(g for g, members in GROUPS.items() if code in members)
        rows.append({"dataset_id": "ensembl_rs12913832_phase3", "population": code,
                     "superpopulation": superpop, "locus": "rs12913832",
                     "allele_1": "A", "allele_2": "G",
                     "hom_1": result["AA"], "het": result["AG"], "hom_2": result["GG"],
                     "n": n, "observation_status": status, "source_id": source_id})
    total = sum(r["n"] for r in rows)
    if total != 2504:
        raise ValueError(f"Component populations total {total}, expected 2504")
    return rows


def orientation_check(rows):
    """Resolve which allele is blue-associated from the population pattern itself.

    This does not assume a database allele label. Northern European panels are
    known to carry the blue-associated allele at high frequency and African and
    East Asian panels at very low frequency; the check asserts that ordering and
    fails if a future data revision reverses it.
    """
    freq = {r["population"]: (2 * r["hom_2"] + r["het"]) / (2 * r["n"]) for r in rows}
    northern = ["FIN", "GBR", "CEU"]
    low = ["YRI", "ESN", "LWK", "CHB", "JPT"]
    if not all(freq[p] > 0.70 for p in northern):
        raise ValueError("Allele 2 is not at high frequency in northern European panels")
    if not all(freq[p] < 0.05 for p in low):
        raise ValueError("Allele 2 is not at low frequency in African and East Asian panels")
    return {"blue_associated_allele": "G", "brown_associated_allele": "A",
            "resolved_by": "population frequency ordering, not the database allele label",
            "northern_european_frequency": {p: round(freq[p], 6) for p in northern},
            "african_east_asian_frequency": {p: round(freq[p], 6) for p in low}}


def write_snapshot(rows, payload, release, destination):
    destination.mkdir(parents=True)
    csv_path = destination / "eye_color_genotypes.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    digest = hashlib.sha256(csv_path.read_bytes()).hexdigest()
    provenance = {
        "dataset_id": "ensembl_rs12913832_phase3",
        "variant": "rs12913832", "gene_context": "HERC2 intron 86, regulatory for OCA2",
        "assembly": "GRCh38", "chromosome": CHROM, "position": POSITION, "strand": 1,
        "url": URL, "ensembl_release": release,
        "retrieved_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "component_populations": sum(len(v) for v in GROUPS.values()),
        "individuals": sum(r["n"] for r in rows),
        "csv_sha256": digest,
        "orientation": orientation_check(rows),
        "observation_note": ("Genotype calls only. The panel carries no eye-colour "
                             "observations, so nothing here validates a phenotype prediction."),
        "reuse_note": ("Aggregate counts from a public reference panel. No individual-level "
                       "data. Superpopulation and ALL rows overlap components and are excluded."),
    }
    (destination / "eye_color_genotypes.provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    return provenance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--live", action="store_true", help="query the public API")
    parser.add_argument("--output-dir", type=Path, help="NEW directory for a refreshed snapshot")
    args = parser.parse_args()

    frozen_csv = ROOT / "data" / "eye_color_genotypes.csv"
    frozen_meta = ROOT / "data" / "eye_color_genotypes.provenance.json"

    if not args.live:
        if not frozen_csv.exists():
            raise SystemExit("No frozen snapshot yet. Run with --live --output-dir <new dir>.")
        meta = json.loads(frozen_meta.read_text(encoding="utf-8"))
        digest = hashlib.sha256(frozen_csv.read_bytes()).hexdigest()
        if digest != meta["csv_sha256"]:
            raise SystemExit("Frozen CSV does not match its recorded digest")
        rows = list(csv.DictReader(frozen_csv.open(encoding="utf-8")))
        for row in rows:
            for key in ("hom_1", "het", "hom_2", "n"):
                row[key] = int(row[key])
        orientation_check(rows)
        print(f"Offline validation passed: {len(rows)} populations, "
              f"{sum(r['n'] for r in rows)} individuals, digest {digest[:16]}")
        return

    with urllib.request.urlopen(RELEASE_URL, timeout=90) as handle:
        release = json.load(handle)["releases"][0]
    with urllib.request.urlopen(URL, timeout=180) as handle:
        payload = json.load(handle)
    rows = extract_rows(payload, "ensembl_rs12913832")

    if args.output_dir is None:
        raise SystemExit("--live requires --output-dir naming a NEW directory")
    if args.output_dir.exists():
        raise SystemExit(f"Refusing to overwrite existing directory {args.output_dir}")
    provenance = write_snapshot(rows, payload, release, args.output_dir)
    print(json.dumps({k: provenance[k] for k in
                      ("individuals", "component_populations", "ensembl_release", "csv_sha256")},
                     indent=2))
    print(f"Written to {args.output_dir}. Move into data/ only after review.")


if __name__ == "__main__":
    main()
