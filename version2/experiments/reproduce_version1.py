"""Reproduce computations visible in the immutable 54-page Version 1 PDF.

This is a Python translation, not an execution of MATLAB. Legacy defects are
preserved deliberately and documented in research/version1_audit.md. The audit
does not import the corrected Version 2 engine, so it is an independent baseline.
Only Python's standard library is required. Run from any working directory.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "version2" / "results"
SOURCE = ROOT / "Bioinformatics-Arshyia Mehran.pdf"
ABO = ("AA", "AO", "BB", "BO", "AB", "OO")
RH = ("DD", "Dd", "dd")
JOINT = tuple(f"{a}/{r}" for a in ABO for r in RH)
SICKLE = ("AA", "AS", "SS")
M_SICKLE = ((F(1), F(1, 2), F(1, 4)),
            (F(0), F(1, 2), F(1, 2)),
            (F(0), F(0), F(1, 4)))
M_ABO_DISPLAYED = tuple(tuple(F(x) for x in row.split()) for row in (
    "0 0 .25 0 0 0", ".5 0 .5 0 .25 1", "0 0 0 0 0 0",
    "0 0 0 .5 .25 0", "0 0 0 0 .25 0", ".5 1 .25 .5 .25 0"))
ABO_TOP6_PAIRS = (("AO", "OO"), ("OO", "OO"), ("AO", "AO"),
                  ("BO", "OO"), ("AO", "BO"), ("AA", "OO"))

# Literal transcription of all three six-column blocks on PDF pp. 44-45.
_JOINT_BLOCKS = (
    ("1 .25 0 .25 .0625 0", "0 .5 0 0 .125 0", "0 .25 1 0 .0625 .25",
     "0 0 0 .5 .125 0", "0 0 0 0 .25 0", "0 0 0 0 .125 .5",
     "0 0 0 0 0 0", "0 0 0 0 0 0", "0 0 0 0 0 0",
     "0 0 0 0 0 0", "0 0 0 0 0 0", "0 0 0 0 0 0",
     "0 0 0 0 0 0", "0 0 0 0 0 0", "0 0 0 0 0 0",
     "0 0 0 .25 .0625 0", "0 0 0 0 .125 0", "0 0 0 0 .0625 .25"),
    ("0 0 0 0 0 0", "0 0 0 0 0 0", "0 0 0 0 0 0",
     "0 0 0 0 0 0", "0 0 0 0 0 0", "0 0 0 0 0 0",
     "1 .25 0 .25 .0625 0", "0 .5 0 0 .125 0", "0 .25 1 0 .0625 .25",
     "0 0 0 .5 .125 0", "0 0 0 0 .25 0", "0 0 0 0 .125 .5",
     "0 0 0 0 0 0", "0 0 0 0 0 0", "0 0 0 0 0 0",
     "0 0 0 .25 .0625 0", "0 0 0 0 .125 0", "0 0 0 0 .0625 .25"),
    (".25 .0625 0 0 0 0", "0 .125 0 0 0 0", "0 .0625 .25 0 0 0",
     "0 0 0 0 0 0", "0 0 0 0 0 0", "0 0 0 0 0 0",
     ".25 .0625 0 0 0 0", "0 .125 0 0 0 0", "0 .0625 .25 0 0 0",
     "0 0 0 0 0 0", "0 0 0 0 0 0", "0 0 0 0 0 0",
     ".5 .125 0 0 0 0", "0 .25 0 0 0 0", "0 .125 .5 0 0 0",
     "0 0 0 1 .25 0", "0 0 0 0 .5 0", "0 0 0 0 .25 1"),
)
M_JOINT_DISPLAYED = tuple(tuple(F(v) for b in _JOINT_BLOCKS for v in b[i].split())
                          for i in range(18))

# PDF p. 42: rank order and percentages as printed, including tied ranks.
PAPER_TOP18 = (
    ("AO/DD", "OO/DD", 8.87), ("OO/DD", "OO/Dd", 8.13),
    ("AO/DD", "OO/Dd", 6.38), ("AO/Dd", "OO/DD", 6.38),
    ("OO/DD", "OO/DD", 5.65), ("AO/DD", "AO/Dd", 5.00),
    ("AO/Dd", "OO/Dd", 4.59), ("AO/DD", "AO/DD", 3.48),
    ("OO/Dd", "OO/Dd", 2.93), ("BO/DD", "OO/DD", 2.63),
    ("AO/DD", "BO/DD", 2.06), ("BO/DD", "OO/Dd", 1.89),
    ("BO/Dd", "OO/DD", 1.89), ("AO/Dd", "AO/Dd", 1.80),
    ("AA/DD", "OO/DD", 1.74), ("AO/Dd", "BO/DD", 1.48),
    ("AO/DD", "BO/Dd", 1.48), ("OO/DD", "OO/dd", 1.46),
)


def matvec(matrix, vector):
    return [sum((a * b for a, b in zip(row, vector)), F(0)) for row in matrix]


def punnett(g1: str, g2: str, classes=ABO):
    """Four equally probable allele combinations, matching PDF pp. 39, 48-49."""
    result = [F(0) for _ in classes]
    for a in g1:
        for b in g2:
            child = "".join(sorted((a, b)))
            result[classes.index(child)] += F(1, 4)
    return result


def punnett18(g1: str, g2: str):
    a1, r1 = g1.split("/")
    a2, r2 = g2.split("/")
    pa, pr = punnett(a1, a2), punnett(r1, r2, RH)
    return [a * r for a in pa for r in pr]


def legacy_sickle(dad: str, mom: str, generations: int):
    """PDF pp. 26-27, preserving record-before-removal and no normalization."""
    x = punnett(dad, mom, SICKLE)
    history = [x.copy()]
    for _ in range(2, generations + 1):
        x_next = matvec(M_SICKLE, x)
        history.append(x_next.copy())  # history(:,gen) = X_next;
        x_next[2] = F(0)               # X_next(3) = 0; after recording
        x = x_next
    return history


def legacy_random_mating(dad: str, mom: str, generations: int, joint=False):
    """Literal ordered-pair loop from pp. 37 and 46, using exact rational values."""
    classes = JOINT if joint else ABO
    cross = punnett18 if joint else punnett
    history = [cross(dad, mom)]
    kernels = [[cross(a, b) for b in classes] for a in classes]
    for _ in range(2, generations + 1):
        previous = history[-1]
        new = [F(0) for _ in classes]
        for i in range(len(classes)):
            for j in range(len(classes)):
                pair_weight = previous[i] * previous[j]
                if pair_weight > 0:
                    for k, probability in enumerate(kernels[i][j]):
                        new[k] += pair_weight * probability
        history.append(new)
    return history


def phenotypes(vector, joint=False):
    if not joint:
        return [vector[0] + vector[1], vector[2] + vector[3], vector[4], vector[5]]
    result = [F(0)] * 8
    for g, value in zip(JOINT, vector):
        a, r = g.split("/")
        base = {"AA": 0, "AO": 0, "BB": 2, "BO": 2, "AB": 4, "OO": 6}[a]
        result[base + (r == "dd")] += value
    return result


def pair_weights(classes, frequencies):
    """Unordered pairs: f_i^2 for equal parents, 2 f_i f_j otherwise."""
    pairs = []
    for i, a in enumerate(classes):
        for j in range(i, len(classes)):
            pairs.append((a, classes[j], frequencies[i] * frequencies[j] * (1 if i == j else 2)))
    return sorted(pairs, key=lambda p: (-p[2], p[0], p[1]))


def coverage_audit():
    a, b, o = map(F, (".26", ".077", ".663"))
    exact = [a*a, 2*a*o, b*b, 2*b*o, 2*a*b, o*o]
    variants = {
        "paper_alleles_exact": exact,
        "paper_p32_genotype_table": list(map(F, (".0676", ".3448", ".0059", ".1021", ".04", ".4396"))),
        "paper_p32_pair_formula_inputs_raw": list(map(F, (".068", ".345", ".006", ".102", ".040", ".440"))),
    }
    rounded = variants["paper_p32_pair_formula_inputs_raw"]
    variants["paper_p32_pair_formula_inputs_normalized"] = [x/sum(rounded) for x in rounded]
    qd = math.sqrt(.07)  # 7% Rh-negative is a phenotype frequency, q_d^2.
    rh = [(1-qd)**2, 2*qd*(1-qd), qd**2]
    summary, rows = {}, []
    for name, freq in variants.items():
        pairs = pair_weights(ABO, freq)
        joint_freq = [float(x)*r for x in freq for r in rh]
        joint_pairs = pair_weights(JOINT, joint_freq)
        covered = sum(p[2] for p in pairs[:6])
        covered_joint = sum(p[2] for p in joint_pairs[:18])
        summary[name] = {
            "genotype_frequencies": [float(x) for x in freq],
            "genotype_mass": float(sum(freq)),
            "total_unordered_pair_mass": float(sum(p[2] for p in pairs)),
            "top6_percent_raw_mass": float(100*covered),
            "top6_percent_of_total_pair_mass": float(100*covered/sum(freq)**2),
            "top18_joint_percent_raw_mass": 100*covered_joint,
            "top18_joint_percent_of_total_pair_mass": 100*covered_joint/float(sum(freq))**2,
            "top18_joint_sum_of_individually_rounded_percentages": sum(round(100*p[2], 2) for p in joint_pairs[:18]),
            "top6": [{"parents": [p[0], p[1]], "percent": float(100*p[2])} for p in pairs[:6]],
            "top18_joint": [{"parents": [p[0], p[1]], "percent": 100*p[2]} for p in joint_pairs[:18]],
        }
        for system, ps in (("ABO", pairs), ("ABO+Rh", joint_pairs)):
            for rank, (p1, p2, w) in enumerate(ps, 1):
                rows.append({"input_variant": name, "system": system, "rank": rank,
                             "parent1": p1, "parent2": p2, "percent": float(w*100)})
    joint_f = {g: float(x)*r for g, (x, r) in zip(JOINT, ((x, r) for x in exact for r in rh))}
    discrepancies = []
    for rank, (p1, p2, printed) in enumerate(PAPER_TOP18, 1):
        exact_percent = 100*joint_f[p1]*joint_f[p2]*(1 if p1 == p2 else 2)
        discrepancies.append({"rank": rank, "parents": [p1, p2], "printed_percent": printed,
                              "recomputed_percent": exact_percent,
                              "difference_percentage_points": printed-exact_percent})
    return {"abo_alleles": {"A": float(a), "B": float(b), "O": float(o)},
            "rh_negative_phenotype": .07, "rh_d_allele": qd, "rh_genotypes": rh,
            "variants": summary, "paper_top6_printed_sum_percent": 83.62,
            "paper_top18_printed_row_sum_percent": sum(F(str(x[2])) for x in PAPER_TOP18),
            "paper_top18_claimed_total_percent": 67.83,
            "paper_top18_vs_exact": discrepancies,
            "exact_abo_implied_phenotype_percent": [float(100*x) for x in phenotypes(exact)]}, rows


def _json_default(value):
    if isinstance(value, F):
        return float(value)
    raise TypeError(type(value).__name__)


def run():
    OUT.mkdir(parents=True, exist_ok=True)
    sickle = legacy_sickle("AS", "AS", 8)
    unmodified = [sickle[0]]
    for _ in range(7):
        unmodified.append(matvec(M_SICKLE, unmodified[-1]))
    examples = [
        ("abo_p35", 35, "AA", "OO", 3, False,
         [25, 50, 0, 0, 0, 25], [75, 0, 0, 25]),
        ("abo_p36", 36, "AO", "BO", 4, False,
         [6.25, 25, 6.25, 25, 12.5, 25], [31.25, 31.25, 12.5, 25]),
        ("aborh_p49", 49, "AO/Dd", "AO/Dd", 3, True,
         [6.25,12.5,6.25,12.5,25,12.5,0,0,0,0,0,0,0,0,0,6.25,12.5,6.25],
         [56.25,18.75,0,0,0,0,18.75,6.25]),
        ("aborh_p50_top", 50, "AA/dd", "BO/Dd", 7, True,
         [1.56,9.38,14.06,1.56,9.38,14.06,.39,2.34,3.52,.78,4.69,7.03,1.56,9.38,14.06,.39,2.34,3.52],
         [21.88,28.12,8.20,10.55,10.94,14.06,2.73,3.52]),
        ("aborh_p50_bottom", 50, "AB/dd", "OO/Dd", 8, True,
         [.39,2.34,3.52,1.56,9.38,14.06,.39,2.34,3.52,1.56,9.38,14.06,.78,4.69,7.03,1.56,9.38,14.06],
         [13.67,17.58,13.67,17.58,5.47,7.03,10.94,14.06]),
    ]
    records, trajectory_rows = [], []
    for name, page, dad, mom, n, joint, printed_g, printed_p in examples:
        history = legacy_random_mating(dad, mom, n, joint=joint)
        final_g = [float(x*100) for x in history[-1]]
        final_p = [float(x*100) for x in phenotypes(history[-1], joint)]
        record = {"id": name, "page": page, "parents": [dad,mom], "generation": n,
                  "genotype_order": JOINT if joint else ABO, "computed_genotype_percent": final_g,
                  "screenshot_genotype_percent": printed_g, "computed_phenotype_percent": final_p,
                  "screenshot_phenotype_percent": printed_p,
                  "phenotype_order": ["A+","A-","B+","B-","AB+","AB-","O+","O-"] if joint else ["A","B","AB","O"],
                  "matches_screenshot_within_display_rounding": all(abs(a-b)<=.005000001 for a,b in zip(final_g+final_p,printed_g+printed_p)),
                  "history_mass": [float(sum(x)) for x in history],
                  "history_exact_fractions": [[str(v) for v in x] for x in history]}
        records.append(record)
        for generation, distribution in enumerate(history, 1):
            for g, x in zip(JOINT if joint else ABO, distribution):
                trajectory_rows.append({"example": name, "generation": generation, "genotype": g,
                                        "probability": float(x), "exact_fraction": str(x)})
    coverage, pair_rows = coverage_audit()
    self_cols = [punnett18(g, g) for g in JOINT]
    self_mismatch = sum(M_JOINT_DISPLAYED[r][c] != self_cols[c][r] for r in range(18) for c in range(18))
    ranked_cols = [punnett18(a,b) for a,b,_ in PAPER_TOP18]
    ranked_mismatch = sum(M_JOINT_DISPLAYED[r][c] != ranked_cols[c][r] for r in range(18) for c in range(18))
    self_history = [records[2]["history_exact_fractions"][0]]
    self_vector = [F(x) for x in self_history[0]]
    for _ in range(2):
        self_vector = matvec(M_JOINT_DISPLAYED, self_vector)
    bogus_abo = punnett("AA", "OO")
    for _ in range(2):
        bogus_abo = matvec(M_ABO_DISPLAYED, bogus_abo)
    result = {
        "evidence_status": "Executed Python translation of visually inspected MATLAB screenshots; MATLAB not executed.",
        "source": {"filename": SOURCE.name, "pages_read": 54,
                   "sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(), "bytes": SOURCE.stat().st_size},
        "sickle_p28": {"parents": ["AS","AS"], "genotype_order": SICKLE,
                       "source_code_pages": [26,27], "screenshot_generation": 3,
                       "screenshot_percent": [75,18.75,0],
                       "history_exact_fractions": [[str(v) for v in row] for row in sickle],
                       "history_percent": [[float(100*v) for v in row] for row in sickle],
                       "history_mass": [float(sum(row)) for row in sickle],
                       "matrix_without_removal_history_percent": [[float(100*v) for v in row] for row in unmodified],
                       "screenshot_matches": sickle[2] == [F(3,4),F(3,16),F(0)]},
        "screenshot_examples": records, "coverage": coverage,
        "state_counts": {"ABO_genotypes":6,"ABO_ordered_parent_pairs":36,"ABO_unordered_parent_pairs":21,
                         "ABO_distinct_only_unordered_pairs":15,"ABORh_genotypes":18,
                         "ABORh_ordered_parent_pairs":324,"ABORh_unordered_parent_pairs":171},
        "displayed_matrix_checks": {
            "joint_table_vs_self_cross_entry_mismatches_of_324": self_mismatch,
            "joint_table_vs_ranked_top18_entry_mismatches_of_324": ranked_mismatch,
            "joint_table_all_column_sums_one": all(sum(row[c] for row in M_JOINT_DISPLAYED)==1 for c in range(18)),
            "p49_if_displayed_self_cross_matrix_iterated_gen3_percent": [float(100*x) for x in self_vector],
            "p35_if_displayed_abo_matrix_iterated_gen3_percent": [float(100*x) for x in bogus_abo],
        },
    }
    assert result["sickle_p28"]["screenshot_matches"]
    assert all(x["matches_screenshot_within_display_rounding"] for x in records)
    assert self_mismatch == 0
    (OUT / "version1_reproduction.json").write_text(json.dumps(result, indent=2, default=_json_default)+"\n", encoding="utf-8")
    for name, rows in (("version1_pair_coverage.csv", pair_rows), ("version1_trajectories.csv", trajectory_rows)):
        with (OUT/name).open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
            writer.writeheader(); writer.writerows(rows)
    for name, matrix, labels, columns in (
        ("version1_joint_displayed_matrix.csv",M_JOINT_DISPLAYED,JOINT,JOINT),
        ("version1_abo_displayed_matrix.csv",M_ABO_DISPLAYED,ABO,
         [f"{a} x {b}" for a,b in ABO_TOP6_PAIRS]),
    ):
        with (OUT/name).open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(["offspring / columns as printed",*columns])
            for label,row in zip(labels,matrix): writer.writerow([label,*map(float,row)])
    print(json.dumps({"screenshot_examples_matched": len(records)+1,
                      "sickle_gen3_mass": result["sickle_p28"]["history_mass"][2],
                      "coverage": {k:{f:v[f] for f in ("top6_percent_raw_mass","top18_joint_percent_raw_mass")} for k,v in coverage["variants"].items()},
                      "joint_printed_row_sum": coverage["paper_top18_printed_row_sum_percent"],
                      "displayed_matrix_checks": result["displayed_matrix_checks"]},indent=2,default=_json_default))
    return result


if __name__ == "__main__":
    run()
