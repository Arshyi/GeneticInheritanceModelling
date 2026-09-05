# Genetics Version II

A research extension of **Bioinformatics-Arshyia Mehran.pdf**, preserved unchanged in the parent directory. The main result is a reproducible audit and complete probabilistic inheritance engine, with measured dense/CSR/hash/streamed comparisons and explicitly limited biological validation.

**Start with [the unified paper](output/pdf/Genetics_Complete.pdf)** (42 pages, [Markdown source](manuscript/genetics_unified.md)) - Version I and Version II merged into one document: the original mathematics preserved and reworked for exposition, then audited, then extended. The [Version II review manuscript](output/pdf/Genetics_Version_II.pdf) ([source](manuscript/version2_manuscript.md)) remains available as the standalone audit. Scientific methods and limitations are in the manuscript. No public upload or submission has been made.

## Findings

- **The original ABO and ABO+Rh MATLAB listings already enumerate all parental pairs.** Their actual support must not be represented by the prose's 83.62%/67.83% truncation figures.
- Exact paper allele frequencies give **83.4764828157%** top-six ABO pair mass; coarse unnormalized inputs explain **83.62%**. The joint top-18 mass **67.8347334921%** reproduces 67.83%.
- The legacy sickle-cell screenshot is reproducible, but its generation-three percentages sum to **93.75%**. The corrected population model separates selection from mating and normalizes breeding mass.
- At five biallelic loci, dense payload is **57,631,824 bytes**, CSR **4,575,976 bytes**. Median validated updates were **12.560 ms** and **9.535 ms**, respectively, in the retained local run. All measured implementation results agreed within **6.94e-17**.
- Public rs334 calls from **26 populations / 2,504 samples** support an HWE and pooling audit. They do not validate clinical prediction or worldwide prevalence.
- A **200-locus synthetic additive model** has 401 score bins. Its nominal 90% interval covered **90.03%** of 10,000 independent simulated outcomes. It is not a human height predictor.
- **LIMITATION - mutation destroys the sparsity the architecture rests on.** With any positive
  per-transmission mutation probability every structural zero becomes positive. At five loci the
  supported-transition count rises from **381,250 to the full dense 7,203,978**, and CSR becomes
  strictly worse than dense - it pays index overhead on a matrix with no zeros. The measured
  12.59x payload advantage is an advantage over a mutation-free model. This argues for the factored
  representation or a structured low-rank correction; neither is implemented.
- **LIMITATION - age is collinear with its own mediators.** In the proposed expression gate the
  methylation fraction is itself a smooth deterministic function of age, so the design carries two
  nearly collinear columns. The age and methylation weights become individually unidentifiable
  **while the objective still looks healthy**, and no diagnostic computed from the fit reveals it.
  Conditioning must be checked before fitting, and the sampling design must supply within-age
  variation in methylation.

## Run

Python 3.12 was used. On this computer, `reproduce.ps1` defaults to the installed bundled interpreter and `run.py` finds the optional project-local `.deps`. To use another environment:

```powershell
cd C:\Users\DELL\Desktop\GENETICS\version2
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe run.py test -q
.\.venv\Scripts\python.exe run.py reproduce
.\.venv\Scripts\python.exe run.py fetch
.\.venv\Scripts\python.exe run.py science
.\.venv\Scripts\python.exe run.py benchmark
.\.venv\Scripts\python.exe run.py test -q --junitxml=results/tests.xml
.\.venv\Scripts\python.exe run.py manuscript
.\.venv\Scripts\python.exe run.py unified
```

The PDF renderer currently uses Windows Georgia and Arial fonts. Poppler `pdftoppm` was used for final visual review. PDF source extraction used `pypdf`. Tests and computations do not need Word or MATLAB. The package is primarily a Python research implementation, not a MATLAB execution or a C/C++ engine.

The local dependencies installed during this session inherited permissions that the Codex sandbox could not read. Approved execution outside that sandbox ran the tests and experiments successfully; this was an environment access issue, not a missing `pytest` or `scipy` module. A fresh normal-user virtual environment avoids relying on those inherited directory permissions. `run.py` prefers `.deps` when that directory exists; remove that optional directory from a copied workspace or rename it before testing a new environment if you want to guarantee use of the new environment's packages.

Fetch defaults to offline validation against the frozen provenance. Use `python run.py fetch --live` for a read-only refresh comparison. Saving a refresh requires a **new output directory**; the script refuses to overwrite the frozen data. API and source changes may produce a deliberately rejected mismatch.

## Code and evidence

| Location | Contents |
|---|---|
| `genetics/core.py` | Canonical catalogs, local Mendelian factors, dense/CSR/hash kernels, streaming, log queries, population operators |
| `genetics/extensions.py` | Two-locus phase/recombination, integer-score DP, exact HWE test, explicitly factorized population updates |
| `tests/` | Independent allele-copy oracle, legacy regression, numerical/resource tests |
| `experiments/reproduce_version1.py` | Exact rational reconstruction with legacy state order |
| `experiments/benchmark.py` | Isolated workers, raw repeated timings, distinct memory metrics |
| `experiments/science.py` | HWE, synthetic dynamics, MNS, ABO/FUT1, linkage and quantitative-trait checks |
| `experiments/build_unified.py` | Unified Version I + II manuscript and PDF, bound to retained results |
| `manuscript/unified_template.md` | Source of the unified paper; edit here, not the generated file |
| `data/` | Frozen aggregate calls and provenance, no personal patient data |
| `results/` | Machine-readable outputs, tests, raw benchmark workers, metadata |
| `figures/` | PNG and SVG figures |
| `sources/source_ledger.json` | Consolidated reproducible source ledger |
| `research/` | Detailed original audit, architecture, biological evidence, adversarial review |

`benchmark_initial_memory_import_included.*` preserves an earlier run whose CSR memory pass included a lazy library import. It is superseded by **benchmark.json / benchmark.csv**, whose memory pass excludes that import. Do not combine measurements from the two revisions. Numeric-buffer bytes, recursive Python object sizes, traced allocations and absolute process peaks measure different things.

## Scientific contract

No threshold prunes positive transmission branches. Structural zeros are omitted in sparse representations. Conditional log queries preserve probabilities too small for ordinary binary64 output; selected guarded population paths explicitly reject zero underflows. This is not universal log-space population inference. General factor-graph/pedigree inference, adaptive radix trees and priority queues are discussed candidates, not delivered implementations. Full sparse enumeration remains exponential.

Counts are complete within the declared simplified catalogs. Full molecular ABO/Rh variation, all HBB disease genotypes, sex-linked inheritance, large linked haplotypes, measured mutation/selection parameters and externally calibrated human-height prediction remain outside the validated implementation. See the manuscript's staged experiments and limitations before extending the models.

## Review and preservation

Original PDF SHA-256: `EE3CB04CE28D4734669D5822FDDBDF03C176FD5C2F5FF3760CB82A4BD1145ECD`.

The supplied PDF and DOCX were not edited. Integration and publication remain the author's decision. External dataset/reference reuse conditions are recorded in the source ledger; the source archive does not relicense external materials or include the optional scientific dependencies.
