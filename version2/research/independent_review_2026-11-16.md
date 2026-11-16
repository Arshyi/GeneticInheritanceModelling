# Independent review: what went well, what can be improved, even better if

**Date:** 16 November 2026
**Reviewed:** this repository at commit `49f55a2` (version 2.1.0, as deposited on Zenodo):
the README, the 14-page bounds paper (`manuscript/genetics_bounds_paper.md`), the relevant
parts of the unified paper, `genetics/`, `experiments/pgs_catalog.py`, the retained
`results/`, and `research/adversarial_review.md`.
**Reviewer:** Claude, an AI assistant made by Anthropic, at the author's request.

## How this review was done

- The documents and code above were read.
- Five checks were run with [`independent_review_checks.py`](independent_review_checks.py),
  which uses the Python standard library only and imports nothing from this package. Its
  complete output is kept in
  [`independent_review_checks_output.txt`](independent_review_checks_output.txt).
- Section 2 of that script reads the retained `results/pgs_catalog.json`; the other sections
  are independent re-computations.

What was not done: the package's own tests and benchmarks were not rerun (the pinned
scientific stack was not installed), PGS scoring files were not downloaded again, and the
literature named in B5 was located by search but not read in full.

## Summary

1. The counting theorems are correct. A brute-force enumeration reproduces them for twelve
   allele configurations.
2. The central lower bound, Omega(15^n), is true for materialising the kernel but is
   presented as a limit on inheritance modelling. One generation of random mating does not
   need the kernel: going through gamete frequencies gives the same result to 1e-16 and
   handles ten loci in 0.05 seconds, where the kernel would hold 288 billion entries.
3. The polygenic-score headline numbers (1,443 standard deviations of error,
   "worthless past about a thousand variants", Theta(n^(5/2))) are a worst-case bound, not
   measured error. The repository's own retained results show the actual error in the mean
   is 0.01 SD at 1,059 variants and 0.02 SD at 3,892. The fitted exponent 1.4987 is fitted
   to the bound formula's own output.
4. "Tested on 88 published scores" describes 8 downloaded scoring files, 6 of them with an
   assumed allele frequency.
5. The work has no biological result yet, and literature that covers its framing is not
   cited.

## What went well

1. **Version I is audited, not defended.** The exact-rational reconstruction finds that the
   83.62 and 67.83 per cent coverage figures describe truncated tables in the prose rather
   than what the MATLAB executed, and that the sickle-cell generation-three output sums to
   93.75 per cent because `X_next(3) = 0` removes mass without renormalising. Reporting
   that your own earlier work did not lose what you set out to restore is rare and
   creditable.
2. **The combinatorics are right.** Check 1 enumerates every unordered parental pair and
   its distinct children for twelve systems (one to five loci, two to five alleles) and
   matches G, U and (T + D)/2 exactly, including 615 for ABO x Rh and 381,250 for five
   biallelic loci.
3. **Proved and measured are kept apart.** Section 7 of the bounds paper states which
   results transfer and which are one machine's numbers, and the growth-exponent discussion
   correctly argues that an asymptotic claim is entitled only to the top of the measured
   range.
4. **The numerical contract is careful.** No probability threshold prunes rare branches,
   positive results that underflow raise instead of returning zero, and log-space queries
   exist. The adversarial review found and fixed real defects of this kind.
5. **Four representations share one semantic kernel and agree to 6.94e-17**, with raw
   per-worker timings and hardware metadata retained.
6. **Data handling is responsible.** 1000 Genomes counts are frozen with digests and
   provenance, PGS files are not redistributed, licences are separated, and AI assistance
   is disclosed.
7. **The limitations are stated in the papers themselves**: linkage equilibrium, no
   phenotypes, no predictive accuracy, Section 14 as derivation only.

## What can be improved

### A. Main findings

#### A1. The Omega(15^n) wall belongs to one data structure, not to the problem

Theorem 6(a) is correct: writing out the complete kernel costs one value per supported
transition. The paper then concludes that "no engineering removes the exponent" and that
"the architecture moves the wall twice and removes it never". That treats the kernel as
necessary. It is not.

Under independent segregation a parent's gamete depends only on that parent, so one
generation of random mating is: genotype distribution -> gamete (haplotype) distribution ->
product of two independent gametes. Check 3 implements this, applied to a random
*correlated* joint distribution (so it is not assuming linkage equilibrium):

| Loci | Kernel transitions | Gamete route, operations | Largest difference |
|---:|---:|---:|---:|
| 5 | 381,250 | 1,868 | 3.1e-17 |
| 8 | 1,281,640,625 | 90,756 | 2.6e-18 (against the Fourier variant) |
| 10 | 288,330,078,125 | 1,280,676 (0.05 s) | (kernel not run) |

The gamete route costs O(n 2^n + 4^n), and with a three-point Fourier transform on each
locus axis it costs O(n 3^n), which is linear in the size of the state itself. The kernel
answers "what are all children of every pair", which no population update needs.

The same observation weakens the Section 14 result that mutation destroys sparsity: in the
gamete route mutation is a 2 x 2 matrix applied per locus to the gamete distribution, and
nothing that the route relies on becomes dense.

*Suggested change.* State Theorem 6(a) as a bound on materialising the kernel, add the
gamete-route update as a theorem (it is a stronger result than the bound), and revise the
"wall" sentences in the abstract, Section 7, the unified paper and the README.

#### A2. The polygenic-score error figures are a bound, not a measurement

`discretise()` returns `worst_case = n * delta`, and every "error in SD" in Section 6 and in
the README is that bound divided by the score's SD. The run also stores the DP's own mean
and variance, so the actual error can be read from `results/pgs_catalog.json` (check 2):

| Score | Variants | Reported error (SD) | Actual error in the mean (SD) | DP SD / true SD |
|---|---:|---:|---:|---:|
| PGS000010 | 27 | 0.012 | 0.0019 | 1.000 |
| PGS002775 | 1,059 | 2.240 | 0.0103 | 1.001 |
| PGS004196 | 3,892 | 7.884 | 0.0225 | 1.001 |
| PGS004197 | 11,490 | 45.947 | 1.4883 | 1.004 |
| PGS000012 | 49,310 | 624.297 | 34.4877 | 1.089 |
| PGS000337 | 75,028 | 1,443.650 | 307.0824 | 0.779 |

Rounding errors on independent weights mostly cancel, so the typical error grows like
sqrt(n) times the step, not n times the step. The distributions are good to a few hundredths
of an SD up to about 4,000 variants, and break down past roughly 10,000, where the step
exceeds the typical weight and most weights round to zero (at least 69 and 95 per cent of
them for the two largest scores, from the stored bin counts). That breakdown is real and
worth reporting, but it is a different and later failure than "worthless past about a
thousand variants".

Consequences:
- **Theorem 7** is correct as a statement about the bound. The claimed cost Theta(n^(5/2))
  is the cost of guaranteeing that bound with a uniform grid, not the cost of computing the
  distribution accurately.
- **The controlled study's exponent 1.4987** is fitted to `bins_required`, which is the
  bound formula evaluated at each n. Only n = 50 to 400 were run. The fit confirms the
  algebra of the formula, not a measured scaling.
- `accuracy['achieved_error_in_sd']` stores the bound under a name that says "achieved".
  At n = 27, 169 and 540 the accuracy-controlled runs with up to 1.66 million bins have
  actual mean errors of 6e-4, 3e-4 and 5e-5 SD; the fixed 40,001-bin runs already had
  2e-3, 9e-5 and 2e-3 SD.
- **Exact weights need no discretisation.** Check 5 computes the score distribution by
  Fourier inversion of the product of per-variant characteristic functions. Against an
  exact integer-weight DP at 40 loci, quantiles from 0.1 to 99.9 per cent agree to within
  0.0005 SD. On the paper's own synthetic design it handles 25,600 variants, which the
  paper refuses, in about 0.6 seconds with 129 evaluations of the characteristic function,
  and recovers the exact skewness and kurtosis. The paper mentions saddlepoint and normal
  approximations; characteristic-function inversion is exact up to a controllable
  truncation and is the natural comparator.

*Suggested change.* Report the actual moment errors next to the bound, rename the stored
field, restate Result 3, Theorem 7 and Result 4 as results about the uniform-grid bound,
and add a CDF comparison against an exact or Fourier reference so that error is measured
on the distribution itself, not only its mean.

#### A3. "88 published scores" were 8 downloaded files

The catalogue lists 88 scores for the trait; `LADDER` downloads eight, and six of those
use `DECLARED_FREQUENCY = 0.30` because no allele frequency is published. The abstract,
Section 6, the commit message of `b42a82c` and the README say the regimes were tested on 88
scores. Write "eight scores spanning the catalogue's range of 88", and note that the
declared frequency makes those six runs cost measurements only, as the code comment
already says.

#### A4. The README's central picture contradicts the paper

"The idea" presents the 3 x 3 sickle-cell matrix as a Markov chain on genotypes to which
"powers, limiting behaviour, eigenvectors" apply. Its columns are crosses (AA x AA,
AA x AS, AS x AS) and its rows are genotypes, so it does not map a state to a state, and
Proposition 2 of the bounds paper proves that random mating is not linear in the genotype
distribution. A first-time reader of the README learns the framing the papers refute. Mark
"The idea" as the Version I framing and point to Proposition 2.

#### A5. The novelty claim misses the literature that owns this framing

The bounds paper states that a matching pair of upper and lower bounds for the multilocus
transmission kernel is, to the author's knowledge, absent from the literature. Neither
paper cites:
- **Genetic algebras**: Etherington (1939) onward, surveyed by Reed, "Algebraic structure of
  genetic inheritance", Bulletin of the AMS 34 (1997), 107-130. Random mating as a quadratic
  map on the simplex, and its linearisation (Holgate, "Sequences of powers in genetic
  algebras", J. London Math. Soc. 42 (1967), 489-496; later work on "Haldane
  linearisation"), is exactly the subject of Proposition 2.
- **Pedigree likelihood algorithms**: Elston-Stewart (peeling over individuals) and
  Lander-Green (hidden Markov model over loci) are the standard ways the field avoids
  enumerating joint multilocus outcomes.
- **Population-genetic recursions on gamete frequencies** with recombination (Geiringer;
  Bennett), which are the route in A1.

The specific counts in Theorems 4 and 5 may still be new; the framing and the escape from
the exponential are not. Narrow the claim to the counts.

#### A6. There is no biological result yet

Every model is either illustrative or reproduced from a textbook, the eye-colour second
locus is a declared modifier, and no phenotype is predicted. The repository already holds
rs334 genotype counts for 26 populations. Estimating the heterozygote advantage (selection
coefficients) that would maintain each observed HbS frequency under the corrected selection
model, with intervals, and comparing it with malaria endemicity from Piel et al. (already
cited) would give the work a falsifiable biological claim.

### B. Smaller points

1. **The one-generation fixed point is stated too generally.** The README says the
   eye-colour distribution "reaches its Hardy-Weinberg fixed point in one generation and
   never moves again". That holds because the joint distribution was built under linkage
   equilibrium. Check 4 starts from half AABB and half aabb at unlinked loci: the
   association halves each generation and the double heterozygote share goes
   0 -> 0.5 -> 0.3125 -> ... -> 0.25. Say "at each locus" or "for a population already in
   linkage equilibrium".
2. **"On real rs12913832 genotype calls"** is true of the HERC2 locus only; the second
   locus is declared and its frequency is set to 0.5, then calibrated to a cohort figure.
   The summary sentences should say "one real locus and one declared modifier".
3. **The legacy code carries no warning.** `matlab/sickle_cell_model.m` and
   `python/gene_modelling.py` still contain the unrenormalised selection step that the
   audit identifies. A header comment pointing to the audit would stop someone reusing
   them.
4. **The MATLAB is stored twice.** `matlab/` and `version2/legacy/matlab/` are identical
   apart from line endings. Keep one copy, or note why both exist.
5. **Releases are not tagged.** Zenodo holds 2.1.0, but the repository has no tags, so the
   deposited commit cannot be checked out by name. Tag `49f55a2` as `v2.1.0`.
6. **The self-review missed the substantive points.** `research/adversarial_review.md`
   checks numerics thoroughly but did not test the headline claims (A1-A3). A short
   "claims tested against alternatives" section would complete it.
7. **"Linkage equilibrium between variants"** is listed as an assumption in Section 6, but
   for LD-pruned or clumped scores it is close to true and for LDpred-style scores it is
   badly false; the eight scores could be labelled by construction method so readers know
   which runs the assumption affects.
8. **Adding files under `version2/` changes the release manifest.** `verify_release.py`
   hashes every file, so this review's files will appear in `results/artifact_manifest.json`
   the next time it is built.

## Even better if

1. **A new main result:** the gamete-route update (A1) as a theorem with its cost, plus a
   linked-loci version using phased haplotypes and recombination fractions. It turns the
   paper from "where the wall is" into "how to get past it".
2. **Measured error everywhere a bound is quoted**, with the Fourier method as the exact
   reference for score distributions.
3. **One biological estimate from data already in the repository** (A6).
4. **A shorter, sharper journal version** centred on the counts (Theorems 4-5) and the
   gamete-route update, with the genetic-algebra and pedigree literature in related work.
5. **A changelog entry or 2.2.0 deposit** recording which claims changed, so readers of the
   Zenodo 2.1.0 version can find the corrections.

## Reproducing the checks

```text
cd version2
python research/independent_review_checks.py
```

Standard library only. Section 2 reads `results/pgs_catalog.json`; sections 1 and 3 to 5
are independent of the package. The retained output is
`research/independent_review_checks_output.txt`. Synthetic weights in check 5 come from
Python's generator, so they differ from the paper's NumPy draw while following the same
design.

## Limits of this review

- The actual-error figures in A2 are errors in the mean and SD read from stored results;
  tail quantiles of the discretised distributions were not compared directly because the
  PGS files were not re-downloaded.
- The literature in A5 was identified by search and from general knowledge; exact citation
  details should be checked before they are added.
- The package's own tests and benchmarks were not rerun.
