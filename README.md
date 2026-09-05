# Genetic Inheritance Modelling

Simulating how genotype frequencies evolve across generations, for sickle-cell
anaemia and for the ABO × Rh blood-group system, using Punnett-square logic
expressed as a linear-algebra problem.

Implemented independently in **Python** and **MATLAB** so the two can be checked
against each other.

### Write-ups

- [`report/IB-Mathematics-Extended-Essay.pdf`](report/IB-Mathematics-Extended-Essay.pdf)
  (29 pages) — the mathematical derivation. Sets up the problem, develops the matrix
  formulation, and covers the linear algebra it needs (matrix inverses, determinants,
  eigenvalues). This is where the modelling argument is actually made.
- [`report/Bioinformatics-Arshyia-Mehran.pdf`](report/Bioinformatics-Arshyia-Mehran.pdf)
  (54 pages) — the computational write-up: implementation, simulation results and
  discussion.

## Version II — removing the square-matrix constraint

[`version2/`](version2/) is a 2026 research extension that treats everything above as
**Version I**: evidence to be audited, not ground truth. The original PDFs, MATLAB and
Python are untouched.

> **Status:** preprint, author-approved for release, not peer reviewed. Availability, AI-assistance
> disclosure and the scope of claims are stated in the paper's front matter.

**The single unified paper is [`version2/output/pdf/Genetics_Complete.pdf`](version2/output/pdf/Genetics_Complete.pdf)**
(75 pages, [Markdown source](version2/manuscript/genetics_unified.md)). It merges Version I
and Version II into one document, in six parts: the original biology and matrix algebra with
every derivation preserved; sickle cell, then ABO, then ABO x Rh, each given its Version I
treatment and then audited; the combinatorial theorem and the complete kernel that replaces
the square matrix; real-data evidence and an explicit comparison against the original
square-matrix baseline; age, mutation and epigenetic regulation; then evaluation.
**Part V adds eye colour, height and P(n).** Eye colour is built as a two-locus epistatic trait on
real rs12913832 genotype calls for 2,504 individuals across 26 populations, run through the same
kernel with four representations agreeing to 2.78e-17, audited for Hardy-Weinberg, and compared
against two published cohorts. Its density is 0.3086, so the sparse architecture earns nothing at
this size - a useful negative result. Under neutral random mating the distribution reaches its
Hardy-Weinberg fixed point in one generation and never moves again: nothing decays. The declared
model implies 25% of HERC2 GG individuals are not blue where a phenotyped cohort of 5,481 reports
33%, and matching that requires fitting a parameter to it. A published six-SNP model reports AUC
0.96; this work reports none, because it has no phenotypes. Height is then treated derivationally,
and P(n) states the general case.

**Complexity is now derived throughout, with lower bounds.** Each model carries its own cost note,
Section 10 compares derived growth against measured growth, and Section 15 gives output-size lower
bounds: materialising the complete kernel is Omega(15^n) for *any* implementation, returning all
children of a heterozygous pair is Omega(3^n), and the implemented factored query is Theta(n) and
therefore asymptotically optimal. The derived exponent is confirmed at the top of the measured
range - supported transitions rise 14.878x from four to five loci, measured construction time rose
15.01x, 15.09x and 14.68x - while a single exponent fitted across all of n = 1..5 returns 8.5-10.4,
because fixed overhead flattens the small-n slope.

Appendix E reserves a documented slot for further polygenic traits, and **Appendix F reproduces
Version I's five MATLAB programs verbatim** - including the `X_next(3) = 0;` line that Section 4.4
audits, so the finding can be checked against the code rather than a description of it. Version I's
own bibliography is now carried in the ledger and cited at point of use; its two secondary population
summaries are named as the origin of the phenotype inputs behind the coverage figures.

The [Version II review manuscript](version2/output/pdf/Genetics_Version_II.pdf)
([source](version2/manuscript/version2_manuscript.md)) remains as the standalone audit.
The research question:

> Can sparse and dynamically allocated data structures overcome the combinatorial
> limitations of dense transition matrices in probabilistic genetic inheritance
> modelling, while preserving a greater proportion of biologically possible outcomes?

The honest answer turned out to be more interesting than the expected one.

**The audit changed the premise.** An exact-rational reconstruction reproduces six of
Version I's displayed examples — and finds that the ABO and ABO × Rh MATLAB listings
already loop over *every* parental pair (36 and 324 ordered pairs). The 83.62% and
67.83% figures describe the prose's truncated tables, not what the code executed. So
Version II cannot claim to restore coverage those programs never lost. It measures
better representations of the complete operation instead.

**What the numbers actually are.** The paper's exact allele frequencies give
83.4764828157% top-six ABO pair mass; the 83.62% figure comes from coarse inputs
summing to 1.001. The joint top-18 mass is 67.8347334921%. The legacy sickle-cell
generation-three output is reproducible but sums to 93.75% — it zeroes SS without
renormalising, which is the defect the corrected population model separates into
selection, mating and normalisation.

**The combinatorics, proved.** G = 3ⁿ genotypes by induction, U = G(G+1)/2 unordered
parental pairs, and exactly (15ⁿ + 5ⁿ)/2 supported transitions — generalised to
arbitrary allele counts. Density falls, but the nonzero count is still exponential,
so sparsity is not a solution to polygenic enumeration.

**Four representations, measured.** Dense, CSR, hash adjacency and streamed, all
sharing one transmission kernel. At five loci CSR holds 4,575,976 payload bytes
against dense's 57,631,824 (12.59×), with median validated updates of 9.535 ms and
12.560 ms. Every implementation agreed to within 6.94e-17. Raw per-worker timings and
hardware metadata are retained — these are one-machine results, not universal claims.

**Biology, staged.** rs334 calls for 2,504 samples across 26 populations drive an
exact Hardy-Weinberg and pooling audit (pooled expectation 1.874 rare homozygotes vs
7.600 fitted per-population). Then simplified M/N, ABO-FUT1 epistasis, phased
two-locus recombination, and a 200-locus additive score whose nominal 90% interval
covered 90.03% of 10,000 simulated draws. None of this is clinical prediction, and
the manuscript says so repeatedly.

**What it deliberately does not claim.** Tries, DAGs, priority queues and general
factor-graph inference are documented as candidates, not delivered implementations.
No predictive-accuracy claim is made anywhere — population coverage is never called
accuracy. An [adversarial review](version2/research/adversarial_review.md) records
four numerical defects found and fixed, plus the limits that remain open.

**Section 14 — where age enters (derivation only, nothing implemented).** Both versions
have carried the same gap: neither model has a clock. Version I named it — "random mutations
can occur in an individual before they breed." Section 14 separates the two routes by which
age acts, which are routinely conflated. The *germline* route (parental age at conception →
de novo mutation) modifies the transmission kernel and is heritable; the *somatic and
epigenetic* route (an individual's own age → mutation load, DNA methylation, histone
acetylation) modifies whether an inherited genotype is expressed, and is not.

It derives an age-parameterised gamete operator that drops into the existing `mutation()`
interface, an expression gate `pi = sigma(b + w_age·t + w_met·m(t) + w_ac·h(t) + ...)` with
its cross-entropy loss, gradient, convexity proof and step-size bound, and a two-state
methylation chain solved to an exponential approach to equilibrium. Anchored to measured
slopes: 1.51 de novo mutations per year of paternal age, 0.37 maternal, ~40 somatic
mutations per year in adult stem cells.

Two results are worth the section on their own. **Mutation destroys the sparsity the whole
architecture rests on** — with any positive mutation rate every structural zero becomes
positive, the kernel goes from 381,250 to 7,203,978 stored entries at five loci, and CSR
becomes strictly worse than dense. And **age is collinear with its own mediators**, so
`w_age` and `w_met` are individually unidentifiable while the fit still looks fine — the
conditioning has to be checked before fitting, never inferred from the residual afterwards.

Nothing in Part V is implemented, tested or validated, and it says so on every page.

**Two limitations are listed in the paper's Findings rather than buried.** A positive mutation
rate destroys the structural sparsity the architecture depends on: every structural zero becomes
positive, the five-locus kernel goes from 381,250 to the full dense 7,203,978 entries, and CSR
becomes strictly worse than dense. And age is collinear with its own epigenetic mediators, so the
expression gate's weights are individually unidentifiable while the objective still looks healthy
- conditioning has to be checked before fitting, never inferred from the fit afterwards.

```bash
cd version2
python -m pip install -r requirements.txt
python run.py test -q          # 32 tests
python run.py reproduce        # exact-rational Version I reconstruction
python run.py science          # HWE, epistasis, linkage, polygenic score
python run.py benchmark        # dense / CSR / hash / streamed
python run.py unified          # rebuild the unified paper + PDF
```

Version II addresses three limitations listed below: mutation and migration now exist
as explicit operators, and selection is a general fitness weighting rather than
complete reproductive exclusion. Drift, genuine polygenic prediction and an
independent parent-offspring validation cohort remain open.


## The idea

A Punnett square answers a one-generation question: given two parent genotypes, what
are the offspring probabilities? Asking it across *n* generations by hand is
hopeless — the branching is exponential.

The modelling step that makes it tractable is to stop treating a Punnett square as a
diagram and start treating it as a **column of a transition matrix**. Write the
population as a vector of genotype proportions, and one generation of reproduction
becomes a matrix multiplication:

```
X_{n+1} = M · X_n        so        X_n = M^n · X_0
```

An exponential branching process collapses into a matrix power. For sickle cell, with
genotypes ordered (HbA/HbA, HbA/HbS, HbS/HbS):

```
        AA×AA   AA×AS   AS×AS
  AA  [  1.0     0.5     0.25  ]
  AS  [  0.0     0.5     0.50  ]
  SS  [  0.0     0.0     0.25  ]
```

Each column is one Punnett square, and the columns sum to 1 because a cross must
produce *some* offspring. That is a Markov chain, and everything the linear algebra
knows about Markov chains — powers, limiting behaviour, eigenvectors — applies.

## What each model does

### Sickle cell (`python/gene_modelling.py`, `matlab/sickle_cell_model.m`)

Three genotypes, propagated over *n* generations from a chosen parent pair, with
genotype and phenotype trends plotted against generation number.

Modelling assumptions, stated explicitly in the source:

- Alleles HbA and HbS, giving three genotypes.
- A non-malarial population (Canada), so the heterozygote carries no survival
  advantage — which is exactly the assumption that would have to be *dropped* to
  reproduce the observed HbS frequencies in malarial regions.
- **HbS/HbS individuals survive but do not reproduce.** This is a selection term, and
  it is the interesting part: it means affected individuals are counted in each
  generation's output but removed from the breeding pool feeding the next. Generation 1
  is computed directly from the parents; generations 2 onward use the fixed matrix.

### ABO × Rh (`python/abo_rh_model.py`, `matlab/aborh_run.m`, `matlab/bloodgroup18_simulation.m`)

Extends the same machinery to two independently inherited loci at once:

- **ABO** — 6 genotypes (AA, AO, AB, BB, BO, OO) with codominance of A and B and
  recessive O, so 6 genotypes collapse to 4 phenotypes.
- **Rh** — 3 genotypes (DD, Dd, dd), simple dominance, 2 phenotypes.

Combined, 6 × 3 = **18 genotypes**, and therefore an 18 × 18 offspring transition
matrix — [`data/18_18_offspring_punnett_matrix.csv`](data/18_18_offspring_punnett_matrix.csv),
which is the object the whole model is built around. `matlab/abo_simulation.m` and
`matlab/abo_run.m` are the single-locus ABO version that came first.

## Running it

```bash
python python/gene_modelling.py
```

Requires NumPy and matplotlib. The MATLAB scripts prompt for parent genotypes and a
generation count; run them from the `matlab/` directory.

## Known limitations

- **Infinite-population, deterministic.** The model propagates expected proportions,
  so it has no genetic drift. Real small populations fluctuate, and drift is the
  dominant force precisely where these models are most interesting.
- **Random mating assumed throughout.** No assortative mating, no population
  structure, no consanguinity.
- **No mutation and no migration**, so allele frequencies change only through the
  selection term.
- **Loci assumed independent** in the ABO × Rh model — reasonable here, since ABO
  (chromosome 9) and Rh (chromosome 1) genuinely do assort independently, but the
  assumption is doing work and is not free in general.
- Selection is modelled only as complete reproductive exclusion of HbS/HbS. Partial
  fitness costs would be a more realistic and more general formulation.

## Citing and licence

> **Preprint** — author-approved for release, not peer reviewed. A DOI will be minted on
> publication; see [`PUBLISHING.md`](PUBLISHING.md) for the steps and
> [`CITATION.cff`](CITATION.cff) for machine-readable citation metadata.

Dual-licensed, see [`LICENSE`](LICENSE):

- **Written work and figures** — CC BY 4.0. Reuse and adapt freely with attribution.
- **Source code** — MIT.
- **Frozen genotype datasets** in `version2/data/` are aggregate counts from the 1000 Genomes
  phase 3 panel via Ensembl. They keep the terms of their originating projects and are **not**
  relicensed here. No third-party figure is reproduced anywhere in this repository.

Build the full research package (paper, code, data, results, figures) with:

```bash
cd version2 && python run.py package
```

## Provenance

Consolidated in August 2026 from four separate repositories
(`SickleCellAnemiaSimulationOverGenerations`, `ABO-RhBloodroupSimulationOverTime`,
`SimulatingBloodGroupInheritanceOverGenerations`,
`ComputationalBiologyGeneExpressionAcrossGenerationsAlgorithm`), which were four views
of one piece of work. The sickle-cell MATLAB model came first; the Python
reimplementation and the two-locus extension followed.
