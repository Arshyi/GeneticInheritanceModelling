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

## Provenance

Consolidated in August 2026 from four separate repositories
(`SickleCellAnemiaSimulationOverGenerations`, `ABO-RhBloodroupSimulationOverTime`,
`SimulatingBloodGroupInheritanceOverGenerations`,
`ComputationalBiologyGeneExpressionAcrossGenerationsAlgorithm`), which were four views
of one piece of work. The sickle-cell MATLAB model came first; the Python
reimplementation and the two-locus extension followed.
