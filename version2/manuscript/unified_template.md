# Inheritance as a Linear Operator, and What Replaces It

## From Punnett squares to a complete probabilistic kernel: an audited, measured and extended treatment of sickle-cell, ABO and ABO x Rh inheritance

Parts I to VII | Arshyia Mehran | Version I written 2022-2025, audited and extended {{DATE}}

Parts I and II preserve the mathematics of the original 54-page study and are reworked only for exposition. Parts III to V are new: they audit that study against its own source code, replace the square-matrix constraint with a complete rectangular kernel, measure four implementations of it, test the population assumptions against public genotype calls, add eye colour as a two-locus trait on real genotype data, and state the general n-locus problem with upper and lower bounds. Part VI derives, without implementing, an extension in which age enters through mutation and epigenetic regulation.

**This is a preprint. It has not been peer reviewed.** Every quantitative claim is reproducible from the accompanying code and retained results; the limits of each claim are stated where it is made, and Section 27 collects the counterarguments that survive.

<!-- pagebreak -->

# Abstract

Version I of this work asked whether linear algebra could replace repeated Punnett-square calculation for sickle-cell anaemia, ABO blood groups and a simplified ABO x Rh system. It answered yes, and identified a limitation: forcing a parental-pair-to-offspring table into a square, diagonalisable matrix appeared to sacrifice population coverage as more alleles were added, falling from {{PRINTED_TOP6}}% for ABO to {{PRINTED_TOP18}}% for ABO x Rh.

This unified treatment reproduces that work exactly, audits it, and then removes the constraint. An exact-rational reconstruction matches six displayed examples from the original. The audit changes the premise: the original ABO and ABO x Rh MATLAB programs already loop over every parental pair, 36 and 324 ordered pairs respectively. The truncation to a square matrix exists in the prose and displayed tables, not in the executed code. No claim of restored coverage is therefore available, and none is made.

The coverage percentages are reproduced and their provenance separated. Exact propagation of the stated allele frequencies gives {{EXACT_TOP6}}% top-six ABO pair mass; the printed {{PRINTED_TOP6}}% follows from coarser inputs summing to 1.001. The joint top-18 mass is {{EXACT_TOP18}}%. The displayed sickle-cell trajectory is reproducible but does not conserve probability: its generation-three output sums to 93.75%, because the algorithm zeroes the affected class without renormalising.

For n independent biallelic diploid loci, an induction proves G = 3^n unphased genotypes and U = G(G+1)/2 unordered parental pairs, generalised to arbitrary allele counts; the complete kernel is G by U with exactly (15^n + 5^n)/2 supported transitions. Dense arrays, compressed sparse rows, hash adjacency and streamed computation are implemented with identical transmission semantics and benchmarked. {{BENCH_ABSTRACT}}

Public rs334 genotype calls from 2,504 individuals across 26 populations support an exact Hardy-Weinberg and pooling audit. A staged extension covers simplified M/N, ABO-FUT1 epistasis, phased two-locus recombination, and a 200-locus additive score whose nominal 90% interval covered 90.03% of 10,000 simulated outcomes. Eye colour is then built as a two-locus epistatic trait on real rs12913832 genotype calls for the same 2,504 individuals, run through the unchanged kernel with four representations agreeing to {{EYE_AGREEMENT}}, and compared against two published cohorts. The general case is stated as P(n), with output-size lower bounds showing that the complete kernel is Omega(15^n) and that the implemented factored query is Theta(n) and therefore optimal. A final part derives, without implementing, an extension in which parental age enters the transmission kernel through de novo mutation and an individual's age enters the genotype-to-phenotype map through methylation and histone acetylation.

Two findings limit the contribution and are stated rather than buried: a positive mutation rate destroys the structural sparsity the architecture depends on, and age is collinear with its own epigenetic mediators, so the weights of the proposed expression gate are not separately identifiable without a design built to break that collinearity.

Keywords: Mendelian inheritance; transition matrices; combinatorial growth; sparse representation; Hardy-Weinberg equilibrium; model audit; reproducibility; epigenetic regulation.

# Findings

1. **The original programs were already complete.** The ABO and ABO x Rh MATLAB listings evaluate all 36 and all 324 ordered parental pairs through random-mating loops. {{DISPLAYED_MATRIX_CHECK}} Restoring coverage to those programs is not an available claim.

2. **The two coverage percentages are reproducible, and their provenance matters.** Exact allele frequencies give {{EXACT_TOP6}}% and {{EXACT_TOP18}}%; the printed {{PRINTED_TOP6}}% arises from unnormalised inputs. Both are model-based pair probabilities, not measured shares of any population.

3. **The legacy sickle-cell trajectory does not conserve probability.** Its displayed generation-three output sums to 93.75%. The corrected model separates selection from mating and renormalises the breeding population, giving a reciprocal rather than exponential decline in the S allele.

4. **Sparsity is real but still exponential.** The kernel has exactly (15^n + 5^n)/2 nonzeros against G*U total. Density falls; the nonzero count does not.

5. **Measured representation results.** At five loci CSR retains 4,575,976 payload bytes against dense's 57,631,824, with median validated updates of 9.535 ms and 12.560 ms. All measured implementations agreed within 6.94e-17. One machine, one workload.

6. **External data audits assumptions, not predictions.** rs334 calls for 2,504 individuals across 26 populations give a pooled expectation of 1.874 rare homozygotes against 7.600 from separately fitted populations. No predictive accuracy is established anywhere in this work.

7. **LIMITATION - mutation destroys the sparsity the architecture rests on.** With any positive per-transmission mutation probability, every structural zero of the kernel becomes a positive number. At five loci the supported-transition count rises from 381,250 to the full dense 7,203,978, and CSR becomes strictly worse than dense because it pays index overhead on a matrix with no zeros. The measured 12.59-fold payload advantage is an advantage over a mutation-free model. This argues for the factored representation or a structured low-rank correction; neither is implemented. See Section 24.

8. **LIMITATION - age is collinear with its own mediators.** In the proposed expression gate, methylation fraction is itself a smooth deterministic function of age, so the design matrix carries two nearly collinear columns. The age and methylation weights are then individually unidentifiable while the objective value remains healthy and no diagnostic computed from the fit reveals it. Conditioning must be checked before fitting, and the sampling design must supply within-age variation. See Section 28.

9. **Eye colour transfers, and shows the architecture earning nothing.** A two-locus epistatic model on real rs12913832 calls for 2,504 individuals runs through the unchanged kernel, with four representations agreeing to {{EYE_AGREEMENT}}. Its density is 0.3086, so compressed storage is not worth its overhead at this size. Under neutral random mating the phenotype distribution is already at its fixed point: nothing decays.

10. **The declared eye-colour model is wrong in a measurable direction.** It implies 25 per cent of HERC2 GG individuals are not blue; a phenotyped cohort of 5,481 reports 33 per cent. Reproducing that figure requires fitting a modifier frequency to it, which is a calibration and not a validation. A published six-SNP model reports AUC 0.96 for blue and brown; this work reports no AUC, because it has no phenotypes.

11. **The exponent is a lower bound, not an implementation defect.** Any algorithm that materialises the complete kernel must write one value per supported transition, so kernel construction is Omega(15^n) for any implementation. Returning all children of a heterozygous pair is Omega(3^n) by output size. The implemented factored single-child query is Theta(n) and therefore asymptotically optimal, and the CSR population update is tight against the same bound.

12. **The derived growth is confirmed at the top of the measured range.** From four to five loci the supported-transition count rises 14.878-fold; measured construction time rose 15.01, 15.09 and 14.68-fold for dense, CSR and hash. A single exponent fitted across all of n = 1..5 returns 8.5 to 10.4 instead, because fixed overhead flattens the small-n slope, and reporting only that would understate the growth by about a third.

13. **Not delivered.** Tries, decision DAGs, priority queues and general factor-graph inference are documented candidates, not implementations. Part VI is derivation only: no code, no test, no dataset. Height is treated derivationally and no height prediction is made.

# Reader's guide

The intellectual progression of Version I is preserved: explain the biology, introduce the mathematics that represents it, implement the model, then evaluate what the output establishes. Each part introduces its concept, derives it, uses it, and then evaluates it against something external.

**Part I** covers the biology, the matrix algebra Version I needs, and the complete sickle-cell derivation, followed by its audit. **Part II** does the same for ABO and then ABO x Rh. **Part III** proves the combinatorial growth by induction, generalises it, and builds and measures the complete kernel that removes the square-matrix constraint. **Part IV** tests the population assumptions against public data and compares the new representation with Version I's square matrices under an explicit taxonomy of what "coverage" can mean. **Part V** adds eye colour as a two-locus epistatic trait driven by real genotype calls, runs it through the same kernel, compares it against published cohorts, then treats height and finally P(n), the general problem of a trait depending on n genes, with upper and lower bounds. **Part VI** derives the age, mutation and epigenetic extension. **Part VII** evaluates the whole.

Numbered Results are stated once and referred to afterwards by number. Every benchmark figure comes from retained machine-readable output; theoretical extrapolations are labelled separately. Sections marked *Audit* revisit Version I with evidence and are the only places where its conclusions are revised. Complexity is treated wherever a model is introduced rather than in one block: a short note for each of sickle cell, ABO and ABO x Rh, the measured comparison in Section 10, and the general upper and lower bounds for n loci in Section 15.

<!-- pagebreak -->

# PART I - FOUNDATIONS AND THE SICKLE-CELL MODEL

# 1. Introduction and research question

A Punnett square answers a one-generation question: given two parental genotypes, what are the offspring probabilities? School biology stops at the second filial generation because the branching becomes unmanageable by hand. To reach generation n one must enumerate every offspring genotype, pair each with every possible mate, and draw a new square for each combination.

Version I's insight was to stop treating a Punnett square as a diagram and start treating it as a **column of a transition matrix**. Write the population as a vector of genotype proportions, and one generation of reproduction becomes a matrix multiplication:

```equation
X(n) = M X(n-1)      so      X(n) = M^n X(0)
```

An exponentially branching process collapses into a matrix power, and everything linear algebra knows about such operators becomes available. Version I attributes this move to a bioinformatics textbook that records a Punnett grid as a matrix and raises it to a power to obtain the genotype ratios of a later generation. [he_petoukhov_2011] [version1_local]

The original research question was whether linear algebra could prove that non-sex-linked traits such as sickle-cell anaemia and blood group have the potential to decay across generations. That question is answered in Parts I and II. It also produced a limitation which Version I stated plainly: as more alleles and genes are added, maintaining a square, diagonalisable matrix forces the model to represent a smaller share of the population.

The present question follows from it:

> **Which representation answers a specified inheritance query completely, correctly and efficiently under a stated biological model?**

Each qualifying word carries weight. Completeness without a declared catalog is undefined. Efficiency without a declared output is misleading. Numerical agreement without external observation establishes nothing about biological accuracy.

# 2. Biological foundations

## 2.1 Alleles, genotype and phenotype

Humans are diploid at autosomal loci: one copy of each gene is inherited from each parent. Two copies of the same allele make a homozygous genotype; two different alleles make a heterozygous one. The pair of alleles carried is the **genotype**; the observable consequence is the **phenotype**. [mcgonegal_2014_biology]

A dominant allele, conventionally written with a capital letter, masks the effect of a recessive allele. For a disease carried on a recessive allele, only the homozygous recessive individual is affected, while the heterozygote is an unaffected carrier who can still transmit the allele. This asymmetry between what is expressed and what is transmitted is the reason genotype, not phenotype, must be the state variable of the model. [kognity_alleles_2022] Part V returns to it, because expression turns out to depend on more than the allele pair.

## 2.2 Sickle-cell anaemia

Sickle-cell anaemia results from inheriting two copies of the HbS allele. The normal HbA allele carries GAG at the sixth codon of the beta-globin coding sequence; HbS carries GTG. The transcribed messenger codons are therefore GAG and GUG, which translate to glutamic acid and valine respectively. Glutamic acid is hydrophilic; valine is hydrophobic, and the substitution distorts haemoglobin polymerisation and hence red-cell shape, reducing oxygen-carrying capacity and causing vaso-occlusion. [{{CLINVAR_SOURCE}}] [kognity_alleles_2022] Version I illustrates the morphology with a clinical figure and the geographic association with a distribution map; both are third-party images, and this document cites them rather than reproducing them. [froedtert_sickle_2022] [britannica_sickle_2022]

The three genotypes of the simplified two-allele model are:

| Genotype | Status | Note |
|---|---|---|
| HbA HbA | Unaffected | No sickle allele |
| HbA HbS | Carrier | Confers a degree of malaria resistance |
| HbS HbS | Affected | Sickle-cell anaemia |

The heterozygote advantage in malarial regions is the standard explanation for the allele's persistence, and it is the assumption that would have to be dropped to model those regions correctly. [{{SCD_SOURCE}}] Version I supported this with an encyclopaedia distribution map [britannica_sickle_2022]; the primary evidence is a geostatistical map of HbS allele frequency built from assembled survey data, which found the areas of high frequency coincident with historical malaria endemicity. [piel2010sickle]

*Scope note, revised in Part I of this unified treatment.* AA/AS/SS is an educational abstraction, not a catalog of HBB disease. Compound heterozygous forms involving other HBB variants also cause disease. Carrier status does not mean half of circulating red cells are permanently sickled. And a universal assumption that affected individuals do not reproduce is not supported; Version I's own conclusion estimated roughly 10% worldwide reproduction among affected individuals, which is small but not zero. Version II therefore treats reproductive weights as declared model parameters and does not forecast eradication dates.

## 2.3 The six parental crosses

With three genotypes per parent and parental order disregarded, there are exactly six parental pairings. Their offspring distributions follow from the Punnett squares:

| Parental cross | HbA HbA | HbA HbS | HbS HbS |
|---|---|---|---|
| HbA HbA x HbA HbA | 1 | 0 | 0 |
| HbA HbA x HbA HbS | 1/2 | 1/2 | 0 |
| HbA HbA x HbS HbS | 0 | 1 | 0 |
| HbA HbS x HbA HbS | 1/4 | 1/2 | 1/4 |
| HbA HbS x HbS HbS | 0 | 1/2 | 1/2 |
| HbS HbS x HbS HbS | 0 | 0 | 1 |

Each row is a valid conditional distribution over offspring and sums to one.

*Complexity.* With one biallelic locus the whole problem is trivially small: G = 3 genotypes, U = 6 unordered parental pairs, and 10 supported transitions out of 18 dense entries, a density of 0.5556. More than half the table is non-zero, so there is nothing for a sparse representation to compress. Section 10.1 records this as the first row of the cost table, and it is the reason Version I's approach was entirely adequate for this trait.

# 3. The matrix algebra Version I requires

This section states the linear algebra used later, at the level at which Version I developed it. [haese_2014_furthermaths] Nothing outside it is assumed.

**Identity matrix.** The n by n identity I has ones on the main diagonal and zeros elsewhere, and satisfies A A^-1 = I.

**Determinant.** For a 2x2 matrix with entries a, b, c, d, det = ad - bc. For a 3x3 matrix with first row a, b, c, expansion along that row gives det = a(ei - fh) - b(di - fg) + c(dh - eg), where each minor is the determinant of the entries excluding the row and column of its coefficient, and signs alternate.

**Inverse of a 2x2.** Swap the main-diagonal entries, negate the off-diagonal entries, and divide by the determinant.

**Eigenvalues.** For a square matrix M, the values of lambda satisfying det(M - I lambda) = 0. For a 2x2 this expands to lambda^2 - (a+d) lambda + (ad - bc) = 0, solvable by the quadratic formula. An n by n matrix has n eigenvalues counted with multiplicity, and their product is the determinant.

**Eigenvectors and diagonalisation.** An eigenvector X satisfies M X = lambda X, equivalently (M - I lambda) X = 0 for non-zero X, which requires det(M - I lambda) = 0. Assembling the eigenvectors as the columns of S and the eigenvalues on the diagonal of D gives

```equation
M = S D S^-1     and therefore     M^n = S D^n S^-1
```

Raising a diagonal matrix to a power raises each diagonal entry, which is why diagonalisation makes high powers cheap. This is the computational motivation for the entire Version I approach.

# 4. The sickle-cell transition matrix

## 4.1 Restricting to the reproducing crosses

Version I assumes a non-malarial population, so the heterozygote carries no survival advantage, and assumes that HbS HbS individuals survive but do not reproduce. Under that assumption the three crosses involving an affected parent are removed, leaving three parental columns and three offspring rows:

| Offspring | AA x AA | AA x AS | AS x AS |
|---|---|---|---|
| HbA HbA | 1 | 1/2 | 1/4 |
| HbA HbS | 0 | 1/2 | 1/2 |
| HbS HbS | 0 | 0 | 1/4 |

Let A(n), B(n) and C(n) be the fractions of the nth generation with genotypes AA, AS and SS. Reading the table as coefficients gives the linear system

```equation
A(n) = A(n-1) + (1/2) B(n-1) + (1/4) C(n-1)
B(n) =          (1/2) B(n-1) + (1/2) C(n-1)
C(n) =                         (1/4) C(n-1)
```

which is exactly X(n) = M X(n-1) with M the upper-triangular matrix of the table.

## 4.2 Diagonalisation

Because M is upper triangular, det(M - I lambda) = (1 - lambda)(1/2 - lambda)(1/4 - lambda), giving eigenvalues

```equation
lambda_1 = 1,   lambda_2 = 1/2,   lambda_3 = 1/4
```

Solving (M - I lambda) X = 0 for each, and normalising the free parameter to x = 1, gives eigenvectors (1, 0, 0), (1, -1, 0) and (1, -2, 1). Assembling them as the columns of S, and obtaining S^-1 by elementary row operations on the augmented matrix [S | I], yields M = S D S^-1 and hence

```equation
M^n = S D^n S^-1,     D^n = diag(1, (1/2)^n, (1/4)^n)
```

## 4.3 The limit and its correct reading

Since (1/2)^n and (1/4)^n both tend to zero, D^n tends to diag(1, 0, 0) as n grows. Version I read this as the carrier and affected classes ceasing to exist over generations, leaving only HbA HbA.

The algebra is correct. Its interpretation requires care, and this is the first place where the unified treatment revises Version I.

> **RESULT 4.1 - WHAT THE EIGENVALUE DECAY DOES AND DOES NOT SHOW**
>
> The eigenvalues 1, 1/2 and 1/4 are algebraically correct for the matrix as written. D^n is a representation in the eigenvector basis, not a table of genotype probabilities in the original basis; an nth-generation statement must be transformed back through S.
>
> More importantly, the columns of M are labelled by parental *pairs* and its rows by single offspring *genotypes*. Multiplying by M a second time requires a rule converting a distribution over offspring genotypes into a distribution over parental pairs. Having three entries on each side does not supply that rule. Section 12 constructs it.

## 4.4 Audit: the displayed program does not conserve probability

The MATLAB listing displayed in Version I performs one further operation. After recording each generation it sets the SS component to zero before the next multiplication, without renormalising. The listing is reproduced verbatim in Appendix F, where the operation is the single line `X_next(3) = 0;` immediately after the generation is recorded. Reconstructing that algorithm in exact rational arithmetic from the AS x AS cross reproduces the displayed screenshot exactly, and also exposes what it costs:

{{V1_LEGACY_SICKLE_TABLE}}

The displayed generation-three output, AA 75%, AS 18.75%, SS 0%, is reproduced. Its components sum to 93.75%. Those values can therefore be obtained, but they are not a normalised distribution of genotype proportions; reading them as such overstates the AA share and understates everything else. Treating the deficit as surviving mass would require an explicit population-size interpretation and a justified transition law, neither of which the genotype-pair relabelling supplies.

This is not a coding error so much as a modelling one, and it is the reason Part III separates the transmission kernel from the mating rule.

## 4.5 The repaired generational model

Selection and mating are separated. Selection reweights the current genotype frequencies by fitness and **renormalises the breeding population**; mating then draws two parents independently from that normalised distribution.

Starting from the AS x AS offspring distribution (1/4, 1/2, 1/4) under neutral random mating, the distribution is unchanged: Mendelian segregation alone does not remove an allele. With reproductive weights (1, 1, 0), the S allele frequency q updates as

```equation
q' = q / (1 + q),     hence     q_t = 1 / (t + 2)   starting from q = 1/2
```

> **RESULT 4.2 - THE DECLINE IS RECIPROCAL, NOT EXPONENTIAL**
>
> Under complete reproductive exclusion of affected individuals in an otherwise neutral randomly mating population, the deleterious allele frequency declines as 1/t, not geometrically. At the third offspring generation the SS birth probability is 1/16 under this model, whereas the displayed legacy algorithm has already set it to zero.
>
> The practical consequence: Version I's projection that the disease approximately disappears within four generations and carriers within eight - converted there into calendar years using an average generation interval of about 25 years [isogg_generation_2016] - does not follow from a correctly normalised model. Recessive alleles persist in carriers, which is precisely why recessive disease is hard to eliminate.

![Generational dynamics and the normalisation audit](../figures/generation_models.png)

*A note on figures and their licensing.* Version I's decay curves were plotted with an online graphing tool [desmos_2022]; every figure in this document is instead regenerated from retained machine-readable output, so that no plotted value is independent of the recorded results.

Version I's two biological illustrations are third-party images: a clinical figure of red-cell morphology [froedtert_sickle_2022] and an encyclopaedia distribution map [britannica_sickle_2022]. Neither carries a licence permitting redistribution, and this document is published publicly, so neither is reproduced. Both are cited at their point of use, which is the reuse that citation does permit, and the underlying claims are additionally supported by primary literature that can be cited without reproduction [piel2010sickle] [{{SCD_SOURCE}}]. Where a figure was carrying explanatory weight it has been replaced rather than dropped: the geographic argument is now made by the per-population allele frequencies of Sections 11 and 13, which are computed here from the frozen reference-panel snapshots and are therefore free of any reuse restriction.

Version I's screenshots of MATLAB source and program output are superseded by Appendix F, which carries the source verbatim, and by the exact reconstruction of the displayed outputs in Section 6.3.

<!-- pagebreak -->

# PART II - BLOOD GROUPS

# 5. ABO: a three-allele locus

## 5.1 Alleles, weights and genotypes

The ABO locus carries three alleles: I^A, I^B and i. A and B are codominant, and both are dominant to O. [dean_abo_2005] [wikipedia_bloodtype_2024] Version I introduced a weighting scheme to score genotypes, assigning full expressive weight to the antigen-producing alleles and none to the recessive one:

| Allele | Weight w |
|---|---|
| I^A | 1 |
| I^B | 1 |
| i | 0 |

with a genotype score defined as the mean of its two allele weights, (w1 + w2)/2. This gives 1.0 for AA, BB and AB, 0.5 for AO and BO, and 0.0 for OO.

That scheme is worth pausing on, because it is the earliest form of an idea this document returns to twice. It is a **linear score over allele dosage**: a weight per allele, averaged over the two copies. Part III generalises it to an additive polygenic score over many loci, and Part V generalises it again into a weighted score with a bias term and a nonlinear link, which is what allows expression to depend on age and epigenetic state as well as on genotype. Version I's weight column is the seed of all of it.

Six genotypes follow: AA, AO, BB, BO, AB, OO.

## 5.2 Counting parental pairs

> **RESULT 5.1 - PAIR COUNTS FOR THE ABO CATALOG**
>
> For six genotype labels, the number of unordered pairs with repetition allowed is 6 x 7 / 2 = 21. The number 15 counts only pairs of *distinct* genotypes and omits the six equal-genotype pairings.
>
> An equal-genotype pair means two parents who happen to share a genotype. It does not imply self-fertilisation, and it cannot be discarded.

Version I states 15 on page 30 and then correctly enumerates 21 pairings on pages 32 and 33. The unified treatment keeps 21.

## 5.3 From phenotype frequencies to genotype frequencies

Version I takes worldwide phenotype shares of approximately O 44%, A 42%, B 10%, AB 4% from a secondary summary of global blood-group prevalence [medicalchannelasia_abo_2023], and applies the Hardy-Weinberg relations as given in a standard population-genetics text to recover allele frequencies p_A = 0.26, p_B = 0.077, p_O = 0.663, which sum to one. [hartl_clark_2007] [{{HWE_SOURCE}}] Genotype probabilities follow as

```equation
P(AA) = p_A^2,  P(AO) = 2 p_A p_O,  P(BB) = p_B^2,
P(BO) = 2 p_B p_O,  P(AB) = 2 p_A p_B,  P(OO) = p_O^2
```

Under independent random mating, an unordered parental pair has probability x_i^2 when the genotypes are equal and 2 x_i x_j when they differ.

*Identifiability note.* The phenotype shares are drawn from a secondary aggregator rather than a primary survey, and carry no stated sampling frame or uncertainty; Section 11.4 shows why pooling them across populations is consequential rather than cosmetic. Converting phenotype frequencies into allele frequencies additionally assumes the genotype-phenotype map and population equilibrium. These frequencies are therefore useful as declared inputs for reproduction, but they are not independent validation data for a model derived from them. Predicting them back would reuse the answer.

## 5.4 The six leading crosses, and where 83.62% comes from

Ranking the 21 unordered pairs by probability gives:

{{V1_TOP6_TABLE}}

Version I sums the printed percentages to {{PRINTED_TOP6}}% and retains these six crosses to form a 6 x 6 square matrix against the six offspring genotypes. The offspring table is:

| Offspring | AO x OO | OO x OO | AO x AO | BO x OO | AO x BO | AA x OO |
|---|---|---|---|---|---|---|
| AA | 0 | 0 | 0.25 | 0 | 0 | 0 |
| AO | 0.50 | 0 | 0.50 | 0 | 0.25 | 1.00 |
| BB | 0 | 0 | 0 | 0 | 0 | 0 |
| BO | 0 | 0 | 0 | 0.50 | 0.25 | 0 |
| AB | 0 | 0 | 0 | 0 | 0.25 | 0 |
| OO | 0.50 | 1.00 | 0.25 | 0.50 | 0.25 | 0 |

## 5.5 Audit: the percentage depends on which inputs are propagated

Recomputing the retained mass under each candidate set of inputs separates arithmetic from provenance:

{{V1_COVERAGE_VARIANTS}}

The coarse genotype inputs used in the pair formulas sum to 1.001, so the total parental-pair mass is 1.001^2 = 1.002001 and the resulting {{PRINTED_TOP6}}% is not a rigorously normalised probability. Propagating the stated allele frequencies exactly gives {{EXACT_TOP6}}%. The discrepancy is small; the point is that a quoted percentage needs its normalisation stated alongside it.

These are model-based pair probabilities under chosen allele frequencies. They are not directly measured shares of the world's population.

# 6. ABO x Rh: two loci at once

## 6.1 The Rh locus

The RHD gene encodes the D antigen. In the simplified model used here, D is dominant and d recessive, giving genotypes DD, Dd and dd and phenotypes positive, positive and negative. Two Rh-positive parents who are both Dd have a one-quarter chance of a dd child. [{{RH_SOURCE}}] [britannica_rh_2025]

*Scope note.* The full Rh system involves RHD and RHCE variation and many variants not captured by a single dominant/recessive pair. The one-quarter statement is valid for Dd x Dd, not for every pair of positive parents.

Version I takes Rh-negative phenotype frequency 0.07, from a secondary chart of regional blood-type prevalence. [statista_bloodtype_2023]

> **RESULT 6.1 - THE ALLELE FREQUENCY IS THE SQUARE ROOT**
>
> Under the simplified Hardy-Weinberg model, a recessive phenotype frequency of 0.07 implies a d allele frequency of sqrt(0.07) = 0.2645751311, not 0.07. Using the phenotype frequency directly as the allele frequency would misstate every genotype probability that follows.

## 6.2 Combined catalog and pair counts

Six ABO states combine with three Rh states to give 18 combined genotypes, and therefore:

{{STATE_COUNT_TABLE}}

*Complexity.* ABO alone has G = 6, U = 21 and 45 supported transitions out of 126 dense entries, a density of 0.3571. Adding Rh gives G = 18, U = 171 and 615 supported transitions out of 3,078, a density of 0.1998. The density has almost halved for one added locus, which is Result 7.4 becoming visible at the smallest scale the paper examines: this is the point at which compressed storage begins to be worth considering, and the point at which manual tabulation stops being reasonable.

Multiplying the exact ABO and Rh genotype frequencies and ranking the unordered joint pairs gives:

{{V1_TOP18_TABLE}}

The top-18 mass is {{EXACT_TOP18}}%, which rounds to the reported {{PRINTED_TOP18}}%. Independently rounded row percentages sum to 67.84%; rounding rows and rounding their unrounded total need not agree.

## 6.3 Audit: what the displayed 18 x 18 table actually is

This is the most consequential audit finding in the document.

{{DISPLAYED_MATRIX_CHECK}}

In other words, the displayed 18 x 18 table is **each combined genotype crossed with itself**, not the ranked top-18-parental-pair table it is described as. Meanwhile the MATLAB code on the following pages loops over all 324 ordered parental pairs, combining independent ABO and Rh transmissions, and is therefore a complete random-mating calculation for its declared model.

The same holds one locus down: the ABO code loops over all 36 ordered pairs, so iterating the displayed six-column table instead would give a different result.

> **RESULT 6.2 - THE LEGITIMATE COMPARISON IS NARROWER THAN IT LOOKS**
>
> Version I's executed programs already include every pairing within their declared catalogs. A successor can therefore show that the prose's truncation is unnecessary, and can measure storage and inference alternatives to the complete algorithm.
>
> It cannot claim that the implemented ABO and ABO x Rh programs had {{PRINTED_TOP6}}% and {{PRINTED_TOP18}}% support and that a new sparse implementation raised them to 100%. That comparison would misrepresent the baseline.

Six displayed examples were reconstructed in exact rational arithmetic to confirm the reconstruction before the audit was trusted:

{{LEGACY_EXAMPLE_TABLE}}

<!-- pagebreak -->

# PART III - THE COMBINATORIAL THEOREM AND THE COMPLETE KERNEL

# 7. How large is the state space?

Version I's limitation was empirical: the matrices got unmanageable. This section makes the growth exact, because a limitation that is proved can be reasoned about, while one that is merely observed cannot.

## 7.1 Genotypes: a proof by induction

> **RESULT 7.1 - GENOTYPE COUNT FOR n BIALLELIC LOCI**
>
> For n independent biallelic diploid loci, the number of distinct unphased multilocus genotypes is G(n) = 3^n.
>
> *Proof by induction on n.*
>
> *Base case.* For n = 1 the catalog is AA, Aa, aa, so G(1) = 3 = 3^1. For completeness, n = 0 admits exactly one empty tuple, and 3^0 = 1.
>
> *Inductive hypothesis.* Suppose G(k) = 3^k for some k >= 1.
>
> *Inductive step.* Every k-locus genotype tuple can be extended at a new independent locus in exactly three ways, by appending AA, Aa or aa. Distinct extensions of distinct tuples are distinct, and every (k+1)-locus tuple arises from exactly one k-locus prefix with exactly one appended state. The extension map is therefore a bijection onto a set three times the size, so G(k+1) = 3 G(k) = 3 * 3^k = 3^(k+1).
>
> By induction G(n) = 3^n for every non-negative integer n. QED

## 7.2 Generalisation to arbitrary allele counts

> **RESULT 7.2 - GENOTYPE COUNT FOR ARBITRARY ALLELE COUNTS**
>
> At a locus with a alleles there are a homozygous genotypes and a(a-1)/2 heterozygous ones, totalling a(a+1)/2. Across n independent loci with allele counts a_1 ... a_n,
>
> ```equation
> G = product over l of [ a_l (a_l + 1) / 2 ]
> ```
>
> The same induction applies: adding a locus multiplies the catalog by that locus's genotype count.
>
> *Check against Part II.* ABO has a = 3, giving 3*4/2 = 6 genotypes. Rh has a = 2, giving 3. Combined, 6 x 3 = 18, which is exactly the catalog Version I used.

## 7.3 Parental pairs

> **RESULT 7.3 - UNORDERED PARENTAL PAIR COUNT**
>
> Choosing two genotypes from G with repetition allowed and order disregarded gives
>
> ```equation
> U = G + G(G-1)/2 = G(G+1)/2
> ```
>
> For ABO, U = 6*7/2 = 21, confirming Result 5.1. For ABO x Rh, U = 18*19/2 = 171. For n biallelic loci, U = 3^n (3^n + 1)/2, which is asymptotically 9^n/2: **each additional biallelic locus multiplies the parental-pair space by about nine.**

This is the growth Version I ran into. It is not an artefact of MATLAB, of matrices, or of implementation quality.

{{SCALING_TABLE}}

## 7.4 Supported transitions: sparse, and still exponential

Not every (parental pair, child) combination is possible. Counting the possible ones bounds what a sparse representation must store.

> **RESULT 7.4 - SUPPORTED-TRANSITION COUNT**
>
> For one biallelic locus, the numbers of possible children across the nine ordered parental pairs form the rows (1,2,1), (2,3,2), (1,2,1), totalling 15. Support products distribute over sums across independent loci, so the ordered kernel has 15^n nonzero entries.
>
> To remove parent-order duplication, pair each supported transition with its parent-swapped counterpart. A transition is fixed under the swap exactly when the two parental genotypes are identical; the single-locus equal-genotype support counts are 1, 3, 1, totalling 5, hence 5^n across loci. Therefore
>
> ```equation
> N_nonzero = (15^n + 5^n) / 2
> N_dense   = G U = (27^n + 9^n) / 2
> density   = (15^n + 5^n) / (27^n + 9^n)
> ```

> **RESULT 7.5 - ARBITRARY ALLELE COUNTS**
>
> Let H = a(a-1)/2 be the heterozygote count at a locus. Across all one-locus parental genotypes the sum of distinct transmissible allele counts is a + 2H = a^2, so ordered pairs contribute a^4, except that two identical heterozygous parents generate one duplicated unordered heterozygous child; subtracting those H duplicates gives T(a) = a^4 - H. Equal-genotype parent pairs contribute D(a) = a + 3H. The unordered support total across independent loci is one half of the product of T(a_l) plus the product of D(a_l).
>
> *Verified by exhaustive allele-copy enumeration:* 45 nonzeros for the ABO catalog and 615 for ABO x Rh.

## 7.5 Hypotheses that follow, and what would refute them

The counting results license three conjectures about representation. They are stated as hypotheses because only the first is settled here.

**H1 - Sparse storage helps, and its help is bounded.** Density falls like (15/27)^n, so compressed storage should win increasingly on payload as n grows, while the absolute nonzero count still grows like 15^n/2. Prediction: a measurable payload advantage at small n that never becomes an asymptotic rescue. *Tested in Section 10; supported.*

**H2 - Crossover is workload-dependent, not universal.** Dense contiguous arithmetic should beat sparse indexing at small n, where overhead dominates. Prediction: dense is competitive or better at n = 1 or 2. *Tested in Section 10; supported.*

**H3 - Only factorisation escapes the exponent.** Any representation that materialises the joint kernel inherits 15^n/2; only methods that keep per-locus factors separate and answer a smaller question avoid it. Prediction: a factored single-child query should scale polynomially in n while full enumeration cannot. *Tested in Section 10 to n = 1000; supported for that query class only.*

A fourth hypothesis is raised in Part V and is **not** supported: that the sparsity of H1 survives biological realism. It does not, once mutation is admitted.

# 8. The complete kernel

## 8.1 Transmission separated from mating

For a diploid genotype g = (u, v), define the gamete probability t_g(a) as one half the number of copies of allele a in g. For an unordered child genotype o = (a, b) the local transmission rule is t_g(a) t_h(a) when a = b, and t_g(a) t_h(b) + t_g(b) t_h(a) otherwise, the second expression summing the two parental routes to the same heterozygote. [{{LINKAGE_SOURCE}}] Under independent segregation the multilocus conditional probability is the product of local ones.

```equation
K[o, pair(i,j)] = P(child = o | parent1 = i, parent2 = j)
shape(K) = G x U ;    sum over o of K[o, c] = 1
```

> **RESULT 8.1 - THE KERNEL NEED NOT BE SQUARE**
>
> K is a rectangular linear map from a distribution over parental pairs to a distribution over children. It requires neither an inverse nor a diagonalisation to perform that map.
>
> The square-matrix constraint of Version I was therefore never a mathematical necessity. It was a consequence of wanting to raise the operator to a power, which requires a map from a space to itself.

## 8.2 Why random mating is nonlinear

If two parents are drawn independently from the same genotype distribution x, set w_ii = x_i^2 and w_ij = 2 x_i x_j for i < j. The next generation is K w(x), and normalisation follows since the sum of w is (sum x)^2 = 1. But the map from x to the next x is quadratic.

> **RESULT 8.2 - NO FIXED LINEAR MAP REPRESENTS RANDOM MATING**
>
> A population of only AA remains AA; a population of only aa remains aa. Mixing the two equally and allowing random mating gives AA 1/4, Aa 1/2, aa 1/4, whereas averaging the two unmixed outcomes gives no heterozygotes at all.
>
> Hence F((x+y)/2) is not (F(x)+F(y))/2, and no fixed matrix power can represent the closed random-mating process.

## 8.3 When a square matrix is legitimate

A square operator is correct for a lineage whose mates are drawn from a *fixed external* distribution q. Define M_q[o, i] = sum over j of K[o, pair(i,j)] q_j. This maps individual genotypes to individual genotypes, preserves normalisation, and may legitimately be raised to a power while q is held fixed.

This rehabilitates part of Version I. Its sickle-cell matrix is a valid fixed-mate operator under a stated mating assumption; what it is not is a closed random-mating population model. The issue was never that square matrices are forbidden, but that the meaning of the state must match the operator applied to it.

# 9. Implementation

## 9.1 One semantic kernel, four physical representations

The engine encodes canonical unordered allele pairs per locus into a mixed-radix integer, validating dimensions and index bounds. Parent-exchange symmetry halves redundant work without deleting equal-genotype pairs. Local Punnett outcomes aggregate repeated allele-copy paths and are cached in a bounded least-recently-used table.

- **Dense** allocates all G x U binary64 entries: 8GU bytes of payload.
- **CSR** stores only positive transmissions with child-row indices and row pointers, built through compressed columns and converted for repeated updates. [{{HWE_SOURCE}}]
- **Hash adjacency** maps a parental pair to a map of its possible children; its footprint is dominated by Python object overhead and is reported under a different metric.
- **Streamed** generates local factor products on demand and accumulates a child distribution without storing the kernel.

All four use the same catalog order and transmission rule, which is what makes the comparison meaningful.

## 9.2 Factored queries

For fully specified parents and a fully specified child, the log probability is the sum of per-locus log probabilities: O(n) arithmetic and O(n) local storage, with no enumeration. This is probability factorisation, not a general pedigree solver; tree-structured factor graphs admit exact messages while arbitrary loops do not, and general variable elimination can require factors exponential in induced width. [kschischang2001sumproduct] [dechter1999bucket]

If the parental population itself factorises across loci, so does the offspring population. But independent segregation does not imply that an arbitrary population factorises: a population concentrated on AA/BB and aa/bb carries a correlation its marginals cannot retain, and a regression test checks that substituting marginals changes the result.

## 9.3 Resource limits and the numerical contract

Before materialising anything, the implementation computes catalog size, pair count, nonzero count and representation estimates, and a byte budget may refuse the request, returning a structured resource error while leaving factored queries available. A caught allocation failure clears the bounded cache and propagates a structured error rather than an incomplete object.

The scientific contract is explicit. **No probability threshold removes rare transmission branches.** Structural zeros are omitted in sparse representations; positive values are never pruned. Conditional queries expose a log-space interface, and the ordinary probability interface raises rather than returning a silent zero when a mathematically positive result underflows. Population updates use binary64 and explicitly reject detected zero underflows; this is not a general log-space population engine, and Section 24 shows why that gap becomes pressing.

## 9.4 Verification before measurement

The primary oracle enumerates maternal and paternal allele-copy choices independently of the kernel builder, converts each path to a canonical child, and sums exact rational weights. Exhaustive tests cover one-locus biallelic and triallelic systems, four-allele loci, two- and three-locus biallelic systems and ABO x Rh, checking parent-exchange symmetry, normalisation, supported-transition counts and exact Mendelian probabilities. Dense and CSR are compared entry for entry; hash and streamed results are compared against dense for seeded non-uniform inputs. Allocation failure is injected deterministically rather than by exhausting the machine.

Legacy regression tests deliberately preserve the original sickle-cell defect, expecting the displayed result to sum to 17/16 while the repaired model sums to one. A test that silently fixed the legacy output would erase the evidence the audit rests on.

Adversarial review found and fixed four numerical defects during development: underflow in a lazy probability product, an absolute tolerance that distorted extreme exact-HWE tails, underflow when multiplying very small positive fitness weights, and a dense convolution where three shifted additions suffice. An independent exact-integer HWE enumeration checked every genotype triple with total sample size 1 to 30, 5,455 triples, with maximum absolute discrepancy below 7e-15.

{{TEST_STATUS}}

# 10. Measurement

{{BENCH_METHODS}}

{{BENCH_TABLE}}

{{BENCH_RESULTS}}

![Measured representation memory and inference costs](../figures/benchmark_overview.png)

Factored single-child queries, a different and smaller output, scale as predicted by H3:

{{QUERY_TABLE}}

At six loci the dense kernel alone would require 1,551,807,720 payload bytes, and the configured 256 MiB budget refuses that materialisation. Sizes outside the bounded experiment are labelled resource-refused or not measured; they are never presented as completed timings.

These results compare implementations written for this project, mixing Python iteration with compiled NumPy and SciPy operations. They do not isolate data-structure choice from every implementation effect, and they report no MATLAB-to-Python speedup, because MATLAB was never executed.

## 10.1 Cost model and the complexity of each model in this paper

A measured time is a fact about one machine. An asymptotic bound is a fact about the problem, and it is the part that transfers. The cost model counts arithmetic operations on probability values and array entries written, treating one binary64 multiply-add as unit cost and ignoring cache behaviour, which the measured timings in Section 10 capture instead.

Under that model the derived quantities of Section 7 are directly the cost of the corresponding operation. Kernel construction writes one value per supported transition, so it is Theta of the nonzero count. A dense population update touches every entry including the structural zeros, so it is Theta(G U). A CSR update touches only the nonzeros. For the models this paper actually builds:

{{COMPLEXITY_TRAIT_TABLE}}

The density column is the quantity that decides whether compressed storage is worth its indexing overhead, and it falls as (15/27)^n by Result 7.4. At one locus it is above one half and sparse storage is pointless; at five loci it is 0.0529 and sparse storage saves an order of magnitude. Section 15.3 gives the matching lower bounds, which is what turns these upper bounds into statements about the problem rather than about this implementation.

## 10.2 Derived growth against measured growth

An asymptotic claim earns its place only if the machine can be seen approaching it. Fitting T = a b^n to the measured construction medians and comparing b against the derived per-locus work ratio gives:

{{COMPLEXITY_GROWTH_TABLE}}

{{COMPLEXITY_NOTE}}

![Derived growth of the four cost regimes, and measured against derived construction growth](../figures/complexity_growth.png)

> **RESULT 10.1 - THE DERIVED EXPONENT IS CONFIRMED AT THE TOP OF THE MEASURED RANGE, NOT ACROSS IT**
>
> From four loci to five, the supported-transition count rises by a factor of 14.878. Measured construction time over the same step rose by 15.01 for dense, 15.09 for CSR and 14.68 for hash: agreement with the derived ratio to about 1.5 per cent.
>
> A single exponent fitted across all of n = 1 to 5 instead returns 8.5 to 10.4, because fixed per-call overhead dominates the smallest problems and flattens the fitted slope.
>
> Both numbers are correct and they answer different questions. The asymptotic claim is entitled to the top of the range and not to the whole of it, and reporting only the fitted exponent would have understated the growth by roughly a third.

## 10.3 Structures considered and not built

A trie can share prefixes and answer partial-genotype queries, but a complete n-locus ternary trie still has (3^(n+1) - 1)/2 nodes; dynamic allocation reduces unused pointers without reducing catalog size. [sedgewickwayne_tries] [grossi2011compressedtries] A hash map gives expected constant probes for fixed-size keys while hashing a long multilocus key still costs length-dependent work. [sedgewickwayne_hash] Reduced decision graphs merge equal suffix computations, but their size depends on the function and the variable order. [bryant1986bdd] A priority queue can visit likely outcomes first, but ordering an incomplete output does not make it a distribution. [python_heapq]

None of these was built. Building one to demonstrate a predetermined winner would have been the wrong experiment; each solves a different subproblem, and a later benchmark must supply the prefix, reuse or ranked-output workload that justifies the extra machinery.

<!-- pagebreak -->

# PART IV - EVIDENCE AND COMPARISON

# 11. Testing the population assumptions against real data

## 11.1 What Hardy-Weinberg is

At one locus with allele probabilities p and q = 1 - p, independent gametes give genotype probabilities p^2, 2pq, q^2, and under the idealised large-population neutral random-mating model these persist across generations. [{{HWE_SOURCE}}] Selection, migration, non-random mating, drift, sampling and genotyping behaviour can all make an observed sample differ. The algebra alone does not license reading a pooled global phenotype percentage as a universal genotype distribution, which is exactly what Part II's ABO inputs did.

## 11.2 Data and provenance

The empirical test uses the public Ensembl rs334 population-genotype endpoint, filtered to the 26 non-overlapping component populations of the 1000 Genomes phase-3 reference panel: 2,504 individuals, with genomic T/T = 2,367, T/A = 137, A/A = 0. These are observed genotype calls, not frequencies inferred under Hardy-Weinberg. Entries absent from the API are zero-completed only after genotype totals, allele counts and reported frequencies agree, and the snapshot is retained with a hash and fetching instructions. [ensembl_rs334_phase3_20260905] [{{1000G_SOURCE}}]

Superpopulation and ALL records overlap component populations and are used only as separate pooled diagnostics. The variant is on the GRCh38 forward strand; HBB's transcript notation has the complementary orientation. The dataset is not a worldwide probability sample, a newborn disease survey, or a family transmission cohort, and its population labels are source metadata, not a basis for inferring any individual's ancestry or clinical outcome. [igsr_disclaimer_2026]

## 11.3 The exact test

Given counts (nTT, nTA, nAA), estimate q = (2 nAA + nTA)/(2N) and compute N times the Hardy-Weinberg probabilities. These fitted expectations are a compatibility check, not out-of-sample predictions. The exact test conditions on N and the observed count of one allele; feasible heterozygote counts share the parity of the allele count and determine both homozygote counts, and the conditional probability of a table is proportional to N! 2^h divided by nTT! h! nAA!. The two-sided probability-ordering p-value sums feasible tables no more probable than the observed one, computed through log-factorials to avoid an absolute tolerance swamping a small tail. [{{EXACT_SOURCE}}] The 26 component tests receive a Holm adjustment. [holm1979]

{{POPULATION_TABLE}}

![External population counts and fitted Hardy-Weinberg expectations](../figures/population_hwe.png)

## 11.4 Reading

The YRI subset has (78, 30, 0) at N = 108, estimated rare-allele frequency 0.138889, expectations approximately (80.083, 25.833, 2.083) and exact p = 0.213792. The pooled AFR diagnostic has (529, 132, 0), expects (535.590, 118.820, 6.590) and gives p = 0.001590. None of the 26 component tests has a Holm-adjusted p below 0.05.

> **RESULT 11.1 - POOLING CHANGES THE PREDICTION BEFORE ANY BIOLOGY IS INVOKED**
>
> Summing fitted rare-homozygote expectations across component populations gives 7.600, against 1.874 from a single globally pooled allele estimate. The difference follows from the variation of q between groups, since the weighted mean of q^2 is at least the square of the weighted mean.
>
> Consequence for Part II: Version I's use of pooled worldwide phenotype shares as though they described one Hardy-Weinberg population is not a harmless simplification. Pooling changes what the model predicts.

No component population records a single A/A call. That is a fact about this ascertained reference panel, not evidence about disease prevalence, and it does not establish selection. External inheritance prediction would require independent parent-offspring observations with a defined sampling protocol.

# 12. Weights, coverage, and comparison with the square-matrix baseline

## 12.1 Four different things called accuracy

Version I's conclusion moves between structural coverage and predictive accuracy. Separating them is necessary before any comparison can be honest.

| Metric | Definition | What it does not establish |
|---|---|---|
| Genotype catalog coverage | represented labels / G | Population prevalence |
| Parent-pair coverage | represented pairs / U | Equal event probabilities |
| Supported-transition coverage | retained nonzeros / complete nonzeros | Clinical prediction |
| Retained probability mass | sum of retained event probabilities under declared inputs | A world-population estimate |
| Computational agreement | difference from an independent implementation or exact oracle | External biological accuracy |
| Predictive performance | held-out observed outcomes against predicted distributions | Portability to other populations |

The top-six ABO table covers 6/21 = 28.57% of pair *labels* while retaining about {{EXACT_TOP6}}% of probability *mass*. Both numbers are meaningful; they answer different questions, and neither is accuracy.

## 12.2 The comparison

| Aspect | Version I square matrix | Complete kernel (this work) |
|---|---|---|
| ABO object | 6 x 6, six retained crosses | 6 x 21 rectangular kernel |
| ABO x Rh object | 18 x 18, described as top-18 crosses | 18 x 171 rectangular kernel |
| What the 18 x 18 table actually contains | each genotype crossed with itself; 114 of 324 cells differ from the ranked construction | not applicable |
| Executed MATLAB behaviour | all 36 and all 324 ordered pairs already enumerated | same semantics, four representations |
| Parent-pair coverage | 28.57% (ABO), 10.53% (ABO x Rh) of labels in the displayed tables | 100% within the declared catalog |
| Retained probability mass | {{EXACT_TOP6}}% / {{EXACT_TOP18}}% for the displayed tables | 100% within the declared catalog |
| Requires diagonalisability | yes, and hence squareness | no |
| Rare outcomes | discarded to preserve squareness | retained; no threshold prunes positive branches |
| Predictive accuracy | not established | not established |

The last row is the important one. The new representation is complete where the displayed tables were truncated, and it is measured rather than asserted. It is not more *accurate*, because neither version has been tested against held-out inheritance outcomes.

> **RESULT 12.1 - WHAT THE REPLACEMENT ACTUALLY BUYS**
>
> Removing the square-matrix constraint converts a representational compromise into a declared modelling choice. Version I discarded rare parental pairs because the algebra it wanted required a square operator; the complete kernel keeps them because a rectangular map needs no inverse.
>
> This is a genuine gain in completeness and an honest one in scope. It is not a gain in predictive power, and the audit of Section 6.3 means it is not even a gain over what the original code executed.

## 12.3 Staged biological extensions

Before attempting anything polygenic, the engine was exercised on systems whose answers are checkable.

**Simplified M/N.** A two-allele codominant abstraction; MN x MN gives MM 1/4, MN 1/2, NN 1/4. The transmission engine needs no change, which is the point: dominance belongs to the genotype-phenotype map, not to the inheritance probabilities. [{{MNS_SOURCE}}]

**ABO-FUT1 epistasis.** ABO antigen expression depends on the H antigen precursor, so an hh state masks A and B expression. For AO/Hh x BO/Hh the independent factors give A, B and AB each 3/16, with the combined O-or-O-like category at 7/16. Internally the engine distinguishes true OO from hh masking, so genotype information survives when phenotype labels merge. Inspired by Bombay-phenotype biology; not a transfusion-compatibility calculator. [{{BOMBAY_SOURCE}}]

**Linkage and phase.** Two individuals both described as Aa/Bb may carry haplotypes AB/ab or Ab/aB. With recombination fraction r the first transmits AB and ab at (1-r)/2 each and recombinants at r/2; the second reverses this. At r = 0 the gamete supports are disjoint; at r = 1/2 they coincide. [{{LINKAGE_SOURCE}}]

> **RESULT 12.2 - THE UNPHASED CATALOG IS NOT ALWAYS SUFFICIENT**
>
> Two individuals with identical unphased multilocus genotypes can have different gamete distributions. The 3^n catalog of Result 7.1 therefore does not determine transmission for linked loci, and phase must be carried explicitly.
>
> Linkage disequilibrium concerns population associations; meiotic recombination concerns transmission. Both may matter, and neither should be silently replaced by an independence assumption.

![Staged epistasis and linkage examples](../figures/epistasis_linkage.png)

<!-- pagebreak -->

# PART V - POLYGENIC TRAITS: EYE COLOUR, HEIGHT, AND P(n)

# 13. Eye colour: the first trait Version I could not reach

Version I's conclusion named eye colour and height as the traits its method could not handle, because both would demand matrices too large to build. Eye colour is the right place to test that claim, for three reasons. It is genuinely polygenic, so it is not a disguised single-locus problem. One locus nevertheless dominates it, so a small model is not absurd. And real genotype data exist at that locus for the same reference panel already used in Part IV, so the model can be driven by observations rather than by assumed frequencies.

## 13.1 The dominant locus

rs12913832 lies in a conserved region of intron 86 of HERC2, in an element that regulates expression of the neighbouring OCA2 gene. It is the strongest known common determinant of blue versus brown iris colour: in the study that identified it, a single ordinal logistic regression on this one SNP reached R-squared 0.68 with an association LOD of 444. [sturm2008herc2] The mechanism is regulatory rather than coding: the brown-associated allele permits chromatin looping that raises OCA2 expression and hence melanin production, and the blue-associated allele does not.

An R-squared of 0.68 is high for a single variant and nowhere near 1. That gap is the reason a second locus is introduced below, and the reason Section 13.7 can compare this model against something external.

## 13.2 Resolving the allele orientation before trusting any label

Part IV recorded that rs334 is represented on the GRCh38 forward strand while HBB's transcript notation runs the other way. The same hazard is sharper here, because much of the eye-colour literature reports rs12913832 as a C/T variant with C blue-associated, while the reference panel reports it on the forward strand as A/G.

Taking a database allele label at face value would therefore be a coin flip. The orientation is instead resolved from the data itself, by a check that would fail if a future data revision reversed it:

> **RESULT 13.1 - THE BLUE-ASSOCIATED ALLELE IS IDENTIFIED BY POPULATION ORDERING, NOT BY LABEL**
>
> Blue eyes are common in northern Europe and rare in African and East Asian populations. In the frozen snapshot the G allele frequency is {{EYE_ORIENTATION}}, while it is 0.0000 in YRI, ESN, LWK, CHB and JPT.
>
> Therefore G is the blue-associated allele on this strand and A is brown-associated, which is the reverse of the C/T convention once strand is accounted for. The fetch script asserts this ordering and fails closed if it ever breaks.

This is a small check and it protects every number that follows. A model built on the reversed assignment would have predicted blue eyes for the entire African panel and would still have passed every internal consistency test, because nothing internal knows what the labels mean.

## 13.3 The declared two-locus model

Locus one is HERC2 rs12913832, with real genotype calls. Locus two is a **declared modifier locus** standing for the residual common pigmentation variation - OCA2, SLC24A4, SLC45A2, TYR, IRF4 and others - compressed into a single biallelic factor. That compression is the model's principal simplification. It is not defended as biology; it is declared so that the boundary between what is observed and what is assumed stays visible.

The phenotype map is epistatic: the modifier does different work on the two HERC2 backgrounds, which is precisely why two loci are needed rather than one.

| HERC2 genotype | Modifier dd | Modifier Dd | Modifier DD |
|---|---|---|---|
| GG (light background) | blue | blue | intermediate |
| AG | hazel | brown | brown |
| AA | hazel | brown | brown |

Four phenotype classes arise from nine genotypes. The map is a modelling choice; Section 13.7 tests the one consequence of it that external data can reach.

## 13.4 The complete kernel for eye colour

The engine needs no modification. Two biallelic loci give {{EYE_MODEL_COUNTS}}

All four representations were built and run on the same population input:

{{EYE_REPRESENTATION_TABLE}}

The maximum absolute deviation from the dense result across every representation is {{EYE_AGREEMENT}}, which is floating-point agreement rather than approximate agreement.

> **RESULT 13.2 - AT THIS SIZE SPARSITY IS NOT WORTH ITS OVERHEAD**
>
> The eye-colour kernel has density 0.3086: roughly one entry in three is a supported transition. Result 7.4 gives density (15/27)^n, which at n = 2 is 0.3086 exactly.
>
> Compressed storage is therefore near-pointless here, and the measured build and update times confirm it. The architecture of Part III earns nothing on this trait. It earns its keep only where the catalog is large, which for eye colour it is not.

That is a negative result about the representation and a positive one about the method: the same kernel, unchanged, handles a two-locus epistatic pigmentation trait and a one-locus blood-group trait, and the correct choice of storage differs between them.

![Support of the eye-colour kernel and its generational dynamics](../figures/eye_color_kernel.png)

## 13.5 Real population data and a Hardy-Weinberg audit

The frozen snapshot holds rs12913832 genotype calls for 2,504 individuals across the 26 component populations of the 1000 Genomes phase 3 panel, retrieved from the Ensembl REST API, validated against independently supplied allele counts, and stored with a SHA-256 digest. [ensembl_rs12913832_phase3] [{{1000G_SOURCE}}]

Exact Hardy-Weinberg tests, computed as in Section 11.3 and Holm-adjusted across the 26 tests, give:

{{EYE_POPULATION_TABLE}}

At the 0.05 level after adjustment, {{EYE_HWE_RESULT}} departs from Hardy-Weinberg proportions. As in Part IV, that is a statement about the compatibility of a sample with a null, not evidence that the null describes the biology.

![rs12913832 allele frequency across 26 populations and the model's phenotype prediction](../figures/eye_color_populations.png)

The left panel shows the sharpest population structure anywhere in this paper. The blue-associated allele is absent or nearly absent across the African and East Asian panels and reaches 0.909 in Finland. This is a far steeper gradient than rs334 shows, and it is the reason a single pooled worldwide eye-colour frequency would be meaningless in a way that Section 11.4 already established in general terms.

## 13.6 Generational dynamics

Starting from the CEU genotype frequencies joined to the modifier locus under linkage equilibrium, and applying the complete random-mating update repeatedly:

{{EYE_GENERATION_TABLE}}

> **RESULT 13.3 - THE DISTRIBUTION IS ALREADY AT ITS FIXED POINT**
>
> The phenotype probabilities move once, between the observed CEU frequencies and generation one, and then never again. That single step is the Hardy-Weinberg theorem itself: the observed sample is not exactly at Hardy-Weinberg proportions, one round of random mating puts it there, and every subsequent round leaves it alone.
>
> This is the fixed point of Section 11.1 appearing in a two-locus epistatic model. Neutral random mating changes neither allele frequency nor, under linkage equilibrium, the joint genotype distribution.
>
> Nothing about eye colour decays across generations. Version I's framing of inherited traits as decaying was a consequence of its selection assumption for sickle cell, not a property of inheritance, and a trait under no selection makes that explicit.

This is worth stating plainly because it is the cleanest available refutation of the intuition that traits dilute over generations. They do not. Only selection, drift, mutation or migration move allele frequencies, and this model contains none of them.

## 13.7 Comparison against real-world observation

The reference panel is not phenotyped. No eye colour is observed anywhere in the dataset, so no prediction made here can be validated against it, and none is. External comparison must come from published cohorts that measured both genotype and phenotype.

The CanPath cohort genotyped and collected self-reported eye colour for 5,481 individuals of European ancestry, 2,757 of them GG and 2,724 AA or AG. It reports that 33 per cent of GG individuals self-reported an eye colour other than blue, and that 512 of the 2,724 AA or AG individuals reported a colour other than brown or hazel. [abbatangelo2026canpath]

That gives one number this model can be held against, and the comparison is made in two clearly separated ways:

{{EYE_EXTERNAL_TABLE}}

> **RESULT 13.4 - A DECLARED PARAMETER UNDERPREDICTS THE OBSERVED DISCORDANCE**
>
> Under the phenotype map of Section 13.3, an individual with HERC2 GG fails to be blue exactly when the modifier locus is DD. With the modifier allele frequency declared at 0.5 and Hardy-Weinberg proportions assumed, that probability is 0.25.
>
> The observed figure is 0.33. The declared model is therefore wrong in a measurable direction: it does not leave enough room for GG individuals to be something other than blue.
>
> Setting the modifier frequency to 0.574456 reproduces 0.33 exactly. That is a **calibration against the observation**, not a validation by it. One parameter fitted to one number cannot be evidence for the model that contains it, and no accuracy follows from the fit.

The second comparison is less comfortable and more informative. A published six-SNP model, IrisPlex, reports areas under the receiver operating characteristic curve of 0.96 for blue and 0.96 for brown. [walsh2011irisplex] This work reports no AUC at all, because it has no phenotypes to compute one against.

> **RESULT 13.5 - THE GAP TO A REAL PREDICTOR IS NOT A GAP IN REPRESENTATION**
>
> The complete kernel of Part III represents every biologically possible outcome of the declared eye-colour model exactly, with structural coverage of 100 per cent and floating-point agreement across four implementations.
>
> It nonetheless predicts nothing, because prediction requires phenotype observations, fitted effect sizes and held-out evaluation, none of which a transmission kernel supplies. The distance between this model and IrisPlex is a distance in data and calibration, not in data structures.

That is the same distinction Section 12.1 draws between coverage and accuracy, arriving here with a concrete external benchmark attached. It is the most direct evidence in this paper that improving the representation does not improve the science on its own.

## 13.8 What the eye-colour model establishes

It establishes that the architecture transfers: a two-locus epistatic trait with real genotype data runs through the unchanged kernel and agrees across four representations to {{EYE_AGREEMENT}}. It establishes that at this catalog size sparse storage is not worth its overhead, which is a useful negative result. It establishes that the trait is at its Hardy-Weinberg fixed point under neutral random mating, and that the allele shows extreme population structure.

It does not establish any eye-colour prediction, any prevalence, or any accuracy. The one external number it can be compared against, it initially missed in a stated direction and then reproduced only by fitting a parameter to it.

# 14. Height, and traits with too many loci to enumerate

## 14.1 Why the kernel cannot be used at all

For eye colour the catalog had nine entries. Height is associated with thousands of common variants: a saturated map of common height-associated variation analysed approximately 5.4 million people and reported 12,111 associated independent SNPs, with materially different out-of-sample variance explained across ancestry groups. [yengo2022height]

Result 7.1 gives the catalog size directly. At n = 12,111 biallelic loci, G = 3^12111, a number with more than five thousand digits. There is no representation of that catalog, sparse or otherwise, and no amount of engineering changes it. The complete-kernel architecture of Part III is not slow for height; it is inapplicable.

## 14.2 Scores instead of catalogs

The escape is to stop enumerating genotypes and start propagating a distribution over a **summary statistic**. Take the additive model S = sum over l of w_l D_l, with dosage D_l in {0, 1, 2} and non-negative integer weights. This is Version I's ABO allele-weight scheme of Section 5.1, extended from one locus to many: a weight per allele, summed over loci instead of averaged over two copies.

Independent dosage factors let the distribution be built by dynamic programming, one locus at a time:

```equation
f_(l+1)(s) = sum over d in {0,1,2} of P(D_(l+1) = d) f_l(s - w_(l+1) d)
f_0(0) = 1
```

Each locus contributes three shifted additions into an array indexed by score, not a new dimension of table. For n loci of unit weight the support has at most 2n + 1 bins.

> **RESULT 14.1 - THE SCORE DISTRIBUTION IS POLYNOMIAL WHERE THE CATALOG IS EXPONENTIAL**
>
> Let B be the number of reachable score bins. Adding one locus costs three shifted multiply-add passes over the current support, so the total work is Theta(sum over l of B_l), bounded by Theta(n B). For unit weights B = 2n + 1 and the whole computation is **Theta(n^2) time and Theta(n) space**.
>
> The corresponding genotype catalog has 3^n entries. At n = 200 that is the difference between an array of 401 numbers and a catalog of about 2.66 x 10^95.
>
> The saving is not free: the output is a distribution over the score, not over genotypes. Any question that genuinely requires the joint genotype distribution is not answered by this method, and Section 15 makes that trade explicit.

## 14.3 Calibration under a known simulator

Because no height phenotypes are available, the score machinery is validated where the truth is stipulated rather than assumed.

{{POLYGENIC_SUMMARY}}

A Wilson 95 per cent interval for the simulated coverage is [89.427%, 90.602%]. [nist_wilson] The independent simulation used binomial sampling plus Gaussian draws rather than resampling the computed array, and no seed was searched for a passing result. [scipy_kstest]

![Synthetic polygenic score and calibrated outcome interval](../figures/polygenic_distribution.png)

This tests the algorithm, the discretisation and the interval construction. It tests nothing about human height.

## 14.4 Why this is not a height predictor

The score is not measured in centimetres, has no fitted human effect sizes and contains no environmental covariates; its residual distribution is known because it was stipulated. Substituting real effect sizes would break three of the assumptions the calibration depends on at once: effects are not integers, loci in linkage disequilibrium are not independent, and the residual is not Gaussian with known variance.

> **RESULT 14.2 - WHAT A REAL HEIGHT MODEL WOULD ADDITIONALLY REQUIRE**
>
> Beyond the machinery demonstrated here: a licensed effect-size source and a declared target population; harmonisation of variant alleles and genome builds; a treatment of linkage disequilibrium and, where relevant, phase; discretisation error bounds for non-integer effects; an explicit environmental and residual model; held-out phenotype evaluation; and reported portability across ancestry groups, which Yengo and colleagues show varies materially. [yengo2022height]
>
> None of these is supplied. Deterministic descendant-height claims would exceed the evidence by a wide margin, and interval claims would exceed it by a narrower one.

# 15. P(n): a trait dependent on n genes

The three traits of this paper sit at three points on one axis. Sickle cell has one locus; eye colour has two in this model; height has thousands. This section states the general case, which is what the question "how does this scale" actually means.

## 15.1 The general object

Let a trait depend on n independent loci with allele counts a_1 to a_n, and let P(n) denote the problem of answering a specified query about the offspring of a specified parental pair, or about a population, under the declared model. The relevant quantities are already derived:

```equation
G = product over l of [ a_l (a_l + 1) / 2 ]        (Result 7.2)
U = G (G + 1) / 2                                  (Result 7.3)
nonzeros = one half of [ product T(a_l) + product D(a_l) ]   (Result 7.5)
```

and for the biallelic case these collapse to G = 3^n, U ~ 9^n/2 and nonzeros = (15^n + 5^n)/2.

{{PN_TABLE}}

## 15.2 Four regimes, and which query each answers

P(n) is not one problem. Its difficulty depends entirely on what output is demanded, and the four regimes below are separated by that and by nothing else.

| Regime | Query answered | Time | Space | Feasible to |
|---|---|---|---|---|
| Complete kernel | The full parental-pair to child map, reusable | Theta(15^n) | Theta(15^n) | n = 5 here; n = 6 refused by budget |
| One full cross | All children of one specified pair | Theta(3^n) worst case | Theta(3^n) output | limited by output size |
| Factored query | Probability of one fully specified child | Theta(n) | Theta(n) | n = 1000 measured |
| Score distribution | Distribution of an additive summary | Theta(n B) | Theta(B) | n = 200 measured, B = 401 |

The first two are exponential and the last two are not. The distinction is not the data structure; it is whether the question asked has an exponentially large answer.

## 15.3 Lower bounds: why the exponent is not an implementation defect

Upper bounds describe this implementation. Lower bounds describe every possible implementation, and they are what justify the claim that the exponential cannot be engineered away.

> **RESULT 15.1 - OUTPUT-SIZE LOWER BOUNDS FOR P(n)**
>
> *(a) Complete kernel.* Any algorithm that materialises the complete transmission kernel must write one value per supported transition. By Result 7.4 there are (15^n + 5^n)/2 of them, so any such algorithm takes **Omega(15^n)** time and space. No representation avoids this while remaining complete, because the bound counts distinct mathematical values, not the encoding used to store them.
>
> *(b) One full cross.* Take both parents heterozygous at every locus. Each locus admits three child genotypes with positive probability, so the child distribution has 3^n entries, all positive. Writing them takes **Omega(3^n)**. This holds regardless of representation, and it is why "return all possible children" is not a reasonable request at large n.
>
> *(c) One specified child.* The probability of a fully specified child depends on the parental genotypes at every one of the n loci; changing any single locus changes the answer. Any correct algorithm must therefore read Omega(n) input. The implemented factored query runs in O(n), so it is **Theta(n) and asymptotically optimal**.
>
> *(d) Score distribution.* The output has B bins and each must be written, giving Omega(B); the dynamic programme achieves O(n B). For unit weights B = 2n + 1, so the method is within a factor of n of the output-size bound.

Parts (a) and (b) are the honest reading of Version I's difficulty. It did not fail because MATLAB was slow or because matrices were the wrong tool. It ran into a lower bound.

Parts (c) and (d) are the honest reading of what Version II adds. Two of the four regimes are provably optimal or near-optimal, and both work by refusing to produce an exponentially large object.

## 15.4 Complexity of the operations actually implemented

| Operation | Upper bound | Lower bound | Notes |
|---|---|---|---|
| Catalog size G | Theta(3^n) | Theta(3^n) | Result 7.1, exact |
| Unordered pairs U | Theta(9^n / 2) | Theta(9^n / 2) | Result 7.3, exact |
| Kernel construction | O(15^n) | Omega(15^n) | Result 15.1(a); tight |
| Dense materialisation | Theta(27^n / 2) | Omega(27^n / 2) | stores structural zeros too |
| One population update, dense | Theta(G U) = Theta(27^n / 2) | Omega(15^n) | dense pays for zeros |
| One population update, CSR | Theta(nnz) = Theta(15^n) | Omega(15^n) | tight |
| Factored single-child query | O(n) | Omega(n) | Result 15.1(c); optimal |
| Additive score distribution | O(n B) | Omega(B) | Result 14.1 |
| Local cross at one locus | O(a^2) | Omega(a) | a alleles at that locus |

Two entries deserve emphasis. The CSR population update is **asymptotically tight**: it does work proportional to the number of mathematically supported transitions, and Result 15.1(a) says nothing complete can do less. The dense update is a factor of (27/15)^n = 1.8^n worse in operations, which is exactly the density ratio of Result 7.4 and the reason the measured payload gap widens with n.

## 15.5 Protocol for adding a further trait

The scaffold accepts a new trait when, and only when, the following are stated in advance. Appendix E gives the same list in the form of a checklist.

1. The allele catalog per locus and the resulting G by Result 7.2, with the feasibility check of Section 9.3 run at that G, and the regime of Section 15.2 identified before any code is written.
2. The genotype-to-phenotype map, including epistasis, kept separate from the transmission rule as in Section 13.3.
3. Whether loci are linked, and if so the recombination fractions and the phase convention (Result 12.2).
4. The source of any effect sizes, with licence and population of origin recorded in the ledger.
5. A held-out evaluation target chosen before fitting, and the population it generalises to, distinguished from coverage per Section 12.1.
6. A conditioning check on the design, performed before any coefficient is reported (Result 28.1).
7. An independent oracle for at least one small case, so the implementation can be checked without trusting itself.

Eye colour has now been added under this protocol and satisfies items 1 to 3 and 7; it does not satisfy 4 to 6, and Section 13.8 says so. Height satisfies none of 4 to 6 and is therefore present only as the derivational treatment of Section 14. A trait that cannot supply items 4 to 6 may still be added as a derivational section, clearly marked, but may not report fitted coefficients or accuracy.

<!-- pagebreak -->

# PART VI - TIME, MUTATION, AND WHAT COMES NEXT

# 16. Where age enters

Everything in Parts I to V is timeless. The kernel depends on two parental genotypes and nothing else; the population operators advance a generation index, not an age. Two parents produce the same offspring distribution at twenty as at fifty, and an individual's genotype is fixed from conception.

Version I identified this gap explicitly, noting that random mutation can occur in an individual before they reproduce and that the model cannot represent it. [version1_local] The gap survived into the complete kernel unchanged: better representation did nothing about it.

> **OBSERVATION 16.1 - TWO ROUTES, ONLY ONE OF WHICH TOUCHES THE KERNEL**
>
> **Germline.** Parental age at conception changes the mutation content of the transmitted gametes. This alters what the next generation inherits, modifies the kernel itself, and is heritable.
>
> **Somatic and epigenetic.** An individual's own age changes which cells carry which mutations and changes the methylation and histone-modification state of their chromatin. This alters whether an inherited genotype is *expressed*, modifies the genotype-to-phenotype map, and is with narrow exceptions not transmitted.

The distinction is load-bearing. A model that adds one "age" term acting on both routes produces an operator that is neither a valid kernel nor a valid penetrance function, and fitting will not reveal the error, because both routes push outcomes in the same direction with age.

**Status of Sections 17 to 29.** Derivation only. No code, no test, no dataset, no figure in this project corresponds to any equation below.

# 17. Notation for Part VI

| Symbol | Meaning |
|---|---|
| t | Age of an individual, in years |
| a_f, a_m | Paternal and maternal age at conception, in years |
| Lambda | Expected genome-wide de novo mutation count in one transmission |
| mu_l | Expected de novo count inside locus l in one transmission |
| eps_l | Probability locus l carries at least one de novo mutation in one transmission |
| Q | Row-stochastic substitution matrix over a locus's alleles, given a mutation |
| M(a) | Age-dependent row-stochastic gamete operator |
| lambda_l(t) | Somatic mutation intensity at locus l at age t |
| m_l(t) | Methylation fraction at locus l at age t, in [0,1] |
| h_l(t) | Standardised activating acetylation signal at locus l at age t |
| alpha, beta | Methylation and demethylation rates |
| z_l(t) | Linear expression score (log-odds scale) |
| w, b_l | Weight vector and locus bias |
| pi_l(t) | Probability locus l is transcriptionally competent at age t |
| J(w) | Mean negative log-likelihood objective |
| eta | Gradient step size |

# 18. Germline mutation and the transmission kernel

## 18.1 What is measured

Sequencing 78 Icelandic parent-offspring trios gave an average de novo mutation rate of 1.20e-8 per nucleotide per generation at mean paternal age 29.7, with the count rising by roughly two mutations per year of paternal age. [kong2012paternal] The larger follow-up, 1,548 trios and 108,778 high-quality de novo mutations, averaged 70.3 per trio and separated the parental contributions: 1.51 additional mutations per year of paternal age against 0.37 per year of maternal age. [jonsson2017trios]

These are regression slopes from one country's cohort, not per-locus probabilities for an arbitrary gene in an arbitrary population.

## 18.2 From a genome-wide count to a per-locus probability

```equation
Lambda(a_f, a_m) = Lambda_0 + beta_f (a_f - a_f0) + beta_m (a_m - a_m0)
beta_f = 1.51 per year,   beta_m = 0.37 per year     [anchored, not fitted here]
```

With G_bp the callable genome length, L_l the locus length and kappa_l a locus enrichment factor absorbing non-uniform mutability,

```equation
mu_l(a_f, a_m) = Lambda(a_f, a_m) * (L_l / G_bp) * kappa_l
eps_l(a_f, a_m) = 1 - exp( -mu_l(a_f, a_m) )
```

Setting kappa_l = 1 asserts uniform mutability, which is known to be false and is retained only as a declared null. For a gene-sized locus eps_l is of order 1e-8 to 1e-6 per transmission.

## 18.3 The age-parameterised gamete operator

The implementation already exposes a validated row-stochastic gamete mutation operator. The extension is to make that matrix a function of parental age, not to add a mechanism.

> **RESULT 18.1 - AGE-PARAMETERISED GAMETE MUTATION**
>
> For a locus with allele alphabet size k,
>
> ```equation
> M_ab(a_f, a_m) = (1 - eps_l) delta_ab + eps_l Q_ab
> ```
>
> and the transmitted gamete law becomes t~_g(a) = sum over b of t_g(b) M_ba, with t_g the Mendelian law of Section 8.1.

> **RESULT 18.2 - NORMALISATION IS PRESERVED**
>
> Each row of M sums to (1 - eps_l) + eps_l * 1 = 1. Since t_g is a distribution and M is row-stochastic, the sum over a of t~_g(a) equals the sum over b of t_g(b) times the row sums of M, which is 1.
>
> Therefore the modified kernel remains a conditional distribution over children, and every normalisation test in the existing suite remains a valid check of it. The allele-copy oracle would then compare against exact Mendelian probabilities composed with a known stochastic matrix, which is equally checkable by enumeration.

# 19. The cost: sparsity does not survive

> **RESULT 19.1 - A POSITIVE MUTATION RATE DESTROYS STRUCTURAL SPARSITY**
>
> Result 7.4 gives (15^n + 5^n)/2 nonzeros against G*U total, and that falling density is the entire justification for compressed sparse storage.
>
> If eps_l > 0 and Q has full support, every allele is reachable from every allele in one transmission. Every structural zero becomes a positive number of order eps_l per mutated locus, and the nonzero count rises to G*U exactly.

{{SPARSITY_MUTATION_TABLE}}

At five loci that is a rise from 381,250 to 7,203,978 stored entries, and CSR becomes strictly worse than dense: it pays index overhead on a matrix with no zeros. The measured 12.59-fold payload advantage of Section 10 is an advantage over a mutation-free model. Introduce biologically realistic mutation and the advantage does not shrink; it inverts.

This falsifies hypothesis H4 of Section 7.5 before it was ever tested, and it is the single most consequential finding in Part V. Three responses exist, and choosing between them is an experiment rather than a preference:

1. **Threshold.** Restores sparsity, but violates the contract of Section 9.3 forbidding the pruning of positive branches. Rejected here.
2. **Factor.** Keeps per-locus operators small and never multiplies them out. This is Section 9.2's argument, and it becomes considerably stronger under mutation.
3. **Structure.** Keep the Mendelian part sparse and carry the mutation part as a low-rank or Kronecker correction, exploiting the fact that M is the identity plus a rank-structured perturbation of size eps_l.

The third is the interesting one, and it is not implemented.

# 20. Somatic accumulation within a lifetime

Somatic mutation does not change what an individual transmits; it changes the individual. Adult stem cells of liver, colon and small intestine accumulate roughly 40 novel mutations per year, at broadly similar rates across those tissues despite very different cancer incidence. [blokzijl2016somatic]

Model somatic mutation at locus l as an inhomogeneous Poisson process with intensity lambda_l(t):

```equation
m_l(t) = integral from 0 to t of lambda_l(s) ds
P(locus l unmutated at age t) = exp( -m_l(t) )
```

> **RESULT 20.1 - CONSTANT INTENSITY GIVES EXPONENTIAL SURVIVAL, NOT LINEAR**
>
> With lambda_l constant, m_l(t) = lambda_l t and the probability of at least one somatic mutation by age t is 1 - exp(-lambda_l t).
>
> For lambda_l t much less than 1 this is approximately lambda_l t, which is why observed counts look linear in age. The linear appearance is the small-argument regime of an exponential, and the two diverge at large lambda_l t.

A single hit is often insufficient. The classical multistage argument supposes k independent rare steps and derives an incidence rising as a power of age; fitting several non-endocrine carcinomas gave approximately a sixth-power dependence, with the log-log slope read as the number of stages minus one. [armitage1954multistage]

> **RESULT 20.2 - MULTISTAGE HAZARD**
>
> If each of k stages occurs independently at small constant rate, the probability all k are complete by age t is of order (c t)^k / k!, so the hazard is
>
> ```equation
> h(t) = c^k t^(k-1) / (k-1)!
> ```
>
> and a log-log plot of hazard against age has slope k - 1, which makes k estimable from incidence curves alone.

This is a hazard for a multi-hit somatic process and must not be substituted into the transmission kernel. Its role is age-of-onset once a genotype has been inherited. The general form, a baseline function of time multiplied by a covariate term, is proportional-hazards regression: h_l(t | x) = h_0(t) exp(gamma^T x). [cox1972hazard]

# 21. Epigenetic state as a reversible process

DNA methylation acts in a context-dependent way at promoters, gene bodies and regulatory elements, with promoter-island methylation associated with transcriptional repression, and the relationship is explicitly not a simple switch. [jones2012methylation] Histone modifications including acetylation regulate chromatin as a responsive scaffold with mark-specific transcriptional consequences. [bannister2011histone]

Model a CpG site as a two-state continuous-time Markov chain with methylation rate alpha and demethylation rate beta:

```equation
dP_M/dt = alpha (1 - P_M) - beta P_M = alpha - (alpha + beta) P_M
```

> **RESULT 21.1 - EXPONENTIAL APPROACH TO EQUILIBRIUM**
>
> ```equation
> P_M(t) = alpha/(alpha+beta) + [ P_M(0) - alpha/(alpha+beta) ] exp( -(alpha+beta) t )
> ```
>
> The methylation fraction approaches alpha/(alpha+beta) monotonically with time constant 1/(alpha+beta).
>
> *Verification.* At t = 0 the bracket cancels the equilibrium term, returning P_M(0); as t grows the exponential vanishes, leaving the equilibrium; differentiating recovers the original equation.

Averaging across a locus's sites gives the observable m_l(t). That methylation carries reproducible age information is established: a multi-tissue age predictor was built from approximately 8,000 samples across 51 healthy tissues and cell types [horvath2013clock], and a quantitative ageing model over more than 450,000 CpG markers in whole blood from 656 individuals aged 19 to 101 measured individual differences in methylome ageing rate [hannum2013aging]. Neither licenses the claim that methylation *causes* ageing.

> **OBSERVATION 21.2 - MUTATION IS ABSORBING, METHYLATION IS ERGODIC**
>
> A mutation fixed in a cell lineage does not revert at a comparable rate: the process is effectively absorbing and P(unmutated) decreases monotonically towards zero. Methylation has strictly positive rates in both directions: the chain is ergodic and settles strictly inside (0,1).
>
> Consequence: an intervention can in principle reverse epigenetic silencing and cannot reverse a mutation. A model representing both with one operator discards the only actionable difference between them.

# 22. The expression gate: weights and bias

The three routes now combine into one scalar per locus. This is the construction Section 5.1's weight column has been building towards.

> **RESULT 22.1 - THE EXPRESSION GATE**
>
> ```equation
> z_l(t) = b_l + w_age t + w_met m_l(t) + w_ac h_l(t) + w_dos g_l + sum over c of w_c x_c
> pi_l(t) = sigma( z_l(t) ) = 1 / (1 + exp(-z_l(t)))
> ```
>
> pi_l(t) is the probability that locus l is transcriptionally competent at age t. The bias b_l is the log-odds of expression for a reference individual at reference age, reference epigenetic state and reference dosage: the locus's baseline propensity to be expressed, carrying everything the covariates do not explain.

This is a generalised linear model with a binomial response and a logit link. [nelder1972glm] The score is nonlinear in age, through m_l(t) and h_l(t) which are themselves nonlinear by Result 21.1, but it is **linear in the parameters**, and that single property determines everything in Section 23.

Two weights have signs constrained in advance, which converts them from free parameters into a falsification test:

| Weight | Predicted sign | Basis | If fitting returns the other sign |
|---|---|---|---|
| w_met | negative | Promoter-island methylation associated with repression [jones2012methylation] | The model is wrong, the sites are gene-body not promoter, or the annotation is misassigned |
| w_ac | positive | Activating acetylation marks open chromatin [bannister2011histone] | The same three candidates, in the same order |
| w_age | unconstrained | Residual age effect after m and h are accounted for | A large value indicates the epigenetic covariates do not capture the age dependence |

The third row is the diagnostic. If m and h genuinely mediate age, conditioning on them should drive w_age towards zero; a large surviving w_age says the mediation story is incomplete.

> **RESULT 22.2 - AGE-DEPENDENT EFFECTIVE PENETRANCE**
>
> If the phenotype requires functional product from at least one allele, and each inherited allele is independently competent with probability pi_l(t), then for c_l functional alleles by inheritance
>
> ```equation
> P(expressed | c_l functional alleles, age t) = 1 - (1 - pi_l(t))^(c_l)
> ```
>
> Penetrance is therefore a function of age and epigenetic state, not a constant of the genotype. A heterozygote whose single functional allele is silenced behaves at that locus and age as though it carried none.

That consequence needs restraint. It follows from the stated model. It is not a claim that any carrier of any variant will develop any condition, and nothing here is evidence for such a claim.

# 23. Estimation: objective, gradient, convexity, step size

Given samples (t_i, m_i, h_i, g_i, y_i) with binary expression label y_i, write x_i for the feature vector and pi_i = sigma(w^T x_i). The objective is the mean negative log-likelihood, that is, cross-entropy:

```equation
J(w) = -(1/M) sum over i of [ y_i log pi_i + (1 - y_i) log(1 - pi_i) ]
```

Squared error is not used: under a binary response it is non-convex after composition with the logistic link, and cross-entropy is the log-likelihood of the assumed response distribution. [nelder1972glm]

> **RESULT 23.1 - GRADIENT**
>
> Using sigma'(z) = sigma(z)(1 - sigma(z)), for one sample:
>
> ```equation
> (1)  d/dw [ -y log pi - (1-y) log(1-pi) ] = [ -(y/pi) + (1-y)/(1-pi) ] dpi/dw
> (2)  dpi/dw = pi (1 - pi) x
> (3)  = [ -(y/pi) + (1-y)/(1-pi) ] pi (1-pi) x
> (4)  = [ -y(1-pi) + (1-y) pi ] x
> (5)  = ( pi - y ) x
> ```
>
> Averaging: grad J(w) = (1/M) X^T (pi - y).
>
> The gradient is the feature matrix applied to the residuals. The factors pi(1-pi) cancel exactly at step (4), which is what makes the logit link canonical for this response.

> **RESULT 23.2 - CONVEXITY**
>
> With S = diag(pi_i (1 - pi_i)), the Hessian is H = (1/M) X^T S X. For any v,
>
> ```equation
> v^T H v = (1/M) sum over i of pi_i (1 - pi_i) (x_i^T v)^2 >= 0
> ```
>
> so H is positive semi-definite and J is convex: every stationary point is a global minimum, unique when X has full column rank. Unlike the linear case H depends on w, so there is no closed-form minimiser and estimation is necessarily iterative.

> **RESULT 23.3 - STEP-SIZE BOUND**
>
> Since pi(1-pi) attains its maximum 1/4 at pi = 1/2, a single sample's Hessian satisfies lambda_max <= ||x||^2 / 4. Gradient descent on a convex function with L-Lipschitz gradient converges for 0 < eta < 2/L, so
>
> ```equation
> 0 < eta < 8 / ||x||^2
> ```
>
> This is derived before any computation and is falsifiable in the ordinary way: step sizes below it should not diverge, and step sizes far above it should.

An online form, consuming one individual at a time, follows with M = 1: w <- w - eta (pi_i - y_i) x_i. Adding a ridge penalty makes J strongly convex and bounds the conditioning, at the cost of biasing estimates towards zero.

# 24. Numerical consequences for the existing engine

With eps_l of order 1e-8, previously structural zeros become entries of order 1e-8 per mutated locus, and products across many loci fall below binary64 range quickly. Section 9.3 records that population updates refuse detected zero underflows rather than returning silent zeros, and that log-space population inference is not implemented.

Enabling mutation makes that refusal path routine rather than exceptional. **A log-domain population update is therefore a prerequisite for the extension of Section 18, not an optional refinement.** Combined with Result 19.1, the practical conclusion is that admitting mutation requires both a different numerical domain and a different storage strategy, and that neither is a small change.

# 25. Predictions and failure modes

| Proposition | Derived from | Would be falsified by |
|---|---|---|
| Germline de novo count rises linearly with paternal age at roughly 1.5 per year | Section 18.1 [jonsson2017trios] | A trio cohort with a materially different or nonlinear slope |
| A positive mutation rate makes the kernel structurally dense | Result 19.1 | Nothing; it is a counting argument. Its practical severity is measurable |
| CSR loses its payload advantage once mutation is enabled | Result 19.1 with Section 10 | A benchmark where CSR still wins, meaning eps was effectively zero |
| Methylation approaches equilibrium exponentially, not linearly | Result 21.1 | Methylation tracking a straight line across a wide age range with no curvature |
| w_met is negative at promoter islands | Section 22 [jones2012methylation] | A fitted positive coefficient on correctly annotated promoter sites |
| Conditioning on m and h shrinks w_age towards zero | Section 22 | A large surviving age coefficient, indicating unmodelled mediation |
| Age-of-onset hazard has log-log slope k-1 | Result 20.2 [armitage1954multistage] | An incidence curve inconsistent with any integer k |

<!-- pagebreak -->

# PART VII - EVALUATION

# 26. What this work supports

The corrected mathematical objects represent every declared Mendelian outcome without requiring the transmission kernel to be square. Exact-rational reconstruction explains the original displayed outputs and identifies divergences between prose, tables and software. Independent enumeration verifies the small-system implementation and the support formulas. Matched local benchmarks measure how the chosen dense, CSR, hash and streamed implementations behave on one stated machine. Factored queries and score dynamic programming avoid materialising outputs the scientific question does not require.

The contribution combines audit, derivation, implementation and evaluation. It does not claim that sparse matrices, hash maps or factorisation are new, and it claims no priority for the counting identities beyond the derivations given here.

# 27. Counterarguments retained

1. The original ABO and Rh software already computes all pairings; a coverage-restoration headline would misrepresent the baseline.
2. Dense arrays can be faster at small problems, because contiguous compiled operations outweigh sparse overhead.
3. Sparse storage remains exponential; success at five loci does not extrapolate to thousands.
4. Factored-query speed concerns a different output from full enumeration.
5. Population factorisation can be wrong even when per-locus transmission factorises correctly.
6. The rs334 panel is an ascertained reference resource; its calls are neither diagnoses nor an unbiased sample of humanity.
7. Exact Hardy-Weinberg p-values assess a narrow null conditional on sampled allele counts and identify no cause of departure.
8. A simulated 90% interval with known parameters does not demonstrate clinical calibration.
9. Byte budgets and caught allocation exceptions are safeguards, not proof of operating-system-level memory safety.
10. Part V is unimplemented throughout, and its two named limitations are unresolved.

# 28. Identifiability: the hardest limitation

This limitation is stated here rather than in Part V's body because it governs whether the extension can ever be estimated, not merely whether it is currently built.

By Result 21.1, m_l(t) is a deterministic function of age up to noise. The design matrix therefore carries a column t and a column m_l(t) that is a smooth monotone transform of it. Over a narrow age range the exponential is close to linear, so the two columns are close to collinear.

> **RESULT 28.1 - COLLINEARITY MAKES THE WEIGHTS UNIDENTIFIABLE WHILE THE FIT LOOKS HEALTHY**
>
> Near-collinearity makes X^T S X near-singular. The estimator covariance, which scales with the inverse of that matrix, becomes enormous in the direction of the offending combination, so w_age and w_met are individually meaningless while their sum remains well determined.
>
> Crucially the objective value does not deteriorate. The model fits, and the coefficients that were the reason for building it are noise. **No diagnostic computed from the fit alone distinguishes this case from a well-conditioned one.**

The defences are structural and must be chosen before fitting:

- **Check conditioning before estimating.** The smallest eigenvalue of X^T S X, or the condition number, is computable in advance and decides whether the fit can answer the question at all.
- **Require within-age variation in methylation.** If every 40-year-old in the sample has the same m_l, the coefficient is unidentified at any sample size. The sampling design must break the collinearity; enlarging it will not.
- **Adjust for cell composition.** Bulk-tissue methylation changes with age partly because the cell-type mixture changes with age. Unadjusted, that confound is attributed to the locus.
- **Do not read causation from the sign.** Transcriptional state can drive methylation as well as follow it, and a regression orients no arrow.

# 29. Boundaries of the implementation

The core catalog assumes autosomal diploidy and unordered allele pairs. It does not handle sex-linked dosage, imprinting, aneuploidy, copy-number variation, somatic mosaicism, penetrance uncertainty, arbitrary pedigrees or large linked haplotype blocks. The linkage extension covers two phased biallelic loci with a supplied recombination fraction. Mutation is a supplied gamete transition process whose rates are not estimated. The score dynamic programme supports non-negative integer effects and independent dosage factors.

These boundaries are explicit because a data structure cannot compensate for a missing biological state variable. A fuller molecular catalog can be inserted only when its inheritance and phenotype rules are specified and independently checked, and unknown parameters should be modelled as uncertain quantities rather than replaced by unexplained noise.

# 30. Conclusion

Version I asked whether linear algebra could replace repeated Punnett squares, and answered correctly that it could. It then identified the cost: forcing the operator to be square, so that it could be diagonalised and raised to a power, appeared to force the model to represent less of the population as more loci were added.

This unified treatment establishes three things about that cost. First, it was smaller than stated: the executed programs already enumerated every parental pair, and the truncation lived in the prose and displayed tables. Second, it was never necessary: the transmission kernel is a rectangular map requiring no inverse, and the square form was a consequence of wanting matrix powers rather than a property of inheritance. Third, removing it is not free: sparse storage remains exponential, dense arrays remain competitive at small n, and the population update is nonlinear whenever mating is closed.

The strongest contribution is not the sparse kernel but the separation it forced: transmission from mating, catalog coverage from probability mass, probability mass from predictive accuracy, and germline from somatic change. Version I's conclusion moved between those quantities. Keeping them apart is what makes the remaining claims defensible and the remaining gaps visible.

Two gaps are worth naming last, because they bound what should be attempted next. Admitting mutation destroys the structural sparsity the architecture depends on, which means the honest successor to this work is factored or structured, not sparse. And age is collinear with its own epigenetic mediators, which means the expression gate of Part VI cannot be estimated from a convenience sample at any size, and needs a design built specifically to break that collinearity.

The next empirical priority is unchanged and unglamorous: an independent parent-offspring transmission dataset with a preregistered evaluation target. Everything in Part VI should wait for it.

<!-- pagebreak -->

# Availability, disclosure and scope

**Code and data.** All source code, tests, frozen datasets with their provenance and digests, raw benchmark output, figures, the source ledger and the build that generates this document are in the accompanying repository. Every number printed here is bound at build time to a retained machine-readable result; none is typed by hand. The reproduction sequence is in Appendix B.

**External data.** Two frozen snapshots of public aggregate allele and genotype counts are used: rs334 and rs12913832, each for the same 2,504 individuals across 26 populations of the 1000 Genomes phase 3 reference panel, retrieved through the Ensembl REST API and stored with SHA-256 digests. [ensembl_rs334_phase3_20260905] [ensembl_rs12913832_phase3] [{{1000G_SOURCE}}] These are aggregate counts. No individual-level genotype and no personal data of any kind is used, held or distributed by this work.

**No human subjects.** No participants were recruited, no samples were collected and no ethical approval was required. All phenotype figures quoted from cohorts are published summary statistics, cited at their point of use.

**AI assistance.** Computation, figure generation, manuscript preparation and adversarial review were assisted by AI agents. Scientific authorship, interpretation and responsibility for every claim are the author's.

**Scope of claims.** This work establishes representational completeness, numerical agreement between four implementations, measured performance on one machine, and asymptotic bounds. It establishes **no predictive accuracy for any trait**, because no held-out phenotype data was available for any model built here. Part VI is derivation only: no code, no test, no dataset corresponds to it. Nothing in this document is medical, diagnostic or reproductive advice.

**Correspondence.** arshyiamehran@gmail.com

<!-- pagebreak -->

# Appendix A. Operational definitions and small examples

The callable interface accepts allele counts and mixed-radix genotype codes. For one biallelic locus the canonical order is (0,0), (0,1), (1,1), giving codes 0, 1, 2. At a triallelic locus it is (0,0), (0,1), (0,2), (1,1), (1,2), (2,2). This canonical order differs from Version I's presentation order; the reproduction module retains the legacy order for comparison.

```python
from genetics import InheritanceModel
model = InheritanceModel((2,))
model.cross(1, 1)                  # {0: 0.25, 1: 0.5, 2: 0.25}
model.kernel('csr').shape          # (3, 6)
model.next_generation([.25,.5,.25])

large = InheritanceModel((2,) * 1000)
parent = large.encode((1,) * 1000)
large.log_probability(parent, parent, 0)
# Finite log probability; no 3^1000 child table is constructed.
```

Population inputs must be non-negative and normalised within the specified tolerance. The next-generation method treats parents as independently sampled from that full distribution. A supplied kernel must match the catalog's G by U shape. Factored population updates require explicit per-locus inputs and the linkage-equilibrium assumption; they do not infer independence from an arbitrary joint vector.

# Appendix B. Reproduction

{{RUN_COMMANDS_UNIFIED}}

{{ARTIFACT_TABLE}}

The checks above concern delivered evidence. They do not assert that MATLAB was executed, that a clinical study was performed, that the software was benchmarked on more than one machine, or that this preprint has been peer reviewed.

# Appendix C. Summary of derived results

| Result | Statement | Section |
|---|---|---|
| 4.1 | Eigenvalue decay is basis-dependent; pair-to-genotype relabelling supplies no generational operator | 4.3 |
| 4.2 | Complete reproductive exclusion gives reciprocal, not exponential, allele decline | 4.5 |
| 5.1 | ABO has 21 unordered parental pairs; 15 omits equal-genotype pairs | 5.2 |
| 6.1 | Recessive phenotype frequency 0.07 implies allele frequency sqrt(0.07) | 6.1 |
| 6.2 | The executed programs were already complete; no coverage-restoration claim is available | 6.3 |
| 7.1 | G(n) = 3^n by induction | 7.1 |
| 7.2 | G = product of a(a+1)/2 over loci | 7.2 |
| 7.3 | U = G(G+1)/2, asymptotically 9^n/2 | 7.3 |
| 7.4 | Nonzeros = (15^n + 5^n)/2 | 7.4 |
| 7.5 | Arbitrary-allele support: half the product of T(a) plus the product of D(a) | 7.4 |
| 8.1 | The kernel is rectangular and needs no inverse | 8.1 |
| 8.2 | No fixed linear map represents closed random mating | 8.2 |
| 10.1 | The derived exponent is confirmed at the top of the measured range, not across it | 10.2 |
| 11.1 | Pooling changes the prediction before any biology is invoked | 11.4 |
| 12.1 | Removing squareness buys completeness, not accuracy | 12.2 |
| 12.2 | The unphased catalog does not determine transmission for linked loci | 12.3 |
| 13.1 | The blue-associated allele is identified by population ordering, not by label | 13.2 |
| 13.2 | At the eye-colour catalog size sparsity is not worth its overhead | 13.4 |
| 13.3 | The eye-colour distribution is already at its Hardy-Weinberg fixed point | 13.6 |
| 13.4 | A declared modifier frequency underpredicts the observed GG discordance | 13.7 |
| 13.5 | The gap to a real predictor is not a gap in representation | 13.7 |
| 14.1 | The score distribution is polynomial where the catalog is exponential | 14.2 |
| 14.2 | What a real height model would additionally require | 14.4 |
| 15.1 | Output-size lower bounds for P(n): Omega(15^n), Omega(3^n), Theta(n), Omega(B) | 15.3 |
| 16.1 | Germline and somatic routes are different objects | 16 |
| 18.1 | Age-parameterised gamete mutation operator | 18.3 |
| 18.2 | Normalisation is preserved | 18.3 |
| 19.1 | A positive mutation rate destroys structural sparsity | 19 |
| 20.1 | Constant intensity gives exponential, not linear, survival | 20 |
| 20.2 | Multistage hazard has log-log slope k-1 | 20 |
| 21.1 | Methylation approaches equilibrium exponentially | 21 |
| 21.2 | Mutation is absorbing; methylation is ergodic | 21 |
| 22.1 | The expression gate | 22 |
| 22.2 | Age-dependent effective penetrance | 22 |
| 23.1 | Gradient is the feature matrix applied to residuals | 23 |
| 23.2 | The cross-entropy objective is convex | 23 |
| 23.3 | Step-size bound eta < 8/||x||^2 | 23 |
| 28.1 | Collinearity makes weights unidentifiable while the fit looks healthy | 28 |

# Appendix D. Version I artefacts preserved

The original 54-page study and its DOCX source are retained unedited alongside this work. Its SHA-256 digest before the extension was EE3CB04CE28D4734669D5822FDDBDF03C176FD5C2F5FF3760CB82A4BD1145ECD. [version1_local]

Preserved without alteration: the six sickle-cell Punnett squares and the six-column parental table; the reduced three-column table and its linear system; the eigenvalues 1, 1/2, 1/4 with their eigenvectors and the diagonalisation M = S D S^-1; the ABO allele-weight scheme and genotype scores; the Hardy-Weinberg derivation of the ABO allele frequencies and the enumeration of all 21 pairings; the ranked top-six table and the 6 x 6 offspring matrix; the Rh model and the ranked top-18 joint table; and the 18 x 18 displayed matrix.

Version I's own bibliography is carried into the source ledger in full, with each entry marked by the role it played and by what it does not establish; its two secondary population summaries are identified as the origin of the phenotype inputs that drive the coverage figures. Its five MATLAB programs are reproduced verbatim in Appendix F.

Revised, with the original retained for comparison: the interpretation of eigenvalue decay (Result 4.1); the normalisation of the generational algorithm (Section 4.4); the pair count 15 versus 21 (Result 5.1); the provenance of the coverage percentages (Section 5.5); and the description of the 18 x 18 table (Section 6.3).

# Appendix E. Protocol for adding a trait

Section 15.5 states this protocol in the body; it is repeated here as a checklist. An addition is accepted when it supplies, in order and before any result is reported:

1. **Catalog.** Allele counts per locus and the resulting G by Result 7.2, with the feasibility check of Section 9.3 run at that G.
2. **Map.** The genotype-to-phenotype function, including epistasis and any masking, written separately from the transmission rule so that Result 12.2's distinction survives.
3. **Linkage.** Whether loci are linked; if so, recombination fractions with their source, and the phase convention.
4. **Effects.** The origin of any effect sizes, their licence, and the population they were estimated in, recorded in the source ledger with a `does_not_establish` field.
5. **Evaluation.** A held-out target chosen before fitting, the population it generalises to, and the metric, distinguished from coverage per Section 12.1.
6. **Conditioning.** The smallest eigenvalue or condition number of the design, reported before any coefficient, per Result 26.1.
7. **Oracle.** An independent enumeration or closed form for at least one small case, so the implementation can be checked without trusting itself.

A trait that cannot supply items 4 to 6 may still be added as a *derivational* section, clearly marked as Part VI is, but may not report fitted coefficients or accuracy. Eye colour was added under this protocol and satisfies items 1 to 3 and 7 only; height satisfies none of items 4 to 6 and is therefore derivational. Both say so at their own point of use.

# Appendix F. Version I source listings

The five MATLAB programs of Version I, reproduced verbatim from `legacy/matlab/` with no edit of any kind. They are the evidence on which Sections 4.4, 5.5 and 6.3 rest, and they are printed here so that the audit can be checked against the code rather than against a description of it. SHA-256 digests of the embedded text are recorded in `results/unified_manifest.json`. Lines longer than the page width are wrapped by the typesetter; the files themselves are unmodified. A Python reimplementation of the same models is held separately in `python/` and was not used for the reconstruction, which is an independent exact-rational translation described in Section 6.3.

{{MATLAB_APPENDIX}}

# References

{{BIBLIOGRAPHY}}
