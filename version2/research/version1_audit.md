# Version 1: source audit and independent numerical reproduction

## Scope and evidence status

The complete 54-page **Bioinformatics-Arshyia Mehran.pdf** was read before the reproduction implementation was written. Text was extracted from every page; every page was rendered and visually inspected, including the embedded mathematical derivations, MATLAB listings and six execution examples. All page references below are the PDF's printed page numbers, which match its physical page numbers. The original was not edited.

- Original size: 3,282,612 bytes.
- Original SHA-256: `ee3cb04ce28d4734669d5822fddbdf03c176fd5c2f5ff3760cb82a4bd1145ecd`.
- Audit script: `version2/experiments/reproduce_version1.py`.
- Machine-readable evidence: `version2/results/version1_reproduction.json`.
- Full pair rankings and trajectories: `version1_pair_coverage.csv` and `version1_trajectories.csv` in the same results directory.
- Literal numeric table transcriptions: `version1_abo_displayed_matrix.csv` and `version1_joint_displayed_matrix.csv`.
- Verification: **8 audit tests passed** using Python's unittest framework. All six screenshot examples agree with the independent Python computation within their printed precision.

This is an **executed Python translation of displayed MATLAB mathematics**, not a MATLAB execution, MATLAB performance measurement, or biological validation. Exact rational arithmetic was used for the inheritance examples; binary64 arithmetic was used where the Rh calculation requires a square root. The interactive input prompts and plot styling were not reproduced. Code was transcribed from images rather than supplied `.m` files. Poppler reported missing display fonts for Symbol and ArialUnicode; some genotype labels in the page 24 graphic are blank, but the numerical curves, matrices and execution examples needed for this audit are legible.

## What the paper establishes, and what the extension should retain

The paper's strongest progression is from a concrete inheritance question to Punnett enumeration, a mathematical representation, executable simulations and an evaluation of limitations. The six sickle-cell parental crosses (pp. 12-14) and six selected ABO crosses (pp. 34-35) give the expected one-generation Mendelian distributions. The eigenvalues and eigenvectors of the specific triangular sickle-cell matrix are computed correctly (pp. 16-23). More significantly, the actual ABO and ABO+Rh MATLAB listings already implement full random mating rather than the truncated square-matrix procedure described in the prose. That sound part of the executable work should be credited explicitly instead of treating the entire first version as the same incorrect model.

The extension can preserve the order **theory and biological definitions -> derivation -> implementation -> execution -> evaluation**, while distinguishing conditional inheritance, a mating rule, a population update and an observation map. Numerical agreement with a screenshot establishes reproducibility of that computation; it does not establish that its modeling assumptions represent a real population or a family's future partners.

## Sickle-cell computation: correct first cross, incompatible subsequent update

Pages 14-15 construct

```text
              parental-pair columns
              AA x AA   AA x AS   AS x AS
offspring AA      1        1/2       1/4
offspring AS      0        1/2       1/2
offspring SS      0         0        1/4
```

Every column is a valid offspring distribution for its stated parental pair. A vector of weights over those **pairs**, however, is different from a vector of frequencies over the individual genotypes **AA, AS, SS**. Page 15 identifies the latter with the former when it writes `x_next = M x`. Having three entries in both spaces does not make their meanings interchangeable. A mating map from individual frequencies to pair frequencies is missing.

Pages 26-27 add a second operation that also matters: the program first records `M*x`, then sets its SS entry to zero, and carries the remaining entries forward **without normalization**. It does not remove SS from the initial generation before the first matrix multiplication. The displayed matrix comment additionally lists AA x AA, AS x AS, SS x SS, which does not match the matrix's actual second and third columns.

For the page 28 example, both parents are `HbAHbS`, and the requested generation is 3. The exact reproduction is:

| Generation | AA (%) | AS (%) | SS (%) | Sum (%) |
|---|---:|---:|---:|---:|
| 1 | 25 | 50 | 25 | 100 |
| 2 | 56.25 | 37.5 | 6.25 | 100 |
| 3 | 75 | 18.75 | 0 | 93.75 |
| 4 | 84.375 | 9.375 | 0 | 93.75 |

The first three rows reproduce the screenshot, including its omitted 6.25% mass at generation 3. These are not normalized genotype frequencies after that point. If the matrix were simply applied twice to generation 1, without the removal step, generation 3 would instead be **76.5625%, 21.875%, 1.5625%**, summing to 100%. Thus the page 28 assertion that its output is verified by raising the same matrix to a power does not describe the implemented algorithm.

The derivation on pages 20-23 obtains a valid diagonalization for the numeric matrix:

```text
S = [[1, 1, 1], [0, -1, -2], [0, 0, 1]],  S^-1 = S,
D = diag(1, 1/2, 1/4),  M^n = S D^n S^-1.
```

But the entries of `D^n` are in the eigenvector coordinate system, not the original genotype coordinates. The page 23 table labels the diagonal entries directly as parental-pair to offspring probabilities; its nonfirst columns consequently do not sum to one. The limit of `D^n` is `diag(1,0,0)`, as shown on page 25, while the limit of **M^n** is `[[1,1,1],[0,0,0],[0,0,0]]`. These are different matrices. Algebraic decay inside this chosen operator does not prove decay of an allele in a population with an unspecified mating rule.

The finite-date elimination claims on page 52 also do not follow from asymptotic decay: `(1/4)^4 = (1/2)^8 = 1/256`, or **0.390625%**, rather than zero. Population size and a stochastic extinction model would be needed to translate a small frequency into the absence of individuals. The geographic assumptions also move from Africa (p. 4) to Canada (p. 25), and the zero-reproduction assumption is later discussed as potentially inappropriate (p. 52). Those are alternative scenarios requiring explicit parameterization, not interchangeable empirical facts.

## ABO computation: the listing already uses the complete random-mating model

Page 30 says there are 15 nonredundant parental genotype pairs among six genotypes. There are **21** unordered pairs when equal-genotype parents are allowed: `6*(6+1)/2`. Fifteen counts only the unequal pairs, `6*(6-1)/2`. The paper's own table on pages 32-33 correctly enumerates 21, including the six diagonal pairs. There are 36 ordered pairs.

The page 35 six-by-six table has offspring rows `AA, AO, BB, BO, AB, OO` and columns `AO x OO, OO x OO, AO x AO, BO x OO, AO x BO, AA x OO`. Its entries are correct conditional crosses. The row and column state types differ, so this table alone is not an iteratable genotype transition matrix. Square dimensions neither supply a mating rule nor guarantee diagonalizability.

The actual MATLAB code on pages 37-39 does something better. It initializes generation 1 with the user-specified Punnett cross, then loops over **all 36 ordered genotype pairs** with weight `x(i)*x(j)` and accumulates their offspring distributions. It does not use the page 35 six-column table or the world-frequency coverage calculation. This is a quadratic population update, not a fixed square matrix power.

| PDF example | Parents | Requested generation | Reproduced genotype percentages, order AA/AO/BB/BO/AB/OO | Phenotypes A/B/AB/O (%) |
|---|---|---:|---|---|
| p. 35 | AA x OO | 3 | 25, 50, 0, 0, 0, 25 | 75, 0, 0, 25 |
| p. 36 | AO x BO | 4 | 6.25, 25, 6.25, 25, 12.5, 25 | 31.25, 31.25, 12.5, 25 |

Both screenshots are reproduced exactly. Their genotype distributions become constant from generation 2 in these idealized examples. As a direct contradiction test, iterating the displayed six-column table from the p. 35 generation-1 vector gives **100% AO** at generation 3, not the screenshot's 25/50/25 distribution. This proves that the table and the execution example are different algorithms.

## Coverage claims and rounding provenance

For an assumed genotype-frequency vector `f`, random mating gives unordered pair weights `f_i^2` on the diagonal and `2*f_i*f_j` off the diagonal. Their sum is `(sum f)^2`. Coverage below is the probability mass of selected **parental-pair events under those assumptions**, not an observed percentage of all humans.

The paper contains three distinct precision levels for its ABO inputs. They yield different results:

| Input used | Sum of genotype entries | ABO top-six pair mass (%) | ABO+Rh top-18 pair mass (%) |
|---|---:|---:|---:|
| Exact p. 32 allele values: A=.26, B=.077, O=.663 | 1 | **83.4764828157** | **67.8347334921** |
| P. 32 genotype table: .0676, .3448, .0059, .1021, .04, .4396 | 1 | 83.489176 | 67.8455062501 |
| Pp. 32-33 pair-formula inputs: .068, .345, .006, .102, .04, .44 | **1.001** | **83.6205** | 67.9390313612 |
| Same coarser inputs normalized to sum one | 1 | 83.4535095274 | 67.8033568442 |

The paper's **83.62%** is reproducible as the rounded sum using its coarser pair-formula inputs. Those inputs have total pair mass **1.002001**, or 100.2001%, so the result is not exactly the coverage of a normalized probability model. Using its exact stated allele values gives **83.48%** to two decimals. The six selected pairs remain the same in every variant; the discrepancy is numerical provenance rather than a changed ranking.

For Rh, the paper's 7% negative figure on page 41 is a phenotype frequency in its own recessive model. The appropriate inference under Hardy-Weinberg assumptions is:

```text
q_d = sqrt(.07) = .2645751311064591
P(DD) = .5408497377870819
P(Dd) = .3891502622129182
P(dd) = .07
```

Combining these values with the exact ABO allele inputs, and assuming population independence of the two loci, reproduces **67.8347334921%**, hence the paper's **67.83%**, for the top 18 of 171 unordered joint-genotype pairs. Every row of the page 42 table matches this computation within 0.005 percentage points. Its individually rounded displayed rows sum to **67.84%**, while the total rounded only once is 67.83%; this is an ordinary rounding discrepancy. Minor order changes within exact ties do not change coverage.

The exact ABO allele values imply phenotype percentages **A=41.236, B=10.8031, AB=4.004, O=43.9569**. They do not exactly reproduce the p. 31 source percentages A=42, B=10, AB=4, O=44. Their derivation or fitting criterion is not supplied. Population pooling, Hardy-Weinberg equilibrium, independence between loci and the claimed global source frequencies must be evaluated separately; the arithmetic above does not validate them.

## ABO+Rh: self-cross table, random-mating code

The joint genotype ordering is ABO-major and Rh-minor:

```text
AA/DD, AA/Dd, AA/dd, AO/DD, AO/Dd, AO/dd,
BB/DD, BB/Dd, BB/dd, BO/DD, BO/Dd, BO/dd,
AB/DD, AB/Dd, AB/dd, OO/DD, OO/Dd, OO/dd.
```

The 18-by-18 table on pages 43-45 is described in prose as the top 18 parental pairs. Its actual columns are single genotype labels. A literal transcription of all **324 entries** agrees exactly with `K[:,g] = Punnett(g,g)`: each column represents two parents with the same joint genotype. In the printed column order, it differs in **114 entries** from the offspring kernels of the ranked pairs on page 42. It is therefore a self-cross operator, not the advertised top-18-pair table.

The MATLAB listing on pages 46-49 does not multiply that table. Despite its initial comment saying 'self-cross transition matrix', the nested loop uses all **324 ordered joint-genotype pairs**. `punnett18` splits ABO and Rh, computes both local crosses and forms their Kronecker product. This is a full random-mating update under conditional independent segregation. The genotype-to-phenotype summation correctly treats `dd` as negative, whereas the prose on p. 40 incorrectly says 'Dd makes Rh negative'. The 25% negative-offspring statement on p. 40 applies to Dd x Dd, not every pair of Rh-positive parents.

| PDF example | Parents | Generation | Reproduced phenotype percentages, order A+/A-/B+/B-/AB+/AB-/O+/O- |
|---|---|---:|---|
| p. 49 | AO/Dd x AO/Dd | 3 | 56.25, 18.75, 0, 0, 0, 0, 18.75, 6.25 |
| p. 50, upper | AA/dd x BO/Dd | 7 | 21.875, 28.125, 8.203125, 10.546875, 10.9375, 14.0625, 2.734375, 3.515625 |
| p. 50, lower | AB/dd x OO/Dd | 8 | 13.671875, 17.578125, 13.671875, 17.578125, 5.46875, 7.03125, 10.9375, 14.0625 |

All 18 genotype entries and all eight phenotype entries match the three screenshots within their two-decimal display rounding. The full exact vectors and trajectories are retained in JSON and CSV. The p. 49 distribution is already stationary at generation 1. The two p. 50 examples become stationary at generation 2. Iterating the displayed self-cross matrix would give only **1.5625% AO/Dd** at generation 3 in the p. 49 example, whereas the screenshot and full random-mating implementation give **25%**. Therefore the p. 49 matrix-power verification claim is inconsistent with the actual code and displayed matrix.

## Further revisions indicated by the complete-paper review

- **Mathematical notation, pp. 7-9 and 17-22:** retain the valid determinant and eigenvector examples, but correct `(M-lambda I)X`, the matrix-vector factor order. Eigenvalues are roots of the characteristic polynomial; their existence is not explained merely by counting diagonal entries. After row reduction of `[S|I]`, the inverse is on the right. `Mv=lambda*v`, not a scalar eigenvalue alone. The final numerical eigenvectors remain correct.
- **Biological exposition, pp. 5-6 and 29-30:** distinguish the number of alleles carried by a diploid individual from the number of alleles in a population. The arbitrary 0/1 allele score collapses AA, BB and AB to the same scalar, so it cannot substitute for a genotype representation. The score is not used in the displayed ABO code.
- **Empirical and clinical statements, pp. 10-11, 14, 25 and 51-53:** the strand terminology, carrier-cell explanation, malaria mechanism, survival/fertility assumptions and uncited 'about 10%' reproductive figure require primary-source verification. They are not established by the calculations and should not be silently carried into the extension as validated inputs.
- **Scaling, pp. 30, 39-43 and 53:** a rectangular inheritance kernel is valid. There is no mathematical requirement to discard parental pairs until it is square. With G genotype states there are G(G+1)/2 unordered pairs, quadratic in G; genotype counts themselves multiply across loci. Exact sparsity and factorization should address growth without probability pruning.
- **Evaluation, pp. 28, 39, 49 and 51-53:** separate reproduction, internal invariants, comparison against independent allele calculations, sensitivity to mating/selection assumptions, external population evidence and performance measurements. A mismatch between one realized child and a probabilistic prediction is not by itself model failure. A clinical counseling application would require evidence and scope beyond these idealized simulations.

## Reproduction commands

From the workspace root, using the bundled runtime that was executed for this audit:

```powershell
& 'C:\Users\DELL\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -X utf8 'version2/experiments/reproduce_version1.py'
& 'C:\Users\DELL\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -X utf8 -m unittest discover -s 'version2/tests' -p test_version1.py -v
```

The audit script has no third-party dependency. Tests intentionally preserve the sickle-cell legacy defect so that a future correction cannot be mistaken for an exact reproduction of Version 1.
