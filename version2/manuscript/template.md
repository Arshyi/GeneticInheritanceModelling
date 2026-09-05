# Complete Inheritance Without a Square-Matrix Constraint

## A computational and biological audit, followed by sparse and factored implementations

Version II research candidate | Extension of Arshyia Mehran's Version I | {{DATE}}

Prepared for author review. Original PDF and DOCX preserved. Computations, manuscript preparation, and adversarial checks were assisted by AI agents; scientific authorship, interpretation, and submission require the author's review. This document reports a completed set of local experiments, not a claim that the wider research programme or clinical validation is finished.

<!-- pagebreak -->

# Abstract

Version I explored whether linear algebra could replace repeated Punnett-square calculations for sickle-cell inheritance, ABO blood groups, and a simplified ABO+Rh system. Its discussion identified an apparent loss of population coverage when parental-pair tables were reduced to square matrices. We re-examine that premise before proposing an alternative architecture. An exact-rational Python reconstruction matches six displayed MATLAB examples. The audit finds that the ABO and ABO+Rh programs already evaluate every parental pair through nonlinear random-mating loops, while the prose and displayed tables describe different constructions. Thus, restoring coverage to those executed algorithms is not a valid claim for the extension.

We distinguish genotype catalogs, parental-pair catalogs, supported transmissions, retained probability mass, and predictive performance. For n biallelic diploid loci, an induction proves G=3^n possible unphased catalog entries and U=G(G+1)/2 unordered parental pairs. Under independent Mendelian segregation, the full rectangular kernel has G by U entries and exactly (15^n+5^n)/2 nonzero transitions. Dense arrays, compressed sparse rows, hash-based adjacency maps, and streamed computation are implemented with the same transmission semantics and benchmarked on one local machine. {{BENCH_ABSTRACT}}

External rs334 genotype-call counts from 2,504 reference-panel participants in 26 populations provide an audit of Hardy-Weinberg assumptions, not clinical predictive validation. A pooled model expects 1.874 rare homozygotes, whereas separately fitted population models expect 7.600 in total; none are recorded in the selected snapshot. A staged extension demonstrates simplified M/N inheritance, ABO-FUT1 epistasis, phase-dependent recombination, explicit population operators, and a 200-locus synthetic score distribution. The latter uses 401 bins and achieves 90.03% coverage for a nominal 90% interval in 10,000 independent simulated draws. The contribution is a reproducible separation of biological assumptions, mathematical objects, numerical guarantees, and workload-specific engineering results.

Keywords: Mendelian inheritance; sparse matrices; probability factorization; combinatorial growth; Hardy-Weinberg equilibrium; reproducibility; model audit.

# Reader's guide and research question

The intellectual progression remains that of Version I: first explain a biological question, then introduce the mathematics that represents it, implement the model, and evaluate what the outputs establish. The present question is: **which representation answers a specified inheritance query completely, correctly, and efficiently under a stated biological model?** The qualifying words matter. Completeness without a declared catalog is undefined; efficiency without a declared output is misleading; numerical agreement without external observations does not establish biological predictive accuracy.

The extension proceeds through four linked arguments. First, the original paper must be reconstructed as it exists, including discrepancies between its prose, tables, and software. Second, the parent-to-child mapping must be separated from the rule that constructs the next generation's parents. Third, storage and inference should exploit actual mathematical structure rather than discard low-probability outcomes. Fourth, increasingly realistic biological assumptions must be introduced only when their meaning can be tested.

Sections 1-3 audit Version I. Sections 4-7 derive the complete representation and evaluate the architecture. Sections 8-10 report population and trait experiments. Sections 11-13 discuss falsification, reproducibility, and remaining work. The appendices give executable examples, detailed results, and a source bibliography. Every reported benchmark number comes from retained machine-readable output. Theoretical extrapolations are identified separately.

# 1. Version I as evidence, rather than an unquestioned baseline

## 1.1 Preservation and reconstruction

The supplied 54-page PDF was read in full. Embedded code and result images were rendered and inspected because text extraction alone omits important MATLAB content. Its SHA-256 digest before the extension was EE3CB04CE28D4734669D5822FDDBDF03C176FD5C2F5FF3760CB82A4BD1145ECD. Both original files remain at the project root. All new work resides in the version2 directory. The source PDF provides the local source of record for the following audit; page references refer to its printed pages.

The reconstruction uses Python's exact rational arithmetic to preserve the probabilities in the screenshots. It is a faithful translation of the displayed algorithms, not an execution of MATLAB. MATLAB performance, graphics behavior, and executable compatibility are consequently unmeasured. Agreement with six examples supports the reconstruction of those examples; it does not verify every undocumented input or establish biological correctness.

Three different objects occur in Version I: a table of child outcomes conditional on a parental pair; a square matrix assembled from selected columns or same-genotype pairings; and a population algorithm that repeatedly combines current genotype frequencies. Their dimensions can look similar, but their state meanings differ. The first audit task is to keep those meanings separate.

## 1.2 Sickle-cell example: a valid local cross, an invalid generational interpretation

The simplified two-allele catalog is AA, AS, SS. An AS x AS cross gives probabilities 1/4, 1/2, 1/4. Version I's page 14 retains AA x AA, AA x AS, and AS x AS as three parental columns. With child rows AA, AS, SS, its matrix is:

```equation
        [1   1/2   1/4]
M  =    [0   1/2   1/2]
        [0    0    1/4]
```

Each column is a valid conditional offspring distribution. However, the input labels are parental pairs and the output labels are single child genotypes. Multiplying by M again requires a new rule that converts child genotypes into parental-pair probabilities. Merely having three entries on each side does not supply that rule.

The eigenvalues 1, 1/2, and 1/4 are algebraically correct. The diagonal matrix of eigenvalues is a representation in an eigenvector basis, not a table of genotype probabilities in the original basis. An nth-generation calculation must transform back to the original basis and must first justify the generational operator. Algebraic decay inside a chosen matrix does not establish disappearance of an allele in a population.

The displayed MATLAB program adds another operation: after recording a generation, it sets the SS component to zero before the following multiplication, without renormalizing. Starting from the AS x AS cross, the displayed generation-three output is AA=75%, AS=18.75%, SS=0%. Its sum is 93.75%. Those values can be reproduced, but they are not a normalized distribution of genotype proportions. Treating them as absolute surviving mass would require an explicit population-size interpretation and a justified transition law, neither of which is provided by the genotype-pair relabeling.

## 1.3 Biological scope of the sickle-cell simplification

AA/AS/SS is an educational two-allele model, not a complete catalog of HBB-related disease. Sickle-cell disease includes compound heterozygous forms involving other HBB variants. Trait status does not mean that half of a person's circulating red cells are permanently sickled. Nor is a universal assumption of zero reproduction among affected individuals supported by this computation. Version II therefore uses reproductive weights only as explicitly chosen model parameters and avoids forecasting disease eradication dates. [{{SCD_SOURCE}}]

The rs334 population experiment later in this paper uses genomic T/A variant-call labels. A reference T at this site does not establish that the rest of HBB is normal. This distinction prevents a single-variant catalog from being mistaken for a clinical diagnostic classification. [{{CLINVAR_SOURCE}}]

# 2. The ABO and ABO+Rh audit changes the contribution

## 2.1 Counting parental pairs correctly

The simplified ABO catalog has six genotypes: AA, AO, BB, BO, AB, OO. [dean_abo_2005] For six genotype labels, the number of unordered pairs with repetition is 6 x 7 / 2 = 21. The number 15 counts only distinct-genotype pairs and omits the six equal-genotype pairings. Version I states 15 on page 30, then correctly lists 21 pairings on pages 32-33. An equal-genotype pair means two parents who share a genotype; it does not imply self-fertilization.

For simplified ABO+Rh, six ABO states combine with DD, Dd, dd to give 18 catalog entries. The corresponding unordered pair count is 18 x 19 / 2 = 171, and the ordered pair count is 324. Rh positive/negative is treated here as a simplified D-presence model. The full Rh system includes RHD/RHCE variation and variants not captured by a single dominant/recessive pair. The statement that two Rh-positive parents have a one-quarter chance of a negative child is valid for Dd x Dd under this model, not for every positive-parent combination. [{{RH_SOURCE}}]

## 2.2 Reproducing the two coverage figures

Version I gives allele frequencies pA=0.26, pB=0.077, pO=0.663. These sum to one. Hardy-Weinberg genotype probabilities are pA^2, 2pA pO, pB^2, 2pB pO, 2pA pB, and pO^2 in its genotype order. To rank unordered parental pairs under independent random mating, equal-genotype pairs receive probability x_i^2 and distinct-genotype pairs receive 2x_i x_j.

The six leading ABO pairs are AO x OO, OO x OO, AO x AO, BO x OO, AO x BO, and AA x OO. The retained probability mass depends on which rounded inputs are used:

| Calculation | Retained top-six pair mass (%) | Interpretation |
|---|---:|---|
| Stated allele frequencies, propagated exactly | 83.4764828157 | Coherent normalized model |
| Four-decimal genotype table on page 32 | 83.489176 | Rounding variant; table sum requires checking |
| Coarser inputs used in the pair formulas | 83.6205 | Raw inputs sum to 1.001 |
| Same coarse inputs, normalized first | 83.4535095274 | Coherent normalized rounding variant |
| Printed six percentages added | 83.62 | Reproduction of displayed arithmetic |

The coarse genotype inputs yield total parental-pair mass (1.001)^2=1.002001. Consequently, the raw 83.6205% is not a rigorously normalized population probability. The discrepancy is small, but it illustrates why provenance and normalization must accompany any quoted percentage.

For Rh-negative phenotype frequency 0.07, the simplified HWE model gives d frequency sqrt(0.07), not 0.07. Multiplying the exact ABO and Rh genotype frequencies and selecting the top 18 unordered pairs yields 67.8347334921%, which rounds to the reported 67.83%. Independently rounded row percentages sum to 67.84%; rounding rows and rounding their unrounded total need not agree. These are model-based pair probabilities under the chosen allele frequencies and assumptions. They are not directly measured percentages of the world's population.

## 2.3 The displayed tables and executed algorithms are different

The six-column ABO table on page 35 contains the selected pairings. However, the MATLAB code on pages 37-39 loops over all six-by-six ordered parental genotypes and weights every cross by the current frequencies. Its output is a complete nonlinear random-mating calculation for its declared model. Its AA x OO generation-three example is AA=25%, AO=50%, OO=25%. Iterating the displayed six-column table instead would give a different result.

The 18-by-18 table on pages 44-45 has another interpretation. All 324 displayed entries match crossing each listed combined genotype with the same genotype. It is not the top-18-ranked-parent-pair table: 114 cells differ from that construction. The code on pages 46-49 again loops over all 324 ordered pairs, combining independent ABO and Rh transmissions. Its AO/Dd x AO/Dd example is stationary from its first offspring generation under the stipulated random-mating model.

This finding changes the legitimate comparison. Version II can show that the prose's truncation is unnecessary, and it can measure storage and inference alternatives to the complete algorithm. It cannot claim that the implemented ABO/Rh programs had 83.62% or 67.83% support and that a new sparse implementation raised those values to 100%. Within their simplified catalogs, the original executed loops already include all pairings.

{{LEGACY_EXAMPLE_TABLE}}

# 3. From local inheritance to a population process

## 3.1 The transmission kernel

For a diploid genotype g=(u,v), define its gamete probability t_g(a) as one-half the number of copies of allele a in g. For an unordered child genotype o=(a,b), the local transmission rule is t_g(a)t_h(a) when a=b, and t_g(a)t_h(b)+t_g(b)t_h(a) when a differs from b. The second expression sums the two parental paths to the same heterozygous genotype.

Under independent segregation at the modeled loci, the full conditional child probability is the product of its local transmission probabilities. This is conditional independence given the parental states. It is a separate assumption from independence of allele frequencies or genotype states in a population.

```equation
K[o, pair(i,j)] = P(child=o | parent1=i, parent2=j)
shape(K) = G x U;    sum_o K[o,c] = 1
```

K is a rectangular linear map from a distribution over parental pairs to a distribution over children. It needs neither an inverse nor a diagonalization to perform that map. A rectangular dense baseline is therefore perfectly valid, although its storage may be expensive.

## 3.2 Random mating is generally nonlinear

If two parents are sampled independently from the same genotype distribution x, define w_ii=x_i^2 and w_ij=2x_i x_j for i<j. The next generation is K w(x). Normalization follows from sum w=(sum x)^2=1, but the mapping from x to the next x is generally quadratic.

A one-locus counterexample makes the issue concrete. A population consisting only of AA remains AA; one consisting only of aa remains aa. Mixing those populations equally and then allowing random mating produces AA=1/4, Aa=1/2, aa=1/4. The average of the two unmixed outcomes instead has no heterozygotes. Thus F((x+y)/2) differs from (F(x)+F(y))/2: a fixed linear map cannot generally represent that random-mating operation.

## 3.3 When a square matrix is legitimate

A square operator is appropriate for a lineage whose mates are drawn from a fixed external distribution q. Define M_q[o,i]=sum_j K[o,pair(i,j)]q_j, with each parent role treated consistently. This operator maps a distribution of individual genotypes back to individual genotypes and preserves normalization. If q stays fixed, M_q raised to a power describes the stipulated process.

For example, with mates always AA, an aa lineage produces only Aa children. Its following generation depends on crossing those Aa descendants with AA again. That is a clear biological mating assumption attached to a valid linear operator. It is distinct from a closed random-mating population. Version II implements both fixed-mate operators and complete random-mating updates to demonstrate that the issue is state meaning, not a prohibition on square matrices.

# 4. Deriving the combinatorial growth

## 4.1 Genotypes: a proof by induction

At one biallelic diploid locus, the unphased catalog is AA, Aa, aa. For zero loci there is one empty tuple, so G(0)=1. Assume G(k)=3^k. Every k-locus tuple can be extended in exactly three ways at a new locus. The extensions are distinct, and every (k+1)-locus tuple has one such prefix. Hence G(k+1)=3G(k)=3^(k+1). By induction, G(n)=3^n for every nonnegative integer n.

For a locus with a allowed alleles, there are a homozygotes and a(a-1)/2 heterozygotes, totaling a(a+1)/2. Adding each locus multiplies the number of catalog tuples by its local count. The same induction gives:

```equation
G = product_l [a_l(a_l+1)/2]
U = G + G(G-1)/2 = G(G+1)/2
```

The catalog count does not assert that every tuple is currently present, equally frequent, viable, or appropriate for a particular biological system. It counts the declared unphased Cartesian representation. Population structure changes probabilities; biological constraints may change the admissible subset; linked inheritance may require additional phase information.

{{SCALING_TABLE}}

## 4.2 Supported transmissions: sparse is still exponential

For one biallelic locus, the numbers of possible children for the nine ordered parental pairs form rows (1,2,1), (2,3,2), and (1,2,1). Their total is 15. For independent loci, support products distribute over sums, so the ordered kernel contains 15^n nonzero entries.

To remove parent-order duplication, pair each supported transition with its parent-swapped counterpart. A transition remains fixed under swapping exactly when its complete parent genotypes are equal. The single-locus equal-genotype support counts are 1,3,1, totaling 5; their multilocus total is 5^n. Therefore:

```equation
N_nonzero = (15^n + 5^n)/2
N_dense   = G U = (27^n + 9^n)/2
density   = (15^n + 5^n)/(27^n + 9^n)
```

The decreasing density makes compressed sparse storage attractive, but the nonzero count still grows exponentially. A sparse matrix does not solve arbitrary polygenic enumeration. This is a mathematical limit of explicit supported-transition storage, not a weakness peculiar to a programming language.

The arbitrary-allele version also has a compact derivation. Let H=a(a-1)/2. Across all one-locus parental genotypes, the sum of distinct transmissible allele counts is a+2H=a^2. The product of parental allele-support counts gives a^4 over ordered pairs, except that equal heterozygous parents generate one duplicated unordered heterozygous child. Subtracting those H duplicates gives T(a)=a^4-H. Equal-genotype parent pairs contribute D(a)=a+3H. Thus the independent-locus unordered support total is one-half of product_l T(a_l) plus product_l D(a_l). Exhaustive allele-copy enumeration verifies the formula for the tested catalogs, including ABO (45 nonzeros) and ABO+Rh (615 nonzeros).

## 4.3 Structural coverage, probability coverage, and predictive performance

| Metric | Numerator / denominator or operation | What it does not establish |
|---|---|---|
| Genotype catalog coverage | Represented labels / G | Population prevalence |
| Parent-pair coverage | Represented pairs / U | Equal event probabilities |
| Supported-transition coverage | Retained nonzero entries / complete nonzeros | Clinical prediction |
| Retained probability mass | Sum of retained event probabilities under declared inputs | A world population estimate |
| Computational agreement | Difference from an independent implementation or exact oracle | External biological accuracy |
| Predictive performance | Held-out observed outcomes compared with predicted distributions | Guaranteed portability |

The top-six ABO table includes 6/21=28.57% of pair labels while retaining about 83.48% of probability mass under the exact paper frequencies. The top-18 joint table would include 18/171=10.53% of pair labels while retaining about 67.83% of that model's mass. These metrics are both meaningful, but they answer different questions. For complete kernels and complete original random-mating loops, pair and transmission coverage are 100% within the declared catalog. None of these values establishes predictive accuracy.

# 5. Choosing and implementing data structures

## 5.1 One semantic kernel, several physical representations

The engine uses canonical unordered allele pairs at each locus. A mixed-radix integer encodes the tuple of local genotype indices. Encoding and decoding validate dimensions and index bounds. Parent-exchange symmetry halves redundant complete-pair work without deleting equal-genotype parents. Local Punnett outcomes aggregate repeated allele-copy paths into one child entry and are cached in a bounded least-recently-used table.

The dense baseline allocates all G by U numeric entries. The CSR representation stores only positive transmission values, child-row indexing, and row pointers. Construction first emits columns into compressed-column buffers, then converts to CSR for repeated population updates. Hash adjacency stores a mapping from a parental pair to the mapping of its possible children. The streamed path generates local-factor products as needed and accumulates a complete child distribution without storing the full kernel. All four use the same catalog order and transmission rule. [scipy_csr_array]

Dense storage costs 8GU bytes of binary64 payload. CSR costs value bytes plus index bytes per nonzero and row-pointer bytes; actual integer widths are measured. Hash memory includes substantial Python object overhead and is reported using a different footprint measurement from numeric buffer bytes. Streamed inference avoids persistent kernel storage but repeats enumeration. Construction and inference are therefore separate experimental questions.

## 5.2 Why a trie is a candidate, not a promised winner

A trie can share prefixes and support partial-genotype queries. A radix structure can compress unbranched paths. These mechanisms are useful for dictionary workloads, but a complete n-locus ternary trie still has (3^(n+1)-1)/2 nodes. Dynamic child allocation reduces unused pointer storage; it does not eliminate exponential catalog size. A hash map supports expected constant table probes for fixed-size encoded keys, while building or hashing a long multilocus key still incurs a length-dependent cost. [sedgewickwayne_tries] [sedgewickwayne_hash] [grossi2011compressedtries]

Reduced decision graphs can merge equal suffix computations, but their size depends on the function and variable order. A priority queue can visit likely outcomes first, but ordering an incomplete output does not make it a complete probability distribution. General DAG, adaptive-radix, and priority-queue implementations were therefore not built merely to demonstrate a predetermined winner. Their possible role is documented; a later benchmark must supply the prefix, reuse, or ranked-output workload that justifies their extra machinery. [bryant1986bdd] [python_heapq]

## 5.3 Factored queries change what must be materialized

For specified parents, local conditional distributions can remain separate. A fully specified child's log probability is the sum of its per-locus log probabilities. With bounded allele counts this takes O(n) arithmetic operations and O(n) local-factor storage. It does not enumerate all children. The underlying integer key can grow with n, so arbitrary-precision encoding adds costs not captured by a constant-word model.

This interface is an instance of probability factorization, not a general pedigree inference solver. Sum-product methods use factor structure to marginalize and reuse intermediate computations. Tree-structured factor graphs support exact messages; arbitrary loops do not become exact automatically. General variable elimination can require factors exponential in induced width. The implemented independent-factor query should therefore not be described as solving every genetic inference problem in linear time. [kschischang2001sumproduct] [dechter1999bucket]

If the parental population itself factorizes across loci, the offspring population also factorizes under independent local transmission. Each locus can be updated independently and its result retained as a local factor. However, independent segregation does not imply that an arbitrary input population factorizes. A two-locus population concentrated on AA/BB and aa/bb has a correlation that its marginals alone cannot retain. An explicit regression test checks that silently replacing this population by independent marginals changes the result.

## 5.4 Dynamic allocation, failure, and numerical limits

Before materialization, the implementation computes catalog size, pair count, nonzero count, and conservative representation estimates. A user-specified byte budget can reject a requested materialization. It returns a resource-limit error and retains the possibility of a factored query; it does not silently truncate outcomes. A caught allocation failure clears the bounded local cache and propagates a structured error without returning an incomplete object. Repeated blind allocation attempts are not a recovery strategy because no available resource has changed.

Budget checks are not guarantees of free physical memory or process-wide peak usage. Python objects, allocator fragmentation, native libraries, temporary arrays, and unrelated processes affect real memory consumption. Successful numerical allocation is not proof that a later operation cannot fail. The experiments consequently report storage and measured allocation quantities with their specific definitions.

No probability threshold removes rare transmission branches. Conditional queries provide a log-space interface, and the ordinary probability interface raises if a mathematically positive result underflows. Lazy enumeration has both ordinary and log-space versions. Full population updates use binary64 arithmetic and explicitly reject detected zero underflows; they are not a general arbitrary-precision population engine. Ordinary floating-point rounding still applies. The distinction between a structural zero and an unrepresentably small positive value is part of the API contract.

# 6. Verification before benchmarking

The primary truth oracle enumerates maternal and paternal allele-copy choices independently of the kernel builder. It converts each path into a canonical child and sums exact rational weights. Exhaustive tests cover one-locus biallelic and triallelic systems, four-allele loci, two- and three-locus biallelic systems, and simplified ABO+Rh. They test parent-exchange symmetry, normalization, supported-transition counts, and exact Mendelian probabilities.

Dense and CSR arrays are compared entry for entry. Hash and streamed population results are compared with the same dense result for seeded nonuniform inputs. Separate tests verify allele conservation under neutral random mating, HWE fixed points, encoding round trips, fixed-mate semantics, phase sensitivity, dynamic-programming score probabilities, and resource-budget refusals. Allocation failure is injected deterministically rather than by exhausting the computer.

Legacy regression tests intentionally preserve the original sickle-cell defect to keep reconstruction separate from repair. The original displayed result is expected to sum to 15/16, while the new normalized model is expected to sum to one. A test that silently repairs the legacy output would erase the very evidence being audited.

Adversarial review found three numerical defects during development: silent underflow in a lazy probability product, an absolute tolerance that distorted extreme exact-HWE tails, and underflow when multiplying very small positive fitness weights. The fixes add log-space conditional output, compare HWE event likelihoods in log space, and normalize selection using log weights. An independent exact-integer HWE enumeration checked every genotype triple with total sample size 1 through 30 (5,455 triples), with maximum absolute discrepancy below 7 x 10^-15. Extreme cases were also checked; this is algorithm verification, not evidence that HWE describes the sample.

{{TEST_STATUS}}

# 7. Benchmark methods and results

## 7.1 Matched workloads

The benchmark compares complete dense, CSR, hash-adjacency, and streamed representations for one through five biallelic loci. Each uses the same supported inheritance mapping and the same normalized population input. Construction is separated from inference. A pure numeric matrix-vector product is distinguished from a validated complete random-mating update, which also forms parental-pair probabilities and performs safeguards. Factored single-child queries at larger n are a separate workload with a smaller requested output.

Repeated measurements use a monotonic high-resolution timer. Warm-up, cache policy, repetitions, and the actual software/hardware environment are recorded with the raw timings. One BLAS thread is used. Separate memory passes avoid conflating tracing overhead with uninstrumented timings. Means or medians cannot turn a one-machine experiment into a universal hardware result; the retained samples allow variation to be examined. [python_timeit]

{{BENCH_METHODS}}

## 7.2 Full-kernel results

{{BENCH_TABLE}}

![Measured representation memory and inference costs](../figures/benchmark_overview.png)

{{BENCH_RESULTS}}

These results compare implementations written for this project, including Python iteration and compiled NumPy/SciPy operations. They do not isolate data structure choice from every implementation effect and do not report MATLAB-to-Python speedup. A carefully optimized C or C++ implementation, different processor, different cache behavior, different input distribution, or repeated-query workload could shift the crossover.

## 7.3 Factored queries and explicit limits

{{QUERY_TABLE}}

A complete table and a single requested child's probability are different outputs. The factored method is valuable precisely when the scientific question does not require every joint outcome. A pair of parents heterozygous at every biallelic locus can have 3^n supported child genotypes. Returning all of those states explicitly has exponential output cost even if the distribution has a compact factored description.

At six loci, the dense full kernel alone requires 1,551,807,720 payload bytes. The configured 256 MiB budget rejects that dense materialization. Sparse storage has a smaller theoretical footprint but still grows rapidly. Sizes outside the bounded experiment are labeled either resource-refused or not measured; they are never shown as completed timings.

# 8. Population assumptions tested against external observations

## 8.1 Why Hardy-Weinberg equilibrium is a model, not a universal decoder

At one locus with allele probabilities p and q=1-p, independent gametes produce genotype probabilities p^2, 2pq, q^2. Under the idealized large-population neutral random-mating model, allele frequencies remain constant and genotype frequencies have this form. Selection, migration, nonrandom mating, drift, sampling, and genotyping behavior can make an observed sample differ from the idealized process. The algebra alone does not license interpreting a pooled global phenotype percentage as a universally applicable genotype distribution. [{{HWE_SOURCE}}]

The ABO frequencies in Version I illustrate a separate identifiability issue: converting phenotype frequencies into allele frequencies already introduces assumptions about the genotype-phenotype map and population equilibrium. They are useful as declared synthetic inputs for a reproduction, but are not independent validation data for the model derived from them. Predicting those same frequencies after fitting to them would reuse the answer.

## 8.2 Data provenance and labeling

The empirical exercise uses the public Ensembl rs334 population-genotype endpoint, filtered to the 26 nonoverlapping component populations of the 1000 Genomes phase-3 reference panel. The preserved table has 2,504 individuals, with genomic T/T=2,367, T/A=137, and A/A=0. The data are observed genotype-call counts, not frequencies inferred using HWE. Entries absent from the API are zero-completed only after genotype totals, allele counts, and reported frequencies agree. The snapshot and validation rules are retained with a hash and fetching instructions. [ensembl_rs334_phase3_20260905] [{{1000G_SOURCE}}]

Superpopulation and ALL records are not appended as additional individuals. They overlap component populations and are used only for separate pooled diagnostics. The variant is represented on the GRCh38 genomic forward strand; HBB's transcript notation has the complementary orientation. The dataset is not a worldwide probability sample, a newborn disease survey, or a family transmission cohort. Its population labels are retained as source metadata rather than used to infer an individual's ancestry or clinical outcome.

## 8.3 Test derivation and numerical implementation

For observed genotype counts (nTT,nTA,nAA), estimate q=(2nAA+nTA)/(2N) and calculate N times the HWE probabilities. These fitted expected counts are a compatibility check, not out-of-sample predictions.

The exact test conditions on N and the observed count of one allele. Feasible heterozygote counts have the same parity as the allele count and determine both homozygote counts. The conditional probability of a table is proportional to N! 2^h divided by nTT! h! nAA!. The two-sided probability-ordering p-value sums feasible tables with probability no greater than the observed table. Log-factorials and log-likelihood comparisons avoid an absolute tolerance overwhelming a small tail. The observed table probability alone is not the two-sided p-value. [{{EXACT_SOURCE}}]

The 26 component-population tests receive a Holm adjustment. [holm1979] Pooled ALL and AFR diagnostics are described separately and are exploratory, not selectively promoted as evidence of a causal mechanism. Monomorphic samples provide little information about HWE departures; a p-value of one there does not validate the broader biology.

{{POPULATION_TABLE}}

![External population counts and fitted HWE expectations](../figures/population_hwe.png)

## 8.4 Interpretation and falsification

The YRI subset contains (78,30,0) calls at N=108; its estimated rare-allele frequency is 0.138889. HWE expects approximately (80.083,25.833,2.083), with an exact p-value of 0.213792. The pooled AFR diagnostic contains (529,132,0), expects (535.590,118.820,6.590), and gives p=0.001590. The ALL pool gives p=0.260974. None of the 26 component tests has a Holm-adjusted p-value below 0.05.

Across component populations, the sum of fitted rare-homozygote expectations is 7.600, compared with 1.874 from one globally pooled allele estimate. The difference follows from the variation of q between groups: the weighted mean of q squared exceeds or equals the square of the weighted mean. Pooling therefore changes what the model predicts even before biological causes are considered. In this dataset it can also make a pooled compatibility result look less discrepant than a selected regional pool. These observations do not prove selection, absence of disease, or a universal population law.

The external data establish that a compact real-data audit can reveal consequential modeling choices. They do not establish predictive accuracy for future descendants. A transmission-validation study would require independent parental and offspring genotypes with a clearly defined sampling protocol. Height or disease prediction would additionally require measured phenotypes, fitted parameters, held-out calibration, and population-specific performance analysis.

# 9. Progressive biological extensions

## 9.1 Neutral inheritance, selection, migration, and mutation

The first repaired generational experiment starts with the AS x AS child distribution (1/4,1/2,1/4). Under neutral random mating it remains unchanged: Mendelian segregation does not by itself drive the S allele to extinction. With reproductive weights (1,1,0), selection first normalizes the breeding population, after which random mating produces the next generation. SS births remain possible because heterozygous parents can transmit S.

For an HWE population under those illustrative weights, the S allele update is q'=q/(1+q). Starting at q=1/2, q_t=1/(t+2) when t=0 denotes the initial offspring population. This is a reciprocal decline, not the general exponential-decay claim of Version I. At the third offspring generation, SS birth probability is 1/16 under this model, while the displayed legacy algorithm has already set it to zero.

![Generational dynamics and normalization audit](../figures/generation_models.png)

Selection is implemented by weighting genotype frequencies and normalizing reproductive mass. Migration mixes local and incoming genotype distributions with an explicit mixture fraction. Mutation acts on a gamete allele distribution through a nonnegative row-stochastic transition matrix. The combined synthetic experiment fixes the sequence selection -> migration -> gamete formation -> mutation -> random fertilization. Its rates and fitness weights are declared example values rather than inferred biological constants. An environmental residual is introduced later in phenotype space and is not confused with mutation or selection.

These operators demonstrate how Version I's acknowledged omissions can enter a coherent probabilistic model. They do not show that arbitrary parameter choices describe Canada, Africa, or the world. A future empirical analysis would estimate age-, environment-, and population-dependent parameters and assess sensitivity to their uncertainty.

## 9.2 Another Mendelian system: simplified M/N

The M/N coding example introduces a named two-allele codominant blood-group abstraction. MN x MN produces MM=1/4, MN=1/2, NN=1/4. The local transmission engine needs no change: what differs from a recessive-trait example is the mapping from genotype to phenotype. The full MNS blood-group system is more complex than this two-allele exercise. [{{MNS_SOURCE}}]

This example is deliberately small. Its purpose is to show that a reusable transmission law should not encode the phenotype's dominance relationships into its inheritance probabilities. It passes the same exhaustive allele-copy oracle used for the initial sickle-cell catalog. No external MNS prediction cohort was used.

## 9.3 A multigenic phenotype: simplified ABO-FUT1 epistasis

ABO antigens depend on precursor biology involving the H antigen. A simplified H/h model can therefore demonstrate a genuine interaction: the hh state masks the usual ABO A/B expression. This is inspired by Bombay-phenotype biology, with substantial molecular and serological detail deliberately omitted. It is not a transfusion-compatibility calculator. [{{BOMBAY_SOURCE}}]

For AO/Hh x BO/Hh, the independent inheritance factors give A, B, AB each with probability 3/16. The combined O or O-like category has probability 7/16. Internally the engine distinguishes ordinary ABO OO from hh masking, so genotype information is retained when phenotype labels merge. The calculation has 18 individual catalog entries and 171 unordered parental pairs, but the phenotype map is no longer a simple concatenation of independent locus labels.

This example bridges Mendelian inheritance and multigenic expression without immediately attempting human height. Its result is exact within the stipulated simplified model and is checked against allele-copy enumeration. It has no independent clinical phenotype validation.

## 9.4 Linkage requires phase, not only more columns

Two individuals described as Aa/Bb can carry haplotypes AB/ab or Ab/aB. With recombination fraction r, the first phase transmits AB and ab at (1-r)/2 each and the recombinant haplotypes at r/2 each. The second phase reverses those probabilities. At r=0 the gamete supports are disjoint; at r=1/2 they coincide. These formulas define the implemented two-locus phased model. [{{LINKAGE_SOURCE}}]

The experiment varies r from 0 to 0.5 and compares the phase-specific gamete distributions. Their total variation falls from one to zero. This directly falsifies the proposition that the unphased 3^n catalog alone is sufficient for every inheritance calculation. Linkage disequilibrium concerns population associations, whereas meiotic recombination concerns transmission. Both may matter, and neither should be silently replaced by an independence assumption.

![Staged epistasis and linkage examples](../figures/epistasis_linkage.png)

# 10. Quantitative traits without an astronomical genome table

## 10.1 The synthetic score model

After the simpler systems validate, the extension considers an additive score S=sum_l w_l D_l, where D_l is dosage 0, 1, or 2 and weights are nonnegative integer score units. Independent dosage factors allow the score distribution to be computed by dynamic programming. At each locus, each current score probability contributes to three shifted scores. Equal totals merge into the same bin.

```equation
f_(l+1)(s) = sum_(d=0,1,2) P(D_(l+1)=d) f_l(s-w_(l+1)d)
f_0(0)=1
```

The implementation uses three shifted additions rather than constructing an enormous genotype table. For n equal-weight loci, the final support has at most 2n+1 bins. The cumulative work grows with the sum of intermediate score supports. Widely spaced weights can make the score range large, so a bin budget is enforced. Arbitrary real-valued effects require a separate discretization, continuous approximation, or sampling strategy; the current exact-integer model does not pretend to solve that problem.

## 10.2 Calibration under a known simulator

The test uses 200 loci, each with dosage probabilities (1/4,1/2,1/4), and all weights equal to one. There are 3^200 possible genotype tuples but only 401 score bins. The resulting score mean is 200, variance 100, and probability-array payload 3,208 bytes. The distribution agrees with an independently evaluated Binomial(400,1/2) probability law.

An independent Gaussian residual with standard deviation 10 produces a continuous synthetic outcome Y=S+E. The model's CDF is the weighted mixture of Gaussian CDFs over the score bins. Numerical inversion gives a central 90% interval [176.7379,223.2621]. Independent simulation uses binomial sampling plus Gaussian draws, not resampling the computed score array.

With seed 20260905 and 10,000 draws, 90.03% fall inside the interval. A Wilson 95% interval for this simulated coverage is [89.427%,90.602%]. [nist_wilson] The probability-integral-transform Kolmogorov-Smirnov statistic is 0.00532. [scipy_kstest] These results support consistency and calibration under a fully stipulated simulator; no selection of a favorable seed or repeated search for a passing interval was performed.

![Synthetic polygenic score and calibrated outcome interval](../figures/polygenic_distribution.png)

## 10.3 Why this is not a height predictor

The synthetic score is not measured in centimeters, has no fitted human effect sizes, and contains no real environmental covariates. The residual distribution is known by construction. The demonstrated calibration cannot be transferred to human height by changing a label.

Large-scale height research illustrates the needed distinction. Yengo and colleagues analyzed approximately 5.4 million people and reported 12,111 associated independent SNPs, with materially different out-of-sample variance explained across ancestry groups. That primary result motivates LD-aware models, independent validation, and population-specific calibration; it is not a performance claim about this repository. [yengo2022height]

A defensible height extension would choose a licensed effect-size source and a target population, harmonize variant alleles and genome builds, account for LD and phase where needed, specify environmental and residual assumptions, and evaluate held-out phenotype distributions. Parameter uncertainty and portability would need to enter the reported intervals. Without those inputs and validations, deterministic descendant-height claims would exceed the evidence.

# 11. Adversarial interpretation and limitations

## 11.1 Claims this study can support

The corrected mathematical objects can represent every declared Mendelian outcome without making the transmission kernel square. Exact-rational reconstruction explains the original displayed outputs and identifies divergences between prose, tables, and software. Independent enumeration verifies the small-system implementation and its support formulas. Matched local benchmarks measure how the chosen Python dense, CSR, hash, and streamed implementations behave on the stated machine. Factored queries and score dynamic programming avoid materializing outputs that the scientific question does not require.

The contribution combines audit, derivation, implementation, and evaluation. It does not claim that sparse matrices, hash maps, or factorization are newly invented, and it does not establish priority for the counting identities beyond the derivation supplied here. The useful advance for this project is a coherent representation and a falsifiable comparison attached to reproducible evidence.

## 11.2 Counterarguments retained in the conclusion

First, the original ABO/Rh software already computes all pairings. A coverage-restoration headline would misrepresent the baseline. Second, dense arrays can be faster for small problems because compiled contiguous operations outweigh sparse overhead. Third, sparse storage remains exponential; its success at five loci does not extrapolate to thousands. Fourth, factored-query speed concerns a different output from full enumeration. Fifth, population factorization can be wrong even when per-locus transmission factors correctly.

Sixth, the rs334 panel is an ascertained reference resource. Its calls are neither diagnoses nor an unbiased sample of humanity. Seventh, exact HWE p-values assess a narrow null conditional on sampled allele counts and do not identify a cause of departure. Eighth, a simulated 90% interval with known parameters does not demonstrate clinical calibration. Ninth, byte budgets and caught allocation exceptions are safeguards, not a proof of operating-system-level memory safety under every failure.

## 11.3 Boundaries of the implementation

The core catalog assumes autosomal diploidy and unordered allele pairs. It does not directly handle sex-linked dosage, imprinting, aneuploidy, copy-number variation, somatic mosaicism, penetrance uncertainty, arbitrary pedigrees, or large linked haplotype blocks. The linkage extension covers two phased biallelic loci with a supplied recombination fraction. Mutation is a supplied gamete transition process; its rates are not estimated. The score DP supports nonnegative integer effects and independent dosage factors.

These boundaries are explicit because adding a data structure cannot compensate for a missing biological state variable. A more complete molecular catalog can be inserted only when its inheritance and phenotype rules are specified and independently checked. Unknown parameters should be modeled as uncertain quantities, not replaced with unexplained random noise.

# 12. Reproducibility and evidence inventory

The package retains source code, regression tests, exact legacy outputs, population counts and provenance, benchmark raw samples and metadata, figures in PNG and SVG, an editable manuscript, a rendered PDF, and machine-readable source ledgers. The launcher accepts test, reproduce, science, benchmark, manuscript, and fetch operations. Dependencies are recorded explicitly. The original PDF's digest is checked again in the final verification record.

The data-fetch script supports refreshing to a separate destination. A later API response can differ from the preserved snapshot, so a new fetch should not silently replace the data underlying these results. Component-population membership, sample totals, allele-count consistency, and variant orientation are part of the validation contract. Derived summaries preserve their input hashes.

The benchmark metadata distinguishes predicted storage, retained payload, traced allocations, and timing. The source ledger distinguishes primary papers, authoritative biological references, official implementation documentation, the local user-provided Version I, and new derivations. The bibliography can be rebuilt from that ledger. Chat discussions supplied design motivation but are not treated as scientific authority.

{{REPRO_COMMANDS}}

# 13. Conclusion and next research decisions

The square-matrix restriction is unnecessary for complete inheritance representation, but the original software audit prevents an exaggerated result: ABO and ABO+Rh were already evaluated through complete random-mating loops. Version II's strongest contribution is the separation of transmission, mating, population assumptions, and computational representation, followed by measured comparisons on equivalent tasks.

The results support a layered architecture: a small exact oracle; bounded dense and CSR kernels for complete enumeration; dynamic maps and streaming for flexible access; and local factors or score distributions when the output permits compression. No single data structure wins every workload, and none removes the need to specify biology correctly.

Before integration or submission, the author should review the changed interpretation of Version I, the intended scope of the computational contribution, the source ledger, and the validation limits. The next empirical priority is an independent parent-offspring transmission dataset and a preregistered evaluation target. A real polygenic phenotype study should follow only after appropriate effects, linkage information, outcomes, and calibration cohorts are available. The original paper remains intact until the author chooses how to present this extension alongside it.

<!-- pagebreak -->

# Appendix A. Operational definitions and small examples

The basic callable interface accepts allele counts and mixed-radix genotype codes. For one biallelic locus the order is (0,0), (0,1), (1,1), giving codes 0,1,2. At a triallelic locus the order is (0,0),(0,1),(0,2),(1,1),(1,2),(2,2). This canonical order differs from Version I's presentation order; the reproduction module explicitly retains the legacy order for comparison.

```python
from genetics import InheritanceModel
model = InheritanceModel((2,))
model.cross(1, 1)                  # {0: 0.25, 1: 0.5, 2: 0.25}
model.kernel('csr').shape         # (3, 6)
model.next_generation([.25,.5,.25])

large = InheritanceModel((2,) * 1000)
parent = large.encode((1,) * 1000)
large.log_probability(parent, parent, 0)
# Finite log probability; no 3^1000 child table is constructed.
```

For population inputs, frequencies must be nonnegative and normalized within the specified numeric tolerance. The next-generation method treats parents as independently sampled from that full input distribution. A supplied kernel must match the catalog's G by U shape. Factored population updates require explicit per-locus inputs and the linkage-equilibrium assumption; they do not infer independence from an arbitrary joint vector.

# Appendix B. Reproduction checklist and artifact manifest

{{ARTIFACT_TABLE}}

The checks above concern delivered evidence. They do not assert that MATLAB was executed, that a clinical study was performed, that the software was benchmarked on multiple computers, or that the manuscript has been peer reviewed. A future revision should retain these distinctions when additional evidence becomes available.

# References

{{BIBLIOGRAPHY}}
