# Complete Inheritance Without a Square-Matrix Constraint

## A computational and biological audit, followed by sparse and factored implementations

Version II research candidate | Extension of Arshyia Mehran's Version I | 05 September 2026

Prepared for author review. Original PDF and DOCX preserved. Computations, manuscript preparation, and adversarial checks were assisted by AI agents; scientific authorship, interpretation, and submission require the author's review. This document reports a completed set of local experiments, not a claim that the wider research programme or clinical validation is finished.

<!-- pagebreak -->

# Abstract

Version I explored whether linear algebra could replace repeated Punnett-square calculations for sickle-cell inheritance, ABO blood groups, and a simplified ABO+Rh system. Its discussion identified an apparent loss of population coverage when parental-pair tables were reduced to square matrices. We re-examine that premise before proposing an alternative architecture. An exact-rational Python reconstruction matches six displayed MATLAB examples. The audit finds that the ABO and ABO+Rh programs already evaluate every parental pair through nonlinear random-mating loops, while the prose and displayed tables describe different constructions. Thus, restoring coverage to those executed algorithms is not a valid claim for the extension.

We distinguish genotype catalogs, parental-pair catalogs, supported transmissions, retained probability mass, and predictive performance. For n biallelic diploid loci, an induction proves G=3^n possible unphased catalog entries and U=G(G+1)/2 unordered parental pairs. Under independent Mendelian segregation, the full rectangular kernel has G by U entries and exactly (15^n+5^n)/2 nonzero transitions. Dense arrays, compressed sparse rows, hash-based adjacency maps, and streamed computation are implemented with the same transmission semantics and benchmarked on one local machine. At five loci, CSR stores 12.59 times less numeric payload than dense storage; the corresponding validated update timings are 12.560 ms (dense) and 9.535 ms (CSR). These are workload- and implementation-specific measurements.

External rs334 genotype-call counts from 2,504 reference-panel participants in 26 populations provide an audit of Hardy-Weinberg assumptions, not clinical predictive validation. A pooled model expects 1.874 rare homozygotes, whereas separately fitted population models expect 7.600 in total; none are recorded in the selected snapshot. A staged extension demonstrates simplified M/N inheritance, ABO-FUT1 epistasis, phase-dependent recombination, explicit population operators, and a 200-locus synthetic score distribution. The latter uses 401 bins and achieves 90.03% coverage for a nominal 90% interval in 10,000 independent simulated draws. The contribution is a reproducible separation of biological assumptions, mathematical objects, numerical guarantees, and workload-specific engineering results.

Keywords: Mendelian inheritance; sparse matrices; probability factorization; combinatorial growth; Hardy-Weinberg equilibrium; reproducibility; model audit.

# Reader's guide and research question

The intellectual progression remains that of Version I: first explain a biological question, then introduce the mathematics that represents it, implement the model, and evaluate what the outputs establish. The present question is: **which representation answers a specified inheritance query completely, correctly, and efficiently under a stated biological model?** The qualifying words matter. Completeness without a declared catalog is undefined; efficiency without a declared output is misleading; numerical agreement without external observations does not establish biological predictive accuracy.

The extension proceeds through four linked arguments. First, the original paper must be reconstructed as it exists, including discrepancies between its prose, tables, and software. Second, the parent-to-child mapping must be separated from the rule that constructs the next generation's parents. Third, storage and inference should exploit actual mathematical structure rather than discard low-probability outcomes. Fourth, increasingly realistic biological assumptions must be introduced only when their meaning can be tested.

Sections 1-3 audit Version I. Sections 4-7 derive the complete representation and evaluate the architecture. Sections 8-10 report population and trait experiments. Sections 11-13 discuss falsification, reproducibility, and remaining work. Section 14 is different in kind: it derives, but does not implement, an extension in which age enters the model through germline mutation, somatic accumulation, and epigenetic regulation, and states what would have to be measured before any of it could be believed. The appendices give executable examples, detailed results, and a source bibliography. Every reported benchmark number comes from retained machine-readable output. Theoretical extrapolations are identified separately.

# 1. Version I as evidence, rather than an unquestioned baseline

## 1.1 Preservation and reconstruction

The supplied 54-page PDF was read in full. Embedded code and result images were rendered and inspected because text extraction alone omits important MATLAB content. [1] Its SHA-256 digest before the extension was EE3CB04CE28D4734669D5822FDDBDF03C176FD5C2F5FF3760CB82A4BD1145ECD. Both original files remain at the project root. All new work resides in the version2 directory. The source PDF provides the local source of record for the following audit; page references refer to its printed pages.

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

AA/AS/SS is an educational two-allele model, not a complete catalog of HBB-related disease. Sickle-cell disease includes compound heterozygous forms involving other HBB variants. Trait status does not mean that half of a person's circulating red cells are permanently sickled. Nor is a universal assumption of zero reproduction among affected individuals supported by this computation. Version II therefore uses reproductive weights only as explicitly chosen model parameters and avoids forecasting disease eradication dates. [2]

The rs334 population experiment later in this paper uses genomic T/A variant-call labels. A reference T at this site does not establish that the rest of HBB is normal. This distinction prevents a single-variant catalog from being mistaken for a clinical diagnostic classification. [3]

# 2. The ABO and ABO+Rh audit changes the contribution

## 2.1 Counting parental pairs correctly

The simplified ABO catalog has six genotypes: AA, AO, BB, BO, AB, OO. [4] For six genotype labels, the number of unordered pairs with repetition is 6 x 7 / 2 = 21. The number 15 counts only distinct-genotype pairs and omits the six equal-genotype pairings. Version I states 15 on page 30, then correctly lists 21 pairings on pages 32-33. An equal-genotype pair means two parents who share a genotype; it does not imply self-fertilization.

For simplified ABO+Rh, six ABO states combine with DD, Dd, dd to give 18 catalog entries. The corresponding unordered pair count is 18 x 19 / 2 = 171, and the ordered pair count is 324. Rh positive/negative is treated here as a simplified D-presence model. The full Rh system includes RHD/RHCE variation and variants not captured by a single dominant/recessive pair. The statement that two Rh-positive parents have a one-quarter chance of a negative child is valid for Dd x Dd under this model, not for every positive-parent combination. [5]

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

| Original page | Displayed input and generation | Reconstruction |
| --- | --- | --- |
| 28 | AS x AS, generation 3 | AA 75%, AS 18.75%, SS 0%; total93.75% |
| 35 | AA x OO, generation 3 | AA25%, AO50%, OO25% |
| 36 | AO x BO, generation 4 | AA6.25%, AO25%, BB6.25%, BO25%, AB12.5%, OO25% |
| 49 | AO/Dd x AO/Dd, generation 3 | A+56.25%, A-18.75%, O+18.75%, O-6.25% |
| 50 | AA/dd x BO/Dd, generation 7 | All18 displayed genotype percentages match after rounding |
| 50 | AB/dd x OO/Dd, generation 8 | All18 displayed genotype percentages match after rounding |

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

| Biallelic loci n | Genotype states G | Unordered pairs U | Supported transitions |
| --- | --- | --- | --- |
| 1 | 3 | 6 | 10 |
| 2 | 9 | 45 | 125 |
| 3 | 27 | 378 | 1,750 |
| 4 | 81 | 3,321 | 25,625 |
| 5 | 243 | 29,646 | 381,250 |
| 6 | 729 | 266,085 | 5,703,125 |
| 10 | 59,049 | 1,743,421,725 | 288,330,078,125 |

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

The dense baseline allocates all G by U numeric entries. The CSR representation stores only positive transmission values, child-row indexing, and row pointers. Construction first emits columns into compressed-column buffers, then converts to CSR for repeated population updates. Hash adjacency stores a mapping from a parental pair to the mapping of its possible children. The streamed path generates local-factor products as needed and accumulates a complete child distribution without storing the full kernel. All four use the same catalog order and transmission rule. [6]

Dense storage costs 8GU bytes of binary64 payload. CSR costs value bytes plus index bytes per nonzero and row-pointer bytes; actual integer widths are measured. Hash memory includes substantial Python object overhead and is reported using a different footprint measurement from numeric buffer bytes. Streamed inference avoids persistent kernel storage but repeats enumeration. Construction and inference are therefore separate experimental questions.

## 5.2 Why a trie is a candidate, not a promised winner

A trie can share prefixes and support partial-genotype queries. A radix structure can compress unbranched paths. These mechanisms are useful for dictionary workloads, but a complete n-locus ternary trie still has (3^(n+1)-1)/2 nodes. Dynamic child allocation reduces unused pointer storage; it does not eliminate exponential catalog size. A hash map supports expected constant table probes for fixed-size encoded keys, while building or hashing a long multilocus key still incurs a length-dependent cost. [7] [8] [9]

Reduced decision graphs can merge equal suffix computations, but their size depends on the function and variable order. A priority queue can visit likely outcomes first, but ordering an incomplete output does not make it a complete probability distribution. General DAG, adaptive-radix, and priority-queue implementations were therefore not built merely to demonstrate a predetermined winner. Their possible role is documented; a later benchmark must supply the prefix, reuse, or ranked-output workload that justifies their extra machinery. [10] [11]

## 5.3 Factored queries change what must be materialized

For specified parents, local conditional distributions can remain separate. A fully specified child's log probability is the sum of its per-locus log probabilities. With bounded allele counts this takes O(n) arithmetic operations and O(n) local-factor storage. It does not enumerate all children. The underlying integer key can grow with n, so arbitrary-precision encoding adds costs not captured by a constant-word model.

This interface is an instance of probability factorization, not a general pedigree inference solver. Sum-product methods use factor structure to marginalize and reuse intermediate computations. Tree-structured factor graphs support exact messages; arbitrary loops do not become exact automatically. General variable elimination can require factors exponential in induced width. The implemented independent-factor query should therefore not be described as solving every genetic inference problem in linear time. [12] [13]

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

The final automated verification recorded 32 passing tests, zero failures, and zero errors. The retained JUnit XML and test log identify the executed suite. All 20 representation-by-dimension benchmark comparisons pass an absolute agreement threshold of 10^-12, with the measured maximum difference reported below.

# 7. Benchmark methods and results

## 7.1 Matched workloads

The benchmark compares complete dense, CSR, hash-adjacency, and streamed representations for one through five biallelic loci. Each uses the same supported inheritance mapping and the same normalized population input. Construction is separated from inference. A pure numeric matrix-vector product is distinguished from a validated complete random-mating update, which also forms parental-pair probabilities and performs safeguards. Factored single-child queries at larger n are a separate workload with a smaller requested output.

Repeated measurements use a monotonic high-resolution timer. Warm-up, cache policy, repetitions, and the actual software/hardware environment are recorded with the raw timings. One BLAS thread is used. Separate memory passes avoid conflating tracing overhead with uninstrumented timings. Means or medians cannot turn a one-machine experiment into a universal hardware result; the retained samples allow variation to be examined. [14]

The retained run uses three construction repetitions and seven inference repetitions per materialized representation and dimension. Each construction sample and each inference or memory configuration runs in a fresh worker. Imports and model setup are outside timed construction. One untimed construction warms allocator/library behavior, followed by clearing the local-cross cache before each timed construction. Inference uses an untimed warm-up and warm local tables. The traced-memory pass begins after scientific-library imports and covers construction plus one update. The complete recorded environment is in benchmark.json; the summary below reproduces its key fields.

| Environment item | Recorded value |
| --- | --- |
| python | 3.12.14 (main, Aug 25 2026, 14:01:42) [MSC v.1944 64 bit (AMD64)] |
| platform | Windows-11-10.0.26200-SP0 |
| processor | Intel64 Family 6 Model 151 Stepping 2, GenuineIntel |
| logical_cpu_count | 24 |
| numpy | 2.5.2 |
| scipy | 1.18.1 |

## 7.2 Full-kernel results

| Loci | Method | Retained bytes* | Build median (ms) | Update median (ms) |
| --- | --- | --- | --- | --- |
| 1 | dense | 144 | 0.092 | 0.032 |
| 1 | csr | 136 | 0.179 | 0.037 |
| 1 | hash | 2,356 | 0.077 | 0.022 |
| 1 | streamed kernel | 0 | No full kernel | 0.039 |
| 2 | dense | 3,240 | 0.323 | 0.045 |
| 2 | csr | 1,540 | 0.507 | 0.050 |
| 2 | hash | 18,756 | 0.333 | 0.067 |
| 2 | streamed kernel | 0 | No full kernel | 0.313 |
| 3 | dense | 81,648 | 3.501 | 0.155 |
| 3 | csr | 21,112 | 3.671 | 0.176 |
| 3 | hash | 191,684 | 3.551 | 0.651 |
| 3 | streamed kernel | 0 | No full kernel | 4.012 |
| 4 | dense | 2,152,008 | 49.604 | 1.129 |
| 4 | csr | 307,828 | 51.231 | 1.026 |
| 4 | hash | 2,260,956 | 50.376 | 8.837 |
| 4 | streamed kernel | 0 | No full kernel | 58.544 |
| 5 | dense | 57,631,824 | 744.577 | 12.560 |
| 5 | csr | 4,575,976 | 772.977 | 9.535 |
| 5 | hash | 29,427,332 | 739.609 | 119.853 |
| 5 | streamed kernel | 0 | No full kernel | 831.938 |

*Dense/CSR: retained numeric buffers only. Hash: recursive Python-object footprint, a different metric. Streamed: no retained full-kernel buffer; this is not zero process memory. The figure uses the separate traced-allocation pass for all methods.

![Measured representation memory and inference costs](../figures/benchmark_overview.png)

At five loci, the full catalog contains G=243 genotypes and U=29,646 unordered pairs. The kernel has 381,250 nonzeros. Dense numeric payload is 57,631,824 bytes; CSR payload is 4,575,976 bytes, a 12.59-fold reduction. The ratio of dense to CSR median validated-update time is 1.32. Construction medians are 0.7446 s and 0.7730 s respectively. The maximum absolute discrepancy from the dense result over every measured method and dimension is 6.94e-17. The complete timing samples, minima, maxima, standard deviations, and separate memory records are retained rather than hidden behind these medians.

These results compare implementations written for this project, including Python iteration and compiled NumPy/SciPy operations. They do not isolate data structure choice from every implementation effect and do not report MATLAB-to-Python speedup. A carefully optimized C or C++ implementation, different processor, different cache behavior, different input distribution, or repeated-query workload could shift the crossover.

## 7.3 Factored queries and explicit limits

| Loci | Median query (microseconds) | Traced peak (bytes) | Absolute log error |
| --- | --- | --- | --- |
| 10 | 11.803 | 12,425 | 0 |
| 50 | 49.378 | 24,609 | 7.11e-14 |
| 100 | 105.590 | 38,801 | 1.14e-13 |
| 1000 | 1246.460 | 266,069 | 2.61e-11 |

A complete table and a single requested child's probability are different outputs. The factored method is valuable precisely when the scientific question does not require every joint outcome. A pair of parents heterozygous at every biallelic locus can have 3^n supported child genotypes. Returning all of those states explicitly has exponential output cost even if the distribution has a compact factored description.

At six loci, the dense full kernel alone requires 1,551,807,720 payload bytes. The configured 256 MiB budget rejects that dense materialization. Sparse storage has a smaller theoretical footprint but still grows rapidly. Sizes outside the bounded experiment are labeled either resource-refused or not measured; they are never shown as completed timings.

# 8. Population assumptions tested against external observations

## 8.1 Why Hardy-Weinberg equilibrium is a model, not a universal decoder

At one locus with allele probabilities p and q=1-p, independent gametes produce genotype probabilities p^2, 2pq, q^2. Under the idealized large-population neutral random-mating model, allele frequencies remain constant and genotype frequencies have this form. Selection, migration, nonrandom mating, drift, sampling, and genotyping behavior can make an observed sample differ from the idealized process. The algebra alone does not license interpreting a pooled global phenotype percentage as a universally applicable genotype distribution. [15]

The ABO frequencies in Version I illustrate a separate identifiability issue: converting phenotype frequencies into allele frequencies already introduces assumptions about the genotype-phenotype map and population equilibrium. They are useful as declared synthetic inputs for a reproduction, but are not independent validation data for the model derived from them. Predicting those same frequencies after fitting to them would reuse the answer.

## 8.2 Data provenance and labeling

The empirical exercise uses the public Ensembl rs334 population-genotype endpoint, filtered to the 26 nonoverlapping component populations of the 1000 Genomes phase-3 reference panel. The preserved table has 2,504 individuals, with genomic T/T=2,367, T/A=137, and A/A=0. The data are observed genotype-call counts, not frequencies inferred using HWE. Entries absent from the API are zero-completed only after genotype totals, allele counts, and reported frequencies agree. The snapshot and validation rules are retained with a hash and fetching instructions. [16] [17]

Superpopulation and ALL records are not appended as additional individuals. They overlap component populations and are used only for separate pooled diagnostics. The variant is represented on the GRCh38 genomic forward strand; HBB's transcript notation has the complementary orientation. The dataset is not a worldwide probability sample, a newborn disease survey, or a family transmission cohort. Its population labels are retained as source metadata rather than used to infer an individual's ancestry or clinical outcome.

## 8.3 Test derivation and numerical implementation

For observed genotype counts (nTT,nTA,nAA), estimate q=(2nAA+nTA)/(2N) and calculate N times the HWE probabilities. These fitted expected counts are a compatibility check, not out-of-sample predictions.

The exact test conditions on N and the observed count of one allele. Feasible heterozygote counts have the same parity as the allele count and determine both homozygote counts. The conditional probability of a table is proportional to N! 2^h divided by nTT! h! nAA!. The two-sided probability-ordering p-value sums feasible tables with probability no greater than the observed table. Log-factorials and log-likelihood comparisons avoid an absolute tolerance overwhelming a small tail. The observed table probability alone is not the two-sided p-value. [18]

The 26 component-population tests receive a Holm adjustment. [19] Pooled ALL and AFR diagnostics are described separately and are exploratory, not selectively promoted as evidence of a causal mechanism. Monomorphic samples provide little information about HWE departures; a p-value of one there does not validate the broader biology.

| Population / diagnostic | N | Observed T/T, T/A, A/A | Exact HWE p |
| --- | --- | --- | --- |
| ACB | 96 | 87, 9, 0 | 1.000000 |
| ASW | 61 | 59, 2, 0 | 1.000000 |
| CLM | 94 | 92, 2, 0 | 1.000000 |
| ESN | 99 | 75, 24, 0 | 0.349703 |
| GWD | 113 | 87, 26, 0 | 0.354696 |
| LWK | 99 | 79, 20, 0 | 0.591267 |
| MSL | 85 | 64, 21, 0 | 0.349575 |
| PUR | 104 | 101, 3, 0 | 1.000000 |
| YRI | 108 | 78, 30, 0 | 0.213792 |
| AFR (pooled) | 661 | 529, 132, 0 | 0.001590 |
| ALL (pooled) | 2504 | 2367, 137, 0 | 0.260974 |

The 17 monomorphic component populations remain in the CSV and the 26-test adjustment. They are omitted here only to keep the displayed diagnostic table compact.

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

The M/N coding example introduces a named two-allele codominant blood-group abstraction. MN x MN produces MM=1/4, MN=1/2, NN=1/4. The local transmission engine needs no change: what differs from a recessive-trait example is the mapping from genotype to phenotype. The full MNS blood-group system is more complex than this two-allele exercise. [20]

This example is deliberately small. Its purpose is to show that a reusable transmission law should not encode the phenotype's dominance relationships into its inheritance probabilities. It passes the same exhaustive allele-copy oracle used for the initial sickle-cell catalog. No external MNS prediction cohort was used.

## 9.3 A multigenic phenotype: simplified ABO-FUT1 epistasis

ABO antigens depend on precursor biology involving the H antigen. A simplified H/h model can therefore demonstrate a genuine interaction: the hh state masks the usual ABO A/B expression. This is inspired by Bombay-phenotype biology, with substantial molecular and serological detail deliberately omitted. It is not a transfusion-compatibility calculator. [21]

For AO/Hh x BO/Hh, the independent inheritance factors give A, B, AB each with probability 3/16. The combined O or O-like category has probability 7/16. Internally the engine distinguishes ordinary ABO OO from hh masking, so genotype information is retained when phenotype labels merge. The calculation has 18 individual catalog entries and 171 unordered parental pairs, but the phenotype map is no longer a simple concatenation of independent locus labels.

This example bridges Mendelian inheritance and multigenic expression without immediately attempting human height. Its result is exact within the stipulated simplified model and is checked against allele-copy enumeration. It has no independent clinical phenotype validation.

## 9.4 Linkage requires phase, not only more columns

Two individuals described as Aa/Bb can carry haplotypes AB/ab or Ab/aB. With recombination fraction r, the first phase transmits AB and ab at (1-r)/2 each and the recombinant haplotypes at r/2 each. The second phase reverses those probabilities. At r=0 the gamete supports are disjoint; at r=1/2 they coincide. These formulas define the implemented two-locus phased model. [22]

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

With seed 20260905 and 10,000 draws, 90.03% fall inside the interval. A Wilson 95% interval for this simulated coverage is [89.427%,90.602%]. [23] The probability-integral-transform Kolmogorov-Smirnov statistic is 0.00532. [24] These results support consistency and calibration under a fully stipulated simulator; no selection of a favorable seed or repeated search for a passing interval was performed.

![Synthetic polygenic score and calibrated outcome interval](../figures/polygenic_distribution.png)

## 10.3 Why this is not a height predictor

The synthetic score is not measured in centimeters, has no fitted human effect sizes, and contains no real environmental covariates. The residual distribution is known by construction. The demonstrated calibration cannot be transferred to human height by changing a label.

Large-scale height research illustrates the needed distinction. Yengo and colleagues analyzed approximately 5.4 million people and reported 12,111 associated independent SNPs, with materially different out-of-sample variance explained across ancestry groups. That primary result motivates LD-aware models, independent validation, and population-specific calibration; it is not a performance claim about this repository. [25]

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

```text
cd version2
python -m pip install -r requirements.txt
python run.py test -q
python run.py reproduce
python run.py fetch
python run.py science
python run.py benchmark
python run.py test -q --junitxml=results/tests.xml
python run.py manuscript
```

For the bundled Windows Python used here, reproduce.ps1 sets the verified interpreter path. The locally installed .deps directory is optional and is excluded from the source archive.

# 13. Conclusion and next research decisions

The square-matrix restriction is unnecessary for complete inheritance representation, but the original software audit prevents an exaggerated result: ABO and ABO+Rh were already evaluated through complete random-mating loops. Version II's strongest contribution is the separation of transmission, mating, population assumptions, and computational representation, followed by measured comparisons on equivalent tasks.

The results support a layered architecture: a small exact oracle; bounded dense and CSR kernels for complete enumeration; dynamic maps and streaming for flexible access; and local factors or score distributions when the output permits compression. No single data structure wins every workload, and none removes the need to specify biology correctly.

Before integration or submission, the author should review the changed interpretation of Version I, the intended scope of the computational contribution, the source ledger, and the validation limits. The next empirical priority is an independent parent-offspring transmission dataset and a preregistered evaluation target. A real polygenic phenotype study should follow only after appropriate effects, linkage information, outcomes, and calibration cohorts are available. The original paper remains intact until the author chooses how to present this extension alongside it. Section 14 sets out, in derivational form only, the one limitation both versions have carried unchanged: neither model has a clock.

<!-- pagebreak -->

# 14. Age, somatic mutation and epigenetic regulation: a derivational proposal

**Status of this section.** Nothing in it is implemented, tested or validated. Sections 1 to 13 report executed work; this section derives a candidate extension and states what would have to be measured before any of it could be believed. It is placed after the conclusion for that reason. Every equation below is a proposal, and the distinction between a derived consequence of stated assumptions and an empirical claim about human biology is maintained throughout. Where a quantity is anchored to published measurement, the measurement is cited and its scope is recorded; no parameter in this section has been estimated from data by this project.

## 14.1 The problem: the kernel has no clock

Every object constructed in Sections 3 to 10 is timeless. The transmission kernel K[o, pair(i,j)] depends on two parental genotypes and nothing else. The population operators advance a generation index, not an age. Two parents produce the same offspring distribution whether they are twenty or fifty, and an individual's genotype is treated as fixed from conception to death.

Both simplifications are wrong in the same direction, and Version I said so. Its conclusion records that "random mutations can occur in an individual before they breed" and identifies this as a scenario the model cannot represent. [1] That is an accurate self-assessment, and it survives into Version II unchanged: adding sparse storage, factored queries and score dynamic programming did nothing about it. The representation improved; the biology did not.

Age enters this system by two mechanisms that are routinely conflated and must not be.

## 14.2 Two routes, one of which touches the kernel

> **OBSERVATION 14.1 — THE GERMLINE ROUTE AND THE SOMATIC ROUTE ARE DIFFERENT OBJECTS**
>
> **Germline.** Parental age at conception changes the mutation content of the gametes actually transmitted. This alters what the next generation inherits. It modifies the transmission kernel itself, and its effects are heritable.
>
> **Somatic and epigenetic.** An individual's own age changes which cells carry which mutations, and changes the methylation and histone-modification state of their chromatin. This alters whether an inherited genotype is expressed. It modifies the genotype-to-phenotype map, and its effects are, with narrow exceptions, not transmitted.

The distinction is not pedantic. A model that adds a single "age" term and lets it act on both routes at once will produce an operator that is neither a valid transmission kernel nor a valid penetrance function, and no amount of fitting will reveal the error, because both routes make outcomes depend on age in the same direction. Sections 14.4, 14.5 and 14.6 therefore treat them separately, and only 14.4 is permitted to alter K.

## 14.3 Notation for this section

| Symbol | Name | Definition |
|---|---|---|
| t | Age | Age of an individual in years, measured from conception or birth as stated. Continuous. |
| a_f, a_m | Parental ages | Paternal and maternal age in years at the conception of the offspring. |
| L | Locus index set | The modelled loci, indexed by l, as in Section 4. |
| Lambda | Expected germline count | Expected number of de novo mutations genome-wide in one transmission. |
| mu_l | Expected locus count | Expected number of de novo mutations falling inside locus l in one transmission. |
| eps_l | Per-transmission mutation probability | Probability that locus l carries at least one de novo mutation in one transmission. |
| Q | Substitution matrix | Row-stochastic matrix over the allele alphabet of a locus, given that a mutation has occurred. |
| M(a) | Age-dependent gamete operator | Row-stochastic allele transition applied to gametes, parameterised by parental age. |
| lambda_l(t) | Somatic intensity | Instantaneous rate of somatic mutation at locus l at age t. |
| m_l(t) | Methylation fraction | Fraction of the relevant CpG sites of locus l that are methylated at age t, in [0,1]. |
| h_l(t) | Acetylation level | Standardised activating histone-acetylation signal at locus l at age t. |
| alpha, beta | Epigenetic rates | Methylation and demethylation rates in the two-state model of 14.6. |
| z_l(t) | Expression score | Linear predictor of the expression gate. Log-odds scale. |
| w | Weight vector | Coefficients of the expression score. The quantities to be estimated. |
| b_l | Bias | Locus-specific intercept of the expression score. |
| pi_l(t) | Expression probability | Probability that locus l is transcriptionally competent at age t. |
| sigma | Logistic function | sigma(z) = 1/(1+exp(-z)). |
| y_i | Observed label | 1 if locus observed expressed in sample i, 0 if silenced. |
| J(w) | Objective | Mean negative log-likelihood of the observed labels. |
| eta | Step size | Scalar multiplying the gradient in an update. |

## 14.4 The germline route: parental age inside the transmission kernel

### 14.4.1 What is measured

Two large trio studies give the empirical anchor. Sequencing 78 Icelandic parent-offspring trios gave an average de novo mutation rate of 1.20e-8 per nucleotide per generation at a mean paternal age of 29.7 years, with the count rising by roughly two mutations per year of paternal age. [26] The larger follow-up, 1,548 trios and 108,778 high-quality de novo mutations, gave a mean of 70.3 mutations per trio and separated the parental contributions: 1.51 additional mutations per year of paternal age against 0.37 per year of maternal age. [27]

These are population-level regression slopes from one country's cohort. They are the best available anchor and they are not a per-locus mutation probability for an arbitrary gene in an arbitrary population.

### 14.4.2 From a genome-wide count to a per-locus probability

Model the expected genome-wide de novo count as affine in both parental ages, with the slopes above and a reference point where the intercept is calibrated:

```equation
Lambda(a_f, a_m) = Lambda_0 + beta_f (a_f - a_f0) + beta_m (a_m - a_m0)
beta_f = 1.51 per year,  beta_m = 0.37 per year   [anchored, not fitted here]
```

De novo mutations are, to a first approximation and away from clustered events, spread across the genome. Writing G_bp for the callable genome length and L_l for the length of locus l, the expected count falling inside that locus is

```equation
mu_l(a_f, a_m) = Lambda(a_f, a_m) * (L_l / G_bp) * kappa_l
```

where kappa_l is a locus-specific enrichment factor absorbing the fact that mutation rate is not uniform: CpG sites, replication timing and sequence context all matter. Setting kappa_l = 1 asserts uniformity, which is known to be false and is retained only as an explicit null. Treating the occurrence of de novo mutations in a short interval as Poisson,

```equation
eps_l(a_f, a_m) = 1 - exp( -mu_l(a_f, a_m) )
```

For a gene-sized locus this is of order 1e-8 to 1e-6 per transmission. It is small, and Section 14.11 returns to what that smallness does to binary64 arithmetic.

### 14.4.3 The age-dependent gamete operator

The implementation already contains the required interface. `genetics/core.py` exposes `mutation(allele_frequencies, transition)`, which applies a row-stochastic allele transition to gametes before fertilization and validates that the matrix is finite, nonnegative and row-stochastic. The extension is not a new mechanism; it is making that matrix a function of parental age.

> **RESULT 14.2 — AGE-PARAMETERISED GAMETE MUTATION OPERATOR**
>
> For a locus with allele alphabet of size k, define
>
> ```equation
> M_ab(a_f, a_m) = (1 - eps_l) * delta_ab + eps_l * Q_ab
> ```
>
> where Q is row-stochastic over the alphabet and delta is the Kronecker delta. The transmitted gamete allele distribution becomes t~_g(a) = sum_b t_g(b) M_ba, with t_g the ordinary Mendelian gamete law of Section 3.1.

The result is a one-line change in form and a substantial change in meaning: the parent-to-child map now carries two continuous parameters that were previously absent from the entire architecture.

> **RESULT 14.3 — NORMALISATION IS PRESERVED**
>
> Each row of M sums to (1 - eps_l) + eps_l * sum_b Q_ab = (1 - eps_l) + eps_l = 1. Since t_g is a probability distribution over alleles and M is row-stochastic, sum_a t~_g(a) = sum_b t_g(b) sum_a M_ba = sum_b t_g(b) = 1.
>
> Therefore the modified kernel remains a conditional distribution over children, and every normalisation test in `tests/` would still be a valid check of it.

This matters because it means the extension does not invalidate the verification apparatus. The allele-copy oracle of Section 6 compares against exact Mendelian probabilities; with mutation switched on it would compare against exact Mendelian probabilities composed with a known stochastic matrix, which is equally checkable by enumeration.

### 14.4.4 The consequence nobody wants: sparsity dies

> **RESULT 14.4 — A POSITIVE MUTATION RATE DESTROYS STRUCTURAL SPARSITY**
>
> Section 4.2 establishes that the unordered kernel has exactly (15^n + 5^n)/2 nonzero entries against G*U = (27^n + 9^n)/2 total entries, and that this falling density is the entire justification for compressed sparse storage.
>
> If eps_l > 0 and Q has full support, then every allele is reachable from every allele in one transmission. Every structural zero of the kernel becomes a positive number of order eps_l per mutated locus. The nonzero count rises from (15^n + 5^n)/2 to G*U exactly.
>
> At five biallelic loci that is a rise from 381,250 to 7,203,978 stored entries, and CSR becomes strictly worse than dense: it pays index overhead on a matrix with no zeros.

This is the most consequential thing in the section and it was not anticipated when the architecture was chosen. The measured 12.59-fold CSR payload advantage of Section 7.2 is an advantage over a mutation-free model. Introduce biologically realistic mutation and the advantage does not shrink; it inverts.

Three responses are available, and choosing between them is an experiment, not a preference. Thresholding restores sparsity but violates the scientific contract of Section 12, which forbids pruning positive transmission branches. Factored representation avoids the problem entirely, because the per-locus operators stay small and are never multiplied out; this is the same argument as Section 5.3 and it becomes considerably stronger here. Structured storage keeps the Mendelian part sparse and represents the mutation part as a low-rank or Kronecker correction, exploiting the fact that M is the identity plus a rank-structured perturbation of size eps_l. The third is the interesting one and it is not implemented.

## 14.5 The somatic route: accumulation within one lifetime

Somatic mutation does not alter what an individual transmits. It alters the individual. The relevant object is therefore not the kernel but a per-locus state that evolves with age.

Adult stem cells of the liver, colon and small intestine accumulate mutations across life at approximately 40 novel mutations per year, at a broadly similar rate across those tissues despite very different cancer incidence. [28] Approximate linearity in age is the empirical starting point.

Model somatic mutation at locus l as an inhomogeneous Poisson process with intensity lambda_l(t). The expected accumulated count and the probability that the locus is still unmutated at age t are

```equation
m_l(t) = integral_0^t lambda_l(s) ds
P(locus l unmutated at age t) = exp( -m_l(t) )
```

> **RESULT 14.5 — CONSTANT INTENSITY GIVES EXPONENTIAL SURVIVAL, NOT A LINEAR ONE**
>
> If lambda_l(s) = lambda_l is constant, then m_l(t) = lambda_l t and the probability that locus l has acquired at least one somatic mutation by age t is 1 - exp(-lambda_l t).
>
> For lambda_l t << 1 this is approximately lambda_l t, which is why counts look linear in age in the data. The linear appearance is the small-argument regime of an exponential, not evidence that the underlying probability is linear, and the two diverge at large lambda_l t.

A single mutation is often insufficient. The classical multistage argument supposes that a cell must complete k independent rare steps, and derives an incidence rising as a power of age; fitting non-endocrine carcinoma incidence gave approximately a sixth-power dependence, with the log-log slope interpreted as the number of required stages minus one. [29]

> **RESULT 14.6 — MULTISTAGE HAZARD AND ITS AGE POWER**
>
> Let each of k stages occur independently at small constant rate. The probability that all k are complete by age t is of order (c t)^k / k!, so the incidence hazard is
>
> ```equation
> h(t) = d/dt [ (c t)^k / k! ] = c^k t^(k-1) / (k-1)!
> ```
>
> A log-log plot of hazard against age therefore has slope k-1, which is what makes k estimable from incidence curves alone.

This is a hazard for a multi-hit somatic process. It is not a model of Mendelian inheritance and must not be substituted into the transmission kernel. Its place in this proposal is as the correct functional form for age-of-onset, once a genotype has been inherited.

The general form covering both, and the standard vocabulary for it, is a hazard factorised into a baseline function of time and a covariate term, which is exactly the structure of proportional-hazards regression. [30] Writing x for covariates including genotype dosage:

```equation
h_l(t | x) = h_0(t) * exp( gamma^T x )
```

## 14.6 The epigenetic route: a reversible two-state process

Methylation is categorically different from mutation, and the difference is mathematical before it is biological.

DNA methylation acts in a context-dependent way at promoters, gene bodies and regulatory elements, with promoter-island methylation associated with transcriptional repression, and the relationship between methylation and transcription is explicitly not a simple switch. [31] Histone modifications, acetylation among them, regulate chromatin as a responsive scaffold with mark-specific transcriptional consequences. [32]

Model a single CpG site as a two-state continuous-time Markov chain with methylation rate alpha and demethylation rate beta. Writing P_M(t) for the probability the site is methylated:

```equation
dP_M/dt = alpha (1 - P_M) - beta P_M
```

> **RESULT 14.7 — EXPONENTIAL APPROACH TO AN EQUILIBRIUM METHYLATION LEVEL**
>
> Rearranging, dP_M/dt = alpha - (alpha + beta) P_M. This is linear with constant coefficients, so with P_M(0) given,
>
> ```equation
> P_M(t) = alpha/(alpha+beta) + [ P_M(0) - alpha/(alpha+beta) ] * exp( -(alpha+beta) t )
> ```
>
> The methylation fraction approaches the equilibrium alpha/(alpha+beta) with time constant 1/(alpha+beta), monotonically, from whichever side it starts.
>
> Verification: at t=0 the bracket cancels the equilibrium term and returns P_M(0); as t grows the exponential vanishes and leaves the equilibrium; substituting back reproduces the differential equation.

Averaging over the sites of a locus gives the observable m_l(t) of the notation table. Methylation state carries enough reproducible age information to build multi-tissue age predictors: one such predictor was constructed from approximately 8,000 samples spanning 51 healthy tissues and cell types [33], and a quantitative ageing model over more than 450,000 CpG markers in whole blood from 656 individuals aged 19 to 101 measured individual differences in methylome ageing rate and characterised epigenetic drift [34]. That such clocks work is evidence that m_l(t) is a real, measurable, age-structured signal. It is not evidence that methylation causes ageing, and neither study licenses that reading.

> **OBSERVATION 14.8 — MUTATION IS ABSORBING; METHYLATION IS ERGODIC**
>
> A mutation, once fixed in a lineage of cells, does not spontaneously revert at a comparable rate: the process has an effectively absorbing state and P(unmutated) decreases monotonically towards zero.
>
> Methylation has strictly positive rates in both directions: the chain is ergodic, has a stationary distribution, and P_M(t) approaches a level strictly inside (0,1) rather than an absorbing endpoint.
>
> Consequence: an intervention can in principle reverse an epigenetic silencing and cannot reverse a mutation. Any model that represents both with the same operator has thrown away that distinction, and with it the only clinically actionable difference between them.

## 14.7 The expression gate: weights, bias and age-dependent penetrance

The three routes now combine into a single scalar per locus. This is the construction the section exists to propose.

Collect everything known about locus l in individual i at age t into a feature vector, and score it linearly:

> **RESULT 14.9 — THE EXPRESSION GATE**
>
> ```equation
> z_l(t) = b_l + w_age * t + w_met * m_l(t) + w_ac * h_l(t) + w_dos * g_l + sum_c w_c x_c
> pi_l(t) = sigma( z_l(t) ) = 1 / (1 + exp(-z_l(t)))
> ```
>
> pi_l(t) is the probability that locus l is transcriptionally competent in that individual at that age. The bias b_l is the log-odds of expression for a reference individual: reference age, reference epigenetic state, reference genotype dosage. It is the locus's baseline propensity to be expressed, and it is the parameter that carries everything the covariates do not explain.

The linear-then-link structure is a generalised linear model with a binomial response and a logit link, and the estimation machinery for that class is standard. [35] The score is nonlinear in age through m_l(t) and h_l(t), which are themselves nonlinear in t by Result 14.7, but it is linear in the parameters w. That single structural property determines everything in Section 14.8, exactly as it does in any linear-in-parameters model.

The signs of two weights are constrained in advance by biology, which makes them a falsification test rather than free parameters:

| Weight | Predicted sign | Basis | If fitting returns the other sign |
|---|---|---|---|
| w_met | negative | Promoter-island methylation associated with repression [31] | Either the model is wrong, the sites are gene-body rather than promoter, or the annotation is misassigned |
| w_ac | positive | Activating acetylation marks open chromatin [32] | Same three candidate explanations, in the same order |
| w_age | unconstrained | Residual age effect after m and h are accounted for | A large value indicates the epigenetic covariates are not capturing the age dependence |

The last row is the diagnostic one. If m_l(t) and h_l(t) genuinely mediate the effect of age, then conditioning on them should leave w_age near zero. A large surviving w_age says the mediation story is incomplete, and that is worth knowing.

> **RESULT 14.10 — AGE-DEPENDENT EFFECTIVE PENETRANCE**
>
> Let the phenotype require a functional product from at least one allele at locus l, and let each inherited allele be independently competent with probability pi_l(t). For an individual carrying c_l functional alleles by inheritance,
>
> ```equation
> P(phenotype expressed | c_l functional alleles, age t) = 1 - (1 - pi_l(t))^(c_l)
> ```
>
> Penetrance is therefore a function of age and epigenetic state, not a constant of the genotype. A heterozygote whose single functional allele is silenced behaves, at that locus and that age, as though it carried none.

That last sentence is the biologically interesting consequence and the one most in need of restraint. It is a consequence of the stated model. It is not a claim that any particular carrier of any particular variant will develop any particular condition, and this section provides no evidence for such a claim.

## 14.8 Estimation: objective, gradient, convexity, step size

Suppose paired observations are available: for sample i, an age t_i, measured methylation and acetylation, a genotype, and a binary expression label y_i. Write x_i for the assembled feature vector and pi_i = sigma(w^T x_i).

The objective is the mean negative log-likelihood of the labels under the Bernoulli model, that is, cross-entropy:

```equation
J(w) = -(1/M) sum_i [ y_i log pi_i + (1 - y_i) log(1 - pi_i) ]
```

Squared error is not used. Under a binary response it is non-convex in w after composition with the logistic link and it penalises confident correct predictions unhelpfully; cross-entropy is the log-likelihood of the assumed response distribution and is the canonical choice for this link. [35]

> **RESULT 14.11 — GRADIENT OF THE CROSS-ENTROPY OBJECTIVE**
>
> The logistic function satisfies sigma'(z) = sigma(z)(1 - sigma(z)). Differentiating the contribution of one sample with respect to w:
>
> ```equation
> (1)  d/dw [ -y log pi - (1-y) log(1-pi) ]  =  [ -(y/pi) + (1-y)/(1-pi) ] * dpi/dw
> (2)  dpi/dw = pi (1 - pi) x
> (3)  substituting:  [ -(y/pi) + (1-y)/(1-pi) ] * pi(1-pi) x
> (4)  = [ -y(1-pi) + (1-y) pi ] x
> (5)  = ( pi - y ) x
> ```
>
> Summing and averaging,
>
> ```equation
> grad J(w) = (1/M) sum_i ( pi_i - y_i ) x_i = (1/M) X^T ( pi - y )
> ```

The gradient is the feature matrix applied to the residual vector: each weight is adjusted in proportion to how strongly its own feature correlates with the remaining error. The factors pi(1-pi) cancel exactly at step (4), which is the property that makes the logit link canonical for this response.

> **RESULT 14.12 — THE OBJECTIVE IS CONVEX**
>
> Differentiating Result 14.11 once more, with S = diag( pi_i (1 - pi_i) ):
>
> ```equation
> H = grad^2 J(w) = (1/M) X^T S X
> ```
>
> For any vector v, v^T H v = (1/M) (X v)^T S (X v) = (1/M) sum_i pi_i(1-pi_i) (x_i^T v)^2 >= 0, since each pi_i(1-pi_i) > 0. So H is positive semi-definite and J is convex: every stationary point is a global minimum.
>
> The minimum is unique when X has full column rank. Unlike the linear case, H depends on w, so the curvature is not constant and no closed-form minimiser exists; estimation is necessarily iterative.

> **RESULT 14.13 — A STEP-SIZE BOUND FOLLOWS FROM THE CURVATURE BOUND**
>
> The scalar pi(1-pi) attains its maximum 1/4 at pi = 1/2. Hence for a single sample the Hessian satisfies
>
> ```equation
> lambda_max( x x^T * pi(1-pi) ) <= ||x||^2 / 4
> ```
>
> Gradient descent on a convex function with Lipschitz-continuous gradient of constant L converges for step sizes 0 < eta < 2/L. Taking L = ||x||^2/4 gives
>
> ```equation
> 0 < eta < 8 / ||x||^2
> ```
>
> This is a bound derived before any computation, and it is falsifiable in exactly the way the analogous bound for a linear model is: step sizes below it should not diverge and step sizes far above it should.

An online form, consuming one individual at a time and discarding them, follows immediately from Result 14.11 with M = 1:

```equation
w <- w - eta ( pi_i - y_i ) x_i
```

Adding a ridge penalty (lambda/2)||w||^2 makes J strongly convex, bounds the condition number of H away from infinity, and is the standard remedy for the identifiability problem of the next subsection. It biases the estimate towards zero, which is a cost and not a free improvement.

## 14.9 Predictions and their failure modes

| Proposition | Derived from | Would be falsified by |
|---|---|---|
| Germline de novo count rises linearly with paternal age at roughly 1.5 per year | 14.4.1, anchored to [27] | A trio cohort with a materially different or nonlinear slope |
| A positive mutation rate makes the kernel structurally dense | Result 14.4 | Nothing; it is a counting argument. Its practical severity is measurable |
| CSR loses its payload advantage once mutation is enabled | Result 14.4 with Section 7.2 | A benchmark showing CSR still wins, which would mean eps was effectively zero |
| Methylation fraction approaches an equilibrium exponentially, not linearly | Result 14.7 | Methylation tracking a straight line across a wide age range with no curvature |
| w_met is negative at promoter islands | 14.7, from [31] | A fitted positive coefficient on correctly annotated promoter sites |
| Conditioning on m and h shrinks w_age towards zero | 14.7 | A large surviving age coefficient, indicating unmodelled mediation |
| Age-of-onset hazard has log-log slope k-1 for a k-stage process | Result 14.6, from [29] | An incidence curve whose slope is inconsistent with any integer k |

## 14.10 Identifiability: age is the hardest covariate there is

The proposal has one serious statistical weakness, and it is structural rather than incidental.

By Result 14.7, m_l(t) is a deterministic function of age up to noise. The design matrix therefore contains a column t and a column m_l(t) which is a smooth monotone transform of it. Over a narrow age range, exp(-(alpha+beta)t) is close to linear in t, so the two columns are close to collinear.

The consequence is standard and severe. Near-collinearity makes X^T S X near-singular; the estimator's covariance, which scales with the inverse of that matrix, becomes enormous in the direction of the offending combination; and the fitted weights w_age and w_met become individually meaningless while their sum remains well determined. Crucially, the objective value does not deteriorate. The model fits, and the coefficients that were the entire point of building it are noise.

No diagnostic computed from the fit alone distinguishes this case from a well-conditioned one. The defences must be structural and must be chosen before fitting:

- **Check the conditioning of the design before estimating anything.** The smallest eigenvalue of X^T S X, or the condition number, is computable in advance and is the quantity that decides whether the fit can answer the question.
- **Require within-age variation in methylation.** If every 40-year-old in the sample has the same m_l, the coefficient is not identified at any sample size. Sampling must be designed to break the collinearity, not merely enlarged.
- **Adjust for cell composition.** Bulk-tissue methylation changes with age partly because the mixture of cell types changes with age. Unadjusted, that confound is attributed to the locus.
- **Do not read causation from the sign.** Transcriptional state can drive methylation as well as follow it. A regression fixes no arrow, and this design cannot orient one.

## 14.11 Boundaries, numerical hazards, and the minimal next step

**Not implemented.** No code, no test, no dataset, no figure in this repository corresponds to any equation in this section. The artifact manifest of Appendix B lists what has been executed, and none of it is this.

**Not validated.** The parameters beta_f, beta_m and the ~40-per-year somatic rate are quoted from published cohorts, not estimated here, and their transfer to another population is an assumption. No value of w, b_l, alpha or beta has been estimated at all.

**Not clinical.** Result 14.10 describes a model in which penetrance varies with age and epigenetic state. It makes no prediction about any individual, and the restraint stated throughout Sections 8 to 11 applies here with more force, not less, because the model is now more suggestive and less tested.

**Numerically hazardous in a way this project has already met.** With eps_l of order 1e-8, the previously structural zeros become entries of order 1e-8 per mutated locus, and products across many loci fall below binary64 range quickly. The adversarial review already records that population updates refuse detected zero underflows rather than returning silent zeros, and that log-space population inference is not implemented. Enabling mutation would make that refusal path routine rather than exceptional. A log-domain population update is therefore a prerequisite for this extension, not an optional refinement.

**The minimal next step is small and worth doing.** Result 14.2 requires no new mechanism: `mutation(allele_frequencies, transition)` already accepts a validated row-stochastic matrix. Constructing that matrix from Result 14.2 with a stated eps_l, running the existing allele-copy oracle against the composed map, and measuring what Result 14.4 predicts about CSR occupancy on the existing benchmark harness would test the computational half of this section using apparatus that already exists. It would establish nothing about human biology, and it would settle the representation question the rest of the paper is about.

The epigenetic half cannot be approached that way. It needs paired methylation, acetylation, expression and age measurements in matched tissue, with cell-composition adjustment and a conditioning check performed before any coefficient is reported. Until that data is in hand, Sections 14.6 to 14.10 remain a derivation with the correct shape and no evidence.

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

| Artifact | Purpose | Evidence status |
| --- | --- | --- |
| research/version1_audit.md | Page-by-page discrepancy explanation | Full paper and screenshots inspected |
| results/version1_reproduction.json | Six exact legacy examples and coverage variants | Executed rational translation |
| genetics/core.py; extensions.py | Reusable kernel and staged extensions | Automated tests passed |
| results/benchmark.json and benchmark_workers/ | Raw timings, memory, metadata and checks | Executed local CPU run |
| data/observed_genotypes.csv + provenance | External genomic calls with validation | Frozen public snapshot |
| results/population_validation.json | HWE and pooling audit | Executed; no clinical accuracy claim |
| results/polygenic_synthetic.json | Score and interval calibration | 10,000 independent synthetic draws |
| sources/source_ledger.json | Reproducible claim/source bibliography | Primary and authoritative sources |
| manuscript/version2_manuscript.md | Editable evidence-derived candidate | Author review required |
| output/pdf/Genetics_Version_II.pdf | Typeset candidate | Rendered and visually checked before delivery |

The checks above concern delivered evidence. They do not assert that MATLAB was executed, that a clinical study was performed, that the software was benchmarked on multiple computers, or that the manuscript has been peer reviewed. A future revision should retain these distinctions when additional evidence becomes available.

# References

[1] Arshyia Mehran. Linear Algebra in Bioinformatics and Computational Biology: Modeling the Generational Decay of Sickle Cell Anemia and Changes in Generational Blood Types Using MATLAB. Local user-supplied PDF, 54 pages. Accessed 2026-09-05.

[2] M. A. Bender, Katie Carlberg. Sickle Cell Disease. 2025. [https://www.ncbi.nlm.nih.gov/books/NBK1377/](https://www.ncbi.nlm.nih.gov/books/NBK1377/) Accessed 2026-09-05.

[3] National Center for Biotechnology Information. NM_000518.5(HBB):c.20A>T (p.Glu7Val), ClinVar VCV000015333.7. [https://www.ncbi.nlm.nih.gov/clinvar/variation/VCV000015333.7/](https://www.ncbi.nlm.nih.gov/clinvar/variation/VCV000015333.7/) Accessed 2026-09-05.

[4] Laura Dean. The ABO blood group. 2005. [https://www.ncbi.nlm.nih.gov/books/NBK2267/](https://www.ncbi.nlm.nih.gov/books/NBK2267/) Accessed 2026-09-05.

[5] Laura Dean. The Rh blood group. 2005. [https://www.ncbi.nlm.nih.gov/books/NBK2269/](https://www.ncbi.nlm.nih.gov/books/NBK2269/) Accessed 2026-09-05.

[6] SciPy. scipy.sparse.csr_array. [https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_array.html](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_array.html) Accessed 2026-09-05.

[7] Robert Sedgewick, Kevin Wayne. Algorithms, 4th edition: 5.2 Tries. [https://algs4.cs.princeton.edu/52trie/](https://algs4.cs.princeton.edu/52trie/) Accessed 2026-09-05.

[8] Robert Sedgewick, Kevin Wayne. Algorithms, 4th edition: 3.4 Hash Tables. [https://algs4.cs.princeton.edu/34hash/](https://algs4.cs.princeton.edu/34hash/) Accessed 2026-09-05.

[9] Roberto Grossi, Giuseppe Ottaviano. Fast Compressed Tries through Path Decompositions. 2011. [https://arxiv.org/abs/1111.5220](https://arxiv.org/abs/1111.5220) Accessed 2026-09-05.

[10] Randal E. Bryant. Graph-Based Algorithms for Boolean Function Manipulation. 1986. IEEE Transactions on Computers C-35(8):677-691. DOI: 10.1109/TC.1986.1676819. [https://www.cs.cmu.edu/~bryant/pubdir/ieeetc86.pdf](https://www.cs.cmu.edu/~bryant/pubdir/ieeetc86.pdf) Accessed 2026-09-05.

[11] Python Software Foundation. heapq: Heap queue algorithm. [https://docs.python.org/3/library/heapq.html](https://docs.python.org/3/library/heapq.html) Accessed 2026-09-05.

[12] Frank R. Kschischang, Brendan J. Frey, Hans-Andrea Loeliger. Factor Graphs and the Sum-Product Algorithm. 2001. IEEE Transactions on Information Theory 47(2):498-519. DOI: 10.1109/18.910572. [https://haloeliger.github.io/papers/2001FG.pdf](https://haloeliger.github.io/papers/2001FG.pdf) Accessed 2026-09-05.

[13] Rina Dechter. Bucket elimination: A unifying framework for reasoning. 1999. Artificial Intelligence 113(1-2):41-85. DOI: 10.1016/S0004-3702(99)00059-4. [https://www.sciencedirect.com/science/article/pii/S0004370299000594](https://www.sciencedirect.com/science/article/pii/S0004370299000594) Accessed 2026-09-05.

[14] Python Software Foundation. timeit: Measure execution time of small code snippets. [https://docs.python.org/3/library/timeit.html](https://docs.python.org/3/library/timeit.html) Accessed 2026-09-05.

[15] G. H. Hardy. Mendelian Proportions in a Mixed Population. 1908. DOI: 10.1126/science.28.706.49. [https://pubmed.ncbi.nlm.nih.gov/17779291/](https://pubmed.ncbi.nlm.nih.gov/17779291/) Accessed 2026-09-05.

[16] Ensembl, The 1000 Genomes Project Consortium. rs334 population genotype and allele counts, 1000GENOMES:phase_3 subset. [https://rest.ensembl.org/variation/human/rs334?population_genotypes=1;pops=1;content-type=application/json](https://rest.ensembl.org/variation/human/rs334?population_genotypes=1;pops=1;content-type=application/json) Accessed 2026-09-05.

[17] The 1000 Genomes Project Consortium. A global reference for human genetic variation. 2015. DOI: 10.1038/nature15393. [https://www.nature.com/articles/nature15393](https://www.nature.com/articles/nature15393) Accessed 2026-09-05.

[18] Janis E. Wigginton, David J. Cutler, Goncalo R. Abecasis. A Note on Exact Tests of Hardy-Weinberg Equilibrium. 2005. DOI: 10.1086/429864. [https://pmc.ncbi.nlm.nih.gov/articles/PMC1199378/](https://pmc.ncbi.nlm.nih.gov/articles/PMC1199378/) Accessed 2026-09-05.

[19] Sture Holm. A Simple Sequentially Rejective Multiple Test Procedure. 1979. Scandinavian Journal of Statistics 6(2):65-70. [https://www.jstor.org/stable/4615733](https://www.jstor.org/stable/4615733) Accessed 2026-09-05.

[20] Laura Dean. The MNS blood group. 2005. [https://www.ncbi.nlm.nih.gov/books/NBK2274/](https://www.ncbi.nlm.nih.gov/books/NBK2274/) Accessed 2026-09-05.

[21] Laura Dean. The Hh blood group. 2005. [https://www.ncbi.nlm.nih.gov/books/NBK2268/](https://www.ncbi.nlm.nih.gov/books/NBK2268/) Accessed 2026-09-05.

[22] T. A. Brown. Mapping Genomes. 2002. [https://www.ncbi.nlm.nih.gov/books/NBK21116/](https://www.ncbi.nlm.nih.gov/books/NBK21116/) Accessed 2026-09-05.

[23] NIST/SEMATECH. e-Handbook of Statistical Methods: 7.2.4.1 Confidence intervals. [https://itl.nist.gov/div898/handbook/prc/section2/prc241.htm](https://itl.nist.gov/div898/handbook/prc/section2/prc241.htm) Accessed 2026-09-05.

[24] SciPy developers. scipy.stats.kstest. [https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kstest.html](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kstest.html) Accessed 2026-09-05.

[25] Loic Yengo, Sailaja Vedantam, Eirini Marouli, and colleagues. A saturated map of common genetic variants associated with human height. 2022. Nature 610:704-712. DOI: 10.1038/s41586-022-05275-y. [https://www.nature.com/articles/s41586-022-05275-y](https://www.nature.com/articles/s41586-022-05275-y) Accessed 2026-09-05.

[26] Augustine Kong, Michael L. Frigge, Gisli Masson, and colleagues. Rate of de novo mutations and the importance of father's age to disease risk. 2012. Nature 488:471-475. DOI: 10.1038/nature11396. [https://www.nature.com/articles/nature11396](https://www.nature.com/articles/nature11396) Accessed 2026-09-05.

[27] Hakon Jonsson, Patrick Sulem, Birte Kehr, and colleagues. Parental influence on human germline de novo mutations in 1,548 trios from Iceland. 2017. Nature 549:519-522. DOI: 10.1038/nature24018. [https://www.nature.com/articles/nature24018](https://www.nature.com/articles/nature24018) Accessed 2026-09-05.

[28] Francis Blokzijl, Joep de Ligt, Myrthe Jager, and colleagues. Tissue-specific mutation accumulation in human adult stem cells during life. 2016. Nature 538:260-264. DOI: 10.1038/nature19768. [https://www.nature.com/articles/nature19768](https://www.nature.com/articles/nature19768) Accessed 2026-09-05.

[29] Peter Armitage, Richard Doll. The Age Distribution of Cancer and a Multi-stage Theory of Carcinogenesis. 1954. British Journal of Cancer 8(1):1-12. DOI: 10.1038/bjc.1954.1. [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2007940/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2007940/) Accessed 2026-09-05.

[30] David R. Cox. Regression Models and Life-Tables. 1972. Journal of the Royal Statistical Society Series B 34(2):187-220. DOI: 10.1111/j.2517-6161.1972.tb00899.x. [https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/j.2517-6161.1972.tb00899.x](https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/j.2517-6161.1972.tb00899.x) Accessed 2026-09-05.

[31] Peter A. Jones. Functions of DNA methylation: islands, start sites, gene bodies and beyond. 2012. Nature Reviews Genetics 13:484-492. DOI: 10.1038/nrg3230. [https://www.nature.com/articles/nrg3230](https://www.nature.com/articles/nrg3230) Accessed 2026-09-05.

[32] Andrew J. Bannister, Tony Kouzarides. Regulation of chromatin by histone modifications. 2011. Cell Research 21:381-395. DOI: 10.1038/cr.2011.22. [https://www.nature.com/articles/cr201122](https://www.nature.com/articles/cr201122) Accessed 2026-09-05.

[33] Steve Horvath. DNA methylation age of human tissues and cell types. 2013. Genome Biology 14:R115. DOI: 10.1186/gb-2013-14-10-r115. [https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115) Accessed 2026-09-05.

[34] Gregory Hannum, Justin Guinney, Ling Zhao, and colleagues. Genome-wide methylation profiles reveal quantitative views of human aging rates. 2013. Molecular Cell 49(2):359-367. DOI: 10.1016/j.molcel.2012.10.016. [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3780611/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3780611/) Accessed 2026-09-05.

[35] John A. Nelder, Robert W. M. Wedderburn. Generalized Linear Models. 1972. Journal of the Royal Statistical Society Series A 135(3):370-384. DOI: 10.2307/2344614. [https://academic.oup.com/jrsssa/article/135/3/370/7110572](https://academic.oup.com/jrsssa/article/135/3/370/7110572) Accessed 2026-09-05.
