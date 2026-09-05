# Version II: inheritance architecture and mathematical specification

Status: research specification and independent mathematical checks, 2026-09-05. This document does not report unexecuted benchmarks or external predictive validation. Its new counting arguments are derived below; cited sources support the underlying algorithm families. The complete source ledger is `../sources/algorithm_sources.json`.

## 1. The object being represented

The defensible comparison is between implementations of the **same complete inheritance mapping**. The companion Version I audit reports that the original ABO and ABO+Rh MATLAB programs loop through all 36 and 324 ordered parental genotype pairs, respectively. Consequently, Version II must not claim that sparse storage restores outcomes that these loops actually omitted. A dense rectangular inheritance kernel supplies a transparent, reusable baseline for the complete computation already expressed procedurally in those loops. The separately reported coverage percentages require their own audit; they cannot establish the programs' support coverage.

A catalogued unphased diploid genotype at locus l is an unordered pair of allele labels. With a_l allowed alleles, there are

\[
g_l=a_l+\binom{a_l}{2}=\binom{a_l+1}{2}=a_l(a_l+1)/2.
\]

The full Cartesian catalog for n loci has

\[
G=\prod_{l=1}^{n}g_l.
\]

This is a count of labelled catalog entries, not an assertion that all entries occur in a sampled population, have equal probability, are viable, or are statistically independent. Population independence is not needed to count a Cartesian catalog. Constraints can reduce a particular model's admissible subset; phase can require an expanded state description.

For biallelic loci, each local catalog is {AA, Aa, aa}. **Induction:** for zero loci the empty tuple gives G(0)=1=3^0. Assume G(k)=3^k. Appending AA, Aa, or aa at a new locus gives three distinct extensions of each existing tuple; every (k+1)-locus tuple has exactly one k-locus prefix and one final choice. Thus G(k+1)=3G(k)=3^(k+1). Therefore G(n)=3^n for every nonnegative integer n.

With parent exchange symmetry, genotype-pair columns are indexed by 0 <= i <= j < G. There are G equal-genotype pairs plus G(G-1)/2 distinct-genotype pairs:

\[
U=G(G+1)/2.
\]

An equal-genotype pair means two parents with the same genotype; it does not require biological self-fertilization. Excluding such pairs is inappropriate for ordinary human random mating. If parental roles differ (for example sex-specific transmission or an explicit maternal effect), use G_m G_f ordered role-specific combinations instead. Catalog assumptions such as autosomal diploidy must be changed explicitly for sex-linked loci.

## 2. Conditional transmission versus population dynamics

For local genotype g=(u,v), define Mendelian gamete probability

\[
t_g(a)=\{1[a=u]+1[a=v]\}/2.
\]

The offspring probability for unordered genotype o=(a,b), a<=b, is

\[
T_l(o\mid g,h)=
\begin{cases}
t_g(a)t_h(a),&a=b,\\
t_g(a)t_h(b)+t_g(b)t_h(a),&a<b.
\end{cases}
\]

The second term aggregates distinct gamete paths to one heterozygous offspring state; it must not become a duplicate state. Each local distribution is nonnegative and sums to one. For independent segregation across the loci being modelled,

\[
T(o\mid g,h)=\prod_l T_l(o_l\mid g_l,h_l).
\]

Define K[o,c(i,j)]=T(o|i,j). Its shape is **G by U**; each column sums to one. The ordered form has shape G by G^2 and a three-index view T[o,i,j]. Rectangular matrices are standard linear maps between different spaces. There is no mathematical requirement to discard columns to make K square, and no requirement to materialize K merely to evaluate it.

For a normalized distribution x over full parental genotypes, independent random draws of two parents give

\[
w_{ij}(x)=\begin{cases}x_i^2&i=j,\\2x_ix_j&i<j,\end{cases}
\qquad x'=K w(x).
\]

Because sum(w)=(sum x)^2=1, normalization is preserved. The mapping x->x' is generally **quadratic**, although K acts linearly on a supplied pair distribution w. Replacing this by a fixed square matrix power changes the model in general. The counterexample is one biallelic locus: x=(1,0,0) and y=(0,0,1) reproduce themselves, whereas equal mixing before random mating produces (1/4,1/2,1/4), different from (F(x)+F(y))/2=(1/2,0,1/2).

A square linear operator is valid in a different, clearly stated model. For fixed mate genotype distribution q,

\[
M_q[o,i]=\sum_j T(o\mid i,j)q_j,
\qquad x'=M_qx.
\]

M_q is G by G and column-stochastic. If q and the transmission law remain fixed over generations, x_t=M_q^t x_0 is valid for that stipulated mating model. A changing q_t gives a sequence of operators. General assortative mating can instead supply a normalized joint parent distribution pi(i,j); marginal distributions alone then do not determine the offspring distribution. Selection, mutation, survival and migration belong in named additional operators or weighting steps, with a stated order.

## 3. Exact support counts and why sparse is still exponential

For one biallelic locus, count possible offspring for each **ordered** parent pair:

| First parent / second parent | AA | Aa | aa |
|---|---:|---:|---:|
| AA | 1 | 2 | 1 |
| Aa | 2 | 3 | 2 |
| aa | 1 | 2 | 1 |

The entries sum to 15. For independent loci, supports form Cartesian products. Distributing the sums over all parental tuples therefore gives

\[
N_{ordered}(n)=15^n.
\]

This counts supported triples (offspring, first-parent, second-parent), not unique offspring states. Swapping the parents partitions the supported triples into two-element orbits except when the **complete** parent genotypes are identical. For identical parents, the one-locus support counts are 1,3,1 and sum to 5. The full fixed-point count is 5^n. Hence

\[
N_{unordered}(n)=(15^n+5^n)/2.
\]

This is not 10^n: globally unordered full parental genotypes cannot be formed by independently forgetting parental identity at each locus. For example, local swaps can exchange which alleles occur together in the same parent's multilocus genotype.

The dense rectangular kernel stores

\[
G U=(27^n+9^n)/2
\]

numeric slots, giving density

\[
\rho_n=(15^n+5^n)/(27^n+9^n)\sim(5/9)^n.
\]

Density falls, but the sparse support still grows as Theta(15^n). A lower density alone is not evidence of practical scalability.

**Generalization of support counts.** At one locus with a alleles, let H=a(a-1)/2 be the number of heterozygotes and s(g) the number of alleles a parent can transmit (one or two). Sum_g s(g)=a+2H=a^2. For parent genotypes g,h the offspring support size is s(g)s(h), except that equal heterozygotes have one duplicate unordered offspring caused by swapping transmitted alleles. Therefore

\[
T(a)=\sum_{g,h}|supp(T(\cdot\mid g,h))|=a^4-H.
\]

The equal-parent support sum is D(a)=a+3H. For independent loci with arbitrary allele counts,

\[
N_{ordered}=\prod_l T(a_l),\qquad
N_{unordered}=\tfrac12\left[\prod_l T(a_l)+\prod_l D(a_l)\right].
\]

These identities assume every catalogued allele is transmitted under the simple Mendelian law, all offspring genotypes remain in the catalog, parent exchange symmetry, and no mutation or viability removal. A positive mutation process may add support; selection may remove support. Do not reuse the identities after altering the law without rederivation.

An independent explicit check enumerated allele-pair catalogs, then constructed each local support by all four gamete choices and set deduplication. It did not call a kernel implementation. Results:

| Allele counts | G | U | Ordered nonzeros | Equal-parent nonzeros | Unordered nonzeros |
|---|---:|---:|---:|---:|---:|
| (2) | 3 | 6 | 15 | 5 | 10 |
| (2,2) | 9 | 45 | 225 | 25 | 125 |
| (2,2,2) | 27 | 378 | 3,375 | 125 | 1,750 |
| (3), simplified ABO | 6 | 21 | 78 | 12 | 45 |
| (3,2), simplified ABO+Rh | 18 | 171 | 1,170 | 60 | 615 |

The numerical check was executed in Python on 2026-09-05. These are exact enumerated support counts, not speed or biological performance measurements. Simplified ABO and Rh refer to the educational allele catalogs; the counts do not claim to exhaust molecular blood-group variation.

## 4. State and probability coverage are different quantities

Use named denominators, not an unexplained coverage percentage:

- **Genotype catalog coverage:** represented genotype labels / G.
- **Parent-pair coverage:** represented parent-pair labels / U under the declared symmetry.
- **Supported-transition coverage:** retained nonzero triples / complete supported triples.
- **Probability mass retained:** sum of probabilities of retained events under a declared, normalized joint distribution. For removed parent-pair columns this is sum of their supplied w entries; for removed offspring entries it is sum over both columns and retained outcomes of w[c]K[o,c].
- **Observed population coverage:** proportion of a specified external sample matched by a precisely stated event or eligibility criterion. It requires actual sample data and a denominator.
- **Predictive performance:** a held-out comparison against observed outcomes, using metrics appropriate to the target. Exact equality between two implementations is computational agreement, not independent biological predictive accuracy.

Uniformly counting pair columns is not weighting them by a population. Even if every genotype has probability 1/G, unordered equal-genotype and distinct-genotype pairs have respective probabilities 1/G^2 and 2/G^2. A representable state can have zero probability for a current input while remaining supported under other parents. Structural zeros need no stored numeric values. A rare positive probability must not be erased by a sparsity threshold and then called impossible.

## 5. Candidate representations and their actual workloads

Let N be the number of supported transitions and d_r the number stored in row r. These costs concern the declared representation; probability precision and arbitrary-sized integer encodings have their own costs.

| Candidate | Representation and operations | Appropriate role | Failure mode or limitation |
|---|---|---|---|
| Dense G by U array | Theta(GU) numeric storage; O(1) indexed lookup; Theta(GU) matrix-vector application | Small, fully inspectable truth baseline; vectorized compiled arithmetic | Stores structural zeros; complete table grows approximately 27^n/2 |
| CSR sparse array | Theta(N+G) numeric/index storage; sorted-row lookup O(log d_r); matvec Theta(N+G) | Repeated complete updates after construction | Still approximately 15^n/2 entries; changing sparsity is costly |
| CSC sparse array | Theta(N+U) storage; cheap column access | Enumerating a selected pair's offspring | More column pointers for this wide kernel; still exponential |
| Adjacency lists / pair-to-offspring graph | Theta(N+U) entries; traverse one pair's degree; accumulate whole output in O(N+G) with dense result | Transparent streaming construction and selected-pair access | Python node/edge object overhead can dominate compact sparse arrays; inheritance across generations need not be a DAG unless generation is part of the node |
| Hash map | Expected constant table probe for pre-encoded fixed-size keys; O(n) to build/hash an n-locus tuple; Theta(N) entries | Dynamic assembly, duplicate aggregation, irregular state keys | Collisions and load factor matter; worst cases and object memory preclude an unconditional O(1) claim |
| Trie / radix structure | O(n) depth for an n-locus key and shared prefixes; output-sensitive enumeration | Prefix, wildcard and partial-genotype queries | A full ternary catalog trie has (3^(n+1)-1)/2 nodes, so prefix sharing alone does not remove exponential growth |
| Reduced decision DAG | Merge genuinely identical suffix subfunctions, with variable ordering and weights | Repeated logical/conditional structure | Possible exponential graph size, sensitivity to variable order, and potentially expensive reductions; no universal compression theorem for genetics |
| Priority queue | Binary-heap insert/pop O(log h) at frontier size h | Lazy highest-probability outcome enumeration, scheduling or bounded search | It is an ordering device, not a complete probability model; top-k output omits mass unless remaining support remains queryable and omitted mass is tracked |
| DP / bounded memoization | Store and reuse repeated local subproblems; cost follows number and size of distinct cached problems | Local crossing tables, score sums, pedigree subproblems | Caching each full parent pair merely stores another exponential table; unbounded caches can exhaust memory |
| Factor tables / factor graph | Local factors and structured marginalization | Independent-locus queries; dependency-aware extension | Small model storage does not make arbitrary full joint output or dense evidence cheap; complexity follows induced width |

CSR's format and practical operation tradeoffs follow the [SciPy CSR-array documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_array.html). Hashing and trie query semantics follow the authors' [Algorithms hash-table](https://algs4.cs.princeton.edu/34hash/) and [trie](https://algs4.cs.princeton.edu/52trie/) references. A heap is available through [Python heapq](https://docs.python.org/3/library/heapq.html). These references justify mechanisms, not an empirical ranking on this task.

Succinct tries can reduce dictionary overhead using specialized encodings, but their published workload is string indexing, not inheritance probability inference. Treat them as candidates requiring a matched implementation and benchmark. [Grossi and Ottaviano, 2011](https://arxiv.org/abs/1111.5220). The reduced-DAG analogy is motivated by sharing equivalent subfunctions; Bryant's foundational Boolean representation also documents ordering sensitivity and exponential worst cases. A weighted multi-valued inheritance DAG would be a further design, not an implementation of Bryant's binary scheme by assertion. [Bryant, 1986](https://www.cs.cmu.edu/~bryant/pubdir/ieeetc86.pdf).

For a dense binary64 array, raw payload is 8GU bytes. CSR raw buffers use b_v N + b_i N + b_p(G+1), where b_v, b_i, b_p are the actual value, index and pointer item sizes. Use measured dtype sizes rather than assuming 32-bit indices. Include arrays of states/pairs and build-time coordinate buffers separately; Python wrapper, allocator, native-library, and peak process memory are not counted by buffer bytes alone. Preflight the dimensions and a conservative peak-memory estimate before construction. A rejected size must return a clear resource-limit status and offer on-demand inference; it must not silently return a truncated model.

## 6. Why factorization is the scalable path for selected queries

For a specified parent pair, storing each local conditional distribution needs O(n) space when allele counts are bounded. A fully specified offspring probability is the product of n local terms, hence O(n) arithmetic operations. For a partial observation of a subset S, multiply the observed local terms and sum each unobserved factor to one. An arbitrary full offspring distribution can still have 3^n nonzero entries, so explicitly returning it cannot have O(n) total cost.

This is an instance of marginalization of local factors, with sum-product messages reusing intermediate results. Cycle-free factor graphs permit exact marginal computation; loopy message passing is not generally exact merely because the local factors are correct. [Kschischang, Frey and Loeliger, 2001](https://haloeliger.github.io/papers/2001FG.pdf). General variable elimination's largest intermediate factor depends on the elimination order and induced width. For bounded domain size d and induced width w, a conventional bound has exponential dependence on w (commonly O(n d^(w+1)) time with factor-storage dependence of the same exponential order). There is no generic polynomial-in-n inference guarantee for an arbitrary genetic dependency graph. [Dechter, 1999](https://www.sciencedirect.com/science/article/pii/S0004370299000594).

Under **both** factorized parental genotype distribution x(g)=product_l x_l(g_l) and independent local transmission, the random-mating offspring distribution also factorizes:

\[
x'_l(o_l)=\sum_{i_l,j_l}T_l(o_l\mid i_l,j_l)x_l(i_l)x_l(j_l),
\quad x'(o)=\prod_l x'_l(o_l).
\]

For bounded allele counts this gives O(n) local update work and storage if the answer stays factored. This is a mathematically exact baseline for the linkage-equilibrium workload, not a replacement for an arbitrary correlated x. A small dense K and CSR K must also match this method where its assumptions hold.

**Dependency counterexample.** If a two-locus parent distribution places probability 1/2 on AA/BB and 1/2 on aa/bb, its marginals do not identify that association. Independent-locus products would put positive probability on parental AA/bb and aa/BB even though the stated joint population has none. Factorizing x solely because transmission factors is invalid. The correlated joint input must be retained in the factors or handled by an exact complete method.

**Linked-locus state limitation.** The unphased Aa/Bb label can stand for coupling haplotypes AB/ab or repulsion Ab/aB. Under the stipulated two-locus recombination model with fraction r, coupling transmits AB and ab at (1-r)/2 each and Ab and aB at r/2 each; repulsion interchanges those roles. At r=0 the two phases have disjoint gamete supports; at r=1/2 both produce four equiprobable gametes. Thus a linked extension needs phased haplotypes or a posterior over phase, plus explicit recombination factors. Linkage equilibrium of a population and independence of meiotic segregation are separate assumptions. The examples here expose limitations of the model; biological linkage sources belong in the companion biological ledger.

Recommended architecture: one canonical genotype/allele encoding and validated local transmission law; a small exact rational reference; complete dense and compressed sparse materializations for bounded dimensions; then an on-demand local-factor interface. General DAG, trie and priority-queue methods remain candidate extensions unless their promised query workload is implemented and measured. Do not label a product-of-local-tables implementation a general pedigree/factor-graph solver.

## 7. Fair, falsifiable benchmark design

Workload A: build the **same full G by U kernel** in dense and CSR representations, then apply each to the same supplied w. Workload B: calculate complete random-mating updates, including the time and memory of forming w(x). Workload C: query a specified offspring or partial observation for specified parents. Workload D: update a factorized population and retain a factored result. Report A-D separately. Timing a single factored query against construction of a complete dense kernel is a useful end-to-end use-case comparison only if its differing outputs and setup costs are explicit; it is not evidence of faster full-kernel construction.

For each bounded n and allele catalog: verify identical canonical state/pair ordering, shape, support, column sums, and output values before timing. Include homogeneous and mixed allele counts. Add random full joint input populations and selected correlated counterexamples; run the factored population baseline only on inputs it can represent exactly. All methods must retain rare positive outcomes under the same probability law.

Record machine/OS/CPU, RAM, interpreter and package versions, numeric/index dtypes, BLAS provider/thread settings, seeds, exact arguments, warmups, repetitions, and raw timings. Separate build from reuse and report amortization. Use repeated measurements after warmup with a monotonic performance counter; retain raw values and report minimum as a low-interference microbenchmark bound plus median and range as descriptive run variability. Do not infer statistical significance from a few timing repetitions. Python's [timeit documentation](https://docs.python.org/3/library/timeit.html) explains setup exclusion, repetition and lower-bound interpretation.

Record payload bytes separately from peak process memory; do not rename payload size to RAM consumed. If using tracemalloc, label it as traced Python allocation and recognize that it is not a complete native/process memory measurement. Ensure both dense and sparse outputs are actually consumed and checked. Consider independent processes and randomized method order for expensive runs. Avoid exhausting the machine just to demonstrate the predicted asymptote. Document skipped dimensions and allocation refusal as censored feasibility results, not infinite measured runtime. Do not extrapolate plotted theoretical estimates as measured data.

Adversarial checks: repeated allele labels; homozygote gamete duplication; parent swap; equal-genotype pair weight; mutation-created support; invalid or unnormalized input; missing phase; unjustified factorization; float underflow in long products; full-kernel allocation limits; hash-key ambiguity; score-bin truncation; and disagreement between structural support and a thresholded floating representation. Long selected-outcome calculations should expose log probabilities so a nonzero mathematical event is not described as impossible after floating underflow.

## 8. Bounded synthetic additive-score extension

After Mendelian cases pass, an additive-score example can investigate inference without pretending to predict human height. Let D_l be a synthetic allele dosage in {0,1,2}, with conditional probabilities p_l(d) from the validated local inheritance model. Stipulate independent D_l and chosen coefficients beta_l:

\[
S=\sum_l\beta_lD_l,
\quad E[S]=\sum_l\beta_l E[D_l],
\quad Var(S)=\sum_l\beta_l^2 Var(D_l).
\]

For correlated loci add 2 sum_(l<m) beta_l beta_m Cov(D_l,D_m); the product method is then invalid unless dependencies are represented. Effects, allele frequencies, environmental variance, and any phenotype mapping must be labelled synthetic if chosen rather than estimated from independent data.

If beta_l=delta b_l with integer b_l, the score law is the coefficient sequence of

\[
P(z)=\prod_l\left[p_l(0)+p_l(1)z^{b_l}+p_l(2)z^{2b_l}\right].
\]

A shifted array handles negative exponents. Dynamic programming convolves the current distribution with each three-point local distribution. It stores O(L) values for an L-bin score range and needs O(sum_l L_l) work with the three-term sparse local update. This avoids enumerating 3^n genotypes but is only efficient when the score range is bounded: coefficients of exponentially increasing size can still create exponential range. Grouped or FFT convolution is an optional numerical implementation, not automatically faster for three-term inputs; direct/FFT choices are described by [SciPy convolve](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.convolve.html).

Non-grid coefficients require either an exact sparse map (possibly exponentially many distinct sums), discretization with a reported error, sampling, or an approximation. If each beta is rounded within delta/2 and D_l<=2, the worst-case score perturbation is at most n delta. This is a score-axis bound, not a probability or phenotype-accuracy guarantee. Handle deterministic loci, negative coefficients, normalization and bin endpoints explicitly.

An illustrative continuous outcome can be defined by the **assumed model** Y=mu_0+S+epsilon, epsilon independent N(0,sigma_e^2). Conditional on exact score weights w_s, this gives the normal-mixture CDF

\[
F_Y(y)=\sum_s w_s\Phi((y-\mu_0-s)/\sigma_e),\qquad \sigma_e>0.
\]

This mixture is exact for the stipulated residual law; replacing the whole distribution by one normal with mean mu_0+E[S] and variance Var(S)+sigma_e^2 is a separate approximation. A positive variance and many loci alone do not guarantee an accurate normal approximation when a few effects dominate or the distribution is multimodal. For sigma_e=0 use the discrete CDF directly.

Validation should compare convolution with exact enumeration for small n, verify analytic moments, and compare large-n simulation to the stipulated distribution with Monte Carlo uncertainty. Report CDF or quantile errors of any normal approximation against the numerical reference. Simulation-based interval coverage checks implementation under its own generative assumptions; they do not establish calibration on humans. There are no justified height, clinical, ancestry-transfer or real-world predictive claims without suitable independent genotype/phenotype data, estimated effects, data-use permission and external validation.

## 9. Bounded conclusion supported by this analysis

Complete inheritance can be represented without a square-table restriction. Sparse storage removes structural-zero payload and can improve feasibility, but full Mendelian support remains exponential. For independent-locus selected queries and factored populations, retaining local factors avoids constructing the global state table and is exact under explicit assumptions. No theory here establishes that a trie, DAG, hash table or sparse matrix is universally fastest, or that a computationally complete model is biologically complete. Version II's contribution should be framed as a reproducible comparison and correction of definitions, with empirical claims limited to the experiments actually executed.
