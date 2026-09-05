# Adversarial review of the Version II candidate

Review date: 2026-09-05. Scope: read-only inspection of `genetics/core.py`, `genetics/extensions.py`, and `tests/test_core.py`, independent mathematical support counting, and executed counterexamples. This is an engineering and mathematical review, not validation of a clinical predictor. The root author made implementation fixes; this reviewer independently reran the specific numerical counterexamples described below. Root may append later verification and benchmark findings without rewriting this history.

## Found, fixed, and independently rechecked

### Positive lazy outcomes could underflow to structural zero

The original `iter_cross` multiplied ordinary probabilities and described its output as positive support. For 1,000 independent biallelic loci with both parents heterozygous everywhere, the first mathematically positive offspring probability is 4^-1000, below binary64 range. The original iterator emitted `(0, 0.0)` even though `probability` already refused that underflow.

The revised iterator raises `FloatingPointError`, and `iter_log_cross` retains a finite log probability. Recheck: the first log outcome was `(0, -1386.2943611198905)`, agreeing with -1000 log(4). This fixes the tested conditional lazy-query path. It does not by itself protect all population-weight arithmetic; see the outstanding numerical limitation below.

### An absolute tolerance invalidated very small exact HWE p-values

The original `hwe_exact` included conditional probabilities below `observed + 1e-12`. This made approximately 1e-12 a spurious tail-selection floor even when the observed probability was orders of magnitude smaller. It was not merely harmless floating rounding.

The revised code compares conditional log weights against the observed log weight. Verification used **independent exact integer combinatorial weights**, not the implementation's SciPy gamma/log formula. For N individuals, h heterozygotes and b second-allele homozygotes, the integer weight is `comb(N,h) * comb(N-h,b) * 2**h`. Feasible counts keep the observed total allele-copy count fixed. The two-sided exact numerator sums weights no larger than the observed integer weight.

| Genotype counts | Original result | Exact integer reference | Rechecked revised result |
|---|---:|---:|---:|
| (0,100,0) | 1.8729491943077044e-13 | 1.51139082730558e-29 | 1.5113908273056113e-29 |
| (0,1000,0) | 1.790382926773772e-12 | 5.3635552016991766e-300 | 5.3635552016929704e-300 |
| (500,0,500) | 1.790382926773772e-12 | 1.3196690976572097e-301 | 1.3196690976571266e-301 |

The remaining differences are consistent with floating evaluation of log weights. A returned binary64 zero for a still smaller mathematical p-value is a numeric range limit, not exact impossibility. Exposing a log p-value and summing the selected tail through logsumexp before exponentiating would strengthen this API further.

### Selection confused underflowed reproductive mass with zero fitness

The original `select([.5,.5], [5e-324,5e-324])` multiplied before normalization, underflowed both terms to zero, and raised `No reproductive mass remains`. The mathematically normalized selected population is unchanged.

The revised code computes relative selected weights in log space, detects genuine zero support separately, and refuses a positive normalized frequency that underflows. Independent recheck returned `[0.5, 0.5]` for the counterexample.

### Additive DP implementation did not match its intended sparse-update complexity

The original score DP built a zero-filled local kernel of length `2*w+1` and called dense `np.convolve`. Only three local coefficients can be nonzero, so large integer effect spacing imposed avoidable arithmetic. Root replaced it with three shifted additions. This matches the O(sum of intermediate score-array lengths) three-term DP analysis for fixed three-value dosage factors. Correctness remains subject to the regression suite and small enumeration oracle; the earlier test suite already provided a useful independent small-score comparison.

## Additional checks that passed

- The allele-copy enumeration in the tests does not call the implementation's local-cross probability table. Exhaustive small-locus comparison is therefore a meaningful reference for aggregation and probability errors, although it still shares catalog encode/decode utilities.
- Independent support enumeration confirmed biallelic unordered nonzero counts 10, 125 and 1,750 for one, two and three loci. It also confirmed 45 for the three-allele ABO catalog and 615 for the (3,2)-allele ABO+Rh catalog. Proofs appear in `architecture_theory.md`.
- Equal-genotype mating pairs are retained. Their random-mating weights use x_i^2, while distinct-genotype unordered pairs use 2 x_i x_j.
- A factorized population with local frequencies `(.1,.4,.5)` and `(.2,.3,.5)` gave full versus linkage-equilibrium factor-update maximum absolute disagreement 8.326672684688674e-17. This supports agreement for the eligible workload; it says nothing about replacing a correlated joint population by marginals.
- The phase example distinguishes coupling and repulsion at zero recombination and makes their gamete laws agree at recombination fraction 1/2. This is a useful falsification test against an unphased linked-locus implementation.
- The model has explicit size preflight checks and an injected allocation-failure test for dense kernel construction. This is narrower than claiming every allocation path is memory-safe.

## Outstanding at the review snapshot

1. **Population products can silently underflow.** With the accepted binary64 input `[1., 1e-200, 0.]` at one biallelic locus, the small heterozygote-pair mass is about 1e-400 and the aa offspring mass about 2.5e-401. `pair_weights` sets that pair to zero, and the streamed population update returns zero for aa. The first input component has rounded to one, so this example is understood within the API's probability-sum tolerance; choosing that tolerance does not make the positive pair impossible. The full kernel still represents aa as possible, but the population value is below the chosen output format. Add product-underflow refusal or log population APIs, or document explicitly that conditional selected-outcome log support does not extend to all population computations. Do not claim that all rare positive probability values are preserved numerically.

2. **The additive bin limit accepted NaN.** Executed call `additive_pmf([[.25,.5,.25]], [2], max_bins=float('nan'))` returned a five-bin result because comparison with NaN is false. A non-finite limit can disable the intended allocation preflight. Validate a strictly positive integer limit before allocating. Related HWE enumeration and mixture-CDF helpers have no general configurable memory budget; these should be described as bounded experiment helpers unless separately guarded.

3. **Supplied kernels are trusted.** `next_generation` accepts an external numeric kernel without explicitly checking full shape and stochastic semantics, and a dictionary can be incomplete or contain arbitrary values. Internal `model.kernel` outputs are validated through tests, but the public method contract should say it trusts this representation or validate inputs at the boundary. Rechecking every stochastic invariant on every timed matrix-vector call is not required; validation can be separated from repeated application.

4. **Allocation handling is not global.** The explicit try/except conversion to `ResourceLimitError` covers cross and kernel materialization, while initialization, `pair_weights`, population output, additive DP, and HWE arrays contain other allocations. Python `MemoryError` still fails without a valid completed result, but a claim that all these paths return the project's custom resource diagnostic would be too broad. Size budgets are operation estimates, not guarantees about OS-wide RAM availability or concurrent work.

5. **Numeric exactness must be scoped.** Small Mendelian probabilities are exactly representable dyadic values, and small rational enumeration supplies an exact reference. Arbitrarily long products, selection normalization, HWE gamma/log calculations, and score convolutions use floating arithmetic. The phrase "exact method" should distinguish an exact mathematical algorithm from exact arithmetic across arbitrary input sizes.

These items are reproducible behavior or documented contract limits, rather than hypothetical claims of an incorrect ordinary small-case result. Root should mark any later fixes with the corresponding rerun evidence.

## Manuscript claims that should survive falsification

The companion Version I audit reports complete ordered-pair loops in the original ABO (36 pairs) and ABO+Rh (324 pairs) MATLAB programs. Version II should therefore compare optimized representations of the complete existing operation. It should not claim to recover omitted outcomes from those actual loops. The prose coverage percentages are a separate claim to audit. This reviewer did not independently re-extract the entire original PDF; the paper-specific statement relies on the named companion audit, while the combinatorial formulas here were independently derived.

Sparse full-kernel storage remains exponential despite falling density. Factored queries avoid full materialization only for specified outputs under specified conditional-independence assumptions. Trie prefix sharing, a heap, and a DAG each solve different subproblems; none is a demonstrated universal replacement for a complete probability model. General factor-graph inference, pedigree analysis and a weighted decision-DAG engine should remain future work unless actually implemented.

Observed genotype counts can test consistency of a population model or estimate allele frequencies. Reconstructing HWE expectations using allele frequencies estimated from those same counts is not held-out predictive validation. A difference between sites or cohorts is not automatically selection or linkage disequilibrium; sampling and population structure can explain differences. External inheritance prediction requires appropriate independent family/offspring observations.

The synthetic additive-score extension validates distribution algorithms under chosen effects, independence and residual assumptions. It does not calibrate height prediction or demonstrate transfer across real populations. Its normal residual is an assumed component of a simulated generative model; a normal approximation to the whole output is a different approximation that needs separate error reporting.

## Benchmark review pending harness availability

At this review snapshot, the dedicated benchmark harness had not yet been available in `experiments/`. Its core comparability checks should be: identical full rectangular kernel, identical pair weights and outputs, build versus reuse separation, measured payload versus peak memory naming, all raw timing repetitions, documented thread settings and installed versions, and separate workloads for factored selected queries and factored populations. Memory-based skips are predicted feasibility limits, not timed failures. The small original 6 by 6 or 18 by 18 table should not be presented as a complete full-pair dense-kernel baseline when its state semantics differ.

No speedup, peak-RAM comparison, predictive accuracy, or publication-readiness claim is certified by this file until the corresponding evidence has been reviewed.

## Follow-up verification and benchmark-harness review

Later on 2026-09-05, root added population-product underflow checks, positive-integer bin-limit validation, and numeric-kernel shape validation. Targeted independent reruns confirmed:

- `next_generation([1.,1e-200,0.])` now raises `FloatingPointError` for streamed, dense, CSR and hash representations. Outstanding item 1's demonstrated silent-underflow defect is fixed for these tested paths. Log-space population inference is still not implemented; refusal is an honest result when binary64 cannot represent the requested value.
- Additive `max_bins` values NaN, positive infinity, -1, 0 and True each now raise `ValueError`. Outstanding item 2's demonstrated bin-limit bypass is fixed. Its separate note about unbudgeted helper allocations remains applicable.
- A supplied numeric kernel of shape (2,6) for the three-genotype model now raises `ValueError: Kernel shape mismatch`. This fixes the shape-validation portion of item 3; arbitrary external numeric values and dictionary contents are still trusted.

Other population helper paths, notably explicit mutation/migration and fixed-mate matrix construction, should not be described as comprehensively log-domain implementations. This review confirms the concrete repaired counterexamples, not a theorem that every floating operation preserves every possible positive contribution.

The new `experiments/benchmark.py` was then reviewed read-only before execution. No run-blocking comparability defect was found. It uses fresh subprocesses, clears the local-cross LRU before each timed construction, separates memory tracing from timing, enforces one detected numerical-library thread, and uses identical full-support seeded input vectors. Its n=6 preflight distinguishes actual guard refusal from a preregistered n<=5 scope limit. Raw per-worker evidence and source hashes are retained.

The measured complete-update operation includes probability validation, pair-weight formation and an additional reachability multiplication; these are end-to-end implementation timings, not isolated matrix-vector multiplication timings. A later claim about pure CSR-versus-dense matvec speed would require a different measurement. The selected conditional log-query workload is explicitly separate from the complete population workload. The linkage-equilibrium population method is verified mathematically and by a small numerical comparison but is not benchmarked by this harness.

Minor reporting issue communicated to root: if every construction sample for a materialized representation fails, its summary currently labels construction `not_applicable`; `failed_or_missing` would be more precise. The overall report already becomes incomplete when a worker fails, so this does not make a failed overall run pass. The fixed method order and uncontrolled desktop background load remain ordinary limits on fine performance comparisons. Small-n single-update samples may be dominated by Python/API overhead; they should be interpreted as the user-visible implemented operation.

This follow-up approves the methodology for bounded execution, not results that had not yet been produced.

A subsequent cross-review by the Version I audit agent identified that importing `scipy` alone does not eagerly load `scipy.sparse`. The first CSR memory pass could therefore include sparse-module import allocations after tracing started, contrary to an all-imports-excluded scope description. Timed construction already had a warmup, so this specifically affects memory-scope comparability and potentially early threadpool metadata. The agreed correction is to import `scipy.sparse` in every worker before metadata and tracing, then collect final results under a consistent revised harness hash. Verification of that harness revision belongs with the final benchmark evidence.
