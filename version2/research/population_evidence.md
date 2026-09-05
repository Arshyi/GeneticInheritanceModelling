# Biological scope and population evidence for Version II

Research evidence checked on 2026-09-05. The machine-readable companion is `sources/population_sources.json`. This note separates textbook inheritance laws, newly derived examples, published reference-panel calls, and statistical model checks. It provides no clinical prediction validation.

## What the Mendelian examples represent

The six ABO genotypes are complete **within the chosen A/B/O allele vocabulary**. They are not a complete catalogue of molecular ABO alleles. The common phenotype map uses codominant A and B expression and recessive O; weak subgroups and the H-antigen pathway make that map a simplification. The reference to six genotypes does not justify dropping any parental cross. [Dean, ABO chapter](https://www.ncbi.nlm.nih.gov/books/NBK2267/) (`dean_abo_2005`).

Likewise, the three labels DD/Dd/dd represent a teaching model of D-antigen presence/absence. The wider Rh system involves closely linked RHD and RHCE, rearrangements, deletion/pseudogene backgrounds, and weak/partial D variants. A mathematically complete 18-state ABO plus D/d model therefore remains biologically restricted. The 2005 source is used for these stable mechanisms; its antigen inventories and clinical recommendations are not presented as current practice. [Dean, Rh chapter](https://www.ncbi.nlm.nih.gov/books/NBK2269/) (`dean_rh_2005`).

The A/S sickle model is an instructive biallelic HBB example. Sickle cell disease also includes compound heterozygosity with other HBB pathogenic variants; transmission probabilities alone cannot determine clinical severity. In the external dataset, a T call at rs334 means only that the HbS alternate A is absent at that coordinate. It does **not** establish a normal HBB sequence. [Bender and Carlberg, GeneReviews, updated 2025-02-13](https://www.ncbi.nlm.nih.gov/books/NBK1377/) (`bender_carlberg_scd_2025`).

### Additional Mendelian system: GYPA M/N

GYPA M and N are codominant, so the restricted diploid vocabulary is MM, MN, NN. Two MN parents yield MM:MN:NN probabilities 1/4:1/2:1/4 by Mendelian segregation. This is a newly calculated cross, not a measured cohort distribution. It covers the M/N teaching example, not the whole MNS blood-group system: GYPB, linked loci, rearrangements, and additional variants matter to the broader system. [Dean, MNS chapter, molecular-basis table and GYPA section](https://www.ncbi.nlm.nih.gov/books/NBK2274/) (`dean_mns_2005`).

### Epistasis example: ABO plus FUT1

FUT1 activity supplies the erythrocyte H precursor required for A/B antigen formation. Two inactive copies can suppress A/B expression despite the ABO genotype. FUT2 and partial H variants add distinctions outside this two-locus teaching map. H-deficient red cells must not be equated clinically with ordinary blood group O. [Dean, Hh chapter](https://www.ncbi.nlm.nih.gov/books/NBK2268/) (`dean_hh_2005`).

For the **specified synthetic cross** AO/Hh × BO/Hh, assume independent assortment of the two loci, equal segregation, no mutation, and the binary FUT1 activity map. ABO outcomes AB/AO/BO/OO each have probability 1/4; the probability of at least one H allele is 3/4. Thus:

| Distinct model outcome | Probability |
|---|---:|
| A expressed, H present | 3/16 |
| B expressed, H present | 3/16 |
| AB expressed, H present | 3/16 |
| Ordinary ABO O with H present | 3/16 |
| H deficient (hh), regardless of ABO genotype | 4/16 |

If an explicitly coarse observable reports only whether A/B antigens are expressed, its combined “O-like” outcome is 3/16 + 4/16 = 7/16. Retain the separate 3/16 ordinary-O and 4/16 H-deficient states internally. This gives a useful architectural lesson: independence of inheritance factors does not imply independence of phenotype expression. The observation function can depend on several loci.

## Hardy–Weinberg: a model to examine, not assume

For a diploid autosomal biallelic locus, let q denote the allele-2 frequency and p=1-q. Independent random union of gametes from a common allele pool gives genotype probabilities p², 2pq, q². In a deterministic population model with ordinary Mendelian inheritance, unchanged allele frequencies across generations additionally require the modeled absence of directional effects from selection, mutation, migration, and drift. Random mating describes the chosen gamete pool; pooling distinct subpopulations is not automatically random mating. Finite sample genotype counts fluctuate around expectations. Historical attribution is to [Hardy (1908)](https://doi.org/10.1126/science.28.706.49); exact-test assumptions and interpretation follow [Wigginton et al. (2005)](https://pmc.ncbi.nlm.nih.gov/articles/PMC1199378/).

For counts (n11,n12,n22), calculate N=n11+n12+n22 and q=(n12+2n22)/(2N). Expectations fitted from those same counts are N(p²,2pq,q²). Comparing these with observations is a goodness-of-fit check, not independent prediction. HWE can hold approximately even where evolutionary processes operate; a non-significant test neither proves equilibrium nor establishes every assumption. A significant test does not identify its cause. Disease-based ascertainment, population structure, and calling error can matter alongside biology. Use an exact test when expected rare-homozygote counts are small; treat a monomorphic sample as uninformative.

The implemented two-sided probability-ordering test conditions on the allele count r=n12+2n22. For feasible heterozygote counts h with the parity of r,

\[
 n_{22}=(r-h)/2,\quad n_{11}=N-h-n_{22},\qquad
 w(h)=\frac{N!\,2^h}{n_{11}!\,h!\,n_{22}!}.
\]

Normalize the weights, then sum probabilities no greater than the observed table's probability. This definition differs from doubling a one-sided tail, taking only the observed-table mass, or using a mid-p correction. Compute and compare log weights to preserve very small probabilities. In particular, an absolute probability tolerance of 10^-12 would incorrectly inflate extreme p-values; the implementation was corrected after this independent audit.

For linked loci, retain parental haplotypes and an explicit recombination fraction. Under the symmetric two-locus teaching model, AB/ab produces AB and ab with probability (1-r)/2 each and Ab and aB with r/2 each. At r=0 only parental haplotypes occur; r=1/2 gives the independent-assortment probabilities. This is a supplied conditional gamete model, not a fitted human recombination map. Recombination varies along chromosomes and physical distance is not identical to recombination frequency. [Brown, Mapping Genomes (2002), Figures 5.15–5.18](https://www.ncbi.nlm.nih.gov/books/NBK21116/) (`brown_linkage_2002`). Population linkage disequilibrium requires joint haplotype observations; single-locus HWE is insufficient evidence for or against it.

## Selected external data: rs334 in 1000 Genomes phase 3

**Frozen files:** `data/observed_genotypes.csv` and `data/observed_genotypes.provenance.json`. **Source:** Ensembl REST, release 116, fetched 2026-09-05. **Primary dataset publication:** [1000 Genomes Project Consortium (2015), Nature 526:68–74, doi:10.1038/nature15393](https://www.nature.com/articles/nature15393). The project combined sequencing and microarray evidence with genotype inference and phasing for 2,504 individuals from 26 populations. These are reference-panel genotype calls, not independently confirmed clinical diagnoses.

Exact API location: [`GET variation/human/rs334?population_genotypes=1;pops=1`](https://rest.ensembl.org/variation/human/rs334?population_genotypes=1;pops=1;content-type=application/json), selecting the `population_genotypes` and `populations` arrays where `population` begins `1000GENOMES:phase_3:`. [Ensembl endpoint documentation](https://rest.ensembl.org/documentation/info/variation_id) describes the aggregate options. The CSV retains only the 26 component populations; overlapping ALL and AFR/AMR/EAS/EUR/SAS rows remain in the provenance for consistency checks.

The allele mapping is GRCh38 positive-strand chr11:5227002 T>A. Genomic alternate A corresponds to HBB coding c.20A>T (p.Glu7Val; also conventionally Glu6Val/E6V), the HbS variant. This strand distinction prevents swapping sickle and reference alleles. [ClinVar VCV000015333.7, Variant Details/HGVS](https://www.ncbi.nlm.nih.gov/clinvar/variation/VCV000015333.7/) (`clinvar_rs334_15333_v7`). Ensembl lists T/A/C/G across its broader rs334 record; only T and A occur in the selected phase-3 aggregate arrays. The exercise is a restricted biallelic contrast.

### Extraction and verification

1. Canonicalize each phased single-locus call by sorting its two genomic alleles. Single-locus HWE does not need phase.
2. Cross-check every reported genotype frequency against its count/N.
3. Cross-check the separate allele-count array: A=2nAA+nTA and T=2nTT+nTA, with total 2N. Also check every reported allele frequency.
4. Complete absent genotype entries with zero only after these checks. A zero means absent in this call set, not biologically impossible. In particular, AA remains representable by the inheritance engine.
5. Verify component sums equal each provided aggregate genotype distribution and total N=2504. Do not count both component and aggregate rows as independent samples.

Frozen totals: T/T=2367, T/A=137, A/A=0; T alleles=4871, A alleles=137, total=5008. The CSV SHA-256 is `2006b6af4722df8ea0d712c3f67d6cb97634a28584a12ed6a03aab82077d9e4a`. The provenance retains the actual relevant API arrays, mapping, retrieval metadata, and checks. It records a response hash, but does not claim the entire raw HTTP response is saved.

| Selected stratum | Observed TT, TA, AA | Derived q(A) | Fitted HWE expected TT, TA, AA | Exact two-sided p |
|---|---|---:|---|---:|
| YRI | 78, 30, 0 | 0.138888889 | 80.083333, 25.833333, 2.083333 | 0.213792186 |
| LWK | 79, 20, 0 | 0.101010101 | 80.010101, 17.979798, 1.010101 | 0.591266935 |
| AFR aggregate | 529, 132, 0 | 0.099848714 | 535.590015, 118.819970, 6.590015 | 0.001590319 |
| ALL aggregate | 2367, 137, 0 | 0.027356230 | 2368.873902, 133.252196, 1.873902 | 0.260973839 |

These displayed p-values are descriptive and unadjusted for multiple testing; strata overlap. They are not four independent confirmatory tests. At a descriptive 0.05 threshold the AFR aggregate departs from the biallelic HWE proportions while the globally pooled table does not. This shows how the choice of sampling unit can conceal a discrepancy. It neither proves selection nor makes the global panel an appropriate panmictic population.

The sum of within-component expected AA counts is 7.600217699, versus 1.873901757 under global pooling; observed AA calls are zero in both views. Algebraically, if groups have weights w_k and allele frequencies q_k, then Σw_k q_k² = q̄² + Var_w(q). The corresponding within-group expected heterozygosity is lower than the global-pool expectation by 2Var_w(q). This identity concerns **model expectations**; it does not force the observed sample to have a heterozygote deficit. Here, observed heterozygosity exceeds both fitted expectations. Do not label these data a pure demonstration of the classical Wahlund deficit.

### Limits and reuse

The population labels describe a diversity sampling design, not biological partitions with universal frequencies. The [project collection principles](https://www.internationalgenome.org/sample_collection_principles/) specify adult sampling and public genomic use. This snapshot cannot estimate newborn disease incidence, ascertain whether missing AA calls arose from selection/survival, sampling or calling, or estimate a selection coefficient. It contains no age-specific trajectories, mating records, independently measured clinical labels, or parent–offspring outcomes. Single-locus aggregate counts also cannot validate linkage disequilibrium or recombination. ABO/Rh phenotype-frequency tables cannot fill that gap by being relabeled genotypes.

The [IGSR reuse policy](https://www.internationalgenome.org/IGSR_disclaimer/) states that original 1000 Genomes data are available without embargo after final publication and should be cited. Third-party rights and EMBL-EBI terms still apply. The saved exercise uses only public aggregate counts and retrieves no individual-level genetic or clinical records. No public upload or publication was performed.

## Reproducible data retrieval

From the workspace root, the following command validates the retained provenance entirely offline and regenerates the CSV bytes in memory, including their SHA-256. It writes nothing:

```powershell
python version2/experiments/fetch_population_data.py
```

To compare a fresh API retrieval against the frozen counts without changing files:

```powershell
python version2/experiments/fetch_population_data.py --live
```

To preserve a new snapshot for review, use `--live --output-dir version2/data/recheck-YYYYMMDD`. That directory must not already exist. The command refuses replacement of the frozen snapshot or any existing output directory. A changed count comparison returns exit status 2 for review. Network access is needed only for `--live`. Retained data can also be verified with `--from-provenance PATH` (alias `--snapshot PATH`); `--compare PATH` chooses the comparison CSV.

Verification completed: the offline snapshot path reproduced the frozen CSV bytes and hash exactly; the `--snapshot` command succeeded; a fresh `--live` read matched every frozen row; six malformed inputs (count, allele count, unexpected genotype, negative count, missing population, wrong strand) were rejected; an existing output-directory request failed without changing the source CSV. These checks assess extraction integrity, not biological predictive performance.

## Rejected clean-table candidate and adversarial source audit

[McAuley et al. (2010), Blood 116:1663–1668, doi:10.1182/blood-2010-01-265249](https://pmc.ncbi.nlm.nih.gov/articles/PMC3073423/) initially appeared suitable because Table 1 contains aggregate HbS classifications in coastal Kenyan cases and community controls. Direct [Europe PMC XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3073423/fullTextXML) verification found:

- Table 1 (XML id T1) reports controls HbAA1252, HbAS218, SCA10, N1480. The narrative reports N1479; its age strata also sum to1479. These cannot both be silently treated as the same denominator.
- The five Table 1 group totals sum to5920 while its total row says5919. SCA group counts sum to121 while that total row says120. Some displayed percentages do not match row denominators.
- Methods explicitly leave HbS/beta0-thalassemia unresolved, so its SCA column cannot be relabeled molecularly confirmed HbSS.

Controls were sampled among children under five in the study area; the disease groups were hospital-ascertained over differing periods. The article's disease-specific HbAS associations motivate context-dependent selection models, but do not supply lifetime reproductive fitness. The article XML states CC BY-NC3.0US. The table was therefore not used as the clean CSV. As a sensitivity check only, control HWE p=0.859938751 for1252/218/10 versus0.859998205 for inferred1251/218/10; neither resolves the source discrepancy or classification ambiguity.

## Independent numerical audit of exact HWE

After replacing an absolute probability tolerance with a log-weight comparison in `genetics/extensions.py`, the implementation was compared with an independent exact-rational oracle. The oracle uses Python integer factorials and `fractions.Fraction`, enumerates all feasible heterozygote counts, orders tables by exact rational weights, and normalizes only at the end. Across all5455 count triples with N=1…30, the maximum absolute difference was6.8834×10^-15. This is a numerical implementation check, independent of the population model's biological applicability.

| Extreme or observed test counts | Exact-rational reference p | Implementation p |
|---|---:|---:|
| 50,0,50 | 1.114224180581451e-30 | 1.1142241805814612e-30 |
| 100,0,100 | 8.795173943032362e-61 | 8.795173943032691e-61 |
| 500,0,500 | 1.3196690976572097e-301 | 1.3196690976571266e-301 |
| 0,100,0 | 1.511390827305580e-29 | 1.5113908273056113e-29 |
| 529,132,0 | 0.001590318956274777 | 0.001590318956275435 |
| 2367,137,0 | 0.2609738386906512 | 0.2609738386903788 |

For the all-heterozygote case, 2^N/binomial(2N,N) is only the observed-table probability. The two-sided exact p-value additionally includes any opposite-tail tables whose probabilities are no larger. For N=100 these differ (1.3999684e-29 versus1.5113908e-29); using the former as a two-sided test oracle is incorrect.

The independent oracle used for this check is reproducible with only the Python standard library:

```python
from fractions import Fraction
from math import factorial

def rational_hwe(n11, n12, n22):
    n, r = n11 + n12 + n22, n12 + 2 * n22
    weights = []
    for h in range(r % 2, min(r, 2 * n - r) + 1, 2):
        b = (r - h) // 2
        a = n - h - b
        w = Fraction(2**h * factorial(n),
                     factorial(a) * factorial(h) * factorial(b))
        weights.append((h, w))
    observed = next(w for h, w in weights if h == n12)
    return sum((w for h, w in weights if w <= observed), Fraction()) / sum(
        (w for h, w in weights), Fraction())
```

The verification ran with the bundled Python3.12 runtime matching the project wheels. The system's unrelated Python3.14 interpreter cannot load the installed cp312 numerical extensions; the fetch script itself needs only the standard library and ran successfully with the system interpreter.
