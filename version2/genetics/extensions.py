"""Bounded linked-locus and synthetic quantitative-trait extensions."""
from collections import defaultdict
from itertools import product
import math
import numpy as np
from .core import distribution, ResourceLimitError


def linked_gametes(haplotypes, recombination):
    """Two phased biallelic loci; exactly one effective recombination fraction r."""
    if not 0 <= recombination <= .5 or not math.isfinite(recombination):
        raise ValueError('Recombination fraction must be in [0, 0.5]')
    if len(haplotypes) != 2 or any(len(h) != 2 or any(a not in (0,1) for a in h) for h in haplotypes):
        raise ValueError('Need two phased haplotypes of two biallelic loci')
    a,b = map(tuple,haplotypes)
    out = defaultdict(float)
    for h,p in [(a,(1-recombination)/2),(b,(1-recombination)/2),((a[0],b[1]),recombination/2),((b[0],a[1]),recombination/2)]:
        if p: out[h] += p
    return dict(out)


def linked_cross(parent1, parent2, recombination):
    out = defaultdict(float)
    for a,p in linked_gametes(parent1,recombination).items():
        for b,q in linked_gametes(parent2,recombination).items():
            out[tuple(sorted((a,b)))] += p*q
    return dict(out)


def additive_pmf(dosage_probabilities, weights=None, max_bins=1_000_001):
    """Exact convolution for independent dosage factors and nonnegative integer effects.

    Effects are score units, not height coefficients. Distinct sums merge by DP;
    support may still grow excessively for widely spaced effects, so bound it.
    """
    if type(max_bins) is not int or max_bins < 1:
        raise ValueError('max_bins must be a positive integer')
    factors = [distribution(p) for p in dosage_probabilities]
    weights = [1]*len(factors) if weights is None else list(weights)
    if len(weights) != len(factors) or any(type(w) is not int or w<0 for w in weights) or any(len(p)!=3 for p in factors):
        raise ValueError('Need three dosages per locus and nonnegative integer effects')
    bins = 2*sum(weights)+1
    if bins > max_bins:
        raise ResourceLimitError('Additive DP exceeds bin budget')
    out = np.ones(1)
    for p,w in zip(factors,weights):
        if w == 0: continue
        updated = np.zeros(len(out)+2*w)
        for dosage,prob in enumerate(p):
            updated[dosage*w:dosage*w+len(out)] += prob*out
        out = updated
    return out


def mixture_cdf(values, score_pmf, residual_sd):
    from scipy.special import ndtr
    p = distribution(score_pmf)
    if not math.isfinite(residual_sd) or residual_sd <= 0:
        raise ValueError('Positive residual standard deviation required')
    values = np.asarray(values)
    return np.array([np.dot(p,ndtr((v-np.arange(len(p)))/residual_sd)) for v in values.flat]).reshape(values.shape)


def hwe_exact(counts):
    """Conditional two-sided probability-ordering test at a biallelic locus.

    Enumerates feasible heterozygote counts given allele count using log weights
    proportional to N! 2^h / (nAA! h! naa!); no chi-square approximation.
    """
    from scipy.special import gammaln, logsumexp
    if len(counts)!=3 or any(type(v) is not int or v<0 for v in counts) or sum(counts)==0:
        raise ValueError('Need three nonnegative integer genotype counts')
    aa,het,bb = counts
    n, copies = sum(counts), 2*bb+het
    hs = np.arange(copies%2, min(copies,2*n-copies)+1,2)
    bs = (copies-hs)//2
    aas = n-hs-bs
    logs = hs*np.log(2)-gammaln(aas+1)-gammaln(hs+1)-gammaln(bs+1)
    probs = np.exp(logs-logsumexp(logs))
    observed_log = logs[np.where(hs==het)[0][0]]
    # Compare relative likelihoods in log space, not with an absolute probability
    # tolerance that would swamp a very small exact tail probability.
    return float(min(1,probs[logs<=observed_log+1e-10].sum()))


def linkage_equilibrium_update(local_frequencies, allele_counts):
    """Return local offspring factors for an explicitly linkage-equilibrium population.

    No selection/mutation/migration. This input cannot express interlocus correlation.
    It must never replace a correlated joint population with its marginals silently.
    """
    from .core import hwe
    from itertools import combinations_with_replacement
    if len(local_frequencies) != len(allele_counts):
        raise ValueError('Locus count mismatch')
    outputs=[]
    for x,a in zip(local_frequencies,allele_counts):
        x=distribution(x)
        states=list(combinations_with_replacement(range(a),2))
        if len(x)!=len(states): raise ValueError('Genotype count mismatch')
        p=np.zeros(a)
        for weight,g in zip(x,states):
            for allele in g: p[allele] += weight/2
        outputs.append(hwe(p))
    return tuple(outputs)
