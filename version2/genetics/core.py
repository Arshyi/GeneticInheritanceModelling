"""Autosomal, unphased diploid catalog; independent segregation conditional on parents.

No probability pruning. Catalog size is distinct from population support.
Full kernels use offspring rows and unordered parental-pair columns.
Linked loci require phased haplotypes (see extensions.py).
"""
from collections import Counter
from functools import lru_cache
from itertools import combinations_with_replacement, product
import math
import numpy as np


class ResourceLimitError(MemoryError):
    """Operation refused without returning an incomplete distribution."""


def positive_product(*values):
    result=math.prod(values)
    if result == 0 and all(v>0 for v in values):
        raise FloatingPointError('Positive population mass underflows float64; log-space population inference is not implemented')
    return result


def distribution(values):
    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or not len(x) or not np.isfinite(x).all() or (x < 0).any():
        raise ValueError('Expected a finite nonnegative probability vector')
    if not np.isclose(x.sum(), 1, atol=1e-12, rtol=0):
        raise ValueError('Probability vector must sum to one')
    return x


@lru_cache(maxsize=4096)
def local_cross(a, b):
    """Canonical allele-index genotypes -> exact dyadic Mendelian probabilities."""
    counts = Counter(tuple(sorted((x, y))) for x in a for y in b)
    return tuple((g, count / 4) for g, count in sorted(counts.items()))


class InheritanceModel:
    def __init__(self, allele_counts, max_bytes=256 * 1024**2):
        self.allele_counts = tuple(allele_counts)
        if not self.allele_counts or any(type(a) is not int or a < 1 for a in self.allele_counts):
            raise ValueError('Each locus needs a positive integer allele count')
        if type(max_bytes) is not int or max_bytes < 1:
            raise ValueError('max_bytes must be positive')
        self.max_bytes = max_bytes
        self.radices = tuple(a * (a + 1) // 2 for a in self.allele_counts)
        self.G = math.prod(self.radices)
        self.U = self.G * (self.G + 1) // 2
        # Bound local catalogs too: do not attempt a giant single-locus allocation.
        self._guard(sum(self.radices) * 160, 'local genotype catalogs')
        self.loci = tuple(tuple(combinations_with_replacement(range(a), 2)) for a in self.allele_counts)
        self.indices = tuple({g: i for i, g in enumerate(states)} for states in self.loci)

    def _guard(self, required, operation):
        if required > self.max_bytes:
            raise ResourceLimitError(f'{operation}: estimated {required:,} bytes exceeds budget {self.max_bytes:,}; use factored queries')

    @property
    def nnz(self):
        t = math.prod(a**4 - a*(a-1)//2 for a in self.allele_counts)
        d = math.prod(a + 3*a*(a-1)//2 for a in self.allele_counts)
        return (t + d) // 2

    def decode(self, code):
        if not isinstance(code, (int, np.integer)) or not 0 <= code < self.G:
            raise ValueError('Genotype code outside catalog')
        out = []
        for radix in reversed(self.radices):
            code, digit = divmod(int(code), radix)
            out.append(digit)
        return tuple(reversed(out))

    def encode(self, state):
        if len(state) != len(self.radices):
            raise ValueError('Incorrect locus count')
        code = 0
        for digit, radix in zip(state, self.radices):
            if not isinstance(digit, (int, np.integer)) or not 0 <= digit < radix:
                raise ValueError('Genotype index outside locus catalog')
            code = code * radix + int(digit)
        return code

    def factors(self, parent1, parent2):
        p, q = self.decode(parent1), self.decode(parent2)
        return tuple(tuple((self.indices[l][g], prob) for g, prob in local_cross(self.loci[l][i], self.loci[l][j]))
                     for l, (i, j) in enumerate(zip(p, q)))

    def log_probability(self, parent1, parent2, child):
        factors, target = self.factors(parent1, parent2), self.decode(child)
        ans = 0.0
        for f, g in zip(factors, target):
            prob = dict(f).get(g, 0.0)
            if prob == 0:
                return -math.inf
            ans += math.log(prob)
        return ans

    def probability(self, parent1, parent2, child):
        logp = self.log_probability(parent1, parent2, child)
        if logp == -math.inf:
            return 0.0
        p = math.exp(logp)
        if p == 0:
            raise FloatingPointError('Positive outcome underflows float64; use log_probability')
        return p

    def iter_cross(self, parent1, parent2):
        """Lazy positive support; caller controls streaming duration and memory."""
        for outcomes in product(*self.factors(parent1, parent2)):
            p = math.prod(x[1] for x in outcomes)
            if p == 0:
                raise FloatingPointError('Positive outcome underflows float64; use iter_log_cross')
            yield self.encode(tuple(x[0] for x in outcomes)), p

    def iter_log_cross(self, parent1, parent2):
        """Stream positive outcomes in log space; no threshold discards rare states."""
        for outcomes in product(*self.factors(parent1, parent2)):
            yield self.encode(tuple(x[0] for x in outcomes)), sum(math.log(x[1]) for x in outcomes)

    def fixed_mate_operator(self, mate_frequencies):
        q = distribution(mate_frequencies)
        if len(q) != self.G:
            raise ValueError('Mate frequency catalog mismatch')
        self._guard(8*self.G*self.G,'fixed-mate operator')
        out = np.zeros((self.G,self.G))
        for i in range(self.G):
            for j in range(self.G):
                if q[j] == 0: continue
                for child,p in self.iter_cross(i,j): out[child,i] += q[j]*p
        return out

    def cross(self, parent1, parent2):
        support = math.prod(len(f) for f in self.factors(parent1, parent2))
        self._guard(support * (128 + 8 * len(self.radices)), 'cross dictionary')
        try:
            return dict(self.iter_cross(parent1, parent2))
        except MemoryError as exc:
            local_cross.cache_clear()
            raise ResourceLimitError('Cross allocation failed; no partial result returned') from exc

    def iter_columns(self):
        for i in range(self.G):
            for j in range(i, self.G):
                yield i, j, self.iter_cross(i, j)

    def kernel(self, kind='csr'):
        """Build all columns directly, never allocate impossible transmission values.

        Dense is an explicit comparison baseline. Budget estimates allow for temporary
        CSC -> CSR conversion but are not an OS-wide memory availability guarantee.
        """
        if kind not in ('dense', 'csr', 'hash'):
            raise ValueError('kind must be dense, csr, or hash')
        required = {'dense': 8*self.G*self.U,
                    'csr': 2*(16*self.nnz + 8*(self.U+self.G+2)),
                    'hash': 220*self.nnz + 100*self.U}[kind]
        self._guard(required, kind + ' full kernel')
        try:
            if kind == 'dense':
                out = np.zeros((self.G, self.U), dtype=np.float64)
                for col, (_, _, entries) in enumerate(self.iter_columns()):
                    for row, p in entries:
                        out[row, col] = p
                return out
            if kind == 'hash':
                return {(i,j): dict(entries) for i,j,entries in self.iter_columns()}
            from scipy.sparse import csc_matrix
            data = np.empty(self.nnz, dtype=np.float64)
            indices = np.empty(self.nnz, dtype=np.int64)
            indptr = np.empty(self.U+1, dtype=np.int64)
            pos = 0
            for col, (_, _, entries) in enumerate(self.iter_columns()):
                indptr[col] = pos
                for row, p in entries:
                    indices[pos], data[pos] = row, p
                    pos += 1
            indptr[-1] = pos
            assert pos == self.nnz
            return csc_matrix((data, indices, indptr), shape=(self.G,self.U)).tocsr()
        except MemoryError as exc:
            local_cross.cache_clear()
            raise ResourceLimitError('Kernel allocation failed; cached local crosses cleared; use lazy factors') from exc

    def pair_weights(self, frequencies):
        x = distribution(frequencies)
        if len(x) != self.G:
            raise ValueError('Frequency catalog mismatch')
        self._guard(8*self.U, 'pair probabilities')
        return np.fromiter((positive_product(x[i],x[j],(1 if i==j else 2))
                            for i in range(self.G) for j in range(i,self.G)), dtype=float, count=self.U)

    def next_generation(self, frequencies, kernel=None):
        x = distribution(frequencies)
        if len(x) != self.G:
            raise ValueError('Frequency catalog mismatch')
        if kernel is not None and not isinstance(kernel, dict):
            if kernel.shape != (self.G,self.U): raise ValueError('Kernel shape mismatch')
            weights=self.pair_weights(x)
            result=np.asarray(kernel @ weights).reshape(-1)
            reachable=np.asarray(kernel @ (weights>0)).reshape(-1)>0
            if (result[reachable]==0).any():
                raise FloatingPointError('Reachable population state underflows float64')
            return result
        self._guard(8*self.G, 'population result')
        out = np.zeros(self.G)
        for i in range(self.G):
            if x[i] == 0: continue
            for j in range(i,self.G):
                weight = positive_product(x[i],x[j],(1 if i==j else 2))
                if weight == 0: continue
                entries = kernel[(i,j)].items() if kernel is not None else self.iter_cross(i,j)
                for g,p in entries:
                    out[g] += positive_product(weight,p)
        return out


def hwe(allele_frequencies):
    p = distribution(allele_frequencies)
    return np.array([positive_product(p[i],p[j],(1 if i==j else 2)) for i in range(len(p)) for j in range(i,len(p))])


def select(frequencies, fitness):
    x = distribution(frequencies)
    w = np.asarray(fitness, dtype=float)
    if w.shape != x.shape or not np.isfinite(w).all() or (w < 0).any():
        raise ValueError('Invalid nonnegative fitness')
    positive = (x>0)&(w>0)
    if not positive.any():
        raise ValueError('No reproductive mass remains')
    logs = np.log(x[positive])+np.log(w[positive])
    relative = np.exp(logs-logs.max())
    result = np.zeros_like(x)
    result[positive] = relative/relative.sum()
    if (result[positive]==0).any():
        raise FloatingPointError('Positive selected frequency underflows float64')
    return result


def migrate(frequencies, immigrants, fraction):
    x, y = distribution(frequencies), distribution(immigrants)
    if x.shape != y.shape or not math.isfinite(fraction) or not 0 <= fraction <= 1:
        raise ValueError('Invalid migration mixture')
    return (1-fraction)*x + fraction*y


def mutation(allele_frequencies, transition):
    """Row-stochastic allele mutation applied to gametes before fertilization."""
    p = distribution(allele_frequencies)
    m = np.asarray(transition, dtype=float)
    if m.shape != (len(p),len(p)) or not np.isfinite(m).all() or (m < 0).any() or not np.allclose(m.sum(axis=1),1,atol=1e-12,rtol=0):
        raise ValueError('Mutation transition must be finite, nonnegative and row-stochastic')
    return p @ m
