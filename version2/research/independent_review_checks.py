"""Checks behind research/independent_review_2026-11-16.md.

Standard library only, so it runs without the pinned scientific stack:

    python research/independent_review_checks.py        (from version2/)

Nothing in the package is imported or modified. Section 2 reads the retained
results/pgs_catalog.json; every other section is an independent re-computation.

Sections
  1  The counting theorems, re-checked by brute-force enumeration
  2  Polygenic scores: the reported worst-case bound against the error actually stored
  3  One generation of random mating without the parent-pair kernel
  4  A correlated two-locus start does not reach its fixed point in one generation
  5  The additive score distribution by Fourier inversion with exact weights
"""

import cmath
import json
import math
import random
import time
from itertools import combinations_with_replacement, product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


# --------------------------------------------------------------------------- 1
def brute_force_counts(allele_counts):
    """G, U and supported transitions by direct enumeration (no formula used)."""
    loci = [list(combinations_with_replacement(range(a), 2)) for a in allele_counts]
    genotypes = list(product(*loci))
    support = 0
    for i, g in enumerate(genotypes):
        for h in genotypes[i:]:
            size = 1
            for gl, hl in zip(g, h):
                size *= len({tuple(sorted((x, y))) for x in gl for y in hl})
            support += size
    G = len(genotypes)
    return G, G * (G + 1) // 2, support


def formula_counts(allele_counts):
    G = math.prod(a * (a + 1) // 2 for a in allele_counts)
    t = math.prod(a ** 4 - a * (a - 1) // 2 for a in allele_counts)
    d = math.prod(a + 3 * a * (a - 1) // 2 for a in allele_counts)
    return G, G * (G + 1) // 2, (t + d) // 2


def section1():
    print("== 1. Counting theorems against brute-force enumeration")
    print(f"   {'allele counts':>18s} {'G':>6s} {'U':>9s} {'supported (brute force)':>24s} {'supported (formula)':>20s} {'match':>6s}")
    systems = [(2,), (3,), (4,), (5,), (3, 2), (2, 2), (3, 3), (4, 2), (2, 2, 2), (3, 2, 2), (2, 2, 2, 2), (2, 2, 2, 2, 2)]
    for system in systems:
        brute, formula = brute_force_counts(system), formula_counts(system)
        print(f"   {str(system):>18s} {brute[0]:6d} {brute[1]:9,d} {brute[2]:24,d} {formula[2]:20,d} {str(brute == formula):>6s}")


# --------------------------------------------------------------------------- 2
def section2():
    print("\n== 2. Polygenic scores: reported worst-case bound against the stored moments (results/pgs_catalog.json)")
    data = json.loads((ROOT / "results" / "pgs_catalog.json").read_text(encoding="utf-8"))
    ladder = data["ladder"]
    print(f"   scores listed for the trait: {data['catalogue_summary']['scores_for_this_trait']};"
          f" scoring files downloaded and run: {len(ladder)}")
    print(f"   {'score':10s} {'n':>7s} {'frequency':>10s} {'bins':>6s} {'bound (SD)':>11s} {'actual mean error (SD)':>23s}"
          f" {'DP SD / true SD':>16s} {'step / mean|w|':>15s} {'weights rounded to zero':>24s}")
    for rec in ladder:
        s = rec["score_dp"]
        sd = s["sd_analytic"]
        ratio = math.sqrt(max(s["variance_from_dp"], 0.0)) / sd
        step = s["discretisation_step"] / rec["weight_summary"]["mean_abs"]
        source = "published" if "published" in rec["effect_allele_frequency_source"] else "declared"
        integer_sum = (s["bins_predicted"] - 1) // 2          # bins = 2 * sum(integer weights) + 1
        zero_share = max(0.0, 1 - integer_sum / rec["n_used"])  # each surviving weight is at least 1
        zero_text = f"at least {zero_share:.0%}" if zero_share > 0 else "not bounded"
        print(f"   {rec['id']:10s} {rec['n_used']:7d} {source:>10s} {s['bins']:6d} {s['worst_case_error_in_sd']:11.3f}"
              f" {s['mean_absolute_discrepancy'] / sd:23.4f} {ratio:16.4f} {step:15.3f} {zero_text:>24s}")
    print("   accuracy-controlled runs (target: bound at 1% of SD):")
    for rec in ladder:
        acc = rec["accuracy_controlled"]
        line = f"     {rec['id']} n={rec['n_used']:6d}: bins the bound demands {acc['bins_required']:14,d}  {acc['status']}"
        if acc["status"] == "measured":
            line += (f"  -> actual mean error {acc['mean_absolute_discrepancy'] / rec['score_dp']['sd_analytic']:.1e} SD;"
                     f" the fixed 40,001-bin run already had {rec['score_dp']['mean_absolute_discrepancy'] / rec['score_dp']['sd_analytic']:.1e} SD")
        print(line)
    study = data["controlled_scaling"]
    print(f"   controlled study: fitted exponent {study['fitted_bin_exponent']:.4f}, fitted to 'bins_required' for n ="
          f" {[row['n'] for row in study['rows']]}")
    print("   'bins_required' is the bound formula 1 + 2 n sum|w| / (0.01 SD) evaluated for each n; measured rows:"
          f" {[row['n'] for row in study['rows'] if row['status'] == 'measured']}")


# --------------------------------------------------------------------------- 3
LOCAL = {}
for _g in range(3):
    for _h in range(3):
        _a, _b = _g / 2.0, _h / 2.0  # chance each parent transmits the counted allele
        _dist = ((0, (1 - _a) * (1 - _b)), (1, _a * (1 - _b) + (1 - _a) * _b), (2, _a * _b))
        LOCAL[(_g, _h)] = tuple((c, p) for c, p in _dist if p > 0)


def digits(code, n):
    out = [0] * n
    for k in range(n - 1, -1, -1):
        code, out[k] = divmod(code, 3)
    return out


def kernel_update(x, n):
    """The kernel route: every unordered parental pair, every supported child."""
    G = 3 ** n
    dig = [digits(c, n) for c in range(G)]
    out = [0.0] * G
    visited = 0
    for i in range(G):
        if x[i] == 0.0:
            continue
        for j in range(i, G):
            w = x[i] * x[j] * (1.0 if i == j else 2.0)
            if w == 0.0:
                continue
            for combo in product(*[LOCAL[(dig[i][l], dig[j][l])] for l in range(n)]):
                code, p = 0, w
                for c, pc in combo:
                    code = code * 3 + c
                    p *= pc
                out[code] += p
                visited += 1
    return out, visited


def apply_axis(arr, shape, axis, mat, zero=0.0):
    """Apply mat (rows = new size, cols = old size) along one axis of a flat row-major array."""
    old, new = shape[axis], len(mat)
    stride = math.prod(shape[axis + 1:])
    block = stride * old
    nblocks = len(arr) // block
    out = [zero] * (nblocks * new * stride)
    ops = 0
    for b in range(nblocks):
        bo, bn = b * block, b * new * stride
        for r in range(new):
            for c in range(old):
                m = mat[r][c]
                if m == 0:
                    continue
                so, sn = bo + c * stride, bn + r * stride
                for k in range(stride):
                    out[sn + k] += m * arr[so + k]
                ops += stride
    return out, shape[:axis] + (new,) + shape[axis + 1:], ops


T_GAMETE = [[1.0, 0.5, 0.0], [0.0, 0.5, 1.0]]  # genotype dosage -> transmitted allele


def gamete_frequencies(x, n):
    arr, shape, ops = list(x), (3,) * n, 0
    for axis in range(n):
        arr, shape, o = apply_axis(arr, shape, axis, T_GAMETE)
        ops += o
    return arr, ops  # indexed by haplotype, one bit per locus


def embed3(h, n):
    """Haplotype bits read as base-3 digits, so that adding two of them gives the dosage code."""
    code = 0
    for k in range(n - 1, -1, -1):
        code = code * 3 + ((h >> k) & 1)
    return code


def gamete_update_direct(x, n):
    gam, ops = gamete_frequencies(x, n)
    H = 2 ** n
    e = [embed3(h, n) for h in range(H)]
    out = [0.0] * (3 ** n)
    for h in range(H):
        if gam[h] == 0.0:
            continue
        for k in range(H):
            out[e[h] + e[k]] += gam[h] * gam[k]
    return out, ops + H * H


def gamete_update_fourier(x, n):
    gam, ops = gamete_frequencies(x, n)
    arr = [0j] * (3 ** n)
    for h in range(2 ** n):
        arr[embed3(h, n)] = complex(gam[h])
    w = cmath.exp(-2j * math.pi / 3)
    forward = [[w ** (r * c) for c in range(3)] for r in range(3)]
    inverse = [[(w ** (-(r * c))) / 3 for c in range(3)] for r in range(3)]
    shape = (3,) * n
    for axis in range(n):
        arr, shape, o = apply_axis(arr, shape, axis, forward, 0j)
        ops += o
    arr = [v * v for v in arr]
    ops += len(arr)
    for axis in range(n):
        arr, shape, o = apply_axis(arr, shape, axis, inverse, 0j)
        ops += o
    return [v.real for v in arr], ops


def section3():
    print("\n== 3. One generation of random mating: parent-pair kernel against gamete frequencies")
    print("   input: a random joint genotype distribution (loci correlated); unlinked biallelic loci")
    print(f"   {'loci':>4s} {'genotypes':>9s} {'kernel transitions':>19s} {'kernel s':>9s} {'gamete ops':>11s} {'gamete s':>9s}"
          f" {'fourier ops':>12s} {'fourier s':>10s} {'largest difference':>19s}")
    rng = random.Random(20261116)
    for n in range(1, 11):
        G = 3 ** n
        x = [rng.random() ** 3 for _ in range(G)]
        total = sum(x)
        x = [v / total for v in x]
        t0 = time.perf_counter()
        direct, direct_ops = gamete_update_direct(x, n)
        direct_s = time.perf_counter() - t0
        diff = 0.0
        kernel_col = f"{(15 ** n + 5 ** n) // 2:>19,d} {'not run':>9s}"
        if n <= 5:
            t0 = time.perf_counter()
            ker, visited = kernel_update(x, n)
            kernel_s = time.perf_counter() - t0
            assert visited == (15 ** n + 5 ** n) // 2
            diff = max(abs(a - b) for a, b in zip(ker, direct))
            kernel_col = f"{visited:>19,d} {kernel_s:9.3f}"
        fourier_col = f"{'not run':>12s} {'':>10s}"
        if n <= 9:
            t0 = time.perf_counter()
            four, four_ops = gamete_update_fourier(x, n)
            four_s = time.perf_counter() - t0
            diff = max(diff, max(abs(a - b) for a, b in zip(four, direct)))
            fourier_col = f"{four_ops:>12,d} {four_s:10.3f}"
        print(f"   {n:4d} {G:9,d} {kernel_col} {direct_ops:11,d} {direct_s:9.3f} {fourier_col} {diff:19.2e}")


# --------------------------------------------------------------------------- 4
def section4():
    print("\n== 4. Two unlinked loci, population starts as half AABB and half aabb")
    x = [0.0] * 9
    x[0], x[8] = 0.5, 0.5
    print(f"   {'generation':>10s} {'AB gamete':>10s} {'Ab gamete':>10s} {'association D':>14s} {'AaBb share':>11s} {'largest change at next step':>28s}")
    for generation in range(9):
        gam, _ = gamete_frequencies(x, 2)
        nxt, _ = gamete_update_direct(x, 2)
        move = max(abs(a - b) for a, b in zip(nxt, x))
        print(f"   {generation:10d} {gam[0]:10.6f} {gam[1]:10.6f} {gam[0] * gam[3] - gam[1] * gam[2]:14.6f} {x[4]:11.6f} {move:28.6f}")
        x = nxt


# --------------------------------------------------------------------------- 5
def cf_coefficients(weights, freqs, M=256, L=24.0):
    """Fourier coefficients of the standardised score's density on a window of L standard deviations."""
    mu = sum(2 * p * w for w, p in zip(weights, freqs))
    sd = math.sqrt(sum(2 * p * (1 - p) * w * w for w, p in zip(weights, freqs)))
    half, du = M // 2, 2 * math.pi / L
    scaled = [w / sd for w in weights]
    phi = [0j] * (M + 1)  # index half + k  <->  u = k * du
    for k in range(half + 1):
        u = du * k
        value = cmath.exp(-1j * u * mu / sd)
        for w_s, p in zip(scaled, freqs):
            term = (1 - p) + p * cmath.exp(1j * u * w_s)   # one allele copy; a diploid dosage squares it
            value *= term * term
            if value == 0:
                break
        phi[half + k] = value
        phi[half - k] = value.conjugate()
    return {"mu": mu, "sd": sd, "phi": phi, "du": du, "L": L, "M": M, "half": half, "evaluations": half + 1}


def cf_cdf(d, z):
    L, du, half, phi = d["L"], d["du"], d["half"], d["phi"]
    total = (z + L / 2) / L * phi[half].real
    for k in range(1, half):
        u = du * k
        total += 2 * (phi[half + k] * (cmath.exp(-1j * u * z) - cmath.exp(1j * u * L / 2)) / (-1j * u * L)).real
    return total


def cf_density(d):
    L, du, half, M, phi = d["L"], d["du"], d["half"], d["M"], d["phi"]
    zs = [-L / 2 + m * L / M for m in range(M)]
    dens = [(phi[half].real + sum(2 * (phi[half + k] * cmath.exp(-1j * du * k * z)).real for k in range(1, half))) / L for z in zs]
    return zs, dens


def exact_pmf(int_weights, freqs):
    lo, pmf = 0, [1.0]
    for w, p in zip(int_weights, freqs):
        probs = ((1 - p) ** 2, 2 * p * (1 - p), p * p)
        shift0 = min(0, 2 * w)
        new = [0.0] * (len(pmf) + 2 * abs(w))
        for dosage, pr in enumerate(probs):
            s = dosage * w - shift0
            for k, v in enumerate(pmf):
                if v:
                    new[s + k] += pr * v
        lo += shift0
        pmf = new
    return lo, pmf


def bisect(fn, level, lo=-12.0, hi=12.0):
    for _ in range(60):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if fn(mid) < level else (lo, mid)
    return (lo + hi) / 2


def section5():
    print("\n== 5a. Fourier inversion against an exact dynamic programme (40 loci, integer weights, mixed frequencies)")
    rng = random.Random(7)
    iw = [rng.choice([-1, 1]) * rng.randint(20, 300) for _ in range(40)]
    fr = [rng.uniform(0.05, 0.95) for _ in range(40)]
    lo, pmf = exact_pmf(iw, fr)
    d = cf_coefficients([float(w) for w in iw], fr)
    mu, sd = d["mu"], d["sd"]
    cum, acc = [], 0.0
    for v in pmf:
        acc += v
        cum.append(acc)

    def exact_cdf(z):
        s = math.floor(mu + sd * z) - lo
        return 0.0 if s < 0 else cum[min(s, len(cum) - 1)]

    worst = max(abs(cf_cdf(d, k / 10) - exact_cdf(k / 10)) for k in range(-60, 61))
    print(f"   exact support {len(pmf):,} integer scores; SD {sd:.1f} score units; Fourier evaluations {d['evaluations']}")
    print(f"   largest CDF difference for z in [-6, 6]: {worst:.2e}   (largest single atom of the exact pmf: {max(pmf):.2e})")
    print(f"   {'level':>7s} {'exact quantile (SD)':>20s} {'Fourier quantile (SD)':>22s} {'difference (SD)':>16s}")
    for level in (0.001, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 0.999):
        zq_exact = (lo + next(i for i, c in enumerate(cum) if c >= level) - mu) / sd
        zq_cf = bisect(lambda z: cf_cdf(d, z), level)
        print(f"   {level:7.3f} {zq_exact:20.4f} {zq_cf:22.4f} {zq_cf - zq_exact:16.4f}")

    print("\n== 5b. The paper's synthetic design: weights ~ Normal(0, 0.05), effect-allele frequency 0.30, exact weights")
    print(f"   {'loci':>6s} {'bins the bound demands':>23s} {'in the paper':>13s} {'Fourier evals':>14s} {'seconds':>8s}"
          f" {'mass':>12s} {'mean':>10s} {'variance':>10s} {'skewness got / exact':>24s} {'excess kurtosis got / exact':>30s}")
    for n, status in ((400, "measured"), (3200, "refused"), (25600, "refused")):
        rng = random.Random(1000 + n)
        w = [rng.gauss(0.0, 0.05) for _ in range(n)]
        p, q = 0.30, 0.70
        t0 = time.perf_counter()
        d = cf_coefficients(w, [p] * n)
        zs, dens = cf_density(d)
        seconds = time.perf_counter() - t0
        dz = d["L"] / d["M"]
        mass = sum(dens) * dz
        mean = sum(z * f for z, f in zip(zs, dens)) * dz
        var = sum(z * z * f for z, f in zip(zs, dens)) * dz - mean ** 2
        m3 = sum((z - mean) ** 3 * f for z, f in zip(zs, dens)) * dz
        m4 = sum((z - mean) ** 4 * f for z, f in zip(zs, dens)) * dz
        sd = d["sd"]
        k3 = sum(v ** 3 for v in w) * 2 * p * q * (q - p) / sd ** 3
        k4 = sum(v ** 4 for v in w) * 2 * p * q * (1 - 6 * p * q) / sd ** 4
        need = math.ceil(1 + 2 * n * sum(abs(v) for v in w) / (0.01 * sd))
        print(f"   {n:6d} {need:23,d} {status:>13s} {d['evaluations']:14d} {seconds:8.2f} {mass:12.9f} {mean:10.2e} {var:10.7f}"
              f" {m3:11.6f} / {k3:10.6f} {m4 - 3:15.6f} / {k4:12.6f}")
    print("   (weights are drawn with Python's generator, so they differ from the paper's NumPy draw; the design is the same)")


if __name__ == "__main__":
    print("Independent review checks, GeneticInheritanceModelling version2")
    section1()
    section2()
    section3()
    section4()
    section5()
