"""The four regimes of Section 15.2, tested on real published polygenic scores.

Coronary artery disease is the useful case because the PGS Catalog holds many
published scores for the *same* phenotype across four orders of magnitude of
variant count. That makes the comparison internal: nothing changes between the
smallest and largest model except n.

Two things are established here. The complete kernel is impossible for every one
of them, including the smallest, by a counting argument rather than by a failed
run. And the additive score dynamic programme handles all of them, at a cost that
is measured rather than asserted, once the real-valued effect weights are
discretised - which introduces an error this module bounds explicitly.

Scoring files are not redistributed. They are fetched from the PGS Catalog at run
time, and only derived quantities and digests are retained, with the provenance
needed to fetch them again.
"""
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import time
import urllib.request

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from genetics.extensions import additive_pmf

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / 'results'
FIGURES = ROOT / 'figures'
API = 'https://www.pgscatalog.org/rest'
AGENT = {'User-Agent': 'GeneticInheritanceModelling/2.0 (arshyiamehran@gmail.com)'}
TRAIT = 'MONDO_0005010'          # coronary artery disorder
COLORS = ['#136f79', '#d18a25', '#293c61', '#b44949']

# One phenotype, four orders of magnitude. The first two publish effect-allele
# frequencies, so for those the population score distribution is real rather than
# assumed; the rest are run with a declared frequency to measure scaling only.
LADDER = ['PGS000010', 'PGS012581', 'PGS002262', 'PGS002775',
          'PGS004196', 'PGS004197', 'PGS000012', 'PGS000337']

DECLARED_FREQUENCY = 0.30        # used only where the score publishes none
TARGET_BINS = 40_001             # fixed-budget run: bins held constant as n grows
TARGET_ERROR_IN_SD = 0.01        # accuracy-controlled run: bound error at 1% of the score SD
MAX_BINS = 4_000_000             # refuse beyond this rather than run for hours
MAX_WORK = 2_000_000_000         # n * bins ceiling, same reason
BUDGET_BYTES = 256 * 1024 ** 2   # the same materialisation budget as Section 9.3


def get_json(url):
    with urllib.request.urlopen(urllib.request.Request(url, headers=AGENT), timeout=120) as r:
        return json.load(r)


def fetch_score(pgs_id):
    """Metadata plus effect weights and, where published, effect-allele frequencies."""
    meta = get_json(f'{API}/score/{pgs_id}')
    url = meta['ftp_scoring_file']
    with urllib.request.urlopen(urllib.request.Request(url, headers=AGENT), timeout=300) as r:
        raw = r.read()
    text = gzip.decompress(raw).decode('utf-8')
    body = [line for line in text.splitlines() if not line.startswith('#')]
    columns = body[0].split('\t')
    wi = columns.index('effect_weight')
    ai = columns.index('allelefrequency_effect') if 'allelefrequency_effect' in columns else None

    weights, freqs = [], []
    for line in body[1:]:
        parts = line.split('\t')
        try:
            w = float(parts[wi])
        except (ValueError, IndexError):
            continue
        if not math.isfinite(w):
            continue
        if ai is not None:
            try:
                f = float(parts[ai])
            except (ValueError, IndexError):
                f = float('nan')
            if not (math.isfinite(f) and 0.0 < f < 1.0):
                continue
            freqs.append(f)
        weights.append(w)
    if ai is not None and len(freqs) != len(weights):
        raise ValueError(f'{pgs_id}: frequency and weight counts disagree')

    pub = meta['publication']
    return {
        'id': meta['id'], 'name': meta['name'],
        'reported_variants': meta['variants_number'],
        'usable_variants': len(weights),
        'method': meta['method_name'], 'weight_type': meta['weight_type'],
        'genome_build': meta['variants_genomebuild'],
        'publication': {'firstauthor': pub['firstauthor'], 'journal': pub['journal'],
                        'date': pub['date_publication'], 'doi': pub['doi'], 'pmid': pub['PMID']},
        'scoring_file_url': url,
        'scoring_file_sha256': hashlib.sha256(raw).hexdigest(),
        'scoring_file_bytes': len(raw),
        'publishes_effect_allele_frequency': ai is not None,
        'licence_note': meta['license'][:300],
    }, np.asarray(weights, dtype=float), (np.asarray(freqs, dtype=float) if ai is not None else None)


def kernel_costs(n):
    """What the complete-kernel regime would demand. A counting argument, not a run."""
    log10_G = n * math.log10(3)
    log10_U = 2 * log10_G - math.log10(2) if n > 1 else math.log10(6)
    log10_nnz = n * math.log10(15) - math.log10(2)
    log10_dense_bytes = log10_G + log10_U + math.log10(8)
    return {
        'log10_genotypes': log10_G,
        'log10_unordered_pairs': log10_U,
        'log10_supported_transitions': log10_nnz,
        'log10_dense_payload_bytes': log10_dense_bytes,
        'exceeds_budget': log10_dense_bytes > math.log10(BUDGET_BYTES),
        'budget_bytes': BUDGET_BYTES,
    }


def discretise(weights, target_bins):
    """Round real effect weights onto an integer grid, and bound the induced error.

    Each locus contributes weight w and dosage 0, 1 or 2, so the score spans
    2*sum|w|. Choosing a step delta that divides that span into target_bins-1
    intervals fixes both the array length and the worst-case error: rounding each
    weight moves it by at most delta/2, and a dosage of at most 2 doubles that, so
    the total displacement is bounded by n*delta.
    """
    span = 2.0 * float(np.abs(weights).sum())
    delta = span / (target_bins - 1)
    integer = np.rint(np.abs(weights) / delta).astype(np.int64)
    integer = np.maximum(integer, 0)
    worst_case = float(len(weights) * delta)
    return integer, delta, worst_case


def dosage_matrix(freqs, signs):
    """Hardy-Weinberg dosage probabilities, with negative weights folded in.

    A locus with weight -w behaves like the same locus with weight +w and the
    dosage distribution reversed, plus a constant offset of 2w. Folding the sign
    into the dosage vector lets one non-negative-weight dynamic programme handle
    both, which is what the implemented routine accepts.
    """
    p = np.asarray(freqs, dtype=float)
    table = np.column_stack([(1 - p) ** 2, 2 * p * (1 - p), p ** 2])
    flip = signs < 0
    table[flip] = table[flip][:, ::-1]
    return table


def run_dp(weights, freqs, signs, bins):
    """One dynamic-programme run at a given bin budget. Returns timing and moments."""
    integer, delta, worst_case = discretise(weights, bins)
    table = dosage_matrix(freqs, signs)
    start = time.perf_counter()
    pmf = np.asarray(additive_pmf(table.tolist(), integer.tolist(), max_bins=MAX_BINS + 16),
                     dtype=float)
    seconds = time.perf_counter() - start
    support = np.arange(len(pmf))
    offset = float(np.sum(2 * weights[signs < 0]))
    mean_dp = float(np.sum(support * pmf)) * delta + offset
    return {'requested_bins': int(bins), 'bins': int(len(pmf)), 'seconds': seconds,
            'discretisation_step': delta, 'worst_case_absolute_error': worst_case,
            'mean_from_dp': mean_dp, 'pmf_payload_bytes': int(pmf.nbytes)}


def required_bins(weights, sd, target_in_sd):
    """Bins needed to hold the worst-case discretisation error at target_in_sd * SD.

    The bound is n*delta with delta = 2*sum|w| / (B-1), so requiring
    n*delta <= target*SD gives B >= 1 + 2*n*sum|w| / (target*SD).
    """
    n = len(weights)
    if sd <= 0:
        return None
    return int(math.ceil(1 + 2 * n * float(np.abs(weights).sum()) / (target_in_sd * sd)))


def analyse(meta, weights, freqs):
    n = len(weights)
    published = freqs is not None
    if not published:
        freqs = np.full(n, DECLARED_FREQUENCY)

    integer, delta, worst_case = discretise(weights, TARGET_BINS)
    signs = np.sign(weights)
    table = dosage_matrix(freqs, signs)

    # Analytic moments of the true (undiscretised) score, for an independent check.
    mean_true = float(np.sum(2 * freqs * weights))
    var_true = float(np.sum(2 * freqs * (1 - freqs) * weights ** 2))

    bins_predicted = int(2 * integer.sum() + 1)
    start = time.perf_counter()
    pmf = additive_pmf(table.tolist(), integer.tolist(), max_bins=8_000_000)
    seconds = time.perf_counter() - start

    pmf = np.asarray(pmf, dtype=float)
    support = np.arange(len(pmf))
    offset = float(np.sum(2 * weights[signs < 0]))     # constant from the sign fold
    mean_dp = float(np.sum(support * pmf)) * delta + offset
    var_dp = float(np.sum((support * delta) ** 2 * pmf) - (np.sum(support * delta * pmf)) ** 2)
    sd = math.sqrt(var_true) if var_true > 0 else 0.0

    # Second run: hold the error bound at a fixed fraction of the score's own SD,
    # and let the bin count be whatever that demands.
    need = required_bins(weights, sd, TARGET_ERROR_IN_SD)
    accuracy = {'target_error_in_sd': TARGET_ERROR_IN_SD, 'bins_required': need}
    if need is None:
        accuracy['status'] = 'undefined; score has zero variance under these assumptions'
    elif need > MAX_BINS or n * need > MAX_WORK:
        accuracy['status'] = 'refused'
        accuracy['reason'] = (f'needs {need:,} bins and about {n * need:,} element operations, '
                              f'beyond the {MAX_BINS:,} bin and {MAX_WORK:,} work ceilings')
    else:
        accuracy['status'] = 'measured'
        accuracy.update(run_dp(weights, freqs, signs, need))
        accuracy['achieved_error_in_sd'] = accuracy['worst_case_absolute_error'] / sd
        accuracy['mean_absolute_discrepancy'] = abs(accuracy['mean_from_dp'] - mean_true)

    return {
        **meta,
        'n_used': n,
        'effect_allele_frequency_source': 'published with the score' if published
                                          else f'declared uniform {DECLARED_FREQUENCY}',
        'weight_summary': {'min': float(weights.min()), 'max': float(weights.max()),
                           'mean_abs': float(np.abs(weights).mean()),
                           'sum_abs': float(np.abs(weights).sum())},
        'complete_kernel': kernel_costs(n),
        'score_dp': {
            'discretisation_step': delta,
            'bins': int(len(pmf)),
            'bins_predicted': bins_predicted,
            'pmf_payload_bytes': int(pmf.nbytes),
            'seconds': seconds,
            'worst_case_absolute_error': worst_case,
            'mean_analytic': mean_true, 'mean_from_dp': mean_dp,
            'variance_analytic': var_true, 'variance_from_dp': var_dp,
            'sd_analytic': sd,
            'worst_case_error_in_sd': (worst_case / sd) if sd > 0 else None,
            'mean_absolute_discrepancy': abs(mean_dp - mean_true),
        },
        'accuracy_controlled': accuracy,
    }


def scaling_study(seed=20260906):
    """Isolate how the accuracy-controlled bin count grows with n.

    The published scores differ in method and weight scale, so an exponent fitted
    across them is confounded. Here everything is held fixed except n: weights are
    drawn from one distribution with a fixed seed, and the effect-allele frequency
    is constant. The derived expectation is B ~ n^(3/2), because the required bins
    scale as n * sum|w| / SD, with sum|w| growing like n and SD like sqrt(n).
    """
    rng = np.random.default_rng(seed)
    rows = []
    for n in [50, 100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600]:
        weights = rng.normal(0.0, 0.05, size=n)
        freqs = np.full(n, DECLARED_FREQUENCY)
        sd = math.sqrt(float(np.sum(2 * freqs * (1 - freqs) * weights ** 2)))
        need = required_bins(weights, sd, TARGET_ERROR_IN_SD)
        row = {'n': n, 'bins_required': need, 'sd': sd,
               'sum_abs_weight': float(np.abs(weights).sum())}
        if need <= MAX_BINS and n * need <= MAX_WORK:
            row.update(run_dp(weights, freqs, np.sign(weights), need))
            row['status'] = 'measured'
        else:
            row['status'] = 'refused'
        rows.append(row)
    ns = np.array([r['n'] for r in rows], dtype=float)
    bins = np.array([r['bins_required'] for r in rows], dtype=float)
    slope, intercept = np.polyfit(np.log(ns), np.log(bins), 1)
    timed = [r for r in rows if r['status'] == 'measured']
    work_slope = None
    if len(timed) >= 3:
        w = np.array([r['n'] * r['bins'] for r in timed], dtype=float)
        s = np.array([r['seconds'] for r in timed], dtype=float)
        work_slope = float(np.polyfit(np.log(w), np.log(s), 1)[0])
    return {'rows': rows, 'fitted_bin_exponent': float(slope),
            'derived_bin_exponent': 1.5,
            'fitted_time_vs_work_exponent': work_slope,
            'note': ('Weights drawn once per n from Normal(0, 0.05) with a fixed seed and a '
                     'constant effect-allele frequency, so only n varies. The bin exponent is '
                     'derived, not assumed: sum|w| grows like n and the score SD like sqrt(n).')}


def figure(records):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
    ns = [r['n_used'] for r in records]
    # 3^n has tens of thousands of digits at the top of this range, so the y axis
    # carries log10 of the entry count rather than the count itself.
    axes[0].semilogx(ns, [r['complete_kernel']['log10_genotypes'] for r in records],
                     'o-', color=COLORS[2], label='complete catalog, log10 of 3^n')
    axes[0].semilogx(ns, [math.log10(r['score_dp']['bins']) for r in records], 's-',
                     color=COLORS[0], label='score DP, log10 of bins retained')
    axes[0].axhline(math.log10(BUDGET_BYTES / 8), color=COLORS[3], ls='--', lw=1.6,
                    label='256 MiB budget, log10 of binary64 entries')
    axes[0].set(xlabel='Variants in the published score, n',
                ylabel='log10(entries)',
                title='One phenotype: what each regime would store')
    axes[0].legend(fontsize=7.5)
    axes[0].grid(which='both', alpha=.15)

    axes[1].loglog(ns, [r['score_dp']['seconds'] for r in records], 'o-', color=COLORS[0],
                   label='measured DP time')
    reference = [r['score_dp']['seconds'] for r in records]
    scale = reference[0] / ns[0]
    axes[1].loglog(ns, [scale * x for x in ns], '--', color=COLORS[1], lw=1.5,
                   label='linear in n, for reference')
    axes[1].set(xlabel='Variants in the published score, n', ylabel='Seconds',
                title='Measured cost of the score distribution')
    axes[1].legend(fontsize=8)
    axes[1].grid(which='both', alpha=.15)
    fig.tight_layout()
    plt.savefig(FIGURES / 'pgs_regimes.png', dpi=190, bbox_inches='tight')
    plt.savefig(FIGURES / 'pgs_regimes.svg', bbox_inches='tight')
    plt.close()


def run():
    RESULTS.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.titleweight': 'bold', 'figure.facecolor': 'white'})

    catalogue = []
    for pgs_id in get_json(f'{API}/trait/{TRAIT}')['associated_pgs_ids']:
        s = get_json(f'{API}/score/{pgs_id}')
        catalogue.append({'id': s['id'], 'variants': s['variants_number'],
                          'firstauthor': s['publication']['firstauthor'],
                          'date': s['publication']['date_publication']})
        time.sleep(0.1)
    catalogue.sort(key=lambda r: r['variants'])

    records = []
    for pgs_id in LADDER:
        meta, weights, freqs = fetch_score(pgs_id)
        record = analyse(meta, weights, freqs)
        records.append(record)
        acc = record['accuracy_controlled']
        tail = (f"needs {acc['bins_required']:,} bins: {acc['status']}"
                if acc['status'] != 'measured'
                else f"{acc['bins']:,} bins in {acc['seconds']:.2f}s for 1% SD")
        print(f"  {record['id']}  n={record['n_used']:>7,}  "
              f"fixed: {record['score_dp']['bins']:>7,} bins "
              f"{record['score_dp']['seconds']:>6.3f}s "
              f"err={record['score_dp']['worst_case_error_in_sd']:>9.3f} SD  |  {tail}")
        time.sleep(0.2)

    scaling = scaling_study()
    print()
    print(f"  controlled scaling: bins ~ n^{scaling['fitted_bin_exponent']:.3f} "
          f"(derived 1.5)", end='')
    if scaling['fitted_time_vs_work_exponent'] is not None:
        print(f"; time ~ work^{scaling['fitted_time_vs_work_exponent']:.3f} (derived 1.0)")
    else:
        print()
    figure(records)
    out = {
        'trait': {'id': TRAIT, 'label': 'coronary artery disorder'},
        'catalogue_summary': {
            'scores_for_this_trait': len(catalogue),
            'smallest': catalogue[0], 'largest': catalogue[-1],
        },
        'ladder': records,
        'controlled_scaling': scaling,
        'assumptions': [
            'Hardy-Weinberg dosage probabilities at every variant.',
            'Linkage equilibrium between variants. This is false for scores built by LD-aware '
            'methods, and it means the reported variance is not the variance of the real score.',
            'Effect weights are treated as fixed and known; their standard errors are ignored.',
            'Where a score publishes no effect-allele frequency, a single declared frequency is '
            'used and the run measures cost only, not any population quantity.',
        ],
        'not_established': [
            'No individual risk is computed for anyone.',
            'No validation of any published score is attempted or implied.',
            'The distributions are of a discretised score under the assumptions above, not of '
            'observed disease risk.',
        ],
        'redistribution_note': ('Scoring files are fetched from the PGS Catalog at run time and '
                                'are not redistributed here. Only derived quantities, digests and '
                                'provenance are retained.'),
        'completed_at_utc': datetime.now(timezone.utc).isoformat(),
    }
    (RESULTS / 'pgs_catalog.json').write_text(json.dumps(out, indent=2, allow_nan=False) + '\n',
                                              encoding='utf-8')
    print(f"\n  {len(catalogue)} published scores for this trait, "
          f"{catalogue[0]['variants']:,} to {catalogue[-1]['variants']:,} variants")
    return out


if __name__ == '__main__':
    run()
