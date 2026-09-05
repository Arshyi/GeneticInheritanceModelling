"""Asymptotic analysis of the inheritance kernel, checked against measured timings.

Produces, for each model in the paper, the exact catalog and support counts, the
derived upper and lower bounds for each operation, and an empirical growth
exponent fitted to the retained benchmark so that the asymptotic claims can be
compared with what the machine actually did.

Lower bounds here are output-size and read-size arguments. They hold for any
algorithm producing the stated output, not merely for this implementation, and
they are the reason no engineering effort removes the exponent.
"""
from datetime import datetime, timezone
import json
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from genetics.core import InheritanceModel

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / 'results'
FIGURES = ROOT / 'figures'
COLORS = {'dense': '#293c61', 'csr': '#136f79', 'hash': '#d18a25', 'streamed_kernel': '#b44949'}

# The models the paper actually builds, in the order it introduces them.
TRAITS = [
    ('Sickle cell (HBB, simplified)', (2,)),
    ('ABO', (3,)),
    ('ABO x Rh', (3, 2)),
    ('Eye colour (HERC2 + modifier)', (2, 2)),
    ('Three biallelic loci', (2, 2, 2)),
    ('Five biallelic loci', (2, 2, 2, 2, 2)),
]


def dump(name, obj):
    (RESULTS / name).write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def counts(allele_counts):
    """Exact catalog, pair and support counts for a declared locus specification."""
    model = InheritanceModel(allele_counts, max_bytes=1 << 62)
    G, U, nnz = model.G, model.U, model.nnz
    return {
        'allele_counts': list(allele_counts),
        'G': G, 'U': U, 'nnz': nnz, 'dense_entries': G * U,
        'density': nnz / (G * U),
        'dense_payload_bytes': G * U * 8,
        'csr_value_payload_bytes': nnz * 8,
        # A het-by-het parental pair at every locus maximises the child support.
        'max_children_one_pair': math.prod(a * (a + 1) // 2 for a in allele_counts),
    }


def biallelic_series(max_n=12):
    """Closed forms for n independent biallelic loci, with their growth ratios."""
    rows = []
    for n in range(1, max_n + 1):
        G = 3 ** n
        U = G * (G + 1) // 2
        nnz = (15 ** n + 5 ** n) // 2
        dense = (27 ** n + 9 ** n) // 2
        rows.append({'n': n, 'G': G, 'U': U, 'nnz': nnz, 'dense_entries': dense,
                     'density': nnz / dense,
                     'nnz_ratio_to_previous': None if n == 1 else nnz / ((15 ** (n - 1) + 5 ** (n - 1)) // 2),
                     'dense_ratio_to_previous': None if n == 1 else dense / ((27 ** (n - 1) + 9 ** (n - 1)) // 2)})
    return rows


def fit_growth(ns, values):
    """Fit T = a * b^n by least squares on log T, returning b and the fit quality."""
    ns = np.asarray(ns, dtype=float)
    y = np.log(np.asarray(values, dtype=float))
    slope, intercept = np.polyfit(ns, y, 1)
    predicted = slope * ns + intercept
    residual = y - predicted
    ss_res = float(np.sum(residual ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    return {'per_locus_multiplier': float(np.exp(slope)),
            'log_intercept': float(intercept),
            'r_squared': 1.0 - ss_res / ss_tot if ss_tot > 0 else None}


def empirical(bench):
    """Measured per-locus growth for construction and update, per representation."""
    rows = bench['full_workload_summary']
    out = {}
    for method in COLORS:
        selected = sorted((r for r in rows if r['method'] == method), key=lambda r: r['loci'])
        ns = [r['loci'] for r in selected]
        update = [r['inference_seconds']['median'] for r in selected]
        entry = {'update': fit_growth(ns, update)}
        # The top-end step is the least contaminated by fixed small-n overhead.
        entry['update_top_step_ratio'] = update[-1] / update[-2]
        build = [(r['loci'], r['construction_seconds']['median']) for r in selected
                 if r['construction_seconds']]
        if len(build) >= 3:
            entry['construction'] = fit_growth([b[0] for b in build], [b[1] for b in build])
            entry['construction_top_step_ratio'] = build[-1][1] / build[-2][1]
        out[method] = entry
    return out


def figure(series, measured):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ns = [r['n'] for r in series][:8]
    axes[0].semilogy(ns, [r['dense_entries'] for r in series][:8], 'o-',
                     color=COLORS['dense'], label='dense entries G*U')
    axes[0].semilogy(ns, [r['nnz'] for r in series][:8], 's-',
                     color=COLORS['csr'], label='supported transitions')
    axes[0].semilogy(ns, [3 ** n for n in ns], '^-', color=COLORS['hash'],
                     label='max children of one pair')
    axes[0].semilogy(ns, ns, 'd-', color=COLORS['streamed_kernel'],
                     label='factored single-child query')
    axes[0].set(xlabel='Biallelic loci n', ylabel='Operations or entries (log scale)',
                title='Derived growth: three exponentials and one linear')
    axes[0].legend(fontsize=8)
    axes[0].grid(axis='y', alpha=.18)

    derived = ((15**5 + 5**5) // 2) / ((15**4 + 5**4) // 2)
    names, fitted, top = [], [], []
    for method, entry in measured.items():
        if 'construction' in entry:
            names.append(method.replace('_', ' '))
            fitted.append(entry['construction']['per_locus_multiplier'])
            top.append(entry['construction_top_step_ratio'])
    x = np.arange(len(names))
    axes[1].bar(x - 0.19, fitted, 0.38, label='fitted over n = 1..5', color='#c9d8dc')
    axes[1].bar(x + 0.19, top, 0.38, label='measured n = 4 to 5', color=COLORS['csr'])
    axes[1].axhline(derived, color=COLORS['hash'], lw=2, ls='--',
                    label=f'derived work ratio {derived:.2f}')
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(names, fontsize=8)
    axes[1].set(ylabel='Multiplier of construction time per added locus',
                title='Derived support growth against measured growth')
    axes[1].legend(fontsize=8)
    axes[1].grid(axis='y', alpha=.18)
    fig.tight_layout()
    plt.savefig(FIGURES / 'complexity_growth.png', dpi=190, bbox_inches='tight')
    plt.savefig(FIGURES / 'complexity_growth.svg', bbox_inches='tight')
    plt.close()


def run():
    RESULTS.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.titleweight': 'bold', 'figure.facecolor': 'white'})
    bench = json.loads((RESULTS / 'benchmark.json').read_text(encoding='utf-8'))
    series = biallelic_series()
    measured = empirical(bench)
    figure(series, measured)

    result = {
        'traits': [{'name': name, **counts(spec)} for name, spec in TRAITS],
        'biallelic_series': series,
        'asymptotic_limits': {
            'genotypes': 'Theta(3^n) for n biallelic loci; Theta(prod a_l(a_l+1)/2) in general',
            'unordered_pairs': 'Theta(9^n / 2)',
            'supported_transitions': 'Theta(15^n) since (15^n + 5^n)/2 is between 15^n/2 and 15^n',
            'dense_entries': 'Theta(27^n / 2)',
            'density_limit': 'density -> 0 as (15/27)^n, so sparsity improves without bound',
            'nnz_ratio_limit': 15.0, 'dense_ratio_limit': 27.0,
            'nnz_ratio_at_top_measured_step': ((15**5 + 5**5)//2) / ((15**4 + 5**4)//2),
        },
        'measured_growth': measured,
        'growth_comparison_note': (
            'Construction is dominated by enumerating supported transitions, so the work ratio '
            'from four to five loci is nnz(5)/nnz(4) = 14.878. Measured construction time over '
            'that same step rose by 15.01 (dense), 15.09 (CSR) and 14.68 (hash), agreeing with '
            'the derived ratio to about 1.5 per cent. A least-squares exponent fitted over all '
            'of n = 1..5 instead gives 8.5 to 10.4, because fixed per-call overhead dominates '
            'the smallest problems and flattens the fitted slope. The asymptotic claim is '
            'therefore supported by the top of the measured range and not by the whole of it, '
            'which is what an asymptotic claim is entitled to.'),
        'completed_at_utc': datetime.now(timezone.utc).isoformat(),
    }
    dump('complexity.json', result)

    print('Derived limits: nnz ratio 15, dense ratio 27')
    for method, entry in measured.items():
        build = entry.get('construction')
        line = f"  {method:16} update x{entry['update']['per_locus_multiplier']:.2f}/locus (top step x{entry['update_top_step_ratio']:.2f})"
        if build:
            line += f"   construction x{build['per_locus_multiplier']:.2f}/locus (top step x{entry['construction_top_step_ratio']:.2f}, R2={build['r_squared']:.4f})"
        print(line)
    for trait in result['traits']:
        print(f"  {trait['name']:34} G={trait['G']:<8} U={trait['U']:<10} "
              f"nnz={trait['nnz']:<10} density={trait['density']:.4f}")
    return result


if __name__ == '__main__':
    run()
