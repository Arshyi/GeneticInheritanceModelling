"""Eye colour as a two-locus epistatic trait, run through the complete kernel.

Locus 1 is HERC2 rs12913832, for which real genotype calls are available and
frozen in data/. Locus 2 is a declared modifier locus standing for the residual
common pigmentation variation (OCA2, SLC24A4, SLC45A2, TYR, IRF4 and others),
compressed into one biallelic factor. That compression is the model's main
simplification and is not defended as biology; it is declared so that what the
model can and cannot support stays visible.

Nothing here validates an eye-colour prediction. The reference panel carries no
eye-colour observations. External phenotype comparisons use published cohort
percentages and are consistency checks, not held-out validation.
"""
from datetime import datetime, timezone
import csv
import json
from pathlib import Path
import time

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from genetics.core import InheritanceModel, hwe
from genetics.extensions import hwe_exact

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / 'results'
FIGURES = ROOT / 'figures'
DATA = ROOT / 'data'
SEED = 20260905
COLORS = ['#136f79', '#d18a25', '#293c61', '#b44949', '#79706e']

# Locus order: (HERC2, modifier). Both biallelic, so G = 9 and U = 45.
# HERC2 local codes: 0 = AA, 1 = AG, 2 = GG, with G the blue-associated allele.
# Modifier local codes: 0 = dd, 1 = Dd, 2 = DD, with D the darkening allele.
HERC2, MODIFIER = 0, 1
PHENOTYPES = ['blue', 'intermediate', 'hazel', 'brown']

# Declared phenotype map. Epistatic: the modifier acts differently on the two
# HERC2 backgrounds, which is the whole point of using two loci rather than one.
PHENOTYPE_MAP = {
    (2, 0): 'blue', (2, 1): 'blue', (2, 2): 'intermediate',
    (1, 0): 'hazel', (1, 1): 'brown', (1, 2): 'brown',
    (0, 0): 'hazel', (0, 1): 'brown', (0, 2): 'brown',
}

# External published anchors. Neither is used to fit this model.
CANPATH = {
    'cohort': 'Canadian Partnership for Tomorrow\'s Health (CanPath)',
    'individuals': 5481, 'gg_individuals': 2757, 'aa_ag_individuals': 2724,
    'gg_non_blue_fraction': 0.33,
    'aa_ag_non_brown_or_hazel_count': 512,
}
IRISPLEX = {'snps': 6, 'auc_blue': 0.96, 'auc_brown': 0.96, 'cohort_size': 6168}


def dump(name, obj):
    (RESULTS / name).write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def savefig(name):
    plt.savefig(FIGURES / (name + '.png'), dpi=190, bbox_inches='tight')
    plt.savefig(FIGURES / (name + '.svg'), bbox_inches='tight')
    plt.close()


def load_population_data():
    rows = []
    with (DATA / 'eye_color_genotypes.csv').open(encoding='utf-8') as handle:
        for row in csv.DictReader(handle):
            for key in ('hom_1', 'het', 'hom_2', 'n'):
                row[key] = int(row[key])
            rows.append(row)
    meta = json.loads((DATA / 'eye_color_genotypes.provenance.json').read_text(encoding='utf-8'))
    if meta['orientation']['blue_associated_allele'] != 'G':
        raise ValueError('Frozen provenance no longer identifies G as blue-associated')
    return rows, meta


def holm(pvalues):
    order = sorted(range(len(pvalues)), key=lambda i: pvalues[i])
    adjusted = [0.0] * len(pvalues)
    running = 0.0
    for rank, index in enumerate(order):
        running = max(running, (len(pvalues) - rank) * pvalues[index])
        adjusted[index] = min(1.0, running)
    return adjusted


def phenotype_of(model, code):
    herc2, modifier = model.decode(code)
    return PHENOTYPE_MAP[(herc2, modifier)]


def phenotype_vector(model, genotype_frequencies):
    totals = {name: 0.0 for name in PHENOTYPES}
    for code, mass in enumerate(genotype_frequencies):
        totals[phenotype_of(model, code)] += float(mass)
    return totals


def representation_agreement(model, frequencies):
    """The same population update through all four representations."""
    timings, results = {}, {}
    for kind in ('dense', 'csr', 'hash'):
        start = time.perf_counter()
        kernel = model.kernel(kind)
        build = time.perf_counter() - start
        start = time.perf_counter()
        out = model.next_generation(frequencies, kernel=kernel)
        update = time.perf_counter() - start
        timings[kind] = {'build_seconds': build, 'update_seconds': update}
        results[kind] = np.asarray(out, dtype=float)
    start = time.perf_counter()
    results['streamed'] = np.asarray(model.next_generation(frequencies), dtype=float)
    timings['streamed'] = {'build_seconds': None,
                           'update_seconds': time.perf_counter() - start}
    reference = results['dense']
    deviation = {k: float(np.max(np.abs(v - reference))) for k, v in results.items()}
    return timings, deviation, reference


def population_audit(model, rows):
    """Exact HWE per population, plus the model's predicted phenotype split."""
    records, pvalues = [], []
    for row in rows:
        counts = (row['hom_1'], row['het'], row['hom_2'])
        n = row['n']
        blue_allele = (2 * row['hom_2'] + row['het']) / (2 * n)
        p = hwe_exact(counts)
        pvalues.append(p)
        observed = np.array(counts, dtype=float) / n
        records.append({
            'population': row['population'], 'superpopulation': row['superpopulation'],
            'n': n, 'counts_AA_AG_GG': list(counts),
            'blue_allele_frequency': blue_allele,
            'observed_genotype_frequencies': observed.tolist(),
            'hwe_expected': (np.array(hwe([1 - blue_allele, blue_allele])) * n).tolist(),
            'hwe_exact_p': p,
        })
    for record, adjusted in zip(records, holm(pvalues)):
        record['hwe_exact_p_holm'] = adjusted
    return records


def predicted_phenotypes(model, records, modifier_frequency):
    """Join the real HERC2 genotype frequencies to the declared modifier locus.

    Linkage equilibrium between the two loci is assumed. HERC2 and the compressed
    modifier are not one physical locus, so this is an assumption about the
    population, not a consequence of Mendelian transmission.
    """
    modifier = hwe([1 - modifier_frequency, modifier_frequency])
    for record in records:
        herc2 = record['observed_genotype_frequencies']
        joint = np.zeros(model.G)
        for i, hi in enumerate(herc2):
            for j, mj in enumerate(modifier):
                joint[model.encode((i, j))] = hi * mj
        record['joint_genotype_mass'] = float(joint.sum())
        record['predicted_phenotypes'] = phenotype_vector(model, joint)
    return records


def figures(records, model, generations, deviation):
    order = sorted(records, key=lambda r: r['blue_allele_frequency'])
    names = [r['population'] for r in order]
    freqs = [r['blue_allele_frequency'] for r in order]
    groups = sorted({r['superpopulation'] for r in records})
    palette = dict(zip(groups, COLORS))

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    axes[0].bar(names, freqs, color=[palette[r['superpopulation']] for r in order])
    axes[0].set(ylabel='rs12913832 blue-associated allele frequency',
                title='Observed allele frequency, 26 populations')
    axes[0].tick_params(axis='x', rotation=90, labelsize=6.5)
    handles = [plt.Rectangle((0, 0), 1, 1, color=palette[g]) for g in groups]
    axes[0].legend(handles, groups, fontsize=7.5, ncol=2)

    blue = [r['predicted_phenotypes']['blue'] for r in order]
    brown = [r['predicted_phenotypes']['brown'] for r in order]
    axes[1].plot(freqs, blue, 'o', color=COLORS[0], label='predicted blue')
    axes[1].plot(freqs, brown, 's', color=COLORS[2], label='predicted brown')
    axes[1].set(xlabel='blue-associated allele frequency',
                ylabel='model phenotype probability',
                title='Declared two-locus map applied to real genotypes')
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.grid(axis='y', alpha=.18)
    fig.tight_layout()
    savefig('eye_color_populations')

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for name in PHENOTYPES:
        axes[0].plot(range(len(generations)), [g[name] for g in generations],
                     marker='o', label=name)
    axes[0].set(xlabel='Generation', ylabel='Phenotype probability',
                title='Random mating from the CEU genotype frequencies')
    axes[0].legend(fontsize=8)
    axes[0].grid(axis='y', alpha=.18)

    kernel = np.asarray(model.kernel('dense'), dtype=float)
    axes[1].imshow(kernel > 0, aspect='auto', cmap='Blues', interpolation='nearest')
    axes[1].set(xlabel='Unordered parental pair (U = %d)' % model.U,
                ylabel='Child genotype (G = %d)' % model.G,
                title='Support of the complete 9 x 45 kernel')
    fig.tight_layout()
    savefig('eye_color_kernel')


def run():
    RESULTS.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.titleweight': 'bold', 'figure.facecolor': 'white'})
    rows, meta = load_population_data()
    model = InheritanceModel((2, 2))

    # Two ways of setting the modifier frequency, kept apart on purpose.
    declared = 0.5
    declared_discordance = declared ** 2
    calibrated = float(np.sqrt(CANPATH['gg_non_blue_fraction']))

    records = population_audit(model, rows)
    records = predicted_phenotypes(model, records, calibrated)

    ceu = next(r for r in records if r['population'] == 'CEU')
    modifier = hwe([1 - calibrated, calibrated])
    start = np.zeros(model.G)
    for i, hi in enumerate(ceu['observed_genotype_frequencies']):
        for j, mj in enumerate(modifier):
            start[model.encode((i, j))] = hi * mj
    generations, frequencies = [phenotype_vector(model, start)], start
    for _ in range(6):
        frequencies = np.asarray(model.next_generation(frequencies), dtype=float)
        generations.append(phenotype_vector(model, frequencies))

    timings, deviation, _ = representation_agreement(model, start)
    figures(records, model, generations, deviation)

    dense_bytes = model.G * model.U * 8
    csr_bytes = model.nnz * 8
    result = {
        'model': {
            'loci': ['HERC2 rs12913832 (real genotype calls)',
                     'declared modifier locus (compressed residual pigmentation variation)'],
            'genotypes_G': model.G, 'unordered_pairs_U': model.U,
            'supported_transitions_nnz': model.nnz,
            'dense_entries_GU': model.G * model.U,
            'dense_payload_bytes': dense_bytes, 'csr_value_payload_bytes': csr_bytes,
            'phenotype_classes': PHENOTYPES,
            'phenotype_map': {f'{k[0]}_{k[1]}': v for k, v in PHENOTYPE_MAP.items()},
        },
        'representation_agreement': {
            'timings_seconds': timings,
            'max_absolute_deviation_from_dense': deviation,
        },
        'modifier_locus': {
            'declared_frequency': declared,
            'declared_implied_gg_non_blue': declared_discordance,
            'calibrated_frequency_from_canpath': calibrated,
            'calibration_note': ('The calibrated value reproduces the published GG discordance '
                                 'by construction. It is a calibration, not a validation, and no '
                                 'accuracy is claimed from it.'),
        },
        'external_anchors': {'canpath': CANPATH, 'irisplex': IRISPLEX},
        'population_records': records,
        'ceu_generations': generations,
        'data_provenance': {k: meta[k] for k in
                            ('dataset_id', 'variant', 'assembly', 'chromosome', 'position',
                             'strand', 'ensembl_release', 'individuals', 'csv_sha256')},
        'orientation': meta['orientation'],
        'validation_status': ('Genotype calls only; the panel carries no eye-colour '
                              'observations. No predictive accuracy is established.'),
        'completed_at_utc': datetime.now(timezone.utc).isoformat(),
    }
    dump('eye_color.json', result)

    significant = [r['population'] for r in records if r['hwe_exact_p_holm'] < 0.05]
    print(f'Eye colour: G={model.G} U={model.U} nnz={model.nnz} '
          f'dense={model.G*model.U}')
    print(f'  max deviation across representations: '
          f'{max(deviation.values()):.3g}')
    print(f'  populations failing Holm-adjusted HWE at 0.05: {significant or "none"}')
    print(f'  calibrated modifier frequency: {calibrated:.6f}')
    return result


if __name__ == '__main__':
    run()
