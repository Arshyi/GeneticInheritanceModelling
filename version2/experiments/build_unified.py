"""Build the unified Version I + Version II manuscript and its typeset PDF.

Version I supplies the narrative spine and every piece of its mathematics.
Version II supplies the audit, the complete kernel, the measured comparisons,
and the age extension. All numeric bindings come from retained machine-readable
results, exactly as in build_manuscript.py; nothing here is typed by hand.
"""
import hashlib
import json

from experiments.build_manuscript import (
    MAN, RES, ROOT, apply_citations, bindings, read, render, table,
)

# Version I's own MATLAB, in the order the paper introduces the models.
LEGACY_MATLAB = [
    ('sickle_cell_model.m',
     'The sickle-cell model of Part I. Generation one is computed from the parental pair by '
     'Punnett logic; generations two onward apply the 3 x 3 transition matrix of Section 4.1. '
     'The line `X_next(3) = 0;` is the operation audited in Section 4.4: the affected class is '
     'zeroed before the next multiplication and the vector is not renormalised, which is why the '
     'displayed generation-three output sums to 93.75% rather than to one.'),
    ('abo_simulation.m',
     'The single-locus ABO offspring computation. It derives offspring genotype probabilities '
     'from parental genotypes directly rather than from the retained six-column table.'),
    ('abo_run.m',
     'The ABO driver. This is the listing that settles Section 5.5: it loops over the full set of '
     'ordered parental genotype pairs and weights each cross by the current frequencies, so the '
     'executed model is a complete random-mating calculation over all 36 ordered pairs, not the '
     'six retained columns of the displayed table.'),
    ('bloodgroup18_simulation.m',
     'The joint ABO x Rh offspring computation, combining independent transmission at the two '
     'loci to give the 18-entry combined catalog.'),
    ('aborh_run.m',
     'The ABO x Rh driver, and the listing behind Result 6.2. It enumerates all 324 ordered '
     'parental pairs. The displayed 18 x 18 table is therefore not what this program iterates.'),
]


def coverage_tables(legacy):
    """Every reproduction of Version I's two coverage percentages, side by side."""
    cov = legacy['coverage']
    variants = cov['variants']
    label = {
        'paper_alleles_exact': 'Stated allele frequencies, propagated exactly',
        'paper_p32_genotype_table': 'Four-decimal genotype table on page 32',
        'paper_p32_pair_formula_inputs_raw': 'Coarser inputs used in the pair formulas',
        'paper_p32_pair_formula_inputs_normalized': 'Same coarse inputs, normalised first',
    }
    rows = [
        (label[k],
         f"{v['top6_percent_raw_mass']:.10f}",
         f"{v['top18_joint_percent_raw_mass']:.10f}",
         f"{v['genotype_mass']:.6f}")
        for k, v in variants.items()
    ]
    rows.append(('Printed percentages added as displayed',
                 f"{cov['paper_top6_printed_sum_percent']:.2f}",
                 f"{cov['paper_top18_claimed_total_percent']:.2f}",
                 'not applicable'))
    variants_table = table(
        ['Calculation', 'Top-six ABO pair mass (%)', 'Top-18 joint pair mass (%)', 'Input genotype mass'],
        rows)

    exact = variants['paper_alleles_exact']
    top6 = table(['Rank', 'Parental cross', 'Model probability (%)'],
                 [(i, ' x '.join(e['parents']), f"{e['percent']:.6f}")
                  for i, e in enumerate(exact['top6'], 1)])
    top18 = table(['Rank', 'Parental cross', 'Model probability (%)'],
                  [(i, ' x '.join(e['parents']), f"{e['percent']:.6f}")
                   for i, e in enumerate(exact['top18_joint'], 1)])
    return variants_table, top6, top18, exact


def sparsity_table():
    """What a positive mutation rate does to the supported-transition count."""
    rows = []
    for n in range(1, 7):
        G = 3 ** n
        U = G * (G + 1) // 2
        mendelian = (15 ** n + 5 ** n) // 2
        dense = G * U
        rows.append((n, f'{G:,}', f'{U:,}', f'{mendelian:,}', f'{dense:,}',
                     f'{dense / mendelian:.2f}x'))
    return table(['Biallelic loci n', 'Genotypes G', 'Unordered pairs U',
                  'Nonzeros, no mutation', 'Nonzeros, any mutation', 'Growth factor'], rows)


def legacy_sickle_table(legacy):
    """The displayed legacy trajectory and the mass it does not conserve."""
    s = legacy['sickle_p28']
    rows = [(i, ', '.join(f'{x:.4f}' for x in pct), f'{mass:.6f}')
            for i, (pct, mass) in enumerate(zip(s['history_percent'], s['history_mass']))]
    return table(['Generation', 'AA, AS, SS (%) as displayed', 'Total mass'], rows[:5])


def matlab_appendix():
    """Version I's MATLAB, embedded verbatim from legacy/matlab, with digests."""
    blocks, digests = [], {}
    for name, lead in LEGACY_MATLAB:
        path = ROOT / 'legacy' / 'matlab' / name
        source = path.read_text(encoding='utf-8').replace('\r\n', '\n').rstrip('\n')
        digests[name] = {
            'sha256': hashlib.sha256(source.encode('utf-8')).hexdigest(),
            'lines': source.count('\n') + 1,
        }
        blocks.append(f'## {name}\n\n{lead}\n\n```text\n{source}\n```')
    return '\n\n'.join(blocks), digests


def eye_tables(eye):
    """Everything the eye-colour part reports, bound to results/eye_color.json."""
    model = eye['model']
    agreement = eye['representation_agreement']
    rep = table(
        ['Representation', 'Build (ms)', 'Update (ms)', 'Max deviation from dense'],
        [(kind.replace('_', ' '),
          'no full kernel' if t['build_seconds'] is None else f"{t['build_seconds']*1000:.3f}",
          f"{t['update_seconds']*1000:.3f}",
          f"{agreement['max_absolute_deviation_from_dense'][kind]:.3g}")
         for kind, t in agreement['timings_seconds'].items()])

    records = sorted(eye['population_records'], key=lambda r: -r['blue_allele_frequency'])
    pop = table(
        ['Population', 'Group', 'N', 'AA, AG, GG', 'Blue allele', 'Exact HWE p', 'Holm p',
         'Model P(blue)'],
        [(r['population'], r['superpopulation'], r['n'], ', '.join(map(str, r['counts_AA_AG_GG'])),
          f"{r['blue_allele_frequency']:.4f}", f"{r['hwe_exact_p']:.4f}",
          f"{r['hwe_exact_p_holm']:.4f}", f"{r['predicted_phenotypes']['blue']:.4f}")
         for r in records])

    gens = table(
        ['Generation'] + [p.capitalize() for p in eye['model']['phenotype_classes']],
        [(i, *[f"{g[p]:.6f}" for p in eye['model']['phenotype_classes']])
         for i, g in enumerate(eye['ceu_generations'])])

    canpath, iris = eye['external_anchors']['canpath'], eye['external_anchors']['irisplex']
    modifier = eye['modifier_locus']
    external = table(
        ['Quantity', 'This model', 'Published', 'Source'],
        [('Loci used', len(model['loci']), f"{iris['snps']} SNPs (IrisPlex)", '[walsh2011irisplex]'),
         ('P(non-blue | GG), declared modifier at 0.5',
          f"{modifier['declared_implied_gg_non_blue']:.4f}",
          f"{canpath['gg_non_blue_fraction']:.2f}", '[abbatangelo2026canpath]'),
         ('Modifier frequency reproducing that discordance',
          f"{modifier['calibrated_frequency_from_canpath']:.6f}", 'not applicable',
          'calibrated here, not measured'),
         ('Phenotyped individuals', 0, f"{canpath['individuals']:,} (CanPath); {iris['cohort_size']:,} (IrisPlex development)",
          '[abbatangelo2026canpath] [walsh2011irisplex]'),
         ('Reported AUC, blue', 'none; no phenotypes available', f"{iris['auc_blue']:.2f}",
          '[walsh2011irisplex]'),
         ('Reported AUC, brown', 'none; no phenotypes available', f"{iris['auc_brown']:.2f}",
          '[walsh2011irisplex]')])
    return rep, pop, gens, external


def complexity_tables(comp):
    traits = table(
        ['Model', 'Allele counts', 'G', 'U', 'Supported transitions', 'Dense entries', 'Density'],
        [(t['name'], ' x '.join(map(str, t['allele_counts'])), f"{t['G']:,}", f"{t['U']:,}",
          f"{t['nnz']:,}", f"{t['dense_entries']:,}", f"{t['density']:.4f}")
         for t in comp['traits']])
    measured = comp['measured_growth']
    growth = table(
        ['Representation', 'Fitted over n = 1..5', 'Measured n = 4 to 5', 'Derived work ratio'],
        [(method.replace('_', ' '),
          f"{entry['construction']['per_locus_multiplier']:.2f}" if 'construction' in entry else 'no full kernel',
          f"{entry['construction_top_step_ratio']:.2f}" if 'construction' in entry else 'no full kernel',
          f"{comp['asymptotic_limits']['nnz_ratio_at_top_measured_step']:.3f}")
         for method, entry in measured.items()])
    return traits, growth


def main():
    ledger, values = bindings()
    legacy = read('version1_reproduction.json')
    poly = read('polygenic_synthetic.json')
    eye = read('eye_color.json')
    comp = read('complexity.json')
    counts = legacy['state_counts']
    checks = legacy['displayed_matrix_checks']

    variants_table, top6, top18, exact = coverage_tables(legacy)
    matlab_source, matlab_digests = matlab_appendix()
    eye_rep, eye_pop, eye_gens, eye_external = eye_tables(eye)
    comp_traits, comp_growth = complexity_tables(comp)

    values.update({
        'V1_COVERAGE_VARIANTS': variants_table,
        'V1_TOP6_TABLE': top6,
        'V1_TOP18_TABLE': top18,
        'V1_LEGACY_SICKLE_TABLE': legacy_sickle_table(legacy),
        'SPARSITY_MUTATION_TABLE': sparsity_table(),
        'STATE_COUNT_TABLE': table(
            ['Counted object', 'Simplified ABO', 'Simplified ABO x Rh'],
            [('Genotypes in the catalog', counts['ABO_genotypes'], counts['ABORh_genotypes']),
             ('Ordered parental pairs', counts['ABO_ordered_parent_pairs'], counts['ABORh_ordered_parent_pairs']),
             ('Unordered parental pairs', counts['ABO_unordered_parent_pairs'], counts['ABORh_unordered_parent_pairs']),
             ('Distinct-genotype pairs only', counts['ABO_distinct_only_unordered_pairs'], 153),
             ('Square matrix Version I retained', 6, 18)]),
        'DISPLAYED_MATRIX_CHECK': (
            f"All {324 - checks['joint_table_vs_self_cross_entry_mismatches_of_324']:,} of the 324 displayed joint-table entries "
            f"match each combined genotype crossed with itself. Against the ranked top-18 construction instead, "
            f"{checks['joint_table_vs_ranked_top18_entry_mismatches_of_324']} of 324 cells differ. Every displayed column sums to one: "
            f"{str(checks['joint_table_all_column_sums_one']).lower()}."),
        'EXACT_TOP6': f"{exact['top6_percent_raw_mass']:.10f}",
        'EXACT_TOP18': f"{exact['top18_joint_percent_raw_mass']:.10f}",
        'PRINTED_TOP6': f"{legacy['coverage']['paper_top6_printed_sum_percent']:.2f}",
        'PRINTED_TOP18': f"{legacy['coverage']['paper_top18_claimed_total_percent']:.2f}",
        'POLYGENIC_SUMMARY': (
            f"The test uses {poly['n_loci']} loci with dosage probabilities (1/4,1/2,1/4) and unit weights. The genotype catalog has "
            f"3^{poly['n_loci']} entries; the score distribution has {poly['score_bins']} bins and a probability array of "
            f"{poly['pmf_bytes']:,} bytes. Score mean is {poly['score_mean']:.1f} and variance {poly['score_variance']:.1f}. With residual "
            f"standard deviation {poly['residual_sd']:.1f}, seed {poly['seed']}, and {poly['independent_synthetic_draws']:,} independent draws, "
            f"{poly['observed_interval_coverage']*100:.2f}% fall inside the nominal {poly['nominal_interval_mass']*100:.0f}% interval. The "
            f"probability-integral-transform KS statistic is {poly['pit_KS_statistic']:.5f}."),
        'EYE_REPRESENTATION_TABLE': eye_rep,
        'EYE_POPULATION_TABLE': eye_pop,
        'EYE_GENERATION_TABLE': eye_gens,
        'EYE_EXTERNAL_TABLE': eye_external,
        'EYE_MODEL_COUNTS': (
            f"G = {eye['model']['genotypes_G']}, U = {eye['model']['unordered_pairs_U']}, "
            f"{eye['model']['supported_transitions_nnz']} supported transitions out of "
            f"{eye['model']['dense_entries_GU']} dense entries, a density of "
            f"{eye['model']['supported_transitions_nnz']/eye['model']['dense_entries_GU']:.4f}. "
            f"Dense payload is {eye['model']['dense_payload_bytes']:,} bytes against "
            f"{eye['model']['csr_value_payload_bytes']:,} bytes of CSR values."),
        'EYE_AGREEMENT': (
            f"{max(eye['representation_agreement']['max_absolute_deviation_from_dense'].values()):.3g}"),
        'EYE_HWE_RESULT': (
            'none of the 26 component populations' if not [
                r for r in eye['population_records'] if r['hwe_exact_p_holm'] < 0.05]
            else ', '.join(r['population'] for r in eye['population_records']
                           if r['hwe_exact_p_holm'] < 0.05)),
        'EYE_ORIENTATION': (
            f"FIN {eye['orientation']['northern_european_frequency']['FIN']:.4f}, "
            f"GBR {eye['orientation']['northern_european_frequency']['GBR']:.4f}, "
            f"CEU {eye['orientation']['northern_european_frequency']['CEU']:.4f}"),
        'COMPLEXITY_TRAIT_TABLE': comp_traits,
        'COMPLEXITY_GROWTH_TABLE': comp_growth,
        'COMPLEXITY_NOTE': comp['growth_comparison_note'],
        'PN_TABLE': table(
            ['Loci n', 'Genotypes G = 3^n', 'Pairs U', 'Supported transitions',
             'Dense entries', 'Density'],
            [(r['n'], f"{r['G']:,}", f"{r['U']:,}", f"{r['nnz']:,}",
              f"{r['dense_entries']:,}", f"{r['density']:.3g}")
             for r in comp['biallelic_series'] if r['n'] <= 10]),
        'ARTIFACT_TABLE': table(['Artifact', 'Purpose', 'Evidence status'], [
            ('research/version1_audit.md', 'Page-by-page discrepancy explanation', 'Full paper and screenshots inspected'),
            ('results/version1_reproduction.json', 'Six exact legacy examples and coverage variants', 'Executed rational translation'),
            ('genetics/core.py; extensions.py', 'Reusable kernel and staged extensions', '32 automated tests passed'),
            ('results/benchmark.json and benchmark_workers/', 'Raw timings, memory, metadata and checks', 'Executed local CPU run'),
            ('results/complexity.json', 'Derived bounds and measured growth', 'Executed; derived limits stated separately'),
            ('data/observed_genotypes.csv + provenance', 'rs334 calls with validation', 'Frozen public snapshot, digest recorded'),
            ('data/eye_color_genotypes.csv + provenance', 'rs12913832 calls with validation', 'Frozen public snapshot, digest recorded'),
            ('results/population_validation.json', 'HWE and pooling audit', 'Executed; no clinical accuracy claim'),
            ('results/eye_color.json', 'Two-locus eye-colour model and audit', 'Executed; no predictive accuracy claim'),
            ('results/polygenic_synthetic.json', 'Score and interval calibration', '10,000 independent synthetic draws'),
            ('legacy/matlab/', 'Version I source, verbatim', 'Embedded as Appendix F with digests'),
            ('sources/source_ledger.json', 'Reproducible claim/source bibliography', 'Primary and authoritative sources'),
            ('manuscript/genetics_unified.md', 'This document, generated from template and results', 'Rebuilt by run.py unified'),
            ('output/pdf/Genetics_Complete.pdf', 'Typeset preprint', 'Rendered and visually checked')]),
        'MATLAB_APPENDIX': matlab_source,
        'RUN_COMMANDS_UNIFIED': (
            '```text\ncd version2\npython -m pip install -r requirements.txt\n'
            'python run.py test -q\npython run.py reproduce\npython run.py fetch\n'
            'python run.py science\npython run.py benchmark\npython run.py manuscript\n'
            'python run.py unified\n```'),
    })

    text = (MAN / 'unified_template.md').read_text(encoding='utf-8')
    for key, val in values.items():
        text = text.replace('{{' + key + '}}', str(val))
    text = apply_citations(text, ledger, citations_name='unified_citations.json')
    (MAN / 'genetics_unified.md').write_text(text, encoding='utf-8')
    (RES / 'unified_manifest.json').write_text(json.dumps({
        'document': 'manuscript/genetics_unified.md',
        'pdf': 'output/pdf/Genetics_Complete.pdf',
        'version1_pages_absorbed': legacy['source']['pages_read'],
        'version1_sha256': legacy['source']['sha256'],
        'version1_matlab_embedded': matlab_digests,
        'status': 'Preprint; not peer reviewed. Author-approved for public release.',
    }, indent=2) + '\n', encoding='utf-8')
    return text


if __name__ == '__main__':
    render(main(),
           out_name='Genetics_Complete.pdf',
           doc_title='Inheritance as a Linear Operator, and What Replaces It - Parts I to VI',
           running_header='INHERITANCE AS A LINEAR OPERATOR, AND WHAT REPLACES IT',
           cover_tag='GENETICS  /  COMPLETE',
           footer='Preprint - not peer reviewed')
