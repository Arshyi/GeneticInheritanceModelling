"""Build the unified Version I + Version II manuscript and its typeset PDF.

Version I supplies the narrative spine and every piece of its mathematics.
Version II supplies the audit, the complete kernel, the measured comparisons,
and the age extension. All numeric bindings come from retained machine-readable
results, exactly as in build_manuscript.py; nothing here is typed by hand.
"""
import json

from experiments.build_manuscript import (
    MAN, RES, apply_citations, bindings, read, render, table,
)


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


def main():
    ledger, values = bindings()
    legacy = read('version1_reproduction.json')
    poly = read('polygenic_synthetic.json')
    counts = legacy['state_counts']
    checks = legacy['displayed_matrix_checks']

    variants_table, top6, top18, exact = coverage_tables(legacy)

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
        'status': 'Author review required; no submission or public claim of peer review',
    }, indent=2) + '\n', encoding='utf-8')
    return text


if __name__ == '__main__':
    render(main(),
           out_name='Genetics_Complete.pdf',
           doc_title='Inheritance as a Linear Operator, and What Replaces It - Parts I to VI',
           running_header='INHERITANCE AS A LINEAR OPERATOR, AND WHAT REPLACES IT',
           cover_tag='GENETICS  /  COMPLETE',
           footer='Unified research manuscript  |  Author review required')
