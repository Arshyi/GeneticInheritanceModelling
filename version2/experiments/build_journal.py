"""Build the short journal version: the bounds result, standalone.

The long-form document is a monograph. This is the slice of it that stands on its
own without reference to the earlier model that motivated it: the rectangular
kernel, the exact combinatorics, the matching lower bounds, the measured
implementations, and the test against published polygenic scores.

It shares every binding with the long form, so the two cannot disagree about a
number. Editing the numbers means re-running the experiments, not editing prose.
"""
import json

from experiments.build_manuscript import MAN, RES, apply_citations, bindings, read, render
from experiments.build_unified import complexity_tables, pgs_tables


def main():
    ledger, values = bindings()
    comp = read('complexity.json')
    pgs = read('pgs_catalog.json')

    comp_traits, comp_growth = complexity_tables(comp)
    (pgs_summary, pgs_kernel, pgs_ladder, pgs_scaling,
     pgs_scaling_result, pgs_extrapolation) = pgs_tables(pgs)

    values.update({
        'COMPLEXITY_TRAIT_TABLE': comp_traits,
        'COMPLEXITY_GROWTH_TABLE': comp_growth,
        'COMPLEXITY_NOTE': comp['growth_comparison_note'],
        'PGS_TRAIT_SUMMARY': pgs_summary,
        'PGS_KERNEL_TABLE': pgs_kernel,
        'PGS_LADDER_TABLE': pgs_ladder,
        'PGS_SCALING_TABLE': pgs_scaling,
        'PGS_SCALING_RESULT': pgs_scaling_result,
        'PGS_EXTRAPOLATION': pgs_extrapolation,
    })

    text = (MAN / 'journal_template.md').read_text(encoding='utf-8')
    for key, val in values.items():
        text = text.replace('{{' + key + '}}', str(val))
    text = apply_citations(text, ledger, citations_name='journal_citations.json')
    (MAN / 'genetics_bounds_paper.md').write_text(text, encoding='utf-8')
    (RES / 'journal_manifest.json').write_text(json.dumps({
        'document': 'manuscript/genetics_bounds_paper.md',
        'pdf': 'output/pdf/Genetics_Bounds_Paper.pdf',
        'relationship': ('Short-form derivative of the long treatment archived at '
                         'doi:10.5281/zenodo.22401514; shares all numeric bindings with it.'),
        'status': 'Preprint; not peer reviewed. Prepared for journal submission.',
    }, indent=2) + '\n', encoding='utf-8')
    return text


if __name__ == '__main__':
    render(main(),
           out_name='Genetics_Bounds_Paper.pdf',
           doc_title='Exact Bounds for Multilocus Mendelian Transmission',
           running_header='EXACT BOUNDS FOR MULTILOCUS MENDELIAN TRANSMISSION',
           cover_tag='GENETICS  /  BOUNDS',
           footer='Preprint - not peer reviewed')
