"""Record artifact integrity and local PDF layout checks; never imply visual QA from text alone."""
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET
from pypdf import PdfReader
import pdfplumber

ROOT=Path(__file__).resolve().parents[1]
EXPECTED_PDF='ee3cb04ce28d4734669d5822fddbdf03c176fd5c2f5ff3760cb82a4bd1145ecd'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    source=ROOT.parent/'Bioinformatics-Arshyia Mehran.pdf'
    assert sha(source)==EXPECTED_PDF,'Original PDF changed'
    bench=json.loads((ROOT/'results/benchmark.json').read_text())
    assert bench['status']=='passed'
    assert bench['core_sha256']==sha(ROOT/'genetics/core.py'),'Core changed since benchmark'
    assert bench['harness_sha256']==sha(ROOT/'experiments/benchmark.py'),'Harness changed since benchmark'
    xml=ET.parse(ROOT/'results/tests.xml').getroot()
    suites=list(xml.iter('testsuite'))
    errors=sum(int(s.get('errors',0))+int(s.get('failures',0)) for s in suites)
    assert errors==0
    pdf=ROOT/'output/pdf/Genetics_Version_II.pdf';reader=PdfReader(pdf)
    text='\n'.join(p.extract_text() or '' for p in reader.pages)
    for forbidden in ['{{','BENCH_ABSTRACT','BENCH_RESULTS','SCD_SOURCE','\ufffd']:
        assert forbidden not in text,forbidden
    problems=[];page_stats=[]
    with pdfplumber.open(pdf) as document:
        for n,page in enumerate(document.pages,1):
            chars=[c for c in page.chars if c.get('text','').strip()]
            outside=[c for c in chars if c['x0']<40 or c['x1']>page.width-40 or c['top']<18 or c['bottom']>page.height-18]
            if outside:problems.append({'page':n,'characters':''.join(c['text'] for c in outside)})
            page_stats.append({'page':n,'words':len(page.extract_words()),'images':len(page.images),
                               'nonempty_chars':len(chars)})
    assert not problems,problems
    paths=[p for p in ROOT.rglob('*') if p.is_file() and not any(s in p.parts for s in ('.deps','.venv','__pycache__','.pytest_cache','benchmark_workers'))
           and p.suffix not in ('.zip','.pyc') and p.name not in ('artifact_manifest.json','verification.json')]
    manifest={'created_at_utc':datetime.now(timezone.utc).isoformat(),'files':{str(p.relative_to(ROOT)).replace('\\','/'):{'sha256':sha(p),'bytes':p.stat().st_size} for p in sorted(paths)}}
    (ROOT/'results/artifact_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    visual_path=ROOT/'results/visual_review.json'
    visual=json.loads(visual_path.read_text()) if visual_path.exists() else {}
    visual_current=visual.get('pdf_sha256')==sha(pdf) and visual.get('status')=='passed'
    report={'verified_at_utc':datetime.now(timezone.utc).isoformat(),'original_pdf_unchanged':True,
            'original_pdf_sha256':sha(source),'original_docx_sha256':sha(ROOT.parent/'Bioinformatics-Arshyia Mehran.docx'),
            'benchmark_status':bench['status'],'core_and_harness_match_benchmarked_source':True,
            'automated_tests':sum(int(s.get('tests',0)) for s in suites),'errors_or_failures':errors,
            'manuscript_pages':len(reader.pages),'manuscript_sha256':sha(pdf),'pdf_text_boundary_check_passed':True,
            'placeholder_check_passed':True,'page_statistics':page_stats,
            'visual_review_current':visual_current,
            'visual_review':'See visual_review.json for recorded visual inspection; text checks alone do not certify layout.'}
    (ROOT/'results/verification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('page_statistics',)},indent=2))


if __name__=='__main__':main()
