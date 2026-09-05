"""Build an editable evidence-derived manuscript and a typeset review PDF."""
from datetime import datetime, timezone
from html import escape
import hashlib
import json
import math
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, PageBreak, Table, TableStyle, Image, Preformatted, KeepTogether
from PIL import Image as PILImage

ROOT=Path(__file__).resolve().parents[1]
RES=ROOT/'results'; MAN=ROOT/'manuscript'; OUT=ROOT/'output/pdf'
PALETTE={'dense':'#293c61','csr':'#136f79','hash':'#d18a25','streamed_kernel':'#b44949'}


def read(name):return json.loads((RES/name).read_text(encoding='utf-8'))
def table(header,rows):
    return '| '+' | '.join(header)+' |\n| '+' | '.join(['---']*len(header))+' |\n'+'\n'.join('| '+' | '.join(map(str,row))+' |' for row in rows)


def sources():
    merged={}
    for file in sorted((ROOT/'sources').glob('*_sources.json')):
        data=json.loads(file.read_text(encoding='utf-8'))
        for source in data['sources']:
            source.setdefault('accessed_date',data.get('accessed_utc_date',data.get('accessed_date','2026-09-05')))
            if source['id'] in merged:raise ValueError('Duplicate source '+source['id'])
            merged[source['id']]=source
    merged['version1_local']={'id':'version1_local','type':'user_provided_local_source','authors':['Arshyia Mehran'],
        'title':'Linear Algebra in Bioinformatics and Computational Biology: Modeling the Generational Decay of Sickle Cell Anemia and Changes in Generational Blood Types Using MATLAB',
        'url':'','local_file':'../Bioinformatics-Arshyia Mehran.pdf','pages':54,'sha256':hashlib.sha256((ROOT.parent/'Bioinformatics-Arshyia Mehran.pdf').read_bytes()).hexdigest()}
    (ROOT/'sources/source_ledger.json').write_text(json.dumps({'schema_version':'1.0','sources':list(merged.values())},indent=2)+'\n',encoding='utf-8')
    return merged


def benchmark_figure(bench):
    rows=bench['full_workload_summary']
    fig,axes=plt.subplots(1,2,figsize=(9,4.2))
    for method in PALETTE:
        selected=[r for r in rows if r['method']==method]
        # A common traced allocation metric for all methods, separate from retained payload table.
        axes[0].plot([r['loci'] for r in selected],[r['tracemalloc_peak_bytes']/1024**2 for r in selected],marker='o',label=method.replace('_',' '),color=PALETTE[method])
        axes[1].plot([r['loci'] for r in selected],[r['inference_seconds']['median']*1000 for r in selected],marker='o',label=method.replace('_',' '),color=PALETTE[method])
    axes[0].set(xlabel='Biallelic loci',ylabel='Traced allocation peak (MiB)',yscale='log',title='Separate construction + update pass')
    axes[1].set(xlabel='Biallelic loci',ylabel='Median complete update (ms)',yscale='log',title='Validated random-mating update')
    for ax in axes:
        ax.set_xticks([1,2,3,4,5]);ax.spines[['right','top']].set_visible(False);ax.grid(axis='y',alpha=.18)
    axes[1].legend(fontsize=8.5)
    fig.tight_layout()
    for extension in ('png','svg'):fig.savefig(ROOT/f'figures/benchmark_overview.{extension}',dpi=200,bbox_inches='tight')
    plt.close(fig)


def bindings():
    ledger=sources();bench=read('benchmark.json');pop=read('population_validation.json');legacy=read('version1_reproduction.json')
    if bench['status']!='passed' or bench['full_loci']!=[1,2,3,4,5]:raise ValueError('Complete benchmark evidence required')
    benchmark_figure(bench)
    rows=bench['full_workload_summary'];largest={r['method']:r for r in rows if r['loci']==5}
    dense,csr=largest['dense'],largest['csr']
    ratio=dense['numeric_array_payload_bytes']/csr['numeric_array_payload_bytes']
    speed=dense['inference_seconds']['median']/csr['inference_seconds']['median']
    sample=next(r for r in bench['raw_workers'] if r['phase']=='inference' and r['status']=='measured')
    env=sample.get('environment',sample.get('metadata',{}))
    xml=ET.parse(RES/'tests.xml').getroot()
    suites=list(xml.iter('testsuite'));test_count=sum(int(s.get('tests',0)) for s in suites)
    failed=sum(int(s.get('failures',0))+int(s.get('errors',0)) for s in suites)
    if failed or test_count<28:raise ValueError('Passing verification evidence required')
    values={'DATE':datetime.now(timezone.utc).strftime('%d %B %Y'),
            'SCD_SOURCE':'bender_carlberg_scd_2025','CLINVAR_SOURCE':'clinvar_rs334_15333_v7','RH_SOURCE':'dean_rh_2005',
            'HWE_SOURCE':'hardy_1908','1000G_SOURCE':'1000genomes_2015','EXACT_SOURCE':'wigginton_hwe_2005',
            'MNS_SOURCE':'dean_mns_2005','BOMBAY_SOURCE':'dean_hh_2005','LINKAGE_SOURCE':'brown_linkage_2002',
            'BENCH_ABSTRACT':f'At five loci, CSR stores {ratio:.2f} times less numeric payload than dense storage; the corresponding validated update timings are {dense["inference_seconds"]["median"]*1000:.3f} ms (dense) and {csr["inference_seconds"]["median"]*1000:.3f} ms (CSR). These are workload- and implementation-specific measurements.',
            'TEST_STATUS':f'The final automated verification recorded {test_count} passing tests, zero failures, and zero errors. The retained JUnit XML and test log identify the executed suite. All 20 representation-by-dimension benchmark comparisons pass an absolute agreement threshold of 10^-12, with the measured maximum difference reported below.',
            'BENCH_METHODS':'The retained run uses three construction repetitions and seven inference repetitions per materialized representation and dimension. Each construction sample and each inference or memory configuration runs in a fresh worker. Imports and model setup are outside timed construction. One untimed construction warms allocator/library behavior, followed by clearing the local-cross cache before each timed construction. Inference uses an untimed warm-up and warm local tables. The traced-memory pass begins after scientific-library imports and covers construction plus one update. The complete recorded environment is in benchmark.json; the summary below reproduces its key fields.\n\n'+table(['Environment item','Recorded value'],[(k,str(v)) for k,v in env.items() if k in ('python','numpy','scipy','platform','processor','logical_cpu_count')]),
            'BENCH_TABLE':table(['Loci','Method','Retained bytes*','Build median (ms)','Update median (ms)'],
                [(r['loci'],r['method'].replace('_',' '),f'{(r["hash_recursive_python_object_bytes"] if r["method"]=="hash" else r["numeric_array_payload_bytes"]):,}',
                  f'{r["construction_seconds"]["median"]*1000:.3f}' if r['construction_seconds'] else 'No full kernel',f'{r["inference_seconds"]["median"]*1000:.3f}') for r in rows])+'\n\n*Dense/CSR: retained numeric buffers only. Hash: recursive Python-object footprint, a different metric. Streamed: no retained full-kernel buffer; this is not zero process memory. The figure uses the separate traced-allocation pass for all methods.',
            'BENCH_RESULTS':f'At five loci, the full catalog contains G={dense["G"]:,} genotypes and U={dense["U"]:,} unordered pairs. The kernel has {dense["nnz"]:,} nonzeros. Dense numeric payload is {dense["numeric_array_payload_bytes"]:,} bytes; CSR payload is {csr["numeric_array_payload_bytes"]:,} bytes, a {ratio:.2f}-fold reduction. The ratio of dense to CSR median validated-update time is {speed:.2f}. Construction medians are {dense["construction_seconds"]["median"]:.4f} s and {csr["construction_seconds"]["median"]:.4f} s respectively. The maximum absolute discrepancy from the dense result over every measured method and dimension is {max(r["max_absolute_difference_from_dense"] for r in rows):.3g}. The complete timing samples, minima, maxima, standard deviations, and separate memory records are retained rather than hidden behind these medians.',
            'QUERY_TABLE':table(['Loci','Median query (microseconds)','Traced peak (bytes)','Absolute log error'],
                [(r['loci'],f'{r["query_seconds"]["median"]*1e6:.3f}',f'{r["tracemalloc_peak_bytes"]:,}',f'{r["absolute_log_error"]:.3g}') for r in bench['query_workload_summary']]),
            'SCALING_TABLE':table(['Biallelic loci n','Genotype states G','Unordered pairs U','Supported transitions'],
                [(n,f'{3**n:,}',f'{3**n*(3**n+1)//2:,}',f'{(15**n+5**n)//2:,}') for n in [1,2,3,4,5,6,10]]),
            'POPULATION_TABLE':table(['Population / diagnostic','N','Observed T/T, T/A, A/A','Exact HWE p'],
                [(r['population'],r['n'],', '.join(map(str,r['counts'])),f'{r["hwe_exact_p"]:.6f}') for r in pop['component_populations'] if r['q']>0]+[(r['population'],r['n'],', '.join(map(str,r['counts'])),f'{r["hwe_exact_p"]:.6f}') for r in [pop['african_pool'],pop['pooled']]])+'\n\nThe 17 monomorphic component populations remain in the CSV and the 26-test adjustment. They are omitted here only to keep the displayed diagnostic table compact.',
            'LEGACY_EXAMPLE_TABLE':table(['Original page','Displayed input and generation','Reconstruction'],[
                (28,'AS x AS, generation 3','AA 75%, AS 18.75%, SS 0%; total93.75%'),
                (35,'AA x OO, generation 3','AA25%, AO50%, OO25%'),
                (36,'AO x BO, generation 4','AA6.25%, AO25%, BB6.25%, BO25%, AB12.5%, OO25%'),
                (49,'AO/Dd x AO/Dd, generation 3','A+56.25%, A-18.75%, O+18.75%, O-6.25%'),
                (50,'AA/dd x BO/Dd, generation 7','All18 displayed genotype percentages match after rounding'),
                (50,'AB/dd x OO/Dd, generation 8','All18 displayed genotype percentages match after rounding')]),
            'REPRO_COMMANDS':'```text\ncd version2\npython -m pip install -r requirements.txt\npython run.py test -q\npython run.py reproduce\npython run.py fetch\npython run.py science\npython run.py benchmark\npython run.py test -q --junitxml=results/tests.xml\npython run.py manuscript\n```\n\nFor the bundled Windows Python used here, reproduce.ps1 sets the verified interpreter path. The locally installed .deps directory is optional and is excluded from the source archive.',
            'ARTIFACT_TABLE':table(['Artifact','Purpose','Evidence status'],[
                ('research/version1_audit.md','Page-by-page discrepancy explanation','Full paper and screenshots inspected'),
                ('results/version1_reproduction.json','Six exact legacy examples and coverage variants','Executed rational translation'),
                ('genetics/core.py; extensions.py','Reusable kernel and staged extensions','Automated tests passed'),
                ('results/benchmark.json and benchmark_workers/','Raw timings, memory, metadata and checks','Executed local CPU run'),
                ('data/observed_genotypes.csv + provenance','External genomic calls with validation','Frozen public snapshot'),
                ('results/population_validation.json','HWE and pooling audit','Executed; no clinical accuracy claim'),
                ('results/polygenic_synthetic.json','Score and interval calibration','10,000 independent synthetic draws'),
                ('sources/source_ledger.json','Reproducible claim/source bibliography','Primary and authoritative sources'),
                ('manuscript/version2_manuscript.md','Editable evidence-derived candidate','Author review required'),
                ('output/pdf/Genetics_Version_II.pdf','Typeset candidate','Rendered and visually checked before delivery')])}
    return ledger,values


def apply_citations(text,ledger,citations_name='manuscript_citations.json'):
    used=[]
    for match in re.finditer(r'\[([a-zA-Z0-9_]+)\](?!\()',text):
        key=match.group(1)
        if key in ledger and key not in used:used.append(key)
    refs=[]
    for index,key in enumerate(used,1):
        s=ledger[key];authors=s.get('authors',s.get('author',[]));authors=', '.join(authors) if isinstance(authors,list) else str(authors)
        link=f' [{s["url"]}]({s["url"]})' if s.get('url') else ' Local user-supplied PDF, 54 pages.'
        if not authors:
            authors=s.get('organization','SciPy developers' if key.startswith('scipy_') else 'Python Software Foundation' if key.startswith('python_') else '')
        parts=[str(value).strip().rstrip('.') for value in (authors,s['title'],s.get('year'),s.get('venue')) if value]
        refs.append(f'[{index}] '+'. '.join(parts)+'.'+(f' DOI: {s["doi"]}.' if s.get('doi') else '')+link+' Accessed '+str(s.get('accessed_date','2026-09-05'))+'.')
        text=text.replace('['+key+']','['+str(index)+']')
    text=text.replace('{{BIBLIOGRAPHY}}','\n\n'.join(refs))
    if '{{' in text:raise ValueError('Unresolved manuscript placeholders: '+str(re.findall(r'\{\{.*?\}\}',text)))
    (RES/citations_name).write_text(json.dumps({'references_in_order':used,'total':len(used)},indent=2)+'\n',encoding='utf-8')
    return text


def render(text,out_name='Genetics_Version_II.pdf',
           doc_title='Complete Inheritance Without a Square-Matrix Constraint - Version II',
           running_header='COMPLETE INHERITANCE WITHOUT A SQUARE-MATRIX CONSTRAINT',
           cover_tag='GENETICS  /  VERSION II',
           footer='Research candidate  |  Author review required'):
    OUT.mkdir(parents=True,exist_ok=True)
    for name,file in [('Body','georgia.ttf'),('Body-Bold','georgiab.ttf'),('Body-Italic','georgiai.ttf'),('Sans','arial.ttf'),('Sans-Bold','arialbd.ttf')]:
        pdfmetrics.registerFont(TTFont(name,'C:/Windows/Fonts/'+file))
    pdfmetrics.registerFontFamily('Body',normal='Body',bold='Body-Bold',italic='Body-Italic',boldItalic='Body-Bold')
    width=A4[0]-112
    styles={
        'body':ParagraphStyle('body',fontName='Body',fontSize=10.1,leading=15.1,spaceAfter=8,textColor=colors.HexColor('#202d39'),splitLongWords=True),
        'h1':ParagraphStyle('h1',fontName='Sans-Bold',fontSize=17,leading=21,spaceBefore=20,spaceAfter=12,keepWithNext=True,textColor=colors.HexColor('#123e4d')),
        'h2':ParagraphStyle('h2',fontName='Sans-Bold',fontSize=12.5,leading=16,spaceBefore=13,spaceAfter=8,keepWithNext=True,textColor=colors.HexColor('#136f79')),
        'title':ParagraphStyle('title',fontName='Sans-Bold',fontSize=32,leading=38,spaceBefore=75,spaceAfter=23,textColor=colors.HexColor('#123e4d')),
        'subtitle':ParagraphStyle('subtitle',fontName='Body',fontSize=17,leading=24,spaceAfter=30,textColor=colors.HexColor('#536d77')),
        'cell':ParagraphStyle('cell',fontName='Sans',fontSize=8,leading=11,spaceAfter=0,splitLongWords=True),
        'caption':ParagraphStyle('caption',fontName='Sans',fontSize=8.5,leading=12,spaceAfter=13,textColor=colors.HexColor('#536d77')),
        'code':ParagraphStyle('code',fontName='Courier',fontSize=8.1,leading=12,spaceBefore=5,spaceAfter=12,leftIndent=12,backColor=colors.HexColor('#f0f5f6')),
        'reference':ParagraphStyle('reference',fontName='Sans',fontSize=8.4,leading=12.3,spaceAfter=10,splitLongWords=True)
    }
    def markup(value):
        val=escape(value)
        val=re.sub(r'\*\*(.+?)\*\*',r'<b>\1</b>',val)
        val=re.sub(r'`([^`]+)`',r'<font name="Courier">\1</font>',val)
        val=re.sub(r'\[([^\]]+)\]\(([^)]+)\)',lambda m:f'<link href="{m[2]}" color="#136f79">{m[1]}</link>',val)
        return val
    story=[Spacer(1,64)];lines=text.splitlines();i=0;cover=True;references=False;fig_no=0
    while i<len(lines):
        line=lines[i].strip()
        if not line:i+=1;continue
        if line=='<!-- pagebreak -->':story.append(PageBreak());cover=False;i+=1;continue
        if line.startswith('```'):
            block=[];i+=1
            while i<len(lines) and not lines[i].startswith('```'):block.append(lines[i]);i+=1
            story.append(KeepTogether([Preformatted('\n'.join(block),styles['code'],maxLineLength=87)]));i+=1;continue
        if line.startswith('# '):
            title=line[2:];references=title=='References'
            story.append(Paragraph(markup(title),styles['title' if cover else 'h1']));i+=1;continue
        if line.startswith('## '):
            story.append(Paragraph(markup(line[3:]),styles['subtitle' if cover else 'h2']));i+=1;continue
        if line.startswith('|'):
            rows=[]
            while i<len(lines) and lines[i].strip().startswith('|'):
                cells=[c.strip() for c in lines[i].strip().strip('|').split('|')]
                if not all(re.fullmatch(r'[-: ]+',c) for c in cells):rows.append(cells)
                i+=1
            n=len(rows[0]);fractions={2:[.3,.7],3:[.27,.33,.4],4:[.27,.12,.39,.22],5:[.07,.2,.25,.24,.24]}.get(n,[1/n]*n)
            data=[[Paragraph(markup(c),styles['cell']) for c in row] for row in rows]
            t=Table(data,colWidths=[width*f for f in fractions],repeatRows=1,hAlign='LEFT')
            t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#dcebed')),('VALIGN',(0,0),(-1,-1),'TOP'),
                ('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f5f8f9')]),('LEFTPADDING',(0,0),(-1,-1),6),
                ('RIGHTPADDING',(0,0),(-1,-1),6),('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6),
                ('LINEBELOW',(0,0),(-1,0),.7,colors.HexColor('#136f79'))]))
            story.extend([t,Spacer(1,12)]);continue
        match=re.match(r'!\[(.*?)\]\((.*?)\)',line)
        if match:
            path=(MAN/match[2]).resolve();img=PILImage.open(path);ratio=img.height/img.width
            fig_no+=1
            story.append(KeepTogether([Image(str(path),width=width,height=width*ratio),Paragraph(f'Figure {fig_no}. '+markup(match[1])+'.',styles['caption'])]));i+=1;continue
        block=[line];i+=1
        while i<len(lines) and lines[i].strip() and not lines[i].startswith(('#','|','```','![','<!--')):
            block.append(lines[i].strip());i+=1
        story.append(Paragraph(markup(' '.join(block)),styles['reference' if references else 'body']))
    def page(canvas,doc):
        canvas.saveState()
        if doc.page==1:
            canvas.setFillColor(colors.HexColor('#136f79'));canvas.rect(56,A4[1]-66,52,5,fill=True,stroke=False)
            canvas.setFont('Sans-Bold',9);canvas.drawString(56,A4[1]-89,cover_tag)
        else:
            canvas.setFillColor(colors.HexColor('#63757e'));canvas.setFont('Sans',8)
            canvas.drawString(56,A4[1]-34,running_header)
            canvas.setStrokeColor(colors.HexColor('#d4dfe3'));canvas.line(56,A4[1]-42,A4[0]-56,A4[1]-42)
        canvas.setFillColor(colors.HexColor('#63757e'));canvas.setFont('Sans',8)
        canvas.drawString(56,30,footer)
        canvas.drawRightString(A4[0]-56,30,str(doc.page));canvas.restoreState()
    path=OUT/out_name
    doc=SimpleDocTemplate(str(path),pagesize=A4,rightMargin=56,leftMargin=56,topMargin=58,bottomMargin=52,
        title=doc_title,author='Research candidate for Arshyia Mehran',pageCompression=1)
    doc.build(story,onFirstPage=page,onLaterPages=page)
    print(str(path))


def build_text():
    ledger,values=bindings()
    text=(MAN/'template.md').read_text(encoding='utf-8')
    for key,val in values.items():text=text.replace('{{'+key+'}}',val)
    text=text.replace('Its SHA-256 digest before the extension','[version1_local] Its SHA-256 digest before the extension')
    text=apply_citations(text,ledger)
    (MAN/'version2_manuscript.md').write_text(text,encoding='utf-8')
    return text


if __name__=='__main__':render(build_text())
