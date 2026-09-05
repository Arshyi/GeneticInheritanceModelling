"""Bundle the paper and its reproducible evidence. One PDF only; no runtimes."""
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import zipfile

ROOT=Path(__file__).resolve().parents[1]


def main():
    excluded={'.deps','.venv','__pycache__','.pytest_cache'}
    # Only the current paper ships. The Version I PDF/DOCX are the author's and are
    # not redistributed; the superseded Version II PDF is not shipped either.
    withheld={'Genetics_Version_II.pdf'}
    files=[p for p in ROOT.rglob('*') if p.is_file() and not excluded.intersection(p.parts)
           and p.suffix not in ('.zip','.pyc') and p.name!='package_manifest.json'
           and p.name not in withheld]
    output=ROOT/'output/Genetics_Research_Package.zip'
    with zipfile.ZipFile(output,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
        for path in sorted(files):
            archive.write(path,'GENETICS/'+str(path.relative_to(ROOT.parent)).replace('\\','/'))
    with zipfile.ZipFile(output) as archive:
        assert archive.testzip() is None
        assert len(archive.infolist())==len(files)
        names={info.filename for info in archive.infolist()}
        assert 'GENETICS/version2/output/pdf/Genetics_Complete.pdf' in names
        assert not [n for n in names if n.endswith('.pdf') and not n.endswith('Genetics_Complete.pdf')]
        assert not [n for n in names if n.endswith('.docx')]
    result={'created_at_utc':datetime.now(timezone.utc).isoformat(),'archive':str(output.relative_to(ROOT)),
            'files':len(files),'bytes':output.stat().st_size,'sha256':hashlib.sha256(output.read_bytes()).hexdigest(),
            'zip_crc_validation':'passed','single_pdf_only':True,'version1_originals_included':False,'dependencies_included':False}
    (ROOT/'results/package_manifest.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
