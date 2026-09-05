"""Bundle the preserved originals and reproducible research evidence, excluding runtimes."""
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import zipfile

ROOT=Path(__file__).resolve().parents[1]


def main():
    excluded={'.deps','.venv','__pycache__','.pytest_cache'}
    files=[ROOT.parent/'Bioinformatics-Arshyia Mehran.pdf',ROOT.parent/'Bioinformatics-Arshyia Mehran.docx']
    files += [p for p in ROOT.rglob('*') if p.is_file() and not excluded.intersection(p.parts)
              and p.suffix not in ('.zip','.pyc') and p.name!='package_manifest.json']
    output=ROOT/'output/Genetics_Version_II_Research_Package.zip'
    with zipfile.ZipFile(output,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
        for path in sorted(files):
            archive.write(path,'GENETICS/'+str(path.relative_to(ROOT.parent)).replace('\\','/'))
    with zipfile.ZipFile(output) as archive:
        assert archive.testzip() is None
        assert len(archive.infolist())==len(files)
        for original in ('Bioinformatics-Arshyia Mehran.pdf','Bioinformatics-Arshyia Mehran.docx'):
            assert hashlib.sha256(archive.read('GENETICS/'+original)).digest()==hashlib.sha256((ROOT.parent/original).read_bytes()).digest()
    result={'created_at_utc':datetime.now(timezone.utc).isoformat(),'archive':str(output.relative_to(ROOT)),
            'files':len(files),'bytes':output.stat().st_size,'sha256':hashlib.sha256(output.read_bytes()).hexdigest(),
            'zip_crc_validation':'passed','originals_match_live_preserved_files':True,'dependencies_included':False}
    (ROOT/'results/package_manifest.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
