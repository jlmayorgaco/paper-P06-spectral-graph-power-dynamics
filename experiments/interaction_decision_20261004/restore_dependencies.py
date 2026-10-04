"""Restore missing byte-identical sources/data; never overwrite differences."""
from pathlib import Path
import zipfile,hashlib,json
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
def main():
    manifest=json.loads((OUT/'DEPENDENCY_MANIFEST.json').read_text())
    with zipfile.ZipFile(OUT/'DEPENDENCIES.zip') as z:
        for item in z.infolist():
            target=(ROOT/item.filename).resolve()
            if not target.is_relative_to(ROOT.resolve()):raise RuntimeError('path outside root')
            data=z.read(item.filename)
            assert hashlib.sha256(data).hexdigest()==manifest['files'][item.filename]
            if target.exists():
                if target.read_bytes()!=data:raise RuntimeError('Refuse to overwrite changed input: '+str(target))
            else:
                target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
    print('Dependencies verified/restored without overwriting differences')
if __name__=='__main__':main()
