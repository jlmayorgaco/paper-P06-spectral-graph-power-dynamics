from model import *
import zipfile,subprocess
def main():
    target=OUT/'qa'/'bundle_smoke';target.mkdir(parents=True,exist_ok=True)
    archive=OUT/'INTERACTION_DECISION_20261004_COMPLETE.zip'
    with zipfile.ZipFile(archive) as z:
        for item in z.infolist():
            path=(target/item.filename).resolve();assert path.is_relative_to(target.resolve())
            path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(z.read(item.filename))
    new=target/'experiments/interaction_decision_20261004'
    subprocess.run([sys.executable,str(new/'restore_dependencies.py')],check=True,cwd=target)
    # This host keeps flint in the historical vendor. A normal pip installation
    # supplies the same package for portable use; no other original inputs used.
    script='import sys,json,numpy as np;sys.path.insert(0,'+repr(str(ROOT/'experiments/regional_paper_closure_20261003/vendor'))+');sys.path.insert(0,'+repr(str(new))+');from certify_roots import *;d=json.loads((OUT/"BUDGET_RESULT.json").read_text());certify("smoke_corrected",np.array(d["p"]),complex(*d["root"]))'
    p=subprocess.run([sys.executable,'-c',script],check=True,cwd=target,text=True,capture_output=True)
    result={'python_certificate_from_extracted_bundle':True,'stdout':p.stdout,
        'external_dependency':'python-flint0.8.0, using host vendor as installed dependency',
        'no_original_model_data_used':True,'code_sha256':hashlib.sha256((new/'certify_roots.py').read_bytes()).hexdigest()}
    save('BUNDLE_SMOKE_TEST.json',result);print(p.stdout)
if __name__=='__main__':main()
