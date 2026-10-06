"""Reproduce the moment identities in an isolated extracted repository layout."""
from moments import *
import zipfile,subprocess
def main():
    archive=OUT/'NETWORK_MOMENT_CODESIGN_20261004_COMPLETE.zip'
    dest=OUT/'qa/smoke_workspace';dest.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(archive) as z:
        for name in z.namelist():
            target=(dest/name).resolve()
            assert target.is_relative_to(dest.resolve()),name
        z.extractall(dest)
    rel='experiments/network_moment_codesign_20261004'
    manifest=json.loads((dest/rel/'BUNDLE_MANIFEST.json').read_text())
    for name,digest in manifest['files'].items():assert hashlib.sha256((dest/name).read_bytes()).hexdigest()==digest,name
    run=subprocess.run([sys.executable,str(dest/rel/'moments.py')],cwd=dest,capture_output=True,text=True)
    (OUT/'qa/smoke_python.log').write_text(run.stdout+'\n'+run.stderr,encoding='utf8')
    assert run.returncode==0,run.stderr
    got=pd.read_csv(dest/rel/'TABLE_03_MOMENT_IDENTITIES.csv')
    reference=pd.read_csv(OUT/'TABLE_03_MOMENT_IDENTITIES.csv')
    assert got.case.tolist()==reference.case.tolist()
    assert np.allclose(got.kappa,reference.kappa,rtol=1e-10,atol=1e-15)
    assert got.coefficient_relative_error.max()<1e-4
    dump('BUNDLE_SMOKE_TEST.json',{'extracted_file_hashes_pass':True,'moment_script_exit':run.returncode,
        'isolated_repository_layout':True,'moment_case_count':len(got),
        'max_coefficient_relative_error':float(got.coefficient_relative_error.max()),
        'scope':'Python coefficient run, not a repeated expensive nonlinear campaign'})
    print('EXTRACTED BUNDLE SMOKE PASS',len(got),'cases')
if __name__=='__main__':main()
