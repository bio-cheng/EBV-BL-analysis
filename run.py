#!/usr/bin/env python3
"""Portable entry point. Outputs may contain sensitive patient-level data."""
import argparse
import hashlib
import importlib.metadata
import json
import platform
from pathlib import Path
import subprocess
import sys
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'src'))
from analysis import prepare, associations

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-dir', type=Path, required=True)
    p.add_argument('--out-dir', type=Path, required=True)
    p.add_argument('--config', type=Path, default=ROOT/'config/analysis.json')
    p.add_argument('--rscript', default='Rscript')
    p.add_argument('--stage', choices=['validate', 'python', 'mast', 'survival', 'plot', 'all'], default='all')
    p.add_argument('--allow-existing-output', action='store_true')
    a = p.parse_args()
    a.data_dir=a.data_dir.resolve(); a.out_dir=a.out_dir.resolve()
    if a.data_dir == a.out_dir or a.data_dir in a.out_dir.parents or a.out_dir in a.data_dir.parents:
        p.error('Data and output directories must be disjoint; never overwrite inputs.')
    if a.out_dir.exists() and any(a.out_dir.iterdir()) and not a.allow_existing_output:
        p.error('Output is nonempty. Use a new directory or --allow-existing-output explicitly.')
    cfg=json.loads(a.config.read_text())
    a.out_dir.mkdir(parents=True, exist_ok=True)
    inputs=sorted(x for x in a.data_dir.iterdir() if x.is_file())
    manifest={'started_utc':datetime.now(timezone.utc).isoformat(),'status':'running',
              'config':cfg,'python':platform.python_version(),'stage':a.stage,
              'packages':{n:importlib.metadata.version(n) for n in ['numpy','pandas','scipy','statsmodels','matplotlib']},
              'inputs':{x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in inputs},
              'code':{str(x.relative_to(ROOT)):hashlib.sha256(x.read_bytes()).hexdigest() for x in sorted(ROOT.rglob('*'))
                      if x.is_file() and x.suffix in ['.py','.R','.json'] and not any(v in x.parts for v in ['data','results','private','__pycache__'])}}
    mf=a.out_dir/'run_manifest.json'
    mf.write_text(json.dumps(manifest,indent=2))
    try:
        prepare(a.data_dir,a.out_dir,cfg,ROOT/'config')
        if a.stage in ['python','all']:
            associations(a.data_dir,a.out_dir,cfg)
        for stage in ['mast','survival']:
            if a.stage in [stage,'all']:
                command=[a.rscript,str(ROOT/'scripts'/f'{stage}.R'),str(a.data_dir),str(a.out_dir),str(ROOT/'config')]
                proc=subprocess.run(command,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
                (a.out_dir/f'{stage}.log').write_text(proc.stdout)
                if proc.returncode:raise RuntimeError(f'{stage} failed; inspect {stage}.log')
        if a.stage in ['plot','all']:
            from plots import render
            render(a.out_dir)
        manifest['status']='complete'
    except Exception as e:
        manifest['status']='failed';manifest['error']=str(e)
        raise
    finally:
        manifest['finished_utc']=datetime.now(timezone.utc).isoformat()
        mf.write_text(json.dumps(manifest,indent=2))
    print(f'Completed {a.stage}. Outputs are LOCAL and may contain sensitive data: {a.out_dir}')

if __name__=='__main__':main()
