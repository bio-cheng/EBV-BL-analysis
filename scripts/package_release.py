#!/usr/bin/env python3
"""Package only reviewed code/config/docs, never inferred patient data or results."""
import argparse
import hashlib
from pathlib import Path
import re
import zipfile
ROOT=Path(__file__).resolve().parents[1]
TOP={'.gitignore','README.md','LICENSE_NOTICE.md','requirements.txt','run.py'}
DIRS={'src':{'.py'},'scripts':{'.py','.R'},'tests':{'.py'},'docs':{'.md','.json'},'config':{'.json','.tsv'},'.github':{'.yml'}}
CONFIG={'analysis.json','R-packages.tsv','mast_gene_sets.tsv','cd8_gene_sets.tsv'}
def selected():
    result=[]
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file():continue
        rel=p.relative_to(ROOT)
        if '__pycache__' in rel.parts:continue
        keep=(len(rel.parts)==1 and p.name in TOP) or str(rel)=='data/README.md'
        if len(rel.parts)>1 and rel.parts[0] in DIRS and p.suffix in DIRS[rel.parts[0]]:
            keep=True
            if rel.parts[0]=='config' and p.name not in CONFIG:raise ValueError('Unreviewed config file: '+str(rel))
        if not keep:continue
        if p.is_symlink():raise ValueError('Symlinks must not be packaged')
        content=p.read_text()
        patterns=[r'/mnt/sda\d/','/ho' + r'me/[^/]+/',r'(?i)EBV-?\d+\b','-----BE' + 'GIN ' + r'.*PRIVATE KEY-----',r'gh[pousr]_[A-Za-z0-9]{20,}']
        if any(re.search(q,content) for q in patterns):raise ValueError('Potential private content in '+str(rel))
        if p.stat().st_size>2_000_000:raise ValueError('Unexpectedly large source file '+str(rel))
        result.append(p)
    return result
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path);p.add_argument('--check-only',action='store_true');a=p.parse_args()
    files=selected();print(f'Privacy-pattern and allowlist checks passed: {len(files)} files. Manual review remains required.')
    if a.check_only:return
    if not a.output:p.error('--output is required unless --check-only')
    dest=a.output.resolve()
    if dest.exists():raise ValueError('Archive already exists; use a new filename')
    if dest==ROOT or ROOT in dest.parents:raise ValueError('Place release archive outside source tree')
    with zipfile.ZipFile(dest,'x',compression=zipfile.ZIP_DEFLATED) as z:
        for f in files:z.write(f,Path(ROOT.name)/f.relative_to(ROOT))
    print('SHA256 '+hashlib.sha256(dest.read_bytes()).hexdigest())
if __name__=='__main__':main()
