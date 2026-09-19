import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'scripts')]
from analysis import pfs_group,boolean,bh,load_counts,logcpm,prepare
from make_demo import generate

class Contracts(unittest.TestCase):
    def test_censoring_and_boundary(self):
        self.assertEqual(pfs_group(5,1),'Short')
        self.assertEqual(pfs_group(5,0),'Early')
        self.assertEqual(pfs_group(7,0),'Long')
        self.assertEqual(pfs_group(6,1),'Boundary')
        self.assertEqual(pfs_group(np.nan,0),'Missing')
    def test_boolean_rejects_ambiguous(self):
        with self.assertRaises(ValueError):boolean(pd.Series(['yes']))
        self.assertEqual(boolean(pd.Series(['true','0'])).tolist(),[True,False])
    def test_bh_missing(self):
        self.assertTrue(np.isnan(bh([.01,np.nan,.02])[1]))
        np.testing.assert_allclose(bh([.01,.02]),[.02,.02])
    def test_logcpm_full_library(self):
        x=pd.DataFrame({'A':[1,3]})
        np.testing.assert_allclose(logcpm(x).A,np.log2([250001,750001]))
    def test_reject_fractional_counts(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/'counts.tsv';pd.DataFrame({'gene':['A'],'sample':[.5]}).to_csv(p,sep='\t',index=False)
            with self.assertRaises(ValueError):load_counts(p)
    def test_synthetic_preparation(self):
        with tempfile.TemporaryDirectory() as t:
            t=Path(t);generate(t/'inputs');out=t/'outputs';out.mkdir()
            cfg=json.loads((ROOT/'config/analysis.json').read_text())
            c=prepare(t/'inputs',out,cfg,ROOT/'config')
            self.assertEqual(len(c),32)
            self.assertEqual(int(c.advanced.sum()),30)
            self.assertTrue((c.loc[c.mast_de_eligible,'pfs6'].isin(['Short','Long'])).all())
            clr=pd.read_csv(out/'clr.tsv',sep='\t',index_col=0)
            np.testing.assert_allclose(clr.sum(axis=1),0,atol=1e-10)
            raw=pd.read_csv(out/'raw_proportions.tsv',sep='\t',index_col=0)
            np.testing.assert_allclose(raw.sum(axis=1),1)
    def test_mast_gate_is_strict(self):
        with tempfile.TemporaryDirectory() as t:
            t=Path(t);generate(t/'inputs');out=t/'outputs';out.mkdir()
            p=t/'inputs/cell_counts.tsv';n=pd.read_csv(p,sep='\t');n.loc[1,'Myeloid_Mast_TPSAB1']=5;n.to_csv(p,sep='\t',index=False)
            c=prepare(t/'inputs',out,json.loads((ROOT/'config/analysis.json').read_text()),ROOT/'config')
            self.assertFalse(c.loc['SYNTHETIC_001','mast_eligible'])

if __name__=='__main__':unittest.main()
