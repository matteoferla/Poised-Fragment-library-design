"""
After the third pass, the chunks are combined here.
And saved as SMILES.
"""

import bz2
import io
import os
import sys
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, PandasTools
import pandas as pd

# ------------------------------------------------------------------------------

keepers = pd.DataFrame()
dfs = []
info = {}
totals = {}
tier = 0

# ------------------------------------------------------------------------------

tier += 1
columns = ['SMILES', 'identifier',  'HAC', 'boringness',
       'synthon_score', 'pip_common_mean', 'pip_uncommon_mean',
       'combined_Zscore']
filename = f'selected_final/shorlist{tier:0>4}.1M.cxsmiles.bz2'
print(filename)
bfh = bz2.open(filename, 'wt')
bfh.write('\t'.join(columns)+'\n')
#
for group in ('Z1', 'Z08-1','Z05-08', 'Z0-05'):
    print(group)
    for path in Path(f'third_pass/{group}').glob('*.bz2'):
        with Chem.ForwardSDMolSupplier(bz2.open(path, 'rb')) as sdfh:
            mols = list(sdfh)
            n = len(mols)
            info[(group, path.name)] = n
            # 'mol': mol,
            df = pd.DataFrame([{'filename': path.name, 'identifier': mol.GetProp('_Name'),  **mol.GetPropsAsDict()} for mol in mols])
            dfs.append(df)
    # end of chunk reading...
    dfs = pd.concat(dfs, ignore_index=True).sort_values('combined_Zscore', ascending=False).reset_index(drop=True).copy()
    totals[group] = len(dfs)
    maxima = (totals[group] // 1_000_000) * 1_000_000
    print(totals)
    # save
    for i, row in dfs.iterrows():
        if i == maxima:
            keepers = dfs.iloc[i:]
            dfs = [keepers]
            break
        bfh.write('\t'.join(row[columns].astype(str))+'\n')
        if i % 1_000_000 == 0 and i != 0:
            bfh.close()
            tier += 1
            filename = f'selected_final/shorlist{tier:0>4}.1M.cxsmiles.bz2'
            print(filename)
            bfh = bz2.open(filename, 'wt')
            bfh.write('\t'.join(columns)+'\n')