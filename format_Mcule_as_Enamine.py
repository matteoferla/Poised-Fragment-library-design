"""
Mcule is in sdf format. For consistency & size it is converted to a similar format to Enamine REAL (cxsmiles)
with quick info used by library classification present.
Wasteful, but overcomes a lot of coding and Mcule is much smaller.
"""
import bz2, gzip
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from pathlib import Path

# --- INPUTS ---
input = ...  # path to Mcule sdf gzipped file
output = ...  # path to Mcule cxsmiles bz2 file
assert input.endswith('.sdf.gz')
assert output.endswith('.cxsmiles.bz2')
assert Path(input).exists()

# --- PROCESSING ---
headers = ['SMILES', 'Identifier', 'HBA', 'HBD', 'HBonds', 'Rotatable_Bonds', 'MW']
with bz2.open(output, 'wt') as bfh:
    bfh.write('\t'.join(headers)+'\n')
    with gzip.open(input) as gzh:
        with Chem.ForwardSDMolSupplier(gzh) as sdfh:
            for mol in sdfh:
                HBA = AllChem.CalcNumHBA(mol)
                HBD = AllChem.CalcNumHBD(mol)
                info = {'SMILES': Chem.MolToSmiles(mol),
                        'Identifier': mol.GetProp('mcule ID'),
                        'HBA': HBA,
                        'HBD': HBD,
                        'HBonds': HBA+HBA,
                        'Rotatable_Bonds': AllChem.CalcNumRotatableBonds(mol),
                        'MW': AllChem.CalcExactMolWt(mol),
                       }
                bfh.write('\t'.join(map(str, info.values()))+'\n')
print('Done')