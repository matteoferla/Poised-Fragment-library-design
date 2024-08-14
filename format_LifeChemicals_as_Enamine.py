"""
LifeChemicals is in sdf format. For consistency & size it is converted to a similar format to Enamine REAL (cxsmiles)
with quick info used by library classification present.
Wasteful, but overcomes a lot of coding and LifeChemicals is much smaller.
"""

import bz2, gzip
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.SaltRemover import SaltRemover
from typing import Sequence
import gzip, contextlib
from pathlib import Path
from collections import defaultdict


saltremover = SaltRemover()

def add_nitrogen_charges(mol):
    #RDLogger.DisableLog('rdApp.*')
    mol.UpdatePropertyCache(strict=False)
    ps = Chem.DetectChemistryProblems(mol)
    if not ps:
        Chem.SanitizeMol(mol)
        return mol
    for p in ps:
        if p.GetType()=='AtomValenceException':
            at = mol.GetAtomWithIdx(p.GetAtomIdx())
            if at.GetAtomicNum()==7 and at.GetFormalCharge()==0 and at.GetExplicitValence()==4:
                at.SetFormalCharge(1)
    Chem.SanitizeMol(mol)
    #RDLogger.EnableLog('rdApp.*')
    return mol

def get_others(mol: Chem.Mol, in_idxs: Sequence[int]):
    return set([n.GetIdx() for a in map(mol.GetAtomWithIdx, in_idxs[1:]) for n in a.GetNeighbors()]).difference(in_idxs)

# ==================================

path = Path('LC_Stock_HTS_Compounds.sdf')
assert path.exists()
headers = ['SMILES', 'Identifier', 'MW', 'HAC', 'LogP', 'HBA', 'HBD', 'Rotatable_Bonds', 'FSP3', 'TPSA']
out_filename = f'LC_Stock_HTS.cxsmiles.bz2'
print(out_filename, 'started')
bfh = bz2.open(out_filename, 'wt')
bfh.write('\t'.join(headers)+'\n')

with Chem.ForwardSDMolSupplier(path.as_posix(), sanitize=False) as sdfh:
    for mol in sdfh:
        if mol is None:
            continue
        mol = saltremover.StripMol(add_nitrogen_charges(mol))
        if len(AllChem.GetMolFrags(mol)) > 1:
            continue
        Chem.SanitizeMol(mol)
        info = {'SMILES': Chem.MolToSmiles(mol),
                'Identifier': mol.GetProp('IDNUMBER'),
                'MW': AllChem.CalcExactMolWt(mol),
                'HAC': mol.GetIntProp('HAC'),
                'LogP': mol.GetDoubleProp('clogP'),
                'HBA': mol.GetIntProp('Acceptor'),
                'HBD': mol.GetIntProp('Donor'),
                'Rotatable_Bonds': mol.GetIntProp('RotBonds'),
                'FSP3': mol.GetDoubleProp('FSP3'),
                'TPSA': mol.GetDoubleProp('TPSA'),
               }
        bfh.write('\t'.join(map(str, info.values()))+'\n')
bfh.close()