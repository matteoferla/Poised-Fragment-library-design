"""
Mcule is in sdf format. For consistency & size it is converted to a similar format to Enamine REAL (cxsmiles)
with quick info used by library classification present.
Wasteful, but overcomes a lot of coding and Mcule is much smaller.
"""

import bz2, gzip
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Draw, rdMolDescriptors, Descriptors
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.SaltRemover import SaltRemover
from typing import Sequence
import gzip, contextlib
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



path = Path('/tmp') / 'mcule_purchasable_full_240610.sdf.gz'
assert path.exists()
headers = ['SMILES', 'Identifier', 'MW', 'HAC', 'LogP', 'HBA', 'HBD', 'Rotatable_Bonds', 'FSP3', 'TPSA']

i = 0
k = 1_000_000
bfh = None
with gzip.open(path.as_posix()) as gzh:
    with Chem.ForwardSDMolSupplier(gzh, sanitize=False) as sdfh:
        for mol in sdfh:
            if mol is None:
                continue
            if i % k == 0:
                if bfh is not None:
                    bfh.close()
                out_filename = f'Mcule_ultimate_1Mchunk{i // k}.cxsmiles.bz2'
                print(out_filename, 'started')
                bfh = bz2.open(out_filename, 'wt')
                bfh.write('\t'.join(headers)+'\n')
            i += 1
            try:
                mol = saltremover.StripMol(add_nitrogen_charges(mol))
                Chem.SanitizeMol(mol)
                HBA = AllChem.CalcNumHBA(mol)
                HBD = AllChem.CalcNumHBD(mol)
                info = {'SMILES': Chem.MolToSmiles(mol),
                        'Identifier': mol.GetProp('mcule ID'),
                        'MW': AllChem.CalcExactMolWt(mol),
                        'HAC': mol.GetNumHeavyAtoms(),
                        'LogP': Descriptors.MolLogP(mol),
                        'HBA': HBA,
                        'HBD': HBD,
                        'Rotatable_Bonds': AllChem.CalcNumRotatableBonds(mol),
                        'FSP3': rdMolDescriptors.CalcFractionCSP3(mol),
                        'TPSA': Descriptors.TPSA(mol),
                       }
                bfh.write('\t'.join(map(str, info.values()))+'\n')
            except Exception as error:
                print(error.__class__.__name__, error)
bfh.close()