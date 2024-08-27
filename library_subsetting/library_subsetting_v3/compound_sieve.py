__all__ = ['InchiType', 'BadCompound', 'SieveMode', 'CompoundSieve']

import io
import json
import enum
import itertools
from pathlib import Path
from typing import List, Dict, Any, Optional, NewType, Union

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, AllChem, rdDeprotect
from rdkit.Chem.rdfiltercatalog import FilterCatalogParams, FilterCatalog
try:
    import torch
    from .USRCAT_sociability import calc_summed_scores
except ImportError:
    torch = None
    calc_summed_scores = None
from .restrictive_decomposition import RestrictiveDecomposer

InchiType = NewType('InchiType', str)

# pains
_params = FilterCatalogParams()
_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)

class BadCompound(Exception):
    pass

class SieveMode(enum.Enum):
    basic = 0   # based on row information, no RDKit
    substructure = 1  # based on RDKit
    synthon = 2  # advanced

class CompoundSieve:
    """
    This class is intended to classify compounds on whether to keep them.
    Initialisation starts the classifier, calling the instance on a row will return a verdict.
    The row is a CXSMILES block row read by the ``read_cxsmiles_block`` static method.
    which will give the columns:

    * 'SMILES'
    * 'HBonds' (sum of HBA and HBD)
    * 'Rotatable_Bonds'
    * 'MW'

    Various filters are used (see ``cutoffs``) and unwanted groups are checked for.
    The method ``assess`` looks at the keys of in ``cutoffs``,
    which are in the format max_ or min_ followed by the key in the ``verdict``.

    NB. If there is no key in verdict that relative to a cutoff, nothing happens.
    This is because the assessment is continuous.
    A weird side effect is having to enforce 'max_N_protection_groups' = 0,
    _i.e._ no protection groups.

    Here are few examples of cutoffs:

    * `min_hbonds` - minimum number of HBonds
    * `min_synthon_sociability` - see below
    * `min_weighted_robogroups` - minumum number of wanted reaction product moieties (amide, sulfonamide, biaryl etc.)
    * `max_rota_per_da` - stop overly long rotatable bonds
    * `max_N_methylene` - like above but specific for too many CH2
    * `max_N_protection_groups` - default = zero protection groups
    * `max_largest_ring_size=8

    The classification stops if a violation is found (run `enable_analysis_mode` to disable cutoffs).
    The error BadCompound is raised if a violation is found, but is caught.
    """

    dps = rdDeprotect.GetDeprotections()
    unwanted = {'carbamate': Chem.MolFromSmiles('NC(=O)O'),
                'exocyclic ester': Chem.MolFromSmarts('[C!R](=O)[OH0!R]'),
                'exocyclic imine': Chem.MolFromSmarts('[C!R]=[N!R]'),
                'alkane': Chem.MolFromSmarts('[CH2!R]-[CH2!R]-[CH2!R]-[CH2!R]'),
                'hydrazine': Chem.MolFromSmarts('[N,n]-[N!R]'),
                }

    # this is a partial repetition of the rxns in RoboDecomposer!
    wanted = {'amide': Chem.MolFromSmarts('[N,n]-[C!R](=O)'),  # lactam is not okay, but on aza-arene is
              'sulfonamide': Chem.MolFromSmarts('[N,n]-[S!R](=O)(=O)'),
              'biaryl': Chem.MolFromSmarts('a-a'),  # suzuki...
              'secondary amine': Chem.MolFromSmarts('[c,C]-[N!R]-[c,C]'),  # Borch & Buchwald-hartwig?
              'substituted aza': Chem.MolFromSmarts('[NR,n]-[C!R]'),  # Buchwald-hartwig, Chan-Lam etc.
              # can the robot do Williamson esterification?
              # ...
              }

    wanted_weights = {'amide': 1,
                      'sulfonamide': 5,  # boost uncommon
                      'biaryl': 5,  # boost uncommon
                      'secondary amine': 0.3,  # not sure we can do thiss
                      'substituted aza': 0.3,  # ditto
                      }
    cutoffs = dict(
                   # these are medchem pickiness
                   min_N_rings=1,
                   max_N_methylene=6,
                   max_N_protection_groups=0,
                   max_largest_ring_size=8,
                   # these remove the worst quartiles
                   min_hbonds_per_HAC=1 / 5,
                   max_rota_per_HAC=1 / 5,
                   min_synthon_sociability_per_HAC=0.354839,
                   min_weighted_robogroups_per_HAC=0.0838,
                   max_boringness=0.1,
                   )

    # PAINS
    pains_catalog = FilterCatalog(_params)

    def __init__(self,
                 mode: SieveMode = SieveMode.synthon,
                 common_synthons_tally: Optional[Dict[InchiType, int]]=None,
                 common_synthons_usrcats: Optional[Dict[InchiType, list]]=None):
        self.mode = mode
        if self.mode == SieveMode.synthon:
            assert common_synthons_tally is not None, 'common_synthons_tally must be provided'
            assert common_synthons_usrcats is not None, 'common_synthons_usrcats must be provided'
            self.common_synthons_tally = torch.tensor(common_synthons_tally, device='cuda')
            self.common_synthons_usrcats = torch.tensor(common_synthons_usrcats, device='cuda')
            self.dejavu_synthons: Dict[InchiType, int] = {}
            self.nuveau_dejavu_synthons: Dict[InchiType, int] = {}
            self.robodecomposer = RestrictiveDecomposer()

    def enable_analysis_mode(self):
        """
        The cutoffs are disabled, so the values are all run...
        """
        self.cutoffs = {k: {'min': 0, 'max': float('inf')}[k[:3]] for k, v in self.cutoffs.items()}

    def __call__(self, row: pd.Series):
        verdict = {'acceptable': False, 'issue': ''}
        try:
            # ## Basic row info based
            self.calc_row_info(row, verdict)
            self.assess(verdict)
            if self.mode == SieveMode.basic:
                verdict['acceptable'] = True
                return verdict
            # ## Mol based
            mol = Chem.MolFromSmiles(row.SMILES)
            self.calc_mol_info(mol, verdict)
            self.assess(verdict)
            self.calc_boringness(mol, verdict)
            self.assess(verdict)
            self.assess_mol_patterns(mol, verdict)
            if self.mode == SieveMode.substructure:
                verdict['acceptable'] = True
                return verdict
            # ## Synthon based
            self.calc_robogroups(mol, verdict)
            self.assess(verdict)
            self.calc_synthon_info(mol, verdict)
            self.assess(verdict)
        except BadCompound as e:
            verdict['issue'] = str(e)
            return verdict
        except Exception as e:
            verdict['issue'] = f'Uncaught {e.__class__.__name__} exception: {e}'
            return verdict
        else:
            verdict['acceptable'] = True
            return verdict

    def calc_row_info(self, row: pd.Series, verdict: dict):
        verdict['hbonds'] = row.HBonds
        verdict['hbonds_per_HAC'] = row.HBonds / row.HAC
        verdict['rota_per_da'] = row.Rotatable_Bonds / row.MW
        verdict['rota_per_HAC'] = row.Rotatable_Bonds / row.HAC

    def assess(self, verdict: dict):
        for key in self.cutoffs:
            if key[4:] not in verdict:
                continue
            elif key[:3] == 'min' and verdict[key[4:]] < self.cutoffs[key]:
                raise BadCompound(f'{key[4:]} too low')
            elif key[:3] == 'max' and verdict[key[4:]] > self.cutoffs[key]:
                raise BadCompound(f'{key[4:]} too high')

    def calc_mol_info(self, mol: Chem.Mol, verdict: dict):
        # ## Mol based
        verdict['N_rings'] = rdMolDescriptors.CalcNumRings(mol)
        verdict['N_methylene'] = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[CH2X4!R]')))
        verdict['N_ring_atoms'] = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[R]')))
        verdict['largest_ring_size'] = max(map(len, mol.GetRingInfo().AtomRings()))
        verdict['N_protection_groups'] = rdDeprotect.Deprotect(mol, deprotections=self.dps) \
            .GetIntProp('DEPROTECTION_COUNT')

    def assess_mol_patterns(self, mol: Chem.Mol, verdict: dict):
        # ## Matching based
        for name, pattern in self.unwanted.items():
            if Chem.Mol.HasSubstructMatch(mol, pattern):
                raise BadCompound(f'Contains {name}')
        if len(self.pains_catalog.GetMatches(mol)):
            raise BadCompound('PAINS')
    def calc_sociability(self, synthon: Chem.Mol) -> float:
        synthon_inchi = Chem.MolToInchi(synthon)
        if synthon_inchi in self.dejavu_synthons:
            return self.dejavu_synthons[synthon_inchi]
        if synthon_inchi in self.nuveau_dejavu_synthons:
            return self.nuveau_dejavu_synthons[synthon_inchi]
        if synthon is None:
            return -1
        AllChem.EmbedMolecule(synthon)
        if Chem.Mol.GetNumHeavyAtoms(synthon) < 3 or Chem.Mol.GetNumConformers(synthon) == 0:
            return -1
        synthon_usrcat = torch.tensor(rdMolDescriptors.GetUSRCAT(synthon), device='cuda')
        sociability = calc_summed_scores(synthon_usrcat, self.common_synthons_usrcats,
                                         self.common_synthons_tally).tolist()
        self.nuveau_dejavu_synthons[synthon_inchi] = sociability
        return sociability
    def calc_synthon_info(self, mol, verdict):
        synthons: List[Chem.Mol] = self.robodecomposer.decompose(mol)
        verdict['N_synthons'] = len(synthons)
        verdict['synthon_sociability'] = sum(
            [self.calc_sociability(synthon) for synthon in synthons])
        verdict['synthon_sociability_per_HAC'] = verdict['synthon_sociability'] / mol.GetNumHeavyAtoms()

    def calc_robogroups(self, mol: Chem.Mol, verdict: dict):
        # ## Scoring wanted groups
        verdict[f'weighted_robogroups'] = 0
        for name, pattern in self.wanted.items():
            verdict[f'N_{name}'] = len(Chem.Mol.GetSubstructMatches(mol, pattern))
            verdict[f'weighted_robogroups'] += verdict[f'N_{name}'] * self.wanted_weights[name]
        verdict[f'weighted_robogroups_per_HAC'] = verdict[f'weighted_robogroups'] / mol.GetNumHeavyAtoms()

    def calc_n_fused_rings(self, mol):
        ars = mol.GetRingInfo().AtomRings()
        return sum([len(set(fore).intersection(aft)) > 1 for fore, aft in list(itertools.combinations(ars, 2))])

    def calc_boringness(self, mol: Chem.Mol, verdict: dict):
        """
        A big problem is that the top sociable compounds are boring compounds
        Namely, phenyls galore.
        """
        verdict['N_spiro'] = rdMolDescriptors.CalcNumSpiroAtoms(mol)
        verdict['N_bridgehead'] = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
        # an `AliphaticRings` includes heterocycles.
        verdict['N_alicyclics'] = rdMolDescriptors.CalcNumAliphaticRings(mol)
        verdict['N_fused_rings'] = self.calc_n_fused_rings(mol)
        verdict['N_heterocyclics'] = rdMolDescriptors.CalcNumHeterocycles(mol)
        verdict['N_aromatic_carbocylics'] = rdMolDescriptors.CalcNumAromaticCarbocycles(mol)
        # previously calculated: # not a methylene radical but a -CH2- group
        verdict['N_methylene'] = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[CH2X4!R]')))
        # make an arbitrary score of coolness
        cool_keys = ['N_spiro', 'N_bridgehead', 'N_alicyclics', 'N_fused_rings']
        # halfcool_keys = ['N_heterocyclics']
        boring_keys = ['N_aromatic_carbocylics']
        # boringish_keys = ['N_methylene']
        verdict['boringness'] = sum(map(verdict.get, boring_keys)) + \
                                verdict['N_methylene'] / 4 - \
                                sum(map(verdict.get, cool_keys)) - \
                                verdict['N_heterocyclics'] / 2

    def classify_df(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Runs the classification on a DataFrame and returns a DataFrame with the verdicts.

        :param df:
        :return:
        """
        _verdicts: pd.Series = df.apply(self, axis=1)
        verdicts: pd.DataFrame = pd.DataFrame(_verdicts.tolist(), index=_verdicts.index)
        print(f'{round(verdicts.acceptable.value_counts().to_dict().get(True, 0) / len(verdicts) * 100)}% accepted')
        return verdicts
