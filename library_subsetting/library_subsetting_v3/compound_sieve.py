__all__ = ['InchiType', 'BadCompound', 'SieveMode', 'CompoundSieve']

import io
import json
import enum
import itertools, functools
from pathlib import Path
from typing import List, Dict, Any, Optional, NewType, Union, Callable, Tuple

import numpy as np
import pandas as pd
from rdkit import Chem, rdBase
from rdkit.Chem import rdMolDescriptors, AllChem, rdDeprotect
from rdkit.Chem import FilterCatalog, rdfiltercatalog

from .pipiteur import Pipiteur, PIPType
from . import data
from .util import ultranormalize, autopass_fun
try:
    import torch
    from .USRCAT_sociability import calc_summed_scores
except ImportError:
    torch = None
    calc_summed_scores = None
from .restrictive_decomposition import RestrictiveDecomposer

InchiType = NewType('InchiType', str)

# pains
pains_catalogue_params = rdfiltercatalog.FilterCatalogParams()
pains_catalogue_params.AddCatalog(rdfiltercatalog.FilterCatalogParams.FilterCatalogs.PAINS)

class BadCompound(Exception):
    pass

class SieveMode(enum.Enum):
    basic = 0   # based on row information, no RDKit
    substructure = 1  # based on RDKit
    synthon_v2 = 2  # advanced, old v2.
    synthon_v3 = 3  # advanced, v3

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
    * `min_synthon_score` - minumum number of wanted reaction product moieties (amide, sulfonamide, biaryl etc.)
    * `max_rota_per_da` - stop overly long rotatable bonds
    * `max_N_methylene` - like above but specific for too many CH2
    * `max_N_protection_groups` - default = zero protection groups
    * `max_largest_ring_size=8

    The classification stops if a violation is found (run `enable_analysis_mode` to disable cutoffs).
    The error BadCompound is raised if a violation is found, but is caught.

    For diagnostics...

    .. code-block:: python

        sieve = CompoundSieve()
        sieve.exception_to_catch = ()
        sieve.cutoffs = {}
        df = CompoundSieve.prep_df(df,smiles_col = 'SMILES', mol_col='mol')
        sieve(df.iloc[0])

    Due to the fact that each 'calc' method adds to the verdict,
    to call one of these, one needs to pass the verdict dictionary.

    .. code-block:: python

        verdict = {}
        sieve.calc_pip(mol, verdict)
    """

    dps = rdDeprotect.GetDeprotections()
    # this could be written as a FilterCatalog
    unwanted = {'exocyclic carbamate': Chem.MolFromSmarts('[N!R]-C(=O)-O'),
                'exocyclic ester': Chem.MolFromSmarts('[C!R](=O)-[OH0!R]'),
                'exocyclic imine': Chem.MolFromSmarts('[C!R]=[N!R]'),
                'exocyclic alkane': Chem.MolFromSmarts('[CH2!R]-[CH2!R]-[CH2!R]-[CH2!R]'),
                'exocyclic hydrazine': Chem.MolFromSmarts('[N,n]-[N!R]'),
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
                   min_synthon_score_per_HAC=0.138470, # v2 is 0.0838
                   min_weighted_robogroups_per_HAC=0.0838,  # quartile
                   max_boringness=0.,
                   min_combined_Zscore=0. # above the arithmetic mean
                   )

    exception_to_catch = (Exception,    )

    # PAINS
    pains_catalog = rdfiltercatalog.FilterCatalog(pains_catalogue_params)

    def __init__(self,
                 mode: SieveMode = SieveMode.synthon_v3,
                 common_synthons_tally: Optional[Dict[InchiType, int]]=None,
                 common_synthons_usrcats: Optional[Dict[InchiType, list]]=None,
                 screening_filename: Optional[str]=None):
        self.mode = mode
        if self.mode == SieveMode.synthon_v2:
            assert common_synthons_tally is not None, 'common_synthons_tally must be provided'
            assert common_synthons_usrcats is not None, 'common_synthons_usrcats must be provided'
            self.common_synthons_tally = torch.tensor(common_synthons_tally, device='cuda')
            self.common_synthons_usrcats = torch.tensor(common_synthons_usrcats, device='cuda')
            self.dejavu_synthons: Dict[InchiType, int] = {}
            self.nuveau_dejavu_synthons: Dict[InchiType, int] = {}
            self.robodecomposer = RestrictiveDecomposer()
        elif self.mode == SieveMode.synthon_v3:
            self.screening_catalog = self.get_screening_library_catalog(screening_filename)
            self.robodecomposer = RestrictiveDecomposer()
            Pipiteur.fdef = data.read_MolChemicalFeatureFactory('Steph_features.fdef')
            Pipiteur.wanted = sorted(['Donor',
                                      'Acceptor',
                                      'NegIonizable',
                                      'PosIonizable',
                                      'Aromatic',
                                      'Aliphatic', ])
            self.pipiteur = Pipiteur(order=3,
                                min_d=2, max_d=8,
                                resolution=0.5,
                                )
            # PIP freqs
            # todo fix:
            # in cumulative_pip dataset its a tuple of each type
            # in unskew params its colon separated...
            self.pip_freqs: Dict[Tuple[str, str, str], np.array] = data.read_pickle('cumulative_pip_smooth_log.pkl.gz')
            self.likelihood_unskew_funs: Dict[str, Callable] = data.parse_unskew_funs('likelihood_skew_params.json')

    def enable_analysis_mode(self):
        """
        The cutoffs are disabled, so the values are all run...
        """
        self.cutoffs = {k: {'min': 0, 'max': float('inf')}[k[:3]] for k, v in self.cutoffs.items()}

    def classify_df(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Runs the classification on a DataFrame and returns a DataFrame with the verdicts.

        :param df:
        :return:
        """
        with rdBase.BlockLogs():
            _verdicts: pd.Series = df.apply(self, axis=1)
        verdicts: pd.DataFrame = pd.DataFrame(_verdicts.tolist(), index=_verdicts.index)
        print(f'{round(verdicts.acceptable.value_counts().to_dict().get(True, 0) / len(verdicts) * 100)}% accepted')
        return verdicts

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
            if 'mol' not in row.index:
                mol = Chem.MolFromSmiles(row.SMILES)
            else:
                mol = row.mol
            self.calc_mol_info(mol, verdict)
            self.assess(verdict)
            self.calc_boringness(mol, verdict)
            self.assess(verdict)
            self.assess_mol_patterns(mol, verdict)
            if self.mode == SieveMode.substructure:
                verdict['acceptable'] = True
                return verdict
            # ## Synthon based
            if self.mode == SieveMode.synthon_v2:
                self.calc_robogroups(mol, verdict)
                self.assess(verdict)
                self.calc_synthon_info_old(mol, verdict)
                self.assess(verdict)
            elif self.mode == SieveMode.synthon_v3:
                self.calc_outtajail_score(mol, verdict)  # boost for matches to XChem screening library
                self.calc_synthon_info(mol, verdict)
                self.assess(verdict)
                if Chem.Mol.GetNumConformers(mol) == 0:
                    mol = AllChem.AddHs(mol)
                    AllChem.EmbedMolecule(mol)
                self.calc_pip(mol, verdict)
                self.calc_score(mol, verdict)
        except BadCompound as e:
            verdict['issue'] = str(e)
            return verdict
        except self.exception_to_catch as e:
            verdict['issue'] = f'Uncaught {e.__class__.__name__} exception: {e}'
            return verdict
        else:
            verdict['acceptable'] = True
            return verdict

    def calc_row_info(self, row: pd.Series, verdict: dict):
        verdict['hbonds'] = row.HBonds
        verdict['HAC'] = row.HAC
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
        verdict['largest_ring_size'] = max([0, *map(len, mol.GetRingInfo().AtomRings())])
        verdict['N_protection_groups'] = rdDeprotect.Deprotect(mol, deprotections=self.dps) \
            .GetIntProp('DEPROTECTION_COUNT')

    def assess_mol_patterns(self, mol: Chem.Mol, verdict: dict):
        # ## Matching based
        for name, pattern in self.unwanted.items():
            if Chem.Mol.HasSubstructMatch(mol, pattern):
                raise BadCompound(f'Contains {name}')
        if len(self.pains_catalog.GetMatches(mol)):
            raise BadCompound('PAINS')

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
        verdict['N_aromatic_carbocycles'] = rdMolDescriptors.CalcNumAromaticCarbocycles(mol)
        # previously calculated: # not a methylene radical but a -CH2- group
        verdict['N_methylene'] = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[CH2X4!R]')))
        # make an arbitrary score of coolness
        cool_keys = ['N_spiro', 'N_bridgehead', 'N_alicyclics', 'N_fused_rings']
        # halfcool_keys = ['N_heterocyclics']
        boring_keys = ['N_aromatic_carbocycles']
        # boringish_keys = ['N_methylene']
        verdict['boringness'] = sum(map(verdict.get, boring_keys)) + \
                                verdict['N_methylene'] / 4 - \
                                sum(map(verdict.get, cool_keys)) - \
                                verdict['N_heterocyclics'] / 2
        verdict['boringness_per_HAC'] = verdict['boringness'] / verdict['HAC']

    def calc_synthon_info(self, mol: Chem.Mol, verdict: dict):
        # version 3
        synthons: List[Chem.Mol] = [s for s in self.robodecomposer.decompose(mol) if s.GetNumHeavyAtoms() > 2]
        verdict['N_synthons'] = len(synthons)
        verdict['synthon_score'] = self.robodecomposer.synthon_score(mol)
        verdict['synthon_score_per_HAC'] = verdict['synthon_score'] / verdict['HAC']

    # this is ad hoc
    score_weights = {'synthon_score_per_HAC': 1,
                     'hbonds_per_HAC': 0.5,
                     'rota_per_HAC': -0.5,
                     #'N_synthons_per_HAC': 1,
                     'N_spiro_per_HAC': 0.1,
                     'N_bridgehead_per_HAC': 0.1,
                     'N_alicyclics_per_HAC': 0.1,
                     'N_fused_rings_per_HAC': 0.1,
                     'N_aromatic_carbocycles_per_HAC': -0.1,
                     'N_heterocyclics_per_HAC': 0.05,
                     'N_methylene_per_HAC': -0.02,
                     'outtajail_score': +2.0,
                     'pip_common_rms': +1.0,
                     'pip_common_rms': +0.5,
                     }
    # these are from Enamine 1M random sample w/o removals
    ref_means = {'synthon_score_per_HAC': 0.21526508936919203,
                 'hbonds_per_HAC': 0.24447480230871893,
                 'rota_per_HAC': 0.2317342518844832,
                 'N_synthons_per_HAC': 0.12582623253284875,
                 'N_spiro_per_HAC': 0.0032131689670107065,
                 'N_bridgehead_per_HAC': 0.005046401318307808,
                 'N_alicyclics_per_HAC': 0.05905642932651799,
                 'N_fused_rings_per_HAC': 0.015589662661026338,
                 'N_aromatic_carbocycles_per_HAC': 0.01918927338610745,
                 'N_heterocyclics_per_HAC': 0.06979145110309398,
                 'N_methylene_per_HAC': 0.08125398462902535,
                 'pip_common_mean': -0.6916386218771066,
                 'pip_uncommon_mean': -1.6803720478038888,
                 'pip_common_rms': 3.0394492640351154,
                 'pip_uncommon_rms': 0.733971015961427,
                 }
    ref_stds = {'synthon_score_per_HAC': 0.1137501096872908,
                'hbonds_per_HAC': 0.06981618332292346,
                'rota_per_HAC': 0.07809299292460986,
                'N_synthons_per_HAC': 0.03198042716946067,
                'N_spiro_per_HAC': 0.010936756469896591,
                'N_bridgehead_per_HAC': 0.020431219793333164,
                'N_alicyclics_per_HAC': 0.0416689554316131,
                'N_fused_rings_per_HAC': 0.028725197477523886,
                'N_aromatic_carbocycles_per_HAC': 0.02464447282361974,
                'N_heterocyclics_per_HAC': 0.03760917968539562,
                'N_methylene_per_HAC': 0.061085330799282266,
                'pip_common_mean': 0.6568192777260955,
                'pip_uncommon_mean': 0.26532635670645766,
                'pip_common_rms': 2.143211511620346,
                'pip_uncommon_rms': 0.7553837646990731,
                }

    def calc_score(self, mol: Chem.Mol, verdict: dict):
        """
        This is a weighted sum of Zscored normalised values.
        These are skewed, but that is not a terrible thing: is something has a really high value for one of the metrics
        then it's actually gook it is mad high!
        """
        for key in self.score_weights:
            if key not in verdict:
                verdict[key] = verdict[key.replace('_per_HAC', '')] / verdict['HAC']
        # only calc;ed for keys in score_weights
        wzscorify = lambda k: self.score_weights[k] * (verdict[k] - self.ref_means.get(k, 0)) / self.ref_stds.get(k, 1)
        verdict['combined_Zscore'] = sum([wzscorify(k) for k in self.score_weights]) / (len(self.score_weights) ** 0.5)

    def calc_outtajail_score(self, mol: Chem.Mol, verdict: dict):
        """
        The out-of-jail-card score is a shift of the combined Zscore to boost substructures of XChem screening library
        It is the shifted number of atoms of the matched substructure:

        * 0–4 HAC gets +0
        * 5 HAC get +0.1
        * 10 HAC gets +0.6
        * 15 HAC gets +1.1
        * 20 HAC gets +1.6

        The distribution of the max HAC of the matches is dominated by zeros,
        while the non-zero values start at 5 due to the ring requirement for a match.
        """
        med = 10.0  # HAC = denominator for score to the power of k
        lower_HAC = 4.  # 5 is min
        k = 1.  # exponent disabled
        # this is the largest HAC match
        outtajail_value = max([0]+[int(match.GetDescription().split(':')[1]) for match in self.screening_catalog.GetMatches(mol)])
        verdict['outtajail_value'] = outtajail_value
        verdict['outtajail_score'] = (max(lower_HAC, outtajail_value)**k - lower_HAC**k) / med**k

    common_pip_trios = ['Acceptor:Acceptor:Acceptor',
                          'Acceptor:Acceptor:Aromatic',
                          'Acceptor:Acceptor:Donor',
                          'Acceptor:Aromatic:Aromatic',
                          'Acceptor:Aromatic:Donor',
                          'Acceptor:Donor:Donor',
                          'Aromatic:Aromatic:Aromatic',
                          'Aromatic:Aromatic:Donor',
                          'Aromatic:Donor:Donor',
                          'Donor:Donor:Donor']
    uncommon_pip_trios =  ['Acceptor:Acceptor:Aliphatic',
                          'Acceptor:Acceptor:NegIonizable',
                          'Acceptor:Acceptor:PosIonizable',
                          'Acceptor:Aliphatic:Aliphatic',
                          'Acceptor:Aliphatic:Aromatic',
                          'Acceptor:Aliphatic:Donor',
                          'Acceptor:Aliphatic:NegIonizable',
                          'Acceptor:Aliphatic:PosIonizable',
                          'Acceptor:Aromatic:NegIonizable',
                          'Acceptor:Aromatic:PosIonizable',
                          'Acceptor:Donor:NegIonizable',
                          'Acceptor:Donor:PosIonizable',
                          'Acceptor:NegIonizable:NegIonizable',
                          'Acceptor:NegIonizable:PosIonizable',
                          'Acceptor:PosIonizable:PosIonizable',
                          'Aliphatic:Aliphatic:Aliphatic',
                          'Aliphatic:Aliphatic:Aromatic',
                          'Aliphatic:Aliphatic:Donor',
                          'Aliphatic:Aliphatic:NegIonizable',
                          'Aliphatic:Aliphatic:PosIonizable',
                          'Aliphatic:Aromatic:Aromatic',
                          'Aliphatic:Aromatic:Donor',
                          'Aliphatic:Aromatic:NegIonizable',
                          'Aliphatic:Aromatic:PosIonizable',
                          'Aliphatic:Donor:Donor',
                          'Aliphatic:Donor:NegIonizable',
                          'Aliphatic:Donor:PosIonizable',
                          'Aliphatic:NegIonizable:NegIonizable',
                          'Aliphatic:NegIonizable:PosIonizable',
                          'Aliphatic:PosIonizable:PosIonizable',
                          'Aromatic:Aromatic:NegIonizable',
                          'Aromatic:Aromatic:PosIonizable',
                          'Aromatic:Donor:NegIonizable',
                          'Aromatic:Donor:PosIonizable',
                          'Aromatic:NegIonizable:NegIonizable',
                          'Aromatic:NegIonizable:PosIonizable',
                          'Aromatic:PosIonizable:PosIonizable',
                          'Donor:Donor:NegIonizable',
                          'Donor:Donor:PosIonizable',
                          'Donor:NegIonizable:NegIonizable',
                          'Donor:NegIonizable:PosIonizable',
                          'Donor:PosIonizable:PosIonizable',
                          'NegIonizable:NegIonizable:NegIonizable',
                          'NegIonizable:NegIonizable:PosIonizable',
                          'NegIonizable:PosIonizable:PosIonizable',
                          'PosIonizable:PosIonizable:PosIonizable']
    def calc_pip(self, mol: Chem.Mol, verdict: dict):
        pipi: Dict[Tuple[str, str, str], np.array] = self.pipiteur(mol)
        for k, m in pipi.items():
            v = np.min(self.pip_freqs[k][np.where(m != 0)], initial=0)
            joint_key = ':'.join(k)
            verdict[f'pip_{joint_key}_rarest'] = v
            verdict[f'pip_{joint_key}_rarest_normalised'] = self.likelihood_unskew_funs[joint_key](v)
        pip_commons = np.array([verdict[f'pip_{k}_rarest_normalised'] for k in self.common_pip_trios])
        pip_uncommons = np.array([verdict[f'pip_{k}_rarest_normalised'] for k in self.uncommon_pip_trios])
        pips = np.concatenate([pip_commons, pip_uncommons])
        verdict[f'pip_mean'] = np.mean(pips)
        verdict[f'pip_common_mean'] = np.mean(pip_commons)
        verdict[f'pip_uncommon_mean'] = np.mean(pip_uncommons)
        # the values are zero centred. I want them min centred.
        lower_bound = pips.min()
        verdict[f'pip_rms'] = np.mean(np.power(pips - lower_bound, 2))
        verdict[f'pip_common_rms'] = np.mean(np.power(pip_commons - lower_bound, 2))
        verdict[f'pip_uncommon_rms'] = np.mean(np.power(pip_uncommons - lower_bound, 2))

    @staticmethod
    def prep_df(df, smiles_col: str = 'SMILES', mol_col=None):
        """
        Fixes in place a dataframe to make it compatible with ``classify_df``

        :param df:
        :param smiles_col:
        :param mol_col:
        :return:
        """
        df = df.copy()
        if smiles_col != 'SMILES':
            df = df.rename(column={smiles_col: 'SMILES'})
        if mol_col is None:
            df['mol'] = df.SMILES.apply(Chem.MolFromSmiles)
        elif mol_col != 'mol':
            df = df.rename(columns={mol_col: 'mol'}).copy()
        else:
            pass  # all good
        df['HAC'] = df.mol.apply(Chem.Mol.GetNumHeavyAtoms)
        df['HBonds'] = df.mol.apply(rdMolDescriptors.CalcNumHBD) + df.mol.apply(rdMolDescriptors.CalcNumHBD)
        df['Rotatable_Bonds'] = df.mol.apply(rdMolDescriptors.CalcNumRotatableBonds)
        df['MW'] = df.mol.apply(rdMolDescriptors.CalcExactMolWt)
        return df.copy()

    @staticmethod
    def get_screening_library_catalog(filename=None):
        filename = filename or Path(__file__).parent / 'data' / 'screening_mols.smi'
        with open(filename) as fh:
            libsynthons = [Chem.MolFromSmiles(line.strip().split('\t')[0]) for line in fh]
        catalog = FilterCatalog.FilterCatalog()
        for i, refmol in enumerate(libsynthons):
            sm = FilterCatalog.SmartsMatcher(refmol)
            entry = FilterCatalog.FilterCatalogEntry(f'#{i}:{Chem.Mol.GetNumHeavyAtoms(refmol)}', sm)
            catalog.AddEntry(entry)
        return catalog

    # ------------------------ DEPRECATED ------------------------

    def calc_sociability(self, synthon: Chem.Mol) -> float:
        "This is v2 code"
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

    def calc_synthon_info_old(self, mol, verdict):
        "This is v2 code"
        synthons: List[Chem.Mol] = self.robodecomposer.decompose(mol)
        verdict['N_synthons'] = len(synthons)
        verdict['synthon_sociability'] = sum(
            [self.calc_sociability(synthon) for synthon in synthons])
        verdict['synthon_sociability_per_HAC'] = verdict['synthon_sociability'] / verdict['HAC']

    def calc_robogroups(self, mol: Chem.Mol, verdict: dict):
        "This is v2 code"
        # ## Scoring wanted groups
        verdict[f'synthon_score'] = 0
        for name, pattern in self.wanted.items():
            verdict[f'N_{name}'] = len(Chem.Mol.GetSubstructMatches(mol, pattern))
            verdict[f'synthon_score'] += verdict[f'N_{name}'] * self.wanted_weights[name]
        verdict[f'synthon_score_per_HAC'] = verdict[f'synthon_score'] / verdict['HAC']

