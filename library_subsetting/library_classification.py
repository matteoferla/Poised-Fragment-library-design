import io
import json
import itertools
from pathlib import Path
from typing import List, Dict, Any, Optional, NewType, Union

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, AllChem, rdDeprotect
from rdkit.Chem.rdfiltercatalog import FilterCatalogParams, FilterCatalog

InchiType = NewType('InchiType', str)

# pains
_params = FilterCatalogParams()
_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)


def calc_n_fused_rings(mol):
    ars = mol.GetRingInfo().AtomRings()
    return sum([len(set(fore).intersection(aft)) > 1 for fore, aft in list(itertools.combinations(ars, 2))])


class RoboDecomposer:
    """
    Creates synthons for reactions the XChem robot can do.
    amidation, Schotten-Baumann, sulfo Schotten-Baumann, Suzuki, but not Chan-Lam, Williamson, Borch, Huisgen etc.
    It is lactam/sulfam safe, which is mainly why I am writing it backwards from scratch.
    ``synthons = RoboDecomposer.decompose(mol)``

    """

    def __init__(self,
                 simplify_halide=True,  # not a real reaction, but a simplification
                 amide=True,  #HATU + Schotten-Baumann
                 sulfonamide=True,  # sulfonamide Schotten–Baumann
                 biaryl=True,  # Suzuki
                 arylamine=False,  # Chan-Lam
                 ether=False,  # Williamson
                 alkyne=False,  # Sonogashira
                 amine=False,  # Borch
                 triazole=False,  # Huisgen
                 ureido=False,  # isocyanate
                 ):
        # order matter as they will be run in that way
        self.rxns: Dict[str, AllChem.ChemicalReaction] = {}
        # =========================================================
        # ## simplification
        if simplify_halide:
            # not a real reaction, but a simplification
            simplification_rxn = AllChem.ReactionFromSmarts('[Cl,Br,I:1]>>[Cl:1]')
            self.rxns['halo-simplification'] = simplification_rxn
        # =========================================================
        # ## HATU & Schotten-Baumann
        if amide:
            # secondary exocyclic amides only:
            # Does not break lactams
            # amide_hydrolysis_rxn = AllChem.ReactionFromSmarts('[C!R:1](=[O:2])-[NH:3]>>[C!R:1](=[O:2])[O].[N:3]')
            # secondary/tertiary exocyclic amides, but no ureido
            amide_hydrolysis_rxn = AllChem.ReactionFromSmarts(
                '[N:1]-[C!R:2](=[O:3])-[!n&!N:4]>>[N:1].[Cl][C:2](=[O:3])-[*:4]')
            self.rxns['alkyl-amide'] = amide_hydrolysis_rxn
            # this is Schotten-Baumann only
            aromatic_amide_hydrolysis_rxn = AllChem.ReactionFromSmarts(
                '[n:1]-[C!R:2](=[O:3])-[!N:4]>>[nH:1].[Cl][C:2](=[O:3])-[*:4]')
            self.rxns['aryl-amide'] = aromatic_amide_hydrolysis_rxn
        # =========================================================
        # ## Schotten-Baumann
        if sulfonamide:
            # sulfonamide
            sulfonamide_cleavage_rxn = AllChem.ReactionFromSmarts(
                '[N:1]-[S!R:2](=[O:3])(=[O:4])>>[N:1].[Cl]-[S:2](=[O:3])(=[O:4])')
            self.rxns['alkyl-sulfonamide'] = sulfonamide_cleavage_rxn
            aromatic_sulfonamide_cleavage_rxn = AllChem.ReactionFromSmarts(
                '[n:1]-[S!R:2](=[O:3])(=[O:4])>>[nH:1].[Cl]-[S:2](=[O:3])(=[O:4])')
            self.rxns['aryl-sulfonamide'] = aromatic_sulfonamide_cleavage_rxn
        # =========================================================
        # ## Suzuki
        if biaryl:
            biaryl_cleavage_rxn = AllChem.ReactionFromSmarts('[aR:1]-[aR:2]>>[aHR:1].[aHR:2]')
            self.rxns['biaryl'] = biaryl_cleavage_rxn
        # =========================================================
        # ## Chan-Lam
        if arylamine:
            arylamine_cleavage_rxn = AllChem.ReactionFromSmarts('[aR:1]-[N:2]>>[aR:1]-B(-[OH])(-[OH]).[NH:2]')
            self.rxns['tertiary-arylamine'] = biaryl_cleavage_rxn
        # =========================================================
        # ##  Williamson
        if ether:
            ether_cleavage_rxn = AllChem.ReactionFromSmarts('[CX4!R:1]-[O!R:2]>>[C:1].[OH:2]')
            self.rxns['ether'] = ether_cleavage_rxn
            thioether_cleavage_rxn = AllChem.ReactionFromSmarts('[CX4!R:1]-[S!R:2]>>[C:1].[SH:2]')
            self.rxns['thioether'] = thioether_cleavage_rxn
        # =========================================================
        # ##  Borch
        # TODO check if the robot do reductive amination (Borch)?
        if amine:
            amine_cleavage_rxn = AllChem.ReactionFromSmarts('[C!R:1]-[N!R:2]>>[C:1].[N:2]')
            self.rxns['amine'] = amine_cleavage_rxn
        # =========================================================
        # ##  Sonogashira
        if alkyne:
            alkyne_cleavage_rxn = AllChem.ReactionFromSmarts('[C!R:1]#[C:2]>>[C!R:1].[C:2]')
            self.rxns['alkyne'] = alkyne_cleavage_rxn
        # =========================================================
        # ##  Huisgen
        if triazole:
            triazole_cleavage_rxn = AllChem.ReactionFromSmarts('[c:1]1:[c:2]:[n:3]:[n:4]:[nX3:5](-[C:6]):1>>[C:1]#[C:2].[N:3]#[N+:4]-[N-:5]-[C:6]')
            self.rxns['triazole'] = triazole_cleavage_rxn
        # =========================================================
        # ##  Ureidation
        if ureido:
            # this expects a primary on one side...
            # as the order matters, it does 2ndary first
            ureido_cleavage_rxn = AllChem.ReactionFromSmarts('[NH0:1]-[C:2](=[O:3])-[N:4]>>[NH1:1].[C:2](=[O:3])-[N:4]')
            self.rxns['ureido_2ary'] = ureido_cleavage_rxn
            ureido_cleavage_rxn = AllChem.ReactionFromSmarts('[NH1!R:1]-[C:2](=[O:3])-[N:4]>>[NH2:1].[C:2](=[O:3])-[N:4]')
            self.rxns['ureido_1ary'] = ureido_cleavage_rxn
        # =========================================================
        # etc.
        # ...
        # =========================================================
        for rxn in self.rxns.values():
            rxn.Initialize()

    def _recursive_cleave(self, mol: Chem.Mol, rxn: AllChem.ChemicalReaction, rxn_name: str) -> List[Chem.Mol]:
        """
        Returns a list of Chem.Mol
        """
        try:
            AllChem.SanitizeMol(mol)  # , catchErrors=True
            if rxn.IsMoleculeReactant(mol):
                prods = rxn.RunReactant(mol, 0)
                if len(prods) > 1:
                    subprods = [subprod for prod in prods[0] for subprod in self._recursive_cleave(mol=prod,
                                                                                                   rxn=rxn,
                                                                                                   rxn_name=rxn_name)]
                else:
                    subprods = prods[0]
                for prod in subprods:
                    AllChem.SanitizeMol(prod)
                    done_rxn_names = prod.GetProp('Reaction').split() if prod.HasProp('Reaction') else []
                    done_rxn_names.append(rxn_name)
                    prod.SetProp('Reaction', ' '.join(done_rxn_names))
                return subprods
        except Exception as error:
            pass
        return [mol]

    def decompose(self, mol) -> List[Chem.Mol]:
        """
        Returns a list of Chem.Mol (synthons).
        See property ``Reactions`` for the cleaved moiety.
        """
        synthons = [mol]
        for rxn_name, rxn in self.rxns.items():
            synthon_groups = [self._recursive_cleave(mol=synthon, rxn=rxn, rxn_name=rxn_name) for synthon in synthons]
            synthons = [synthon for synthons in synthon_groups for synthon in synthons]
        return synthons

    @classmethod
    def test(cls):
        robodecomposer = cls()
        results = []
        for name, smi in dict(amide=('CCNC(=O)CC', 2),
                              sulfonamide=('CCNS(=O)(=O)CC', 2),
                              diglycine=('NCC(=O)NCC(=O)O', 2),
                              triglycine=('NCC(=O)NCC(=O)NCC(=O)O', 3),
                              lactam=('C1NC(=O)CC1', 1),
                              tertamide=('CCN(C)C(=O)CC', 2),
                              cycloalkyl_amide=('N1(CCCCC1)C(=O)CC', 2),
                              aryl_amide=('n1(cccc1)C(=O)CC', 2),
                              saccharin=('O=C2c1ccccc1S(=O)(=O)N2', 1),
                              urea=('CCNC(=O)NCC', 1),
                              biaryl=('n1(cccc1)-c1(ccccc1)', 2)
                              ).items():
            mol = Chem.MolFromSmiles(smi[0])
            synthons: List[Chem.Mol] = robodecomposer.decompose(mol)
            assert len(synthons) == smi[1], f'{name} failed'
            results.append(dict(name=name, input=mol, output=synthons, expected_products=smi[1]))
        return results


class BadCompound(Exception):
    pass


class Classifier:
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
    Namely:

    * `min_hbonds` - minimum number of HBonds
    * `min_synthon_sociability` - see below
    * `min_weighted_robogroups` - minumum number of wanted reaction product moieties (amide, sulfonamide, biaryl etc.)
    * `max_rota_per_da` - stop overly long rotatable bonds
    * `max_N_methylene` - like above but specific for too many CH2
    * `max_N_protection_groups` - default = zero protection groups
    * `max_largest_ring_size=8, )

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
    cutoffs = dict(min_N_rings=1,
                   max_N_methylene=6,
                   max_N_protection_groups=0,
                   max_largest_ring_size=8,
                   # these are optional effectively
                   min_hbonds_per_HAC=1 / 5,
                   max_rota_per_HAC=1 / 5,
                   min_synthon_sociability_per_HAC=0.354839,
                   min_weighted_robogroups_per_HAC=0.0838,
                   max_boringness=0.1,
                   )

    # PAINS
    pains_catalog = FilterCatalog(_params)

    def __init__(self, sociability: Optional[Union[Dict[InchiType, float], str]] = None):
        if isinstance(sociability, dict):
            self.sociability = sociability
        elif isinstance(sociability, str) and '.json' in sociability:
            self.sociability = json.loads(Path(sociability).read_text())
        else:
            self.sociability = {}
        self.robodecomposer = RoboDecomposer()

    def enable_analysis_mode(self):
        """
        The cutoffs are disabled, so the values are all run...
        """
        self.cutoffs = {k: {'min': 0, 'max': float('inf')}[k[:3]] for k, v in self.cutoffs.items()}

    def __call__(self, row: pd.Series):
        verdict = {'acceptable': False, 'issue': ''}
        try:
            self.calc_row_info(row, verdict)
            self.assess(verdict)
            mol = Chem.MolFromSmiles(row.SMILES)
            self.calc_mol_info(mol, verdict)
            self.assess(verdict)
            self.calc_boringness(mol, verdict)
            self.assess(verdict)
            self.assess_mol_patterns(mol, verdict)
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

    def calc_synthon_info(self, mol, verdict):
        synthons: List[Chem.Mol] = self.robodecomposer.decompose(mol)
        verdict['N_synthons'] = len(synthons)
        # sociability is Dict of inchi to value,
        # where value is N of VCs with synthon times N of synthons with USRCAT > 0.7
        # the sum (not rare synthon penalised by roots)
        verdict['synthon_sociability'] = sum(
            [self.sociability.get(Chem.MolToInchi(synthon), -0.) for synthon in synthons])
        verdict['synthon_sociability_per_HAC'] = verdict['synthon_sociability'] / mol.GetNumHeavyAtoms()

    def calc_robogroups(self, mol: Chem.Mol, verdict: dict):
        # ## Scoring wanted groups
        verdict[f'weighted_robogroups'] = 0
        for name, pattern in self.wanted.items():
            verdict[f'N_{name}'] = len(Chem.Mol.GetSubstructMatches(mol, pattern))
            verdict[f'weighted_robogroups'] += verdict[f'N_{name}'] * self.wanted_weights[name]
        verdict[f'weighted_robogroups_per_HAC'] = verdict[f'weighted_robogroups'] / mol.GetNumHeavyAtoms()


    def calc_boringness(self, mol: Chem.Mol, verdict: dict):
        """
        A big problem is that the top sociable compounds are boring compounds
        Namely, phenyls galore.
        """
        verdict['N_spiro'] = rdMolDescriptors.CalcNumSpiroAtoms(mol)
        verdict['N_bridgehead'] = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
        # an `AliphaticRings` includes heterocycles.
        verdict['N_alicyclics'] = rdMolDescriptors.CalcNumAliphaticRings(mol)
        verdict['N_fused_rings'] = calc_n_fused_rings(mol)
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

    @staticmethod
    def read_cxsmiles_block(content: str, header_info: Dict[str, Any]) -> pd.DataFrame:
        """
        Reads a CXSMILES block and returns a DataFrame.
        Header info is a dict of column names and dtypes.

        content lacks the header row and instead header_info keys will be used.

        For consistency, ``Identifier``, ``SMILES``

        :param content:
        :param header_info:
        :return:
        """
        df = pd.read_csv(io.StringIO(content),
                         delimiter='\t',
                         names=list(header_info.keys()),
                         dtype=header_info).fillna(False)
        # for k in ['no_idea_1', 'no_idea_2', 'no_idea_3', 'no_idea_4', 'no_idea_5', 'no_idea_6']:
        #     df[k] = df[k].astype(bool)
        df['HBonds'] = df.HBA + df.HBD
        return df.copy()

    enamine_header_info = {
        "SMILES": str,  # default is smiles
        "Identifier": str,  # default is id
        "MW": float,
        "HAC": int,
        "LogP": float,  # default is sLogP
        "HBA": int,
        "HBD": int,
        "Rotatable_Bonds": int,  # default is RotBonds
        "FSP3": float,
        "TPSA": float,
        "lead-like": float,
        "350/3_lead-like": float,
        "fragments": float,
        "strict_fragments": float,
        "PPI_modulators": float,
        "natural_product-like": float,
        "mol_type": str,  # default is Type
        "InChIKey": str
    }
