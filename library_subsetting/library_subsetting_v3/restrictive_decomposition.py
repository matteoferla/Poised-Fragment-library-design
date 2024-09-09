__all__ = ['RestrictiveDecomposer', 'InchiType']

import itertools
import yaml
from pathlib import Path
from typing import List, Dict, Any, Optional, NewType, Union, Sequence
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, AllChem, rdDeprotect, Draw

InchiType = NewType('InchiType', str)


class RxnDetails:
    def __init__(self, name: str,
                 smarts: str,
                 weight:float=1.,
                 group:str='none',
                 reaction_name:str=''):
        self.name = name
        self.rxn_smarts = smarts
        self.weight = weight
        self.reaction_name = reaction_name
        self.group = group
        # filled in manually due to faff and for testing only anyway...
        self.valids = []
        self.invalids = []
        # compute AllChem.ChemicalReaction
        self.rxn: AllChem.ChemicalReaction = AllChem.ReactionFromSmarts(smarts)
        self.rxn.Initialize()
        self.tally = 0 # internal purposes

class RestrictiveDecomposer:
    """
    Creates synthons for reactions the XChem robot can do.
    amidation, Schotten-Baumann, sulfo Schotten-Baumann, Suzuki, but not Chan-Lam, Williamson, Borch, Huisgen etc.
    It is lactam/sulfam safe, which is mainly why I am writing it backwards from scratch.
    ``synthons = RoboDecomposer.decompose(mol)``
    """

    def __init__(self,
                yaml_path: Optional[str]=None,
                only_reactions: Optional[Sequence[str]] = None,
                only_groups: Optional[Sequence[str]] = None):
        # order matter as they will be run in that way
        self.rxns: List[RxnDetails] = self.load_reactions(yaml_path=yaml_path,
                                                          only_reactions=only_reactions,
                                                          only_groups=only_groups)


    @staticmethod
    def load_reactions(yaml_path: Optional[str]=None,
                       only_reactions: Optional[Sequence[str]] = None,
                       only_groups: Optional[Sequence[str]] = None) -> List[RxnDetails]:
        """
        Load reactions from a yaml file.
        """
        # ## read data
        if yaml_path is None:
            yaml_path = Path(__file__).parent / 'data' / 'decompositions.yaml'
        with Path(yaml_path).open('r') as fh:
            data = yaml.safe_load(fh)
        # ## Parse data
        rxns:List[RxnDetails] = []
        for key, datum in data.items():
            if only_groups and datum.get('group') not in only_groups:
                continue
            elif only_reactions and key not in only_reactions:
                continue
            else:
                rxn = RxnDetails(name=key,
                                 smarts=datum['smarts'],
                                 weight=datum.get('weight', 1.),
                                 group=datum.get('group', 'none'),
                                 reaction_name=datum.get('reaction', ''))
                rxns.append(rxn)
        return rxns

    def apply_rxn(self, mol: Chem.Mol, rxn_details: RxnDetails) -> List[Chem.Mol]:
        """
        Returns a list of Chem.Mol
        """
        try:
            AllChem.SanitizeMol(mol)  # , catchErrors=True
            if rxn_details.rxn.IsMoleculeReactant(mol):
                prods = rxn_details.rxn.RunReactant(mol, 0)
                if len(prods) > 1:
                    subprods = [subprod for prod in prods[0] for subprod in self.apply_rxn(mol=prod,
                                                                                                   rxn_details=rxn_details)]
                else:
                    subprods = prods[0]
                for prod in subprods:
                    AllChem.SanitizeMol(prod)
                    done_rxn_names = prod.GetProp('Reaction').split() if prod.HasProp('Reaction') else []
                    done_rxn_names.append(rxn_details.name)
                    prod.SetProp('Reaction', ' '.join(done_rxn_names))
                return subprods
        except KeyboardInterrupt as error:
            raise error
        except Exception as error:
            pass
        return [mol]

    def decompose(self, mol: Chem.Mol) -> List[Chem.Mol]:
        """
        Returns a list of Chem.Mol (synthons).
        See property ``Reactions`` for the cleaved moiety.
        """
        synthons: List[Chem.Mol] = [mol]
        rxn: RxnDetails
        for rxn_details in self.rxns:
            synthon_groups = [self.apply_rxn(mol=synthon, rxn_details=rxn_details) for synthon in synthons]
            synthons = [synthon for synthons in synthon_groups for synthon in synthons]
        return synthons

    def tally_groups(self, mol: Chem.Mol) -> List[str]:
        rxn_details: RxnDetails
        moieties = []
        for rxn_details in self.rxns:
            if rxn_details.rxn.IsMoleculeReactant(mol):
                moieties.append(rxn_details.name)
        return moieties

    def synthon_score(self, mol: Chem.Mol) -> float:
        """
        Returns a score based on the weighted synthons.
        """
        return sum(rxn.weight for rxn in self.rxns if rxn.rxn.IsMoleculeReactant(mol))

    @classmethod
    def test(cls):
        robodecomposer = cls()
        results = []
        for name, smi in dict(amide=('CCNC(=O)CC', 2),
                              sulfonamide=('CCNS(=O)(=O)CC', 2),
                              diglycine=('NCC(=O)NCC(=O)O', 2),
                              triglycine=('NCC(=O)NCC(=O)NCC(=O)O', 3),
                              lactam=('C1NC(=O)CC1', 1),
                              tertamide=('CCN(C)C(=O)CC', 3),
                              cycloalkyl_amide=('N1(CCCCC1)C(=O)CC', 2),
                              aryl_amide=('n1(cccc1)C(=O)CC', 2),
                              sec_amine=('CCNCC', 2),
                              saccharin=('O=C2c1ccccc1S(=O)(=O)N2', 1),
                              urea=('CCNC(=O)NCC', 2),
                              biaryl=('n1(cccc1)-c1(ccccc1)', 2)
                              ).items():
            mol = Chem.MolFromSmiles(smi[0])
            synthons: List[Chem.Mol] = robodecomposer.decompose(mol)
            if len(synthons) != smi[1]:
                print(f'{name} failed. Expected {smi[1]} got {len(synthons)}')
                from IPython.display import display
                display(mol)
                display( Draw.MolsToGridImage(synthons, legends=[s.GetProp('Reaction') for s in synthons]) )
            results.append(dict(name=name, input=mol, output=synthons, expected_products=smi[1]))
            return results


