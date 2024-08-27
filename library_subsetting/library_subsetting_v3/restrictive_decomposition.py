__all__ = ['RestrictiveDecomposer', 'InchiType']

import itertools
from typing import List, Dict, Any, Optional, NewType, Union
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, AllChem, rdDeprotect

InchiType = NewType('InchiType', str)


class RestrictiveDecomposer:
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