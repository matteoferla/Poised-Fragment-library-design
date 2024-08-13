"""
Ligity is not open source. It uses an algorithm called PIP.
Here is it's reimplementation.
"""

import itertools
import os
import numpy as np
import numpy.typing as npt
from rdkit import RDConfig, Chem, Geometry
from rdkit.Chem import AllChem
from rdkit.Chem.FeatMaps import FeatMaps
from typing import Tuple, Sequence, Dict, List, Iterable
import warnings

PIPType = Dict[Tuple[str, str, str], npt.NDArray[int]]


class Pipiteur:

    def __init__(self, order: int = 3, min_d: int = 2, max_d: int = 8, resolution: float = 1.):
        self.order = int(order)
        if self.order > 3:
            warnings.warn('Higher order was not tested => not bothered with')
            # for 4th order? I need to read up as there will lots of pairings
        self.n_bins = int((max_d - min_d) / resolution) + 2
        self.bin_edges = np.linspace(min_d, max_d, num=self.n_bins - 1)

    fdef: AllChem.MolChemicalFeatureFactory = AllChem.BuildFeatureFactory(
        os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))
    wanted = ['Acceptor', 'Aromatic', 'Donor', 'NegIonizable', 'PosIonizable']

    @classmethod
    def get_N_feats(self, mol: Chem.Mol) -> int:
        feats: Sequence[AllChem.MolChemicalFeature] = self.fdef.GetFeaturesForMol(mol)
        feats = [feat for feat in feats if feat.GetFamily() in self.wanted]
        return len(feats)

    def __call__(self, mol: Chem.Mol) -> PIPType:
        key_combinations = list(itertools.combinations_with_replacement(self.wanted, r=self.order))
        # returned variable:
        pip: PIPType = {keys: np.zeros([self.n_bins] * self.order, dtype=int) for keys in key_combinations}
        # compute
        feats: Sequence[AllChem.MolChemicalFeature] = self.fdef.GetFeaturesForMol(mol)
        feats = sorted([feat for feat in feats if feat.GetFamily() in self.wanted], key=lambda feat: feat.GetFamily())
        # this is just [(0, 1), (0, 2), (1, 2)] for 3rd order.
        # r=2 for all as it is for distance do not inadvertently correct!
        pairings = list(itertools.combinations(list(range(0, self.order)), r=2))
        combofeats: List[AllChem.MolChemicalFeature]  # noqa len == order
        for combofeats in itertools.combinations(feats, r=self.order):
            keys: Tuple[str, str, str] = tuple([feat.GetFamily() for feat in combofeats])  # noqa
            distances: List[float] = [combofeats[i].GetPos().Distance(combofeats[j].GetPos()) for i, j in pairings]
            if any([d == 0. for d in distances]):
                # one of the atoms has two or more behaviours
                continue
            #print(keys, distances)
            pip[keys][tuple(np.digitize([distances], bins=self.bin_edges)[0])] = 1
        return pip

    def describe_pip(self, pip: PIPType) -> None:
        print('Pharmacophoric shapes: ', self.flatten(pip).sum())
        for k in pip:
            nonzeros = list(zip(*np.nonzero(pip[k])))
            if nonzeros:
                print(k, nonzeros)

    def flatten(self, pip: PIPType) -> np.ndarray:
        return np.hstack(list(pip.values())).flatten()

    def score_pips(self, pip1: PIPType, pip2: PIPType) -> float:
        """
        Tanimoto
        """
        flat_pip1 = (self.flatten(pip1) > 0).astype(int)
        flat_pip2 = (self.flatten(pip2) > 0).astype(int)
        intersection = np.dot(flat_pip1, flat_pip2)
        union = np.sum(flat_pip1) + np.sum(flat_pip2) - intersection
        if union == 0.:
            return float('nan')
        return intersection / union