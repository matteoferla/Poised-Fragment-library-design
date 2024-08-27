"""
This is a variant of ``amicability.py`` but stops at USRCAT.
"""

import bz2
import sys
import time
import functools
from collections import defaultdict
from pathlib import Path
from typing import Dict, NewType, Optional
import numpy as np
import pandas as pd
from library_subsetting.deprecated.library_classification import RoboDecomposer, Classifier
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem, rdMolDescriptors
from scipy import stats
import operator

RDLogger.DisableLog('rdApp.*')

InchiType = NewType('InchiType', str)
classifier = Classifier()

def read_library(filename, header_info: Optional[dict]=None) -> pd.DataFrame:
    with bz2.open(filename, 'rt') as bfh:
        headers = next(bfh).strip().split('\t')
        content = bfh.read()
    if header_info is not None:
        pass
    elif len(headers) != len(Classifier.enamine_header_info):
        # the order ought to be the same... but just in case
        header_info = {h: Classifier.enamine_header_info.get(h, str) for h in headers}
    else:
        header_info = Classifier.enamine_header_info
    print(header_info)
    return Classifier.read_cxsmiles_block(content, header_info=header_info)



def count_synthons(df: pd.DataFrame) -> Dict[InchiType, int]:
    """
    Count the number of synthons in a dataframe.

    :param df: has to have column ``SMILES``
    :return: dict of inchi to counts
    """
    robodecomposer = RoboDecomposer()
    tally = defaultdict(int)

    ex_carbamate = Chem.MolFromSmiles('[Cl]C(=O)O')
    for smiles in df.SMILES:
        try:
            mol = Chem.MolFromSmiles(smiles)
            classifier.assess_mol_patterns(mol, {})
        except Exception as error:
            continue
        for part in robodecomposer.decompose(mol):
            if part.HasSubstructMatch(ex_carbamate):  # this was fixed, but here just in case.
                continue
            tally[Chem.MolToInchi(part)]+=1
    return dict(sorted(tally.items(), key=operator.itemgetter(1), reverse=True))

def add_mol(synthons: pd.DataFrame) -> pd.DataFrame:
    """
    From inchi add ``mol`` column to the dataframe.

    :param synthons: has ``inchi``
    :return:
    """
    synthons['mol'] = synthons.inchi.apply(Chem.MolFromInchi)
    synthons.mol.apply(AllChem.EmbedMolecule)
    return synthons

def get_usrcat(mol):
    if Chem.Mol.GetNumHeavyAtoms(mol) < 3 or Chem.Mol.GetNumConformers(mol) == 0:
        return [0] * 60
    return rdMolDescriptors.GetUSRCAT(mol)

def get_weighted_USRCAT07(query_usrcat, synthons: pd.DataFrame) -> float:
    """
    Get the number of synthons USRCAT > 0.7 scaled by their counts.

    :param query_usrcat:
    :param synthons: has ``USRCAT`` and ``counts`` columns
    :return:
    """
    score = 0
    for i, row in synthons.iterrows():
        target_usrcat = row.USRCAT
        score += row.counts if rdMolDescriptors.GetUSRScore(query_usrcat, target_usrcat) > 0.7 else 0
    return score

def calculate_amicability(tally: Dict[InchiType, float], n_compounds: int) -> pd.DataFrame:
    synthons = pd.DataFrame({'inchi': tally.keys(), 'counts': tally.values()})
    #synthons = pd.DataFrame({'inchi': k, 'counts': c for k, c in tally.items() if c > 5})
    synthons = synthons.loc[synthons.counts > 4].copy()
    add_mol(synthons)  # add `mol` columns (3D embeds the molecules)
    synthons['USRCAT'] = synthons.mol.apply(get_usrcat)
    synthons['weighted_USRCAT0.7'] = synthons.USRCAT.apply(functools.partial(get_weighted_USRCAT07, synthons=synthons))
    synthons['nor_weighted_USRCAT0.7'] = synthons['weighted_USRCAT0.7'] / n_compounds
    return synthons

def get_alpha(array: np.array):
    unique, counts = np.unique(array, return_counts=True)
    # unique and counts as the same shape as these arent np.histogram edges
    # deal with shift:
    xmin = min(unique)
    filtered_counts = counts[unique >= xmin]
    filtered_unique = unique[unique >= xmin]
    # space out vector:
    start = filtered_unique.min()
    stop = filtered_unique.max()
    new_vector = np.arange(start, stop + 1)
    new_counts_vector = np.zeros(new_vector.shape, dtype=int)
    new_counts_vector[filtered_unique - start] = filtered_counts
    # go!
    alpha, _, _ = stats.powerlaw.fit(new_counts_vector, method="MLE")
    return alpha

# ========================================
def main(filename: str):
    tick = time.time()
    # ----------------------------------------
    # ## read
    df: pd.DataFrame = read_library(filename)
    n_compounds = len(df)
    print('Number of compounds:', n_compounds)
    tock = time.time()
    print('Elapsed time:', tock - tick)
    # ----------------------------------------
    # ## tally synthons
    tally: Dict[InchiType, int] = count_synthons(df)
    print('Number of unique synthons:', len(tally))
    nonsingleton_tally = {inchi: count for inchi, count in tally.items() if count > 1}
    print('Number of nonsingleton synthons:', len(nonsingleton_tally))
    tock = time.time()
    print('Elapsed time:', tock - tick)
    # ----------------------------------------
    # ## get the amicability
    synthons = pd.DataFrame({'inchi': tally.keys(), 'counts': tally.values()})
    add_mol(synthons)  # add `mol` columns (3D embeds the molecules)
    synthons['USRCAT'] = synthons.mol.apply(get_usrcat)
    synthons.to_pickle(f'{Path(filename).stem}_synthons.pkl.gz')
    print('written file', f'{Path(filename).stem}_synthons.pkl.gz')
    tock = time.time()
    print('Elapsed time:', tock-tick)
    # ----------------------------------------
    # ## Stats
    # use integer to make automatic binning easier/smoother
    alpha = get_alpha(synthons['counts'].values.astype(int))
    print('Alpha: ', alpha)
    tock = time.time()
    print('Elapsed time:', tock-tick)


if __name__ == '__main__':
    main(sys.argv[1])
