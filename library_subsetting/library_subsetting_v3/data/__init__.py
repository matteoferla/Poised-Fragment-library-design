import importlib.resources
import pickle, gzip
import json
from typing import Optional
from pathlib import Path
from rdkit.Chem import AllChem, RDConfig
from ..util import ultranormalize, autopass_fun
import functools


def read_data_file() -> str:
    with importlib.resources.open_text(__package__, 'data.txt') as file:
        return file.read()

def read_json(filename: str) -> dict:
    with open(filename, 'r') as fh:
        return json.load(fh)

def read_pickle(filename: str) -> dict:
    if filename.endswith('.gz'):
        with gzip.open(filename, 'rb') as gfh:
            return pickle.load(gfh)
    else:
        return pickle.load(open(filename, 'rb'))

def read_MolChemicalFeatureFactory(filename: Optional[str]=None) -> AllChem.MolChemicalFeatureFactory:
    if filename is None:
        path = Path(RDConfig.RDDataDir) / 'BaseFeatures.fdef'
    else:
        path = Path(__file__).parent / filename
    return AllChem.BuildFeatureFactory(path.as_posix())

def parse_unskew_funs(filename='likelihood_skew_params.json'):
    """
    Remove the skew from the log freq of pipi data
    """
    likelihood_skew_params = read_json(filename)
    fundex = {}
    for key in likelihood_skew_params:
        upper_bound = likelihood_skew_params[key]['upper_bound']
        if likelihood_skew_params[key]['autopass']:
            fundex[key] = functools.partial(autopass_fun, bound=upper_bound)
        else:
            fundex[key] = functools.partial(ultranormalize,
                                           skew_shape=likelihood_skew_params[key]['alpha'],
                                           skew_loc=likelihood_skew_params[key]['loc'],
                                           skew_scale=likelihood_skew_params[key]['scale'],
                                           nan_replacement=-2,
                                           nor_bound=2,
                                           flip_sign=-1)
