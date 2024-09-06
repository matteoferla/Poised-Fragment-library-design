import importlib.resources
import pickle, gzip
import json
from typing import Optional
from pathlib import Path
from rdkit.Chem import AllChem, RDConfig


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
