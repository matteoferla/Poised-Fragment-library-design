from typing import Dict, Any
import pandas as pd
import io

class DatasetConverter:
    # this is modded
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
