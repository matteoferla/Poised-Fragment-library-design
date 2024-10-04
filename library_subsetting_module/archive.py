__all__ = ['SimultaneousArchiveReader',
           'HashSimultaneousArchiveReader',
           'HistorySimultaneousArchiveReader',
           'CounterArchiveReader',
           'SequentialArchiveWriter']

"""
Merging large sorted datasets, while avoiding duplicates.

Class ``SimultaneousArchiveReader`` reads multiple files simultaneously and yields the lines in sorted order.
The inputs are CXSMILES files, which are tab-separated files with the SMILES in the first column and a score in the last column.
The files are read in simultaneously and the lines are sorted by the score in the last column.
Three subclasses were implemented to avoid duplicates:

- ``HashSimultaneousArchiveReader`` uses a hash of the SMILES to avoid duplicates.
- ``HistorySimultaneousArchiveReader`` uses a history of the last 10,000 SMILES to avoid duplicates.
- ``CounterArchiveReader`` uses a list of SMILES and a history to avoid duplicates.

The history approaches assume that the scores are the same for duplicates,
which is not actually true, due to conformer generation for pip scores.

Class ``SequentialArchiveWriter`` writes the output files in a rolling manner:
if the file exceeds 1 million lines, a new file is started.

Example usage:

.. code-block:: python

import bz2
import time
import json

inpaths = []
for inpath in Path('third_pass').glob('Z1/*.bz2'):
    inpaths.append(inpath)

outarx = SequentialArchiveWriter(17)
inarx = HistorySimultaneousArchiveReader(inpaths)

for line in inarx:
    outarx.write(line)

print(outarx.current_tally, inarx.duplicate_tally)
del inarx   # close the files. I was too lazy to implement a context manager
del outarx
"""


import numpy as np
import numpy.typing as npt
import bz2
from typing import List, TextIO
from pathlib import Path
import bz2
import io
import os
import sys
import time
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, PandasTools
import pandas as pd

class SimultaneousArchiveReader:
    """
    The files in ``inpaths`` are read simultaneously. The lines are sorted by the score in the last column.

    This is an iterator that yields the next line in the sorted order.

    The method ``assess_duplicate`` is called to determine if a line is a duplicate, if True, the next line is read.
    """

    def __init__(self, inpaths: List[Path]):
        self.handle_names: List[str] = [inpath.stem for inpath in inpaths]
        self.handles: List[TextIO] = [bz2.open(inpath, 'rt') for inpath in inpaths]
        self.current_lines: List[str] = [''] * len(inpaths)
        self.current_scores: List[float] = [-666] * len(inpaths)
        self.duplicate_tally = 0
        for i in range(len(inpaths)):
            self.move_forward(i)

    def __iter__(self):
        return self

    def __next__(self):
        m = max(self.current_scores)
        if m == -666:  # max w/ nan is nan
            raise StopIteration
        i = self.current_scores.index(m)
        line = self.current_lines[i]
        self.current_lines[i] = ''
        self.move_forward(i)
        return line

    def close(self):
        for fh in self.handles:
            fh.close()

    def __del__(self):
        self.close()

    def move_forward(self, i):
        try:
            newline = next(self.handles[i])
        except StopIteration:
            self.current_scores[i] = -666
            print(self.handle_names[i], 'exhausted')
            return None
        parts = newline.strip().split('\t')
        smiles: str = parts[0]
        if smiles == 'SMILES':
            # rogue header!
            return self.move_forward(i)
        elif self.assess_duplicate(smiles):
            return self.move_forward(i)
        else:
            self.current_lines[i] = newline
            self.current_scores[i] = float(parts[-1])
            return None

    def assess_duplicate(self, smiles: str) -> bool:
        # for subclasses to implement
        pass

class HashSimultaneousArchiveReader(SimultaneousArchiveReader):
    def __init__(self, inpaths: List[Path]):
        self.seen: npt.NDArray[np.int64] = np.empty(0)
        self.new: List[int] = []  # to be added to seen
        super().__init__(inpaths)


    def assess_duplicate(self, smiles: str) -> bool:
        hashed = hash(smiles)
        if bool(np.isin(hashed, self.seen)) or hashed in self.new:
            # skip: seen before!
            self.duplicate_tally += 1
            return True
        self.new.append(hashed)
        if len(self.new) > 1e3:  # move new to array
            print('added 1,000 more')
            self.seen = np.concatenate((self.seen, self.new))
            self.new = []
        return False

class HistorySimultaneousArchiveReader(SimultaneousArchiveReader):
    def __init__(self, inpaths: List[Path], history_size: int = 10_000):
        self.history = list(' ' * history_size)
        self.history_i = 0
        super().__init__(inpaths)

    def assess_duplicate(self, smiles: str) -> bool:
        """history deduplication, assuming they score the same!"""
        if smiles in self.history:
            self.duplicate_tally += 1
            return True
        else:
            self.history_i += 1
            if self.history_i >= len(self.history):
                self.history_i = 0
            self.history[self.history_i] = smiles
            return False

class CounterArchiveReader(HistorySimultaneousArchiveReader):
    def __init__(self, inpaths: List[Path], counter_smiles: List[str], history_size: int = 10_000):
        self.counter_smiles = counter_smiles
        super().__init__(inpaths, history_size)

    def assess_duplicate(self, smiles: str) -> bool:
        if not super().assess_duplicate(smiles) and smiles in self.counter_smiles:
            self.duplicate_tally += 1
            return True
        return False

# ------------------------------------------------------------------------------

class SequentialArchiveWriter:

    def __init__(self, i=0, template:str='additional_finalists/shorlist{i:0>4}.1M.cxsmiles.bz2'):
        self.i = 0
        self.template = template
        self.current_fh = None
        self.current_tally = 0
        self.start(i)

    def start(self, i=None):
        if self.current_fh is not None:
            self.current_fh.close()
        self.i += 1
        if i:
            self.i = i
        self.current_fh = bz2.open(self.template.format(i=self.i), 'wt')
        headers = 'SMILES	Identifier	HAC	HBA	HBD	Rotatable_Bonds	boringness	synthon_score	pip_common_mean	pip_uncommon_mean	combined_Zscore\n'
        self.current_fh.write(headers)
        return self.current_fh

    def __del__(self):
        self.current_fh.close()

    def write(self, text: str):
        self.current_fh.write(text.strip() + '\n')
        self.current_tally += 1
        if self.current_tally >= 1e6:
            self.current_tally = 0
            self.start()