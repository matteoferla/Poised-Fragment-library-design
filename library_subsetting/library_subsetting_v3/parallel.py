import os
import traceback
import bz2
import itertools
import json
import os
import sys
import traceback
from pathlib import Path
from types import FunctionType
from typing import List, Dict, Iterable, Optional
from pebble import ProcessPool
import pandas as pd

try:
    pd.set_option('future.no_silent_downcasting', True)
except pd.errors.OptionError:
    pass # version is old



def first_pass_process(chunk: str, filename: str, i: int, headers: List[str], out_filename_template: str, **kwargs):
    output_file = out_filename_template.format(i=i)
    classifier = CompoundSieve(mode=SieveMode.basic)
    # header_info is based off headers, but modified a bit
    df = DatasetConverter.read_cxsmiles_block('\n'.join(chunk), header_info=GPUClassifier.enamine_header_info)
    # ## Process the chunk
    verdicts = classifier.classify_df(df)
    Path(output_file).parent.mkdir(exist_ok=True, parents=True)
    if sum(verdicts.acceptable):
        for key in ['N_synthons', 'synthon_sociability', 'weighted_robogroups', 'boringness']:
            df[key] = verdicts[key]
        cols = ['SMILES', 'Identifier', 'HAC', 'HBA', 'HBD', 'Rotatable_Bonds', 'synthon_sociability', 'N_synthons', 'weighted_robogroups', 'boringness']
        txt = '\t'.join(map(str, cols)) + '\n'
        for idx, row in df.loc[verdicts.acceptable].iterrows():
            txt += '\t'.join([str(row[k]) for k in cols]) + '\n'
        with bz2.open(output_file, 'wt') as fh:
            fh.write(txt)
    # ## wrap up
    info = {'filename': filename, 'output_filename': output_file, 'chunk_idx': i,
            **verdicts.issue.value_counts().to_dict()}
    #write_jsonl(info, summary_cache)
    return info

class ParallelChunker:
    """
    The reason for this is to not overload the system with too many futures.
    As the blocks are big and I don't want to clog up /tmp with too many files.
    """
    max_workers = os.cpu_count() - 1
    exceptions_to_catch = (Exception,)

    def __init__(self,
                 chunk_size=100_000,
                 task_func: FunctionType=print):
        """
        ``filename`` and ``out_filename_template`` are arguments of ``.process_file``
        because they are file specific, ``chunk_size`` and ``task_func`` are not.
        However, ``.results`` property is shared with all.

        :param chunk_size:
        :param task_func:
        """
        self.futures = []
        self.results = []
        self.chunk_size = chunk_size
        self.task_func = task_func

    def resolve(self):
        for future in self.futures:
            try:
                self.results.append(future.result())
            except KeyboardInterrupt as e:
                raise e
            except self.exceptions_to_catch as e:
                tb = '\n'.join(traceback.format_exception(e))
                print(f"Processing failed: {e.__class__.__name__} {e}\n", tb)

    def wait(self):
        """
        Wait for all running futures to complete.

        :return:
        """
        if len(self.futures) >= self.max_workers:
            self.resolve()
            self.futures = []

    def process_file(self, filename: str,
                     **kwargs):
        """

        :param filename:
        :param kwargs: these are passed to self.task_func
        :return: summary pd.DataFrame
        """
        path = Path(filename)
        assert path.exists(), 'file does not exist'
        # process_chunk writes out...
        # ------
        with ProcessPool(max_workers=self.max_workers) as pool:
            with bz2.open(filename, 'rt') as fh:
                headers = next(fh).strip().split('\t')
                for i, chunk in enumerate(self.chunked_iterator(fh, self.chunk_size)):
                    self.wait()
                    chunk_kwargs = {'chunk': chunk,
                                    'filename': filename,
                                    'i': i,
                                    'headers': headers,
                                    **kwargs}
                    future = pool.schedule(self.task_func, kwargs=chunk_kwargs)
                    self.futures.append(future)
        self.wait()
        df = pd.DataFrame(self.results)
        return df

    @staticmethod
    def chunked_iterator(iterable: Iterable, size: int) -> Iterable:
        """Yield successive chunks of a specified size from an iterable."""
        iterator = iter(iterable)
        for first in iterator:
            yield list(itertools.chain([first], itertools.islice(iterator, size - 1)))
