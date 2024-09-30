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

class ParallelChunker:
    """
    The reason for this is to not overload the system with too many futures.
    As the blocks are big and I don't want to clog up /tmp with too many files.
    """
    max_workers = os.cpu_count() - 1
    exceptions_to_catch = (Exception,)

    def __init__(self,
                 chunk_size=100_000,
                 task_func: FunctionType=print,
                 verbose: bool=False):
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
        self.verbose = False

    def resolve(self):
        for future in self.futures:
            try:
                self.results.append(future.result())
            except KeyboardInterrupt as e:
                raise e
            except self.exceptions_to_catch as e:
                tb = '\n'.join(traceback.format_exception(e))
                if self.verbose:
                    print(f"Task failed but caught: {e.__class__.__name__} {e}\n", tb)
                self.results.append({'error': f'{e.__class__.__name__} {e}'})

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
        self.resolve()  # wait for all futures to complete
        self.futures = []
        df = pd.DataFrame(self.results)
        return df

    @staticmethod
    def chunked_iterator(iterable: Iterable, size: int) -> Iterable:
        """Yield successive chunks of a specified size from an iterable."""
        iterator = iter(iterable)
        for first in iterator:
            yield list(itertools.chain([first], itertools.islice(iterator, size - 1)))
