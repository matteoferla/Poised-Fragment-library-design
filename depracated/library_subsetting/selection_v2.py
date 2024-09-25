"""
This script is used to filter a large file of CXSMILES strings
to only include those that are sociable to the RoboDecomposer+Classifier pipeline.

The huge bz2 import file is read in chunks of 100,000 lines and processed in parallel.
Temporarily,  bz2 files with the filtered chunks will be written to /tmp/output/.
The output is written to a bz2 file with the same name as the input file, but with 'selection_' prepended.

seen-synthons.jsonl and seen-synthons_extra.jsonl are inchi caches.

Version 2 does URSCAT score on torch, cf

"""

prior_seen_synthons = 'seen-synthons.jsonl'
newly_seen_synthons = 'seen-synthons_extra.jsonl'
summary_cache       = 'results-backup.jsonl'

import bz2
import itertools
import json
import os
import sys
import traceback
from pathlib import Path

import pandas as pd
from library_classification import RoboDecomposer, InchiType
from library_classification_torch import GPUClassifier
from pebble import ProcessPool
from rdkit import RDLogger
from typing import List, Dict
import numpy as np
import numpy.typing as npt

num_cpus = os.cpu_count()
pd.set_option('future.no_silent_downcasting', True)
RDLogger.DisableLog('rdApp.*')

# ========================================================================================

def write_jsonl(obj, filename):
    with open(filename, 'a') as fh:
        fh.write(json.dumps(obj) + '\n')

def read_jsonl(filename):
    if not os.path.exists(filename):
        return []
    with open(filename) as fh:
        data = []
        for line in fh:
            try:
                data.append(json.loads(line))
            except json.JSONDecodeError as error:
                pass  # burnt line
        return data

def process_chunk(chunk: str, filename: str, i: int, headers: List[str],
                  common_synthons_tally,
                    common_synthons_usrcats
                  ):
    output_file = '/tmp/output/' + Path(filename).stem + f'/filtered_chunk{i}.bz2'
    classifier = GPUClassifier(common_synthons_tally=common_synthons_tally,
                               common_synthons_usrcats=common_synthons_usrcats)
    # header_info is based off headers, but modified a bit
    df = GPUClassifier.read_cxsmiles_block('\n'.join(chunk), header_info=GPUClassifier.enamine_header_info)
    dejavu: Dict[InchiType, int]
    for dejavu in read_jsonl(prior_seen_synthons):
        classifier.dejavu_synthons.update(dejavu)
    for dejavu in read_jsonl(newly_seen_synthons):
        classifier.dejavu_synthons.update(dejavu)
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
    write_jsonl(info, summary_cache)
    write_jsonl(classifier.nuveau_dejavu_synthons, newly_seen_synthons)
    return info


def test_process_chunk(chunk, *args, **kwargs):
    return f"Processed {len(chunk)} lines"


def chunked_iterator(iterable, size):
    """Yield successive chunks of a specified size from an iterable."""
    iterator = iter(iterable)
    for first in iterator:
        yield list(itertools.chain([first], itertools.islice(iterator, size - 1)))


# ========================================================================================
# ## Process the file




class ParallelMaster:
    """
    The reason for this is to not overload the system with too many futures.
    As the blocks are big
    """
    max_workers = num_cpus - 1

    def __init__(self):
        self.futures = []
        self.results = []

    def resolve(self):
        for future in self.futures:
            try:
                self.results.append(future.result())
            except Exception as e:
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

    def process_file(self, filename: str, common_synthons_tally, common_synthons_usrcats):
        """

        :param filename:
        :return: summary pd.DataFrame
        """
        path = Path(filename)
        assert path.exists(), 'file does not exist'
        # process_chunk writes out...
        # ------
        with ProcessPool(max_workers=self.max_workers) as pool:
            with bz2.open(filename, 'rt') as fh:
                headers = next(fh).strip().split('\t')
                for i, chunk in enumerate(chunked_iterator(fh, chunk_size)):
                    self.wait()
                    # test version:
                    # process_chunk = test_process_chunk
                    args = (chunk, filename, i, headers, common_synthons_tally, common_synthons_usrcats)
                    future = pool.schedule(process_chunk, args=args)
                    self.futures.append(future)
        self.wait()
        df = pd.DataFrame(self.results)
        return df

# ========================================================================================
if __name__ == '__main__':
    # make sure all is in order
    RoboDecomposer().test()
    chunk_size = 100_000

    GPUClassifier.cutoffs['min_hbonds_per_HAC'] = 1 / 5  # quartile
    GPUClassifier.cutoffs['max_rota_per_HAC'] = 1 / 5 # quartile (~.22)
    GPUClassifier.cutoffs['min_synthon_sociability_per_HAC'] = 0.354839  # quartile
    GPUClassifier.cutoffs['min_weighted_robogroups_per_HAC'] = 0.0838 # quartile
    GPUClassifier.cutoffs['max_boringness'] = 0
    print(GPUClassifier.cutoffs)

    # 'common_synthons.pkl.gz'
    common_synthons = pd.read_pickle(sys.argv[2])  # .set_index('inchi')
    common_synthons_tally: npt.NDArray[np.int_] = common_synthons.tally.values  # Shape: (N, 1)
    common_synthons_usrcats: List[List[float]] = common_synthons.USRCAT.to_list()  # Shape: (N, 60)

    master = ParallelMaster()
    filename = sys.argv[1]
    df = master.process_file(filename,
                            common_synthons_tally=common_synthons_tally,
                            common_synthons_usrcats=common_synthons_usrcats)
    df.to_csv(f'{Path(filename).stem}_reduction_results.csv')
    # assert len(df), 'No compounds were selected'
    # for some reason the results are not being passed back to the main process
    output_file = f'{Path(filename).stem}_selection.bz2'
    # ## Combine the outputs
    n = 0
    with bz2.open(output_file, 'wt') as output_fh:
        for path in (Path('/tmp/output') / Path(filename).stem).glob('filtered_chunk*.bz2'):
            with bz2.open(path.as_posix(), 'rt') as input_fh:
                if n > 0:
                    next(input_fh)  # skip header line
                for line in input_fh:
                    n += 1
                    output_fh.write(line)

    print(f"Combined output written to {output_file} - {n} lines")
