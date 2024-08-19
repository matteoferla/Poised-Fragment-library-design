"""
> This script is deprecated. Use selection_v2.py instead.

This script is used to filter a large file of CXSMILES strings
to only include those that are sociable to the RoboDecomposer+Classifier pipeline.

The huge bz2 import file is read in chunks of 100,000 lines and processed in parallel.
Temporarily,  bz2 files with the filtered chunks will be written to /tmp/output/.
The output is written to a bz2 file with the same name as the input file, but with 'selection_' prepended.
"""

import bz2
import itertools
import json
import os
import sys
import traceback
from pathlib import Path

import pandas as pd
from library_classification import Classifier, RoboDecomposer
from pebble import ProcessPool
from rdkit import RDLogger
from typing import List

num_cpus = os.cpu_count()
pd.set_option('future.no_silent_downcasting', True)
RDLogger.DisableLog('rdApp.*')
RoboDecomposer().test()

sociability = json.loads(Path('sociability.json').read_text())

chunk_size = 100_000
classifier = Classifier(sociability)


def write_jsonl(obj, filename):
    with open(filename, 'a') as fh:
        fh.write(json.dumps(obj) + '\n')


def process_chunk(chunk:str, filename:str, i:int, headers: List[str]):
    output_file = '/tmp/output/' + Path(filename).stem + f'/filtered_chunk{i}.bz2'
    # header_info is based off headers, but modified a bit
    df = Classifier.read_cxsmiles_block('\n'.join(chunk),header_info=Classifier.enamine_header_info)
    verdicts = classifier.classify_df(df)
    Path(output_file).parent.mkdir(exist_ok=True, parents=True)
    if sum(verdicts.acceptable):
        for key in ['N_synthons', 'synthon_sociability', 'weighted_robogroups']:
            df[key] = verdicts[key]
        txt = '\t'.join(map(str, headers)) + '\n'
        cols = ['SMILES', 'Identifier', 'synthon_sociability', 'N_synthons', 'weighted_robogroups']
        for idx, row in df.loc[verdicts.acceptable].iterrows():
            txt += '\t'.join([str(row[k]) for k in cols]) + '\n'
        with bz2.open(output_file, 'wt') as fh:
            fh.write(txt)
    info = {'filename': filename, 'chunk_idx': i, **verdicts.issue.value_counts().to_dict()}
    write_jsonl(info, 'results-backup.jsonl')
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
max_workers = num_cpus - 1
path = Path(sys.argv[1])
filename = path.as_posix()
print(filename)
assert path.exists(), 'file does not exist'
output_file = f'selection_{path.stem}.bz2'

results = []
def resolve(futures):
    for future in futures:
        try:
            results.append(future.result())
        except Exception as e:
            tb = '\n'.join(traceback.format_exception(e))
            print(f"Processing failed: {e.__class__.__name__} {e}\n {tb}")



with ProcessPool(max_workers=max_workers) as pool:
    futures = []
    with bz2.open(filename, 'rt') as fh:
        headers = next(fh).strip().split('\t')
        for i, chunk in enumerate(chunked_iterator(fh, chunk_size)):
            if len(futures) >= max_workers:
                resolve(futures)
                futures = []
            # test version:
            # process_chunk = test_process_chunk
            future = pool.schedule(process_chunk, args=(chunk, filename, i, headers))
            futures.append(future)

df = pd.DataFrame(results)
df.to_csv(f'{path.stem}_reduction_results.csv')

# ========================================================================================
# ## Combine the outputs
n = 0
with bz2.open(output_file, 'wt') as output_fh:
    for path in Path('/tmp/output/' + path.stem).glob('*.bz2'):
        with bz2.open(path, 'rt') as input_fh:
            if n > 0:
                next(input_fh) # skip header line
            for line in input_fh:
                n += 1
                output_fh.write(line)

print(f"Combined output written to {output_file} - {n} lines")
