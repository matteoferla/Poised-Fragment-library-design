# chunk basic so that more multiple tasks can then be run as opposed to megafiles!

import traceback
from pathlib import Path
from typing import List
from pathlib import Path
import sys, os

WORKINGDIR = '/opt/xchem-fragalysis-2/mferla/library_making'

sys.path.append(f'{WORKINGDIR}/repo/library_subsetting')

from library_subsetting_v3 import ParallelChunker, CompoundSieve, SieveMode, DatasetConverter, sieve_chunk

CompoundSieve.cutoffs['min_hbonds_per_HAC'] = 1 / 5  # quartile
CompoundSieve.cutoffs['max_rota_per_HAC'] = 1 / 5 # quartile (~.22)
CompoundSieve.cutoffs['min_synthon_sociability_per_HAC'] = 0.354839  # quartile
CompoundSieve.cutoffs['min_weighted_robogroups_per_HAC'] = 0.0838 # quartile
CompoundSieve.cutoffs['max_boringness'] = 0

path = Path(sys.argv[1])
assert path.exists(), f'{path} does not exists'

master = ParallelChunker(chunk_size = 10_000_000, task_func=sieve_chunk)
out_filename_template=f'{WORKINGDIR}/first_pass/{path.stem}_chunk{{i}}.bz2'
df = master.process_file(filename=path.as_posix(),
                         out_filename_template=out_filename_template,
                        summary_cache='first_pass_summary.jsonl',
                        mode=SieveMode.basic,
                        )
print(path, len(df))
df.to_csv(f'{WORKINGDIR}/csvs/{path.stem}_reduction_results.csv')