# chunk basic so that more multiple tasks can then be run as opposed to megafiles!

from pathlib import Path
import sys

WORKINGDIR = '/opt/xchem-fragalysis-2/mferla/library_making'

sys.path.append(f'{WORKINGDIR}/repo/library_subsetting_module')

from library_subsetting_module import ParallelChunker, SieveMode, sieve_chunk

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