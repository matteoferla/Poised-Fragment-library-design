# chunk basic so that more multiple tasks can then be run as opposed to megafiles!
import bz2
from pathlib import Path
import sys, os, shutil

WORKINGDIR = '/opt/xchem-fragalysis-2/mferla/library_making'

sys.path.append(f'{WORKINGDIR}/repo/library_subsetting_module')

from library_subsetting_module import ParallelChunker, SieveMode, sieve_chunk

path = Path(sys.argv[1])
assert path.exists(), f'{path} does not exists'
os.makedirs('/tmp/second_pass', exist_ok=True)
master = ParallelChunker(chunk_size = 100_000, task_func=sieve_chunk)
out_filename_template=f'/tmp/second_pass/{path.stem}_chunk{{i}}.bz2'
df = master.process_file(filename=path.as_posix(),
                         out_filename_template=out_filename_template,
                        summary_cache='second_pass_summary.jsonl',
                        mode=SieveMode.substructure,
                        )
out_filename_template=f'{WORKINGDIR}/second_pass/{path.name}' # same name, diff folder
header_added = False
with bz2.open(out_filename_template, 'wt') as out_fh:
    for chunk_path in Path('/tmp/second_pass').glob('*.bz2'):
        with bz2.open(chunk_path, 'rt') as in_fh:
            # skip first line (header)
            if not header_added:
                header = next(in_fh)
                out_fh.write(header)
                header_added = True
            for line in in_fh:
                out_fh.write(line)

shutil.rmtree('/tmp/second_pass')
print(path, len(df))
df.to_csv(f'{WORKINGDIR}/csvs/{path.stem}_2nd_reduction_results.csv')