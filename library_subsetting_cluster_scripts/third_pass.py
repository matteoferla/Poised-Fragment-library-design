# chunk basic so that more multiple tasks can then be run as opposed to megafiles!
import bz2
from pathlib import Path
import sys, os, shutil
import contextlib

WORKINGDIR = '/opt/xchem-fragalysis-2/mferla/library_making'
os.chdir(WORKINGDIR)
sys.path.append(f'{WORKINGDIR}/repo/library_subsetting_module')

from library_subsetting_module import ParallelChunker, SieveMode, sieve_chunk2

sdf_mode = False
path = Path(sys.argv[1])
assert path.exists(), f'{path} does not exists'
os.makedirs('/tmp/third_pass', exist_ok=True)
master = ParallelChunker(chunk_size = 100_000, task_func=sieve_chunk2)
out_filename_template=f'/tmp/third_pass/{path.stem}_{{tier}}_chunk{{i}}.bz2'
df = master.process_file(filename=path.as_posix(),
                         out_filename_template=out_filename_template,
                        summary_cache=f'{WORKINGDIR}/third_pass_summary.jsonl',
                        mode=SieveMode.synthon,
                        )
for tier in ['Zn2-n1', 'Zn1-n05', 'Zn05-0', 'Z0-05', 'Z05-08', 'Z08-1', 'Z1']:
    out_filename=f'{WORKINGDIR}/third_pass/{tier}/{path.name}'
    if sdf_mode:
        out_file = out_filename.replace('.cxsmiles.bz2', '.sdf.bz2')
    os.makedirs(Path(out_filename).parent, exist_ok=True)
    with bz2.open(out_filename, 'wt') as out_fh:
        for chunk_path in Path('/tmp/third_pass').glob(f'*{tier}*.bz2'):
            with contextlib.suppress(Exception):
                with bz2.open(chunk_path, 'rt') as in_fh:
                    for line in in_fh:
                        out_fh.write(line)

shutil.rmtree('/tmp/third_pass')
print(path, len(df))
df.to_csv(f'{WORKINGDIR}/csvs/{path.stem}_3nd_reduction_results.csv')