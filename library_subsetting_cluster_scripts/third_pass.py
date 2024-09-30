# chunk basic so that more multiple tasks can then be run as opposed to megafiles!
import bz2
from pathlib import Path
import sys, os, shutil
import contextlib

WORKINGDIR = '/opt/xchem-fragalysis-2/mferla/library_making'
os.chdir(WORKINGDIR)
sys.path.append(f'{WORKINGDIR}/repo')

from library_subsetting_module import ParallelChunker, SieveMode, sieve_chunk2

sdf_mode = False
infile_path = Path(sys.argv[1])
assert infile_path.exists(), f'{infile_path} does not exists'
os.makedirs('/tmp/third_pass', exist_ok=True)
master = ParallelChunker(chunk_size = 100_000, task_func=sieve_chunk2)
out_filename_template=f'/tmp/third_pass/{infile_path.stem}_{{tier}}_chunk{{i}}.bz2'
df = master.process_file(filename=infile_path.as_posix(),
                         out_filename_template=out_filename_template,
                        summary_cache=f'{WORKINGDIR}/third_pass_summary.jsonl',
                        mode=SieveMode.synthon,
                        )
headers = ['SMILES', 'Identifier', 'HAC', 'HBA', 'HBD', 'Rotatable_Bonds', 'boringness', 'synthon_score',
               'pip_common_mean', 'pip_uncommon_mean', 'combined_Zscore']
for tier in ['Zn2-n1', 'Zn1-n05', 'Zn05-0', 'Z0-05', 'Z05-08', 'Z08-1', 'Z1']:
    out_filename=f'{WORKINGDIR}/third_pass/{tier}/{infile_path.name}'
    if sdf_mode:
        out_file = out_filename.replace('.cxsmiles.bz2', '.sdf.bz2')
    os.makedirs(Path(out_filename).parent, exist_ok=True)
    with bz2.open(out_filename, 'wt') as out_fh:
        out_fh.write('\t'.join(headers) + '\n')
        for chunk_path in Path('/tmp/third_pass').glob(f'{infile_path.stem}_{tier}*.bz2'):
            with contextlib.suppress(Exception):  # why would this happen??
                with bz2.open(chunk_path, 'rt') as in_fh:
                    for line in in_fh:
                        out_fh.write(line)

shutil.rmtree('/tmp/third_pass')
print(infile_path, len(df))
df.to_csv(f'{WORKINGDIR}/csvs/{infile_path.stem}_3nd_reduction_results.csv')