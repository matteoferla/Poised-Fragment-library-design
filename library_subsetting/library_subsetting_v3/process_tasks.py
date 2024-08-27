"""
Multiprocessing task, in Pebble at least, cannot be from __main__.
They need to be imported from a module.
"""

__all__ = ['sieve_chunk', 'test_process_chunk']

from . import CompoundSieve, SieveMode, DatasetConverter, write_jsonl
from typing import List
import bz2
from pathlib import Path

def sieve_chunk(chunk: List[str],
                       filename: str,
                       i: int,
                       out_filename_template: str,
                       summary_cache:str,
                       mode:SieveMode=SieveMode.synthon,
                       **kwargs):
    """
    Example usage:

    ... code-block:: python

    import traceback

        from pathlib import Path
        from typing import List
        from pathlib import Path
        import sys, os

        sys.path.append('repo/library_subsetting')

        from library_subsetting_v3 import ParallelChunker, CompoundSieve, SieveMode, DatasetConverter, first_pass_process

        CompoundSieve.cutoffs['min_hbonds_per_HAC'] = 1 / 5  # quartile
        CompoundSieve.cutoffs['max_rota_per_HAC'] = 1 / 5 # quartile (~.22)
        CompoundSieve.cutoffs['min_synthon_sociability_per_HAC'] = 0.354839  # quartile
        CompoundSieve.cutoffs['min_weighted_robogroups_per_HAC'] = 0.0838 # quartile
        CompoundSieve.cutoffs['max_boringness'] = 0

        WORKINGDIR = '/Users/user/Coding/library_exploration'

        master = ParallelChunker(chunk_size = 10_000_000, task_func=first_pass_process)
        for path in Path('/tmp').glob('Enamine_REAL_HAC_*.cxsmiles.bz2'):
            try:
                df = master.process_file(filename=path.as_posix(),
                                         out_filename_template=f'{WORKINGDIR}/first_pass_chunked/{path.stem}_chunk{{i}}.bz2',
                                        summary_cache='first_pass.jsonl',
                                        mode=SieveMode.basic,
                                        )
                print(path)
                df.to_csv(f'{WORKINGDIR}/{path.stem}_reduction_results.csv')
                del df
            except KeyboardInterrupt as beep:
                raise beep
            except Exception as error:
                print(error.__class__.__name__, error, traceback.format_exception(error))

    :param chunk:
    :param filename:
    :param i:
    :param out_filename_template:
    :param summary_cache:
    :param mode:
    :param kwargs:
    :return:
    """
    output_file = out_filename_template.format(i=i)
    classifier = CompoundSieve(mode=mode)
    # header_info is based off headers, but modified a bit
    df = DatasetConverter.read_cxsmiles_block('\n'.join(chunk), header_info=DatasetConverter.enamine_header_info)
    # ## Process the chunk
    verdicts = classifier.classify_df(df)
    Path(output_file).parent.mkdir(exist_ok=True, parents=True)
    if sum(verdicts.acceptable):
        cols = ['SMILES', 'Identifier', 'HAC', 'HBA', 'HBD', 'Rotatable_Bonds', 'synthon_sociability', 'N_synthons', 'weighted_robogroups', 'boringness']
        txt = '\t'.join(map(str, cols)) + '\n'
        for idx, row in df.loc[verdicts.acceptable].iterrows():
            txt += '\t'.join([str(row.get(k, 0.)) for k in cols]) + '\n'
        with bz2.open(output_file, 'wt') as fh:
            fh.write(txt)
    else:
        print(f"No compounds selected in {filename} chunk {i}", flush=True)
    # ## wrap up
    info = {'filename': filename, 'output_filename': output_file, 'chunk_idx': i,
            **verdicts.issue.value_counts().to_dict()}
    write_jsonl(info, summary_cache)
    return info

def test_process_chunk(chunk, *args, **kwargs):
    return f"Test: received {len(chunk)} lines ({args}, {kwargs})"
