"""
Multiprocessing task, in Pebble at least, cannot be from __main__.
They need to be imported from a module.
"""

__all__ = ['sieve_chunk', 'test_process_chunk', 'sieve_chunk2sdf']

from . import CompoundSieve, SieveMode, DatasetConverter, write_jsonl
from typing import List, Optional
import bz2
from pathlib import Path

def sieve_chunk(chunk: List[str],
                       filename: str,
                       i: int,
                       summary_cache:str,
                       out_filename_template: str,
                       mode:SieveMode=SieveMode.basic,
                       **kwargs):
    """
    The chunk is processed and saved to disk.

    :param chunk: the list of lines (str) to process
    :param filename: the original filename (for record keeping)
    :param i: chunk index (for record keeping and for ``filename_template.format(i=i)``)
    :param summary_cache:
    :param out_filename_template: out filename with {i} placeholder
    :param mode: ``SieveMode.basic``, ``SieveMode.substructure`` or ``SieveMode.synthon``
    :param kwargs: ParallelChunker may pass arguments that are not needed.
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
        cols = ['SMILES', 'Identifier', 'HAC', 'HBA', 'HBD', 'Rotatable_Bonds', 'synthon_sociability', 'N_synthons', 'synthon_score', 'boringness']
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

def sieve_chunk2sdf(chunk: List[str],
                   filename: str,
                   i: int,
                   summary_cache:str,
                   out_filename_template: str,
                   **kwargs):
    """
    The chunk is processed and saved to disk.

    :param chunk: the list of lines (str) to process
    :param filename: the original filename (for record keeping)
    :param i: chunk index (for record keeping and for ``filename_template.format(i=i)``)
    :param summary_cache:
    :param out_filename_template: out filename with {i} and {tier} placeholder
    :param kwargs: ParallelChunker may pass arguments that are not needed.
    :return:
    """
    output_files = {tier: out_filename_template.format(i=i, tier=tier) for tier in ['Z0-05', 'Z05-08', 'Z08-1', 'Z1']}
    classifier = CompoundSieve(mode=SieveMode.synthon_v3, use_row_info=False)
    # header_info is based off headers, but modified a bit
    df = DatasetConverter.read_cxsmiles_block('\n'.join(chunk), header_info=DatasetConverter.enamine_header_info)
    # ## Process the chunk
    verdicts = classifier.classify_df(df)
    Path(out_filename_template).parent.mkdir(exist_ok=True, parents=True)
    if sum(verdicts.acceptable):
        # 'Z0-05', 'Z05-08', 'Z08-1', 'Z1'
        masks = {'Z0-05': (verdicts.combined_Zscore >= 0.) & (verdicts.combined_Zscore < 0.5),
                 'Z05-08': (verdicts.combined_Zscore >= 0.5) & (verdicts.combined_Zscore < 0.8),
                  'Z08-1': (verdicts.combined_Zscore >= 0.8) & (verdicts.combined_Zscore < 1.),
                  'Z1': (verdicts.combined_Zscore >= 1.)
                 }
        for tier, mask in masks.items():
            with bz2.open(output_files[tier], 'wt') as fh:
                for sdfblock in verdicts.loc[verdicts.acceptable & mask].sdfblock:
                    if isinstance(sdfblock, str):
                        fh.write(sdfblock) # the $$$$\n is already in the sdfblock end
    else:
        print(f"No compounds selected in {filename} chunk {i}", flush=True)
    # ## wrap up
    info = {'filename': filename, 'output_filename': out_filename_template.format(i=i, tier='-'), 'chunk_idx': i,
            **verdicts.issue.value_counts().to_dict()}
    write_jsonl(info, summary_cache)
    return info

def test_process_chunk(chunk, *args, **kwargs):
    return f"Test: received {len(chunk)} lines ({args}, {kwargs})"
