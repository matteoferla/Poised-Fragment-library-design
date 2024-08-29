This will be incorportated into docs

```python
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
```

```python

from library_subsetting_v3 import ParallelChunker

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

```