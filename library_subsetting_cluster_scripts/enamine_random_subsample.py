"""
A small snippet to subsample the very large Enamine REAL dataset.
See Enamine downloader in https://github.com/matteoferla/Fragment-hit-follow-up-chemistry
for what this extends
"""

import re, random, bz2
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import List

REAL = [
    'Enamine_REAL_HAC_11_21_666M_CXSMILES.cxsmiles.bz2',
    'Enamine_REAL_HAC_22_23_828M_CXSMILES.cxsmiles.bz2',
    'Enamine_REAL_HAC_24_686M_CXSMILES.cxsmiles.bz2',
    'Enamine_REAL_HAC_25_789M_CXSMILES.cxsmiles.bz2',
    'Enamine_REAL_HAC_26_766M_CXSMILES.cxsmiles.bz2',
    'Enamine_REAL_HAC_27_872M_CXSMILES.cxsmiles.bz2',
    'Enamine_REAL_HAC_28_803M_CXSMILES.cxsmiles.bz2',
    'Enamine_REAL_HAC_29_38_1.3B_Part_1_CXSMILES.cxsmiles.bz2',
    'Enamine_REAL_HAC_29_38_1.3B_Part_2_CXSMILES.cxsmiles.bz2'
    ]

def read_random_lines(fh, num_lines: int) -> list[str]:
    """
    Reservoir sampling
    """
    selected_lines = []
    for line_number, line in enumerate(fh):
        if line_number < num_lines:
            selected_lines.append(line)
        else:
            r = random.randint(0, line_number)
            if r < num_lines:
                selected_lines[r] = line
    return selected_lines

def process_file(file_name: str, total: int, wanted: int) -> List[str]:
    try:
        _mantissa, _exponent = re.match(r'Enamine_REAL_HAC_.*_([\d\.]+)(\w)', Path(file_name).name).groups()
        mantissa = float(_mantissa)
        exponent = 1e6 if _exponent == 'M' else 1e9
        subtotal = int(mantissa * exponent)
        target = int(subtotal / total * wanted) + 1
        print(file_name, target, subtotal, flush=True)
    except Exception as error:
        print(file_name, error.__class__.__name__, str(error), flush=True)
        return []
    path = Path(file_name, flush=True)
    if not path.exists():
        print(f'{path} missing')
        return []
    try:
        with bz2.open(path, 'rt') as fh:
            header = next(fh).split('\t')
            return read_random_lines(fh, target)
    except Exception as error:
        print(f'Error with {file_name}')
        return []

def parallel_process_files(file_names: List[str], total: int, wanted: int) -> List[str]:
    sampled = []
    with ProcessPoolExecutor() as executor:
        future_to_file = {executor.submit(process_file, file_name, total, wanted): file_name for file_name in file_names}
        for future in as_completed(future_to_file):
            sampled.extend(future.result())
    return sampled

if __name__ == '__main__':
    total = 6.75e9
    wanted = 100_000
    paths = [f'/tmp/{f}' for f in REAL]
    sampled = parallel_process_files(paths, total, wanted)

    with bz2.open(f'Enamine_REAL_{wanted}_sampled.bz2', 'wt') as fh:
        for line in sampled:
            fh.write(line)
    print('Done')
