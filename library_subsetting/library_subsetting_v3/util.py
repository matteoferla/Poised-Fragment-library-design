import json, os
from typing import Any, NewType

InchiType = NewType('InchiType', str)

def write_jsonl(obj: Any, filename: str) -> None:
    with open(filename, 'a') as fh:
        fh.write(json.dumps(obj) + '\n')

def read_jsonl(filename: str) -> Any:
    if not os.path.exists(filename):
        return []
    with open(filename) as fh:
        data = []
        for line in fh:
            try:
                data.append(json.loads(line))
            except json.JSONDecodeError as error:
                pass  # burnt line
        return data
