from collections import defaultdict
import subprocess
import sys
import re
import shutil

result = subprocess.run(
    ["nextflow", "log", sys.argv[1], "-f", "process,workdir"],
    capture_output=True
)

process = defaultdict(lambda: set())

for line in result.stdout.decode("utf-8").splitlines():
    items = line.split('\t')
    if len(items) == 2 and re.search(sys.argv[2], items[0]) is not None:
        process[items[0]].add(items[1])
for k, paths in process.items():
    print(k)
    for p in paths:
        shutil.rmtree(p)
