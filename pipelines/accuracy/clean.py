import subprocess
import sys

result = subprocess.run(
    ["nextflow", "log", sys.argv[1], "-f", "process,workdir"],
    capture_output=True
)

process = {}
for line in result.stdout.decode("utf-8").splitlines():
    items = line.split('\t')
    if len(items) == 2:
        print(items)
