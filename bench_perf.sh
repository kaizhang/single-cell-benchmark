work_dir="./work_perf"
nextflow run src/performance/main.nf -w $work_dir -with-trace -resume -qs 1 $@
