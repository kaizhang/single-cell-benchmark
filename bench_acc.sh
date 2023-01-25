
work_dir="./work_acc"
nextflow run src/accuracy/main.nf -resume -w $work_dir $@
