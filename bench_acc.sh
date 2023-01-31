current_dir=`pwd`
pipeline_dir="pipelines/accuracy"
cd $pipeline_dir
nextflow run main.nf -resume $@
