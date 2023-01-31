current_dir=`pwd`
pipeline_dir="pipelines/performance"
cd $pipeline_dir
nextflow run main.nf -with-report $current_dir/perf_report.html -resume -qs 1 $@
