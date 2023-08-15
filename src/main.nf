nextflow.enable.dsl=2

params.outdir = 'results'

include { bench_atac; bench_atac_simulated; } from './ATAC/main.nf'
include { bench_rna } from './RNA/main.nf'
include { bench_hic } from './HiC/main.nf'
include { parseYamlString } from './common/utils.gvy'

datasets = Channel.from(
    parseYamlString(file(params.input_data).text)
)

workflow atac_simulated {
    bench_atac_simulated()
}

workflow {
    if (params.pipelines == null || params.pipelines.contains('ATAC')) {
        bench_atac(datasets)
    }
    
    if (params.pipelines == null || params.pipelines.contains('RNA')) {
        bench_rna(datasets)
    }

    if (params.pipelines == null || params.pipelines.contains('HiC')) {
        bench_hic()
    }
}