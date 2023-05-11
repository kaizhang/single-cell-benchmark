nextflow.enable.dsl=2

params.outdir = 'results'

include { bench_ATAC } from './ATAC/main.nf'
include { bench_RNA } from './RNA/main.nf'
include { bench_HiC } from './HiC/main.nf'

workflow ATAC {
    bench_ATAC()
}

workflow RNA {
    bench_RNA()
}

workflow HiC {
    bench_HiC()
}

workflow {
    ATAC()
    RNA()
    HiC()
}