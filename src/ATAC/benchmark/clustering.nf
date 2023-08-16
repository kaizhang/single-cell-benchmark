nextflow.enable.dsl=2

params.outdir = "results"

include { dim_reduct_cosine as snapatac2 } from '../software/snapatac2.nf'
include { bench_clustering } from '../../common/benchmark.nf' params(resultDir: "${params.outdir}/ATAC/clustering")
include { knn_leiden_exact;
          knn_leiden_hora;
          knn_leiden_exact_weighted;
          kmeans;
        } from '../../common/algorithms/clustering.nf'

include { json } from '../../common/utils.gvy'

workflow bench {
    take: datasets
    main:
        embeddings = snapatac2(datasets)
            | map({ [json(it[0]).data_name, it] })
            | combine(
                datasets | map({ [json(it[0]).data_name, it[1]] }),
                by: 0
            )
            | map ({ [it[1][0], it[1][1], it[2]] })
 
        clusters = knn_leiden_exact(embeddings) | concat(
            knn_leiden_hora(embeddings),
            knn_leiden_exact_weighted(embeddings),
            kmeans(embeddings),
        )

        clusters
            | map({ [json(it[0]).data_name, it] })
            | combine(
                datasets | map({ [json(it[0]).data_name, it[1]] }),
                by: 0
            )
            | map ({ [it[1][0], it[1][1], it[2]] })
            | bench_clustering
}