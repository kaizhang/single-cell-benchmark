nextflow.enable.dsl=2

params.outdir = "results"

include { dim_reduct_cosine as snapatac2 } from '../software/snapatac2.nf'

include { bench_clustering as bench
        } from '../../common/benchmark.nf' params(resultDir: "${params.outdir}/ATAC/clustering")
include { knn_leiden_exact;
          knn_leiden_hora;
          knn_leiden_exact_weighted;
          kmeans;
        } from '../../common/algorithms/clustering.nf'


workflow bench_clustering {
    take: datasets
    main:
        data = datasets | map { [it[0], it[1]] }
        embeddings = snapatac2(data) | combine(data, by: 0)

        knn_leiden_exact(embeddings) | concat(
            knn_leiden_hora(embeddings),
            knn_leiden_exact_weighted(embeddings),
            kmeans(embeddings),
        )
            | combine(data, by: 0)
            | bench
}
