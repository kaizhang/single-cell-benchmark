nextflow.enable.dsl=2

include { download_hg38; download_hg19; } from '../dataset.nf'

include { dim_reduct_snapatac as snapatac; } from '../software/snapatac.nf'
include { dim_reduct_jaccard as snapatac2_jaccard;
          dim_reduct_cosine as snapatac2_cosine;
          dim_reduct_svd as snapatac2_svd;
        } from '../software/snapatac2.nf'
include { dim_reduct_signac_1 as signac_1;
          dim_reduct_signac_2 as signac_2;
          dim_reduct_signac_3 as signac_3;
          dim_reduct_signac_4 as signac_4;
        } from '../software/signac.nf'
include { dim_reduct_archr_1 as archr_1;
          dim_reduct_archr_2 as archr_2;
          dim_reduct_archr_3 as archr_3;
        } from '../software/archr.nf'
include { dim_reduct_peakvi as peakvi } from '../software/peakvi.nf'
include { dim_reduct_scale as scale } from '../software/scale.nf'
include { dim_reduct_pycistopic as cisTopic } from '../software/pycistopic.nf'
include { dim_reduct_scbasset as scBasset } from '../software/scbasset.nf'
include { dim_reduct_episcanpy as epiScanpy } from '../software/episcanpy.nf'

include { benchmark } from '../../../common/benchmark.nf'

workflow bench_dim_reduct {
    take: datasets
    main:
        datasets_ = datasets | map {[it[0], it[1]]}

        snapatac2_jaccard(datasets_) | concat(
            snapatac2_cosine(datasets_),
            snapatac2_svd(datasets_),

            snapatac(datasets_),

            epiScanpy(datasets_),

            signac_1(datasets_),
            signac_2(datasets_),
            signac_3(datasets_),
            signac_4(datasets_),

            archr_1(datasets_),
            archr_2(datasets_),
            archr_3(datasets_),

            cisTopic(datasets_),

            peakvi(datasets_),

            scale(datasets_),

            scBasset(datasets),
        )
            | combine(datasets_, by: 0)
            | benchmark
}
