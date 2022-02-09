nextflow.enable.dsl=2

include { dim_reduct_snapatac_2;
          dim_reduct_snapatac_2_svd;
          dim_reduct_snapatac_2_cosine;
          dim_reduct_snapatac_2_nystrom;
          dim_reduct_snapatac_2_v2_nystrom_full;
          dim_reduct_snapatac_1;
          dim_reduct_snapatac_1_nystrom;
          dim_reduct_archr_1_log_tf_idf;
          dim_reduct_archr_1_logtf_logidf;
          dim_reduct_archr_1_tf_logidf;
          dim_reduct_archr_1_subsample;
        } from './software'
include { import_dataset } from './dataset'
include { benchmark_dim_reduct;
          generate_report;
          plot_report;
        } from './report'

workflow {
    datasets = import_dataset(Channel.fromList([
        tuple("Real", "dataset/Real"),
        tuple("Synthetic", "dataset/Synthetic"),
    ])).flatten()
    
    benchData = datasets.map { tuple(it, 50) }

    nystromBenchData = datasets.flatMap { data -> [
            [0.1, 0.2, 0.5],
            [1, 2, 3, 4, 5],
        ].combinations().collect { x -> 
            d = data.clone()
            d["samplingFraction"] = x[0]
            d["randomSeed"] = x[1]
            tuple(d, 50)
        }
    }

    runResult = dim_reduct_snapatac_2(benchData).concat(
        dim_reduct_archr_1_log_tf_idf(benchData),
        dim_reduct_archr_1_logtf_logidf(benchData),
        dim_reduct_archr_1_tf_logidf(benchData),

        dim_reduct_snapatac_1(benchData),
        dim_reduct_snapatac_2_cosine(benchData),
        dim_reduct_snapatac_2_svd(benchData),

        dim_reduct_snapatac_1_nystrom(nystromBenchData),
        dim_reduct_snapatac_2_nystrom(nystromBenchData),
        dim_reduct_archr_1_subsample(nystromBenchData),
    )

    benchResult = benchmark_dim_reduct(runResult)
        .map { method, data, result, plot -> [
            "datasetType": data.dataType,
            "datasetName": data.name,
            "datasetMatrix": data.anndata,
            "methodName": method,
            "samplingFraction": data.get("samplingFraction", 1),
            "randomSeed": data.randomSeed,
            "resultOutput": result,
            "umap": plot,
        ]}.collect()

    plot_report(generate_report(benchResult))
}