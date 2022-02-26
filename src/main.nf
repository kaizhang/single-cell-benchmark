nextflow.enable.dsl=2

include { dim_reduct_snapatac_2;
          dim_reduct_snapatac_2_svd;
          dim_reduct_snapatac_2_cosine;
          dim_reduct_snapatac_2_nystrom;

          dim_reduct_snapatac_1;
          dim_reduct_snapatac_1_nystrom;

          dim_reduct_archr_1;
          dim_reduct_archr_2;
          dim_reduct_archr_3;
          dim_reduct_archr_subsample;

          dim_reduct_signac_1;
          dim_reduct_signac_2;
          dim_reduct_signac_3;
          dim_reduct_signac_4;

          dim_reduct_cistopic;
        } from './software'
include { import_dataset } from './dataset'
include { benchmark_dim_reduct;
          benchmark_clustering;
        } from './benchmark'
include { gen_dim_reduct_report
          gen_clust_report;
          plot_dim_reduct_report;
          plot_clust_report;
          plot_umap;
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
        dim_reduct_archr_1(benchData),
        //dim_reduct_archr_2(benchData),
        //dim_reduct_archr_3(benchData),

        dim_reduct_signac_1(benchData),
        dim_reduct_signac_2(benchData),
        dim_reduct_signac_3(benchData),
        dim_reduct_signac_4(benchData),

        dim_reduct_snapatac_1(benchData),

        dim_reduct_snapatac_2_cosine(benchData),
        dim_reduct_snapatac_2_svd(benchData),

        dim_reduct_cistopic(benchData),

        dim_reduct_snapatac_1_nystrom(nystromBenchData),
        dim_reduct_snapatac_2_nystrom(nystromBenchData),
        //dim_reduct_archr_subsample(nystromBenchData),
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
    dim_reduct_report = gen_dim_reduct_report(benchResult)
    plot_dim_reduct_report(dim_reduct_report)
    plot_umap(dim_reduct_report)

    clustering = benchmark_clustering(
        runResult.filter { name, data, result -> name == "SnapATAC2" && data.samplingFraction == null }
    ).map { method, data, result -> [
        "datasetType": data.dataType,
        "datasetName": data.name,
        "resultOutput": result,
    ]}.collect()
    plot_clust_report(gen_clust_report(clustering))
}
