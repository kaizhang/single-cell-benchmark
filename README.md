scATAC-benchmark
================

The list of benchmark datasets
------------------------------

| Name                 | #cells | #features | #cell types |
| -------------------- | ------ | --------- | ----------- |
| BoneMarrow_Chen_2019 | 9,600  | 156,311   | 6           |
| 10x_Brain5k          | 2,317  | 155,093   | 10          |
| 10x_PBMC_10k         | 9,631  | 107,194   | 19          |
| Buenrostro_2018      | 2,034  | 237,440   | 10          |
| Chen_NBT_2019        | 9,190  | 241,757   | 22          |
| GSE194122            | 69,249 | 116,490   | 22          |
| Ma_Cell_2020         | 32,231 | 340,341   | 22          |
| Trevino_Cell_2021    | 8,981  | 467,315   | 13          |
| Yao_Nature_2021      | 54,844 | 148,814   | 11          |
| Zhang_Cell_2021_GI   | 3,750  | 1,171,805 | 9           |

NOTE: the benchmark dataset is not uploaded here.

How to prepare benchmark dataset
--------------------------------

Benchmark data is stored in `anndata` format. The `anndata` object should contain:

1. `adata.X`: cell by feature count matrix. The format of feature names should be: `chr:start-end`.
2. `adata.obs["cell_annotation"]`: cluster/cell-type annotation of each barcode.

How to run benchmark
--------------------

### Dimenstion reduction

```
./bench_acc
```

### Running time and memory usage

```
./bench_perf
```

Visualize the resource usage:

```
nextflow log <run name> -f hash,name,status,realtime,peak_rss
```
