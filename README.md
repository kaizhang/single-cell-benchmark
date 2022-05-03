# scATAC-benchmark

NOTE: the benchmark dataset is not uploaded here.

## Dimenstion reduction

```
nextflow run src/accuracy/main.nf -resume
```

## Running time and memory usage

```
nextflow run src/resource_usage/main.nf -with-trace -qs 1 -resume
```

Visualize the resource usage:

```
nextflow log <run name> -f hash,name,status,realtime,peak_rss
```

## How to prepare benchmark dataset

Benchmark data is stored in `anndata` format. The `anndata` object should contain:

1. `adata.X`: cell by feature count matrix. The format of feature names should be: `chr:start-end`.
2. `adata.obs["cell_annotation"]`: cluster/cell-type annotation of each barcode.