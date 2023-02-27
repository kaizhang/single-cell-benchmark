scATAC-benchmark
================

The list of benchmark datasets
------------------------------

https://1drv.ms/x/s!AuDfAP0z9Zsma9Hsr9Ggoz3MtFU?e=0RKKLU

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
