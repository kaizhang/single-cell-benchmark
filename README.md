# Comprehensive benchmark of dimension reduction methods for single-cell data analysis

Software requirements

- [Nextflow](https://www.nextflow.io/)
- Docker or Singularity

The list of scATAC-seq datasets
-------------------------------

https://1drv.ms/x/s!AuDfAP0z9Zsma9Hsr9Ggoz3MtFU?e=0RKKLU

How to prepare benchmark dataset
--------------------------------

Benchmark data is stored in `anndata` format. The `anndata` object should contain:

1. `adata.X`: cell by feature count matrix. The format of feature names should be: `chr:start-end`.
2. `adata.obs["cell_annotation"]`: cluster/cell-type annotation of each barcode.
3. `data.obs['batch']`: batch annotation of each barcode (optional).

The list of scRNA-seq datasets
-------------------------------

The list of scHi-C datasets
-------------------------------

How to run benchmark
--------------------

Use `./bench.sh` within each subdirectory to run benchmark. The benchmark results will be stored in `./results` folder.

You can choose to run specific benchmark by specifying the benchmark name. For example, `./bench.sh -entry ATAC` will run ATAC benchmark only.