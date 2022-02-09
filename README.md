# scATAC-benchmark

NOTE: the benchmark dataset is not uploaded here.

Run the benchmarking pipeline using `nextflow run src/main.nf -resume`.

## How to prepare benchmark dataset

Benchmark data is stored in `anndata` format. The `anndata` object should contain:

1. `adata.X`: cell by feature count matrix. The format of feature names should be: `chr:start-end`.
2. `adata.obs["cell_annotation"]`: cluster/cell-type annotation of each barcode.