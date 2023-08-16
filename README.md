# Comprehensive benchmark of dimension reduction methods for single-cell data analysis

## Software requirements

- [Nextflow](https://www.nextflow.io/)
- Docker or Singularity

## The list of benchmarking datasets

- scATAC-seq: https://1drv.ms/x/s!AuDfAP0z9Zsma9Hsr9Ggoz3MtFU?e=0RKKLU
- scRNA-seq:
- scHi-C:

## How to prepare new benchmark datasets

Benchmark data is stored in `anndata` format. The `anndata` object should contain:

1. `adata.X`: cell by feature count matrix. For scATAC-seq, the format of feature names should be: `chr:start-end`.
2. `adata.obs["cell_annotation"]`: cluster/cell-type annotation of each barcode.
3. `data.obs['batch']`: batch annotation of each barcode (optional).

## How to run benchmark

The pipeline uses docker or singularity containers to run the benchmark. Therefore, you need to install either docker or singularity on your machine.

Use `./bench.sh -profile singularity` or `./bench.sh -profile docker` to run benchmarks.
The benchmark results will be stored in `./results` folder.

See `nextflow.config` for additional configuration options.