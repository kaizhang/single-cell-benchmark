nextflow.enable.dsl=2

process import_dataset {
    cache false
    input:
      tuple val(type), val(dir)
    output:
      val list

    exec:
    list = []
    file(dir).eachDir { item -> 
        list << [
            "dataType": type,
            "name": item.getBaseName(),
            "anndata": item + "/matrix.h5ad",
        ]
    }
}

process import_raw_dataset {
    cache false
    input:
      tuple val(type), val(dir)
    output:
      val list

    exec:
    list = []
    file(dir).eachDir { item -> 
        list << [
            "dataType": type,
            "name": item.getBaseName(),
            "fragments_name_sorted": item + "/fragments_nsrt.tsv.gz",
            "fragments_sorted": item + "/fragments_srt.tsv.gz",
            "gene_annotations": item + "/gene_anno.gff3.gz",
            "metadata": item + "/metadata.tsv.gz",
        ]
    }
}
