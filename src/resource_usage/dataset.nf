nextflow.enable.dsl=2

process import_dataset {
    cache false
    input:
      val(dir)
    output:
      val list

    exec:
    list = []
    file(dir).eachDir { item -> 
        list << [
            "name": item.getBaseName(),
            "fragments_name_sorted": item + "/fragments_nsrt.bed.gz",
            "fragments_sorted": item + "/fragments_srt.tsv.gz",
            "gene_annotations": item + "/gene_anno.gff3.gz",
        ]
    }
}
