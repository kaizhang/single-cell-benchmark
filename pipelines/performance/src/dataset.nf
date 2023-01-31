nextflow.enable.dsl=2

process import_dataset {
    cache false
    input:
      val(dir)
    output:
      val list

    exec:
    list = []
    file("../../" + dir).eachDir { item -> 
        list << tuple(
          item.getBaseName(),
          item + "/fragments_nsrt.bed.gz",
          item + "/fragments_srt.tsv.gz",
          item + "/gene_anno.gff3.gz",
        )
    }
}
