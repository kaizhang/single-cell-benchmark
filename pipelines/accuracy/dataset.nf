nextflow.enable.dsl=2

process import_dataset {
    cache false
    input:
      tuple val(type), val(dir)
    output:
      val list

    exec:
    list = []
    file("../../" + dir).eachDir { item -> 
        list << tuple(type, item.getBaseName(), item + "/matrix.h5ad")
    }
}
