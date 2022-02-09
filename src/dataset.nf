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
            "matrixMarket": item + "/matrix.mtx.gz"
        ]
    }
}