nextflow.enable.dsl=2

import groovy.transform.Immutable
import nextflow.io.ValueObject
import nextflow.util.KryoHelper

@ValueObject
@Immutable(copyWith=true, knownImmutables = ['name', 'anndata', 'matrixMarket'])
class BenchDataset {
    static { 
        // Register this class with the Kryo framework that serializes and deserializes objects
        // that pass through channles. This allows for caching when this object is used.
        KryoHelper.register(BenchDataset)
    }
    String name
    String anndata
    String matrixMarket
}

process import_dataset {
    input:
      val dir
    output:
      val list

    exec:
    list = []
    file(dir).eachDir { item -> 
        list << new BenchDataset(
            name: item.getBaseName(),
            anndata: item + "/matrix.h5ad",
            matrixMarket: item + "/matrix.mtx.gz"
        )
    }
}