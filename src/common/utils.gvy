import org.yaml.snakeyaml.Yaml

def is_included(method, include, exclude) {
    (include == null || include.contains(method)) &&
    (exclude == null || !exclude.contains(method))
}

def json_string(obj) {
    new groovy.json.JsonBuilder(obj).toString()
}

def json(text) {
    new groovy.json.JsonSlurper().parseText(text)
}

def add_meta(metadata, key, val) {
    def obj = json(metadata)
    obj[key] = val
    json_string(obj)
}

def parseYamlString(text) {
    def yaml = new Yaml()
    def obj = yaml.load(text)
    return obj
}

def genBenchId(metadata) {
    def obj = json(metadata)
    def method = obj['method']
    if ('hvg' in obj) {
        method = method + '/hvg_' + obj['hvg'].toString()
    } 
    if ('atac_feats' in obj) {
        method = method + '/atac_feats_' + obj['atac_feats'].toString()
    } 
    if ('scaled' in obj && obj['scaled']) {
        method = method + '/scaled'
    }
    add_meta(metadata, 'bench_id', method)
}