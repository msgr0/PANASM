#!/usr/bin/env nextflow

process PUBLISH {
    publishDir "${params.output}", mode: 'copy', overwrite: true, pattern: "${name}"
    
    input:
    tuple val(meta), path(item)
    
    output:
    tuple val(meta), path(item), val(name)
    
    script:
    name = item.getName()
    """
    """
}
