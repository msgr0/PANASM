#!/usr/bin/env nextflow

process PUBLISH {
    publishDir "${params.output}", mode: 'copy', overwrite: true, pattern: "${item}"
    
    input:
    tuple val(meta), path(item)
    
    output:
    tuple val(meta), path(item)
    
    script:
    """
    """
}
