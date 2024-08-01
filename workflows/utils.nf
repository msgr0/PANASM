#!/usr/bin/env nextflow

process PUBLISH {
    publishDir "${outpath}", mode: 'copy', overwrite: true, pattern: "${item}"
    
    input:
    tuple val(meta), path(item)
    path(outpath)
    
    output:
    tuple val(meta), path(item)

    """
    echo
    """
}
