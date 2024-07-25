#!/usr/bin/env nextflow

process GPLASPAN {

    // errorStrategy 'ignore'
    maxForks 1
    input:
    tuple val(meta), path(gfa), path(pred)

    output:
    tuple val(meta), path(res), emit: res
    // tuple val(meta), path(bins), emit: bins

    script:
    slim_graph = "${meta.id}.slim.gfa"
    trimmed = "${meta.id}.slimmed.gfa"
    gplas_pred = "${meta.id}.pasm.gplas.ml.pred"

    bins    = "results/${meta.id}.pasm_bins.tab"
    results = "results/${meta.id}.pasm_results.tab"
    images = "results/${meta.id}.pasm_plasmidome_network.png"
    res = "${meta.id}.pasm.gplas.pred.tab"

    """
    #!/bin/bash

    python $projectDir/bin/easy-pangenome.py --input ${gfa} --output ${trimmed}

    python $projectDir/bin/mlpl.asmtopan.py --pred ${pred} --graph ${gfa}  --output ${gplas_pred} --thr 1000


    gplas -c predict -i ${trimmed} -P ${gplas_pred} -n ${meta.id}.pasm

    python $projectDir/bin/evaluation/transform_gplas_pred.py --input ${bins} --gfa ${gfa} --output ${res} 

    """
    // python $projectDir/bin/remove_nodes.py --input ${slim_graph} --output ${trimmed} --threshold 2500

}

process GPLASUNI {
    // errorStrategy 'ignore'
    maxForks 1


    input:
    tuple val(meta), path(gfa), path(pred)

    output:
    tuple val(meta), path(res), emit: res

    script:
    ren_gfa = "${meta.id}.u.ren.gfa"
    ren_fasta = "${meta.id}.u.ren.fasta"

    gplas_pred = "${meta.id}.u.gplas.ml.pred"

    bins    = "results/${meta.id}.u_bins.tab"
    results = "results/${meta.id}.u_results.tab"
    images = "results/${meta.id}.u_plasmidome_network.png"

    res = "${meta.id}.u.gplas.pred.tab"


    """
    #!/bin/bash
    python $projectDir/bin/rename_gfa.py --input ${gfa} --output ${ren_gfa} --prefix ""

    awk '/^S/{print ">"\$2; print ""\$3}' ${ren_gfa} > ${ren_fasta}

    python $projectDir/bin/mlpl.asm.py --pred ${pred} --graph ${ren_gfa}  --output ${gplas_pred} --prefix uni

  
    gplas -c predict -i ${ren_gfa} -P ${gplas_pred} -n ${meta.id}.u
    
    python $projectDir/bin/evaluation/transform_gplas_pred.py --input ${bins} --gfa ${ren_gfa} --output ${res} --prefix uni

    """

}

process GPLASSKE {

    // errorStrategy 'ignore'
    maxForks 1

    input:
    tuple val(meta), path(gfa), path(pred)

    output:
    tuple val(meta), path(res), emit: res

    script:
    ren_gfa = "${meta.id}.s.ren.gfa"
    ren_fasta = "${meta.id}.s.ren.fasta"

    gplas_pred = "${meta.id}.s.gplas.pred"

    bins    = "results/${meta.id}.s_bins.tab"
    results = "results/${meta.id}.s_results.tab"
    images = "results/${meta.id}s._plasmidome_network.png"
    res = "${meta.id}.s.gplas.pred.tab"


    // here the script converts SKESA into the naming convention gplas uses: i.e.
    // with numbered contigs 1.2.3.4.... and dp score instead of KC
    // rename gfa (CALL IT CONVERT GFA???)
    """
    #!/bin/bash
    python $projectDir/bin/rename_gfa.py --input ${gfa} --convert --output ${ren_gfa} --prefix ""
    awk '/^S/{print ">"\$2; print ""\$3}' ${ren_gfa} > ${ren_fasta}

    python $projectDir/bin/mlpl.asm.py --pred ${pred} --graph ${ren_gfa}  --output ${gplas_pred} --prefix ske

    gplas -c predict -i ${ren_gfa} -P ${gplas_pred} -n ${meta.id}.s

    python $projectDir/bin/evaluation/transform_gplas_pred.py --input ${bins} --gfa ${ren_gfa} --output ${res} --prefix ske
    """
}

