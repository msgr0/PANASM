#!/usr/bin/env nextflow

process GPLASPAN {

    cache 'lenient'
    errorStrategy 'ignore'
    maxForks 1
    input:
    tuple val(meta), path(gfa), path(pred)

    output:
    tuple val(meta), path(res), emit: res
    // tuple val(meta), path(bins), emit: bins

    script:
    slim_graph = "${meta.id}.slim.gfa"
    species = convert_species(meta.species)
    gplas_pred = "${meta.id}.pasm.gplas.pred"
    pbf_pred = "${meta.id}.pasm.pbf.pred"

    bins    = "results/${meta.id}.pasm_bins.tab"
    results = "results/${meta.id}.pasm_results.tab"
    images = "results/${meta.id}.pasm_plasmidome_network.png"
    res = "${meta.id}.pasm.gplas.pred.tab"

    """
    #!/bin/bash

    python $projectDir/bin/easy-pangenome.py --input ${gfa} --output ${slim_graph}

    python $projectDir/bin/mlpl.asmtopan.py --pred ${pred} --graph ${gfa}  --output ${gplas_pred} --pbf ${pbf_pred}

    gplas -c predict -i ${slim_graph} -P ${gplas_pred} -n ${meta.id}.pasm

    python $projectDir/bin/evaluation/transform_gplas_pred.py --input ${bins} --gfa ${gfa} --output ${res} 

    """

}

process GPLASUNI {
    cache 'lenient'
    errorStrategy 'ignore'
    maxForks 1


    input:
    tuple val(meta), path(gfa), path(pred)

    output:
    tuple val(meta), path(bins), emit: bins

    script:
    ren_gfa = "${meta.name}.u.ren.gfa"
    ren_fasta = "${meta.name}.u.ren.fasta"
    species = convert_species(meta.species)

    gplas_pred = "${meta.name}.u.gplas.pred"
    pbf_pred = "${meta.name}.u.pbf.pred"

    bins    = "results/${meta.name}.u_bins.tab"
    results = "results/${meta.name}.u_results.tab"
    images = "results/${meta.name}.u_plasmidome_network.png"

    res = "${meta.id}.u.gplas.pred.tab"


    """
    #!/bin/bash
    python $projectDir/bin/rename_gfa.py --input ${gfa} --output ${ren_gfa} --prefix ""

    awk '/^S/{print ">"\$2; print ""\$3}' ${ren_gfa} > ${ren_fasta}

    python $projectDir/bin/mlpl.asm.py --pred ${pred} --graph ${ren_gfa}  --output ${gplas_pred} --pbf ${pbf_pred}

  
    gplas -c predict -i ${ren_gfa} -P ${gplas_pred} -n ${meta.name}.u
    
    python $projectDir/bin/evaluation/transform_gplas_pred.py --input ${bins} --gfa ${gfa} --output ${res} 

    """

}

process GPLASSKE {
    cache 'lenient'
    errorStrategy 'ignore'
    maxForks 1

    input:
    tuple val(meta), path(gfa), path(pred)

    output:
    tuple val(meta), path(bins), emit: bins

    script:
    ren_gfa = "${meta.id}.s.ren.gfa"
    ren_fasta = "${meta.id}.s.ren.fasta"

    species = convert_species(meta.species)
    gplas_pred = "${meta.id}.s.gplas.pred"
    pbf_pred = "${meta.id}.s.pbf.pred"

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

    python $projectDir/bin/mlpl.asm.py --pred ${pred} --graph ${ren_gfa}  --output ${gplas_pred} --pbf ${pbf_pred}

    gplas -c predict -i ${ren_gfa} -P ${gplas_pred} -n ${meta.id}.s

    python $projectDir/bin/evaluation/transform_gplas_pred.py --input ${bins} --gfa ${gfa} --output ${res} 
    """
}

