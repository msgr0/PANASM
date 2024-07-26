#!/usr/bin/env nextflow
paramfile = "$projectDir/workflows/pangenome-params.json"
release = '1.1.0'
include {PUBLISH } from "./utils"

process PREPROCESS {
    input:
    tuple val(meta), path(uni_gfa), path(skesa_gfa)

    output:
    tuple val(meta), path(uni_gfa_o), path(ske_gfa_o), emit: assembly_graphs
    tuple val(meta), path(uni_fasta), path(ske_fasta), emit: assembly_fasta
    tuple val(meta), path(mixed_fa), path(mixed_fagz), emit: mixed_fasta

    script:

    uni_orig = "${meta.id}.orig.u.gfa"
    uni_trim = "${meta.id}.trim.u.gfa"
    uni_gfa_o = "${meta.id}.u.gfa"
    uni_fasta = "${meta.id}.u.fasta"

    ske_orig = "${meta.id}.orig.s.gfa"
    ske_trim = "${meta.id}.trim.s.gfa"
    ske_gfa_o = "${meta.id}.s.gfa"
    ske_fasta = "${meta.id}.s.fasta"

    mixed_fa = "${meta.id}.mixed.fasta"
    mixed_fagz = "${meta.id}.mixed.fasta.gz"
    """
    #!/usr/bin/env bash

    bgzip -d -c ${uni_gfa} > ${uni_orig}
    bgzip -d -c ${skesa_gfa} > ${ske_orig}

    if [ ${meta.thr} == 0]; then
        cp ${uni_orig} ${uni_trim}
        cp ${ske_orig} ${ske_trim}
    else
        # python $projectDir/bin/remove_nodes.py -i ${uni_orig} -o ${uni_trim} -t ${meta.thr}
        # python $projectDir/bin/remove_nodes.py -i ${ske_orig} -o ${ske_trim} -t ${meta.thr}
        python $projectDir/bin/cedric/filter_GFA.py ${uni_orig} ${meta.thr} ${uni_trim}
        python $projectDir/bin/cedric/filter_GFA.py ${ske_orig} ${meta.thr} ${ske_trim}
    
    fi

    python $projectDir/bin/rename_gfa.py -i ${uni_trim} -o ${uni_gfa_o} -p 'uni'
    python $projectDir/bin/rename_gfa.py -i ${ske_trim} -o ${ske_gfa_o} -p 'ske'

    awk '/^S/{print ">"\$2; print ""\$3}' ${uni_gfa_o} > ${uni_fasta}
    awk '/^S/{print ">"\$2; print ""\$3}' ${ske_gfa_o} > ${ske_fasta}

    cat ${uni_fasta} ${ske_fasta} > ${mixed_fa}
    bgzip -k -c ${mixed_fa} > ${mixed_fagz}
    """
}

process MAKEPANGENOME {
    maxForks 4
    time '30m'

    input:
    tuple val(meta), path(mixfagz)

    output:
    tuple val(meta), path(pangenome), emit: pangenome

    script:
    haplos = 2
    pangenome = "${meta.id}.pan.gfa"
    """
    nextflow run nf-core/pangenome -r ${release} -profile apptainer -resume --input ${mixfagz} --n_haplotypes ${haplos} --outdir . -params-file ${paramfile} -w ${task.workDir} 
    cp ./FINAL_GFA/${mixfagz}.gfaffix.unchop.Ygs.view.gfa ${pangenome}
    """

}

process AUGMENT {
    
    input:
    tuple val(meta), path(pangenome), path(uni), path(ske)

    output:
    tuple val(meta), path(panassembly), emit: panassembly


    script:
    panassembly = "${meta.id}.pasm.gfa"
    pangenomecl = "${meta.id}.pan.cl.gfa"
    // trimmed = ${meta.id}.pasm.trim.gfa
    // trimmed_fasta = ${meta.id}.pasm.trim.fasta
    // TODO merge the python scripts?? check if reused somewhere
    """
    python $projectDir/bin/gfa_cleaner.py --input ${pangenome} --output ${pangenomecl}
    python $projectDir/bin/pangenome_enancher.py --input ${pangenomecl} --output ${panassembly} --assembly ${uni} ${ske}
    """
}

workflow PANASSEMBLY {
    take:
    input_ch

    main:

    PREPROCESS(input_ch)
    MAKEPANGENOME(PREPROCESS.out.mixed_fasta.map{id, fa, fagz -> [id, fagz]})
    AUGMENT(MAKEPANGENOME.out.pangenome.join(PREPROCESS.out.assembly_graphs))

    PUBLISH(AUGMENT.out.panassembly)

    emit:
    panassembly = AUGMENT.out.panassembly
    mixed_fasta = PREPROCESS.out.mixed_fasta
    assembly_graphs = PREPROCESS.out.assembly_graphs
    assembly_fasta = PREPROCESS.out.assembly_fasta
    
    
}


// workflow {

//     input_ch = Channel.fromFilePairs(
//         "${params.input}/*-S*-{s,u}.gfa.gz", type: 'file'
//     )
//     input_ch = input_ch.map{meta, files -> 
//         def fmeta = [:]
//         fmeta.id = meta[5..-1]
//         fmeta.species = meta[0..3]

//         [fmeta, files]
//     }


//     thresholds = Channel.of('0', '500')
//     input_ch = input_ch.combine(thresholds).map {
//         meta, files, thr ->
//         [meta + [thr:thr, name:"${meta.species}-${meta.id}-${thr}"], files]
//     }
//     input_ch | view

//     PANASSEMBLY(input_ch)
//     PUBLISH(PANASSEMBLY.out.panassembly)

// }