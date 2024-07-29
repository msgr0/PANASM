#!/usr/bin/env nextflow

// include {PUBLISH; PUBLISH as REF} from "./utils.nf"

process NCBI {
    publishDir "${params.input}/", mode: 'copy', overwrite: true
    errorStrategy 'ignore'

    input:
    val (meta)

    output:
    tuple val(meta), path(reference_ren), emit: ref

    script:
    referencegz = "${meta.id}.fna.gz"
    reference = "${meta.id}.fna"
    reference_ren = "${meta.id}.ren.fna"

    id = "${meta.id}".split("-")[0]
    """
    python $projectDir/bin/evaluation/ncbi_link.py --input ${id} --output ${referencegz}
    bgzip -d -c ${referencegz} > ${reference}

    python $projectDir/bin/evaluation/strip_plasmid_fasta.py --input ${reference} --output ${reference_ren}
    """
}

process BLAST {

    input:
    tuple val(meta), path(graph), path(mix), path(uni), path(ske), path(reference)

    output:
    tuple val(meta), path(pan_mix_gt), path(mix_gt), emit: pangt
    tuple val(meta), path(pan_uni_gt), path(uni_gt), emit: unigt
    tuple val(meta), path(pan_ske_gt), path(ske_gt), emit: skegt


    script:
    outmix = meta.id + "." + "mix"
    // mapping_file = output + ".mapping.tsv"
    mix_gt = outmix + ".gt.tsv"
    pan_mix_gt = outmix + ".pan.gt.tsv"

    outuni = meta.id + "." + "uni"
    // mapping_file = output + ".mapping.tsv"
    uni_gt = outuni + ".gt.tsv"
    pan_uni_gt = outuni + ".pan.gt.tsv"

    outske = meta.id + "." + "ske"
    // mapping_file = output + ".mapping.tsv"
    ske_gt = outske + ".gt.tsv"
    pan_ske_gt = outske + ".pan.gt.tsv"
    
    """
    python $projectDir/bin/evaluation/build_truth.py --pangenome ${graph} --assembly ${mix} --reference ${reference} --output ${outmix}
    
    python $projectDir/bin/evaluation/build_truth.py --pangenome ${graph} --assembly ${uni} --reference ${reference} --output ${outuni}

    python $projectDir/bin/evaluation/build_truth.py --pangenome ${graph} --assembly ${ske} --reference ${reference} --output ${outske}
    """

}

workflow DOWNLOAD_GT {
    take:
    id

    main:
    NCBI(id)

    emit:
    reference =  NCBI.out.ref

}

workflow BUILD_GT {
    take:

    input_ch // pangenome, mixed_fasta, assembly_fasta
    reference


    main:
    
    BLAST ( input_ch.join(reference) )

    // PUBLISH(BLAST.out.pangt.map{it, a, b -> [it, a]
    //     }.mix(
    //         BLAST.out.pangt.map{it, a, b -> [it, b]}
    //     ).mix(
    //         BLAST.out.unigt.map{it, a, b -> [it, a]}
    //     ).mix(
    //         BLAST.out.unigt.map{it, a, b -> [it, b]}
    //     ).mix(
    //         BLAST.out.skegt.map{it, a, b -> [it, a]}
    //     ).mix(
    //         BLAST.out.skegt.map{it, a, b -> [it, b]}
    //     )
    // )
    emit:
    pan_mix = BLAST.out.pangt
    pan_uni = BLAST.out.unigt
    pan_ske = BLAST.out.skegt
}
