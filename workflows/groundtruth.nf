#!/usr/bin/env nextflow

include {PUBLISH } from "./utils.nf"

process NCBI {
    errorStrategy 'ignore'

    input:
    val (meta)

    output:
    tuple val(meta), path(reference), emit: ref

    script:
    referencegz = "${meta.id}.fna.gz"
    reference = "${meta.id}.fna"

    id = "${meta.id}".split("-")[0]
    """
    python $projectDir/bin/evaluation/ncbi_link.py --input ${id} --ouput ${referencegz}
    bgzip -d -c ${reference} > ${referencegz}
    """
}

process BLAST {

    input:
    tuple val(meta), val(graph), val(mix), val(uni), val(ske), val(reference)

    output:
    tuple val(meta), val(pan_mix_gt), val(mix_gt), emit: pangt
    tuple val(meta), val(pan_uni_gt), val(uni_gt), emit: unigt
    tuple val(meta), val(pan_ske_gt), val(ske_gt), emit: skegt


    script:
    outmix = meta.id + "." + "mix"
    // mapping_file = output + ".mapping.tsv"
    mix_gt = outmix + ".gt.tsv"
    pan_mix_gt = output + ".pan.gt.tsv"

    outuni = meta.id + "." + "uni"
    // mapping_file = output + ".mapping.tsv"
    uni_gt = outuni + ".gt.tsv"
    pan_uni_gt = output + ".pan.gt.tsv"

    outske = meta.id + "." + "ske"
    // mapping_file = output + ".mapping.tsv"
    ske_gt = outske + ".gt.tsv"
    pan_ske_gt = output + ".pan.gt.tsv"
    
    """
    python $projectDir/bin/evaluation/build_truth.py --input --pangenome ${graph} --assembly ${mix} --reference ${reference} --output ${outmix}
    
    python $projectDir/bin/evaluation/build_truth.py --input --pangenome ${graph} --assembly ${uni} --reference ${reference} --output ${outuni}

        
    python $projectDir/bin/evaluation/build_truth.py --input --pangenome ${graph} --assembly ${ske} --reference ${reference} --output ${outske}
    """

}


workflow BUILD_GT {
    take:

    input_ch // pangenome, mixed_fasta, assembly_fasta


    main:
    id = input_ch.map{id, pan, mixed, asm1, asm2 -> id}

    NCBI(id)

    reference = NCBI.out.ref

    BLAST ( input_ch.join(reference) )

    PUBLISH(BLAST.out.pangt.map{it, a, b -> [it, a]
        }.mix(
            BLAST.out.pangt.map{it, a, b -> [it, b]}
        ).mix(
            BLAST.out.unigt.map{it, a, b -> [it, a]}
        ).mix(
            BLAST.out.unigt.map{it, a, b -> [it, b]}
        ).mix(
            BLAST.out.skegt.map{it, a, b -> [it, a]}
        ).mix(
            BLAST.out.skegt.map{it, a, b -> [it, b]}
        )
    )
    emit:
    pan_mix = BLAST.out.pangt
    pan_uni = BLAST.out.unigt
    pan_ske = BLAST.out.skegt
}
