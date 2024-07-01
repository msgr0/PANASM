#!/usr/bin/env nextflow

// DESCRIPTION:
// this workflows is intended to run 
// 1. PAN_ASSEMBLY on an intended threshold,
include { PANASSEMBLY } from "./workflows/panassembly.nf"
// 2. GENERATE A LIST of suitable samples (GT from NCBI)
include { BUILD_GT    } from "./workflows/groundtruth.nf"
// 3. GENERATE the preliminary prediction with MLPLASMIDS  
include { MLPLASMIDS } from "./workflows/mlplasmids.nf"
// 4. RUN the following tools
//     a. GPLAS
include { GPLASPAN; GPLASUNI; GPLASSKE  } from "./workflows/gplas.nf"
//     b. PBF
//     c. PBF pangenome (only on pangenome graphs)
include { PBFPANSTAR; PBFPAN; PBF as PBFUNI; PBF as PBFSKE } from "./workflows/pbf.nf"

include {PUBLISH } from "./workflows/utils.nf"

// 5. EVALUATE (GT, predition) for Unic, Skesa, Pangenome(U) and Pangenome(S)
// include { EVALUATE } from "./workflows.evaluate.nf"
// 6. PLOT (Unicyler, Pangenome(U)) PLOT (Skesa, Pangenome)
// do this last step by hand-ipynb-rnb-julianb: compute global statistics! compairison plots!


workflow {
    // take input data folder
    // read a pair of files, assembly s and assembly u

    thresholds = Channel.of('0', '150', '250', '500')
    input_ch = Channel.fromFilePairs(
        "${params.input}/*-S*-{s,u}.gfa.gz", type: 'file'
    )

     // s files will be placed first cause of colexigraphical order
    // sort them
    // fmeta.id = actual sample name
    // fmeta.species = actual sample species 4 character encoding
    input_ch = input_ch.map{meta, files -> 
        def fmeta = [:]
        fmeta.id = meta[5..-1]
        fmeta.species = meta[0..3]

        [fmeta, files]
    }
    // adding thresholds
    input_ch = input_ch.combine(thresholds).map {
        meta, files, thr ->
        meta.id = meta.id + "-${thr}"
        [meta + [thr:thr], files]
    }
    // extracting unicycler assembly
    uni_ch = input_ch.map {
        meta, files ->
        uni = files.findAll {it.toString().contains("-u.gfa.gz")}
        [meta, uni[0]]
    }
    //extracting skesa assembly
    ske_ch = input_ch.map {
        meta, files ->
        ske = files.findAll {it.toString().contains("-s.gfa.gz")}
        [meta, ske[0]]
    }
    // getting the actual input channel
    input_ch = uni_ch.join(ske_ch) | view

    input_ch | PANASSEMBLY

    pasm_ch = PANASSEMBLY.out.panassembly
    mixed_fasta_ch = PANASSEMBLY.out.mixed_fasta
    assembly_ch = PANASSEMBLY.out.assembly_graphs
    assembly_fa_ch = PANASSEMBLY.out.assembly_fasta

    
    // todo DE-Couple build-gt and NCBI, ncbi will be ran before buildGT
    // and download the assembly just once per sample instead of once
    // per SAMPLE-THRESHOLD
    build_ch = pasm_ch.join(mixed_fasta_ch.map{id, fa, fagz -> [id, fa]}).join(assembly_fa_ch) | view
    BUILD_GT(build_ch)

    // if no results form BUILDGT, it will hopefully exit here
    MLPLASMIDS(mixed_fasta_ch.map{id, fa, fagz -> [id, fagz]}.join(assembly_fa_ch)) // each sample at a certain threshold will have his prediction done here

    // DEcouple the transformation of the prediction here with a process?
    // mlplsmid --> convert --> gplas-pbf pred.

    (
        GPLASPAN(pasm_ch.join(MLPLASMIDS.out.mixed)) &// meta, pangenome, prediction
        GPLASUNI(assembly_ch.map{id, uni, ske -> [id, uni]}.join(MLPLASMIDS.out.uni)) &
        GPLASSKE(assembly_ch.map{id, uni, ske -> [id, ske]}.join(MLPLASMIDS.out.ske)) &

        PBFPAN(pasm_ch.join(MLPLASMIDS.out.mixed)) &
        PBFPANSTAR(pasm_ch.join(MLPLASMIDS.out.mixed)) &
        PBFUNI(assembly_ch.map{id, uni, ske -> [id, uni]}.join(MLPLASMIDS.out.uni), "u") &
        PBFSKE(assembly_ch.map{id, uni, ske -> [id, ske]}.join(MLPLASMIDS.out.ske), "s")
    ) | PUBLISH


    // (GPLASPAN.out.res &
    // GPLASUNI.out.res &
    // GPLASSKE.out.res &
    // PBFPAN.out.res &
    // PBFPANSTAR.out.res &
    // PBFUNI.out.res &
    // PBFSKE.out.res)




    // WORKFLOW

    // input_ch | PANASSEMBLY

    // build the pangenome with panassembly

    // run mlplasmids

}