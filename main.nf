#!/usr/bin/env nextflow

// DESCRIPTION:
// this workflows is intended to run 
// 1. PAN_ASSEMBLY on an intended threshold,
include { PANASSEMBLY } from "./workflows/panassembly.nf"
// 2. GENERATE A LIST of suitable samples (GT from NCBI)
include { DOWNLOAD_GT } from "./workflows/groundtruth.nf"
include { BUILD_GT    } from "./workflows/groundtruth.nf"
// 3. GENERATE the preliminary prediction with MLPLASMIDS  

include { MLPLASMIDS } from "./workflows/mlplasmids.nf"
// 4. RUN the following tools
//     a. GPLAS
include { GPLASPAN; GPLASUNI; GPLASSKE  } from "./workflows/gplas.nf"
//     b. PBF
//     c. PBF pangenome (only on pangenome graphs)
include { PBFPANSTAR; PBFPAN; PBF as PBFUNI; PBF as PBFSKE } from "./workflows/pbf.nf"

include { EVALUATE as EVAL } from "./workflows/evaluation.nf"


workflow {

    thresholds = Channel.of(params.thr)
    input_ch = Channel.fromFilePairs("${params.input}/*-S*-{s,u}.gfa.gz", type: 'file', checkIfExists: true)

    input_ch = input_ch.map{ meta, filepair -> def fmeta = [:]; fmeta.id = meta[5..-1]; fmeta.species = meta[0..3]
        [fmeta, filepair]
    }
    meta_ch = input_ch.map{meta, files -> meta}

    if (params.local) {
        reference_ch = meta_ch.combine(Channel.fromPath("${params.input}/*.ren.fna", type: 'file', checkIfExists: true))
    }
    else {
        DOWNLOAD_GT(meta_ch)
        reference_ch = DOWNLOAD_GT.out.reference
    }
    reference_ch | view

    input_ch = input_ch.combine(thresholds) | map{ meta, files, thr -> [meta +=[thr:thr], files]} | view
    unicycler_ch = input_ch | map { meta, files -> uni = files.findAll { it.toString().contains("-u.gfa.gz")}; [meta, uni[0]]}
    skesa_ch = input_ch | map { meta, files -> ske = files.findAll { it.toString().contains("-s.gfa.gz")}; [meta, ske[0]]}
    input_ch = unicycler_ch.join(skesa_ch) | view

    output_ch = Channel.empty()

    /*
    Main pipeline based on choices
    */
    reference_ch | view


    input_ch | PANASSEMBLY
    panassembly_ch = PANASSEMBLY.out.panassembly

    mixed_fasta_ch = PANASSEMBLY.out.mixed_fasta
    assembly_ch = PANASSEMBLY.out.assembly_graphs
    assembly_fa_ch = PANASSEMBLY.out.assembly_fasta

    ground_truth_ch = BUILD_GT(panassembly_ch.join(mixed_fasta_ch.map{id, fa, fagz -> [id, fa]}).join(assembly_fa_ch), reference_ch)

    MLPLASMIDS(mixed_fasta_ch.map{id, fa, fagz -> [id, fa]}.join(assembly_fa_ch)) // each sample at a certain threshold will have his prediction done here
        // MLPLASMIDS.out.mixed | PUBLISH

    GPLASPAN(pasm_ch.join(MLPLASMIDS.out.mixed))
    PBFPAN(pasm_ch.join(MLPLASMIDS.out.mixed)) 
    PBFPANSTAR(pasm_ch.join(MLPLASMIDS.out.mixed)) 

    gplas_panu_ch = [GPLASPAN.out.res.join(BUILD_GT.out.pan_uni.map{id, pan, uni -> [id, pan]}), "gplas.pan.uni"]
    gplas_pans_ch = [GPLASPAN.out.res.join(BUILD_GT.out.pan_ske.map{id, pan, ske -> [id, pan]}), "gplas.pan.ske"]
    
    pbf_panu_ch = [PBFPAN.out.res.join(BUILD_GT.out.pan_uni.map{id, pan, uni -> [id, pan]}), "pbf.pan.uni"]
    pbf_pans_ch = [PBFPAN.out.res.join(BUILD_GT.out.pan_ske.map{id, pan, ske -> [id, pan]}), "pbf.pan.ske"]
        
    pbfs_panu_ch = [PBFPANSTAR.out.res.join(BUILD_GT.out.pan_uni.map{id, pan, uni -> [id, pan]}), "pbfstar.pan.uni"]
    pbfs_pans_ch = [PBFPANSTAR.out.res.join(BUILD_GT.out.pan_ske.map{id, pan, ske -> [id, pan]}), "pbfstar.pan.ske"]

    evaluate_ch = gplas_panu_ch.mix(gplas_pans_ch).mix(pbf_panu_ch).mix(pbf_pans_ch).mix(pbfs_panu_ch).mix(pbfs_pans_ch)
    EVAL(evaluate_ch)
        


        
}


// // //     if (params.inpt.pangenome) {
// // //         input_ch | PANASSEMBLY
// // //         pasm_ch = PANASSEMBLY.out.panassembly
// // //         mixed_fasta_ch = PANASSEMBLY.out.mixed_fasta
// // //         assembly_ch = PANASSEMBLY.out.assembly_graphs
// // //         assembly_fa_ch = PANASSEMBLY.out.assembly_fasta

// // //         build_ch = 

// // //         BUILD_GT(build_ch, reference_ch)
// // //         MLPLASMIDS(mixed_fasta_ch.map{id, fa, fagz -> [id, fa]}.join(assembly_fa_ch)) // each sample at a certain threshold will have his prediction done here
// // //         MLPLASMIDS.out.mixed | PUBLISH
// // //     }

// // //     if (params.tool.gplas) {
// // //         GPLASPAN(pasm_ch.join(MLPLASMIDS.out.mixed))
// // //         // out_ch = params.inpt.assembly ? GPLASUNI(assembly_ch.map{id, uni, ske -> [id, uni]}.join(MLPLASMIDS.out.uni)) : Channel.empty()
// // //         // GPLASUNI(assembly_ch.map{id, uni, ske -> [id, uni]}.join(MLPLASMIDS.out.uni)) 
// // //         // GPLASSKE(assembly_ch.map{id, uni, ske -> [id, ske]}.join(MLPLASMIDS.out.ske)) 

// // //         // gppanuni_ch = 
// // //         EVAL1(GPLASPAN.out.res.join(BUILD_GT.out.pan_uni.map{id, pan, uni -> [id, pan]}),
// // //         "gplas.pan.uni")
        
// // //         // gpuni_ch =
// // //         // EVAL2(GPLASUNI.out.res.join(BUILD_GT.out.pan_uni.map{id, pan, uni -> [id, uni]}),
// // //         // "gplas.uni")

// // //         ///////////////////// ske gplas
        
// // //         // gppanske_ch = 
// // //         EVAL3(GPLASPAN.out.res.join(BUILD_GT.out.pan_ske.map{id, pan, ske -> [id, pan]}),
// // //         "gplas.pan.ske")
        
// // //         // gpske_ch =
// // //         // EVAL4(GPLASSKE.out.res.join(BUILD_GT.out.pan_ske.map{id, pan, ske -> [id, ske]}),
// // //         // "gplas.ske")

// // //     }

// // //     if (params.tool.pbf) {
// // //         // PBFPAN(pasm_ch.join(MLPLASMIDS.out.mixed)) 
// // //         // PBFPANSTAR(pasm_ch.join(MLPLASMIDS.out.mixed)) 
// // //         // PBFUNI(assembly_ch.map{id, uni, ske -> [id, uni]}.join(MLPLASMIDS.out.uni), "u") 
// // //         // PBFSKE(assembly_ch.map{id, uni, ske -> [id, ske]}.join(MLPLASMIDS.out.ske), "s")

// // //         // // pbfpanuni_ch =
// // //         // EVAL5(PBFPAN.out.res.join(BUILD_GT.out.pan_uni.map{id, pan, uni -> [id, pan]}),
// // //         // "pbf.pan.uni")
        
// // //         // // pbfstarpanuni_ch =
// // //         // EVAL6(PBFPANSTAR.out.res.join(BUILD_GT.out.pan_uni.map{id, pan, uni -> [id, pan]}),
// // //         // "pbf.panstar.uni")
        
// // //         // // pbfuni_ch =
// // //         // EVAL7(PBFUNI.out.res.join(BUILD_GT.out.pan_uni.map{id, pan, uni -> [id, uni]}),
// // //         // "pbf.uni")

// // //         // ///////////////////// ske pbf
// // //         // // pbfpanske_ch =
// // //         // EVAL8(PBFPAN.out.res.join(BUILD_GT.out.pan_ske.map{id, pan, ske -> [id, pan]}),
// // //         // "pbf.pan.ske")

// // //         // // pbfstarpanske_ch = 
// // //         // EVAL9(PBFPANSTAR.out.res.join(BUILD_GT.out.pan_ske.map{id, pan, ske -> [id, pan]}),
// // //         // "pbf.panstar.ske")

// // //         // // pbfske_ch = 
// // //         // EVAL10(PBFSKE.out.res.join(BUILD_GT.out.pan_ske.map{id, pan, ske -> [id, ske]}),
// // //         // "pbf.ske")

// // //     }



// // //     ///////////////////// uni pbf


// // //     // WORKFLOW

// // //     // input_ch | PANASSEMBLY

// // //     // build the pangenome with panassembly

// // //     // run mlplasmids


// // //     // GPLASPAN.out.res | PUBLISH //.mix(
// // //         // GPLASUNI.out.res.mix(
// // //         //     GPLASSKE.out.res)
// // //     // ) | PUBLISH

// // //     //         .mix(
// // //     //             PBFPAN.out.res.mix(
// // //     //                 PBFPANSTAR.out.res.mix(
// // //     //                     PBFUNI.out.res.mix(
// // //     //                         PBFSKE.out.res
// // //     //                     )
// // //     //                 )
// // //     //             )
// // //     //         )
// // //     //     )
// // //     // )
// // // }

