#!/usr/bin/env nextflow
binfolder = "~/bin"
process PBFPAN {
    // NOT include penalty term, as the assembly
    input:
    tuple val(meta), path(graph), path(pred)

    output:

    tuple val(meta), path(res), emit: res

    script:

    putils = "${binfolder}/pbfpan/code/plasbin_utils.py"
    pflow = "${binfolder}/pbfpan/code/plasbin_flow.py"
    seed_len = 1000
    seed_score = 0.5
    min_pls_len = 1000
    asm = "pasm"

    input_csv = "${meta.id}.${asm}-input.csv"
    gc_content = "${meta.id}.${asm}.gc.tsv"
    bins = "${meta.id}.${asm}.pbf.bins.tsv"
    res = "${meta.id}.${asm}.pbf.pred.tab"
    asm_pred = "null.pred"
    pbf_pred = "${meta.id}.${asm}.pbf.pred"

    """
    python $projectDir/bin/mlpl.asmtopan.py --pred ${pred} --graph ${graph}  --output ${asm_pred} --pbf ${pbf_pred}

    bgzip -k ${graph}
    echo "sample,gfa,pls_score" > "${input_csv}"
    echo "${meta.id}.${asm},${graph}.gz,${pbf_pred}" >> "${input_csv}"
    
    python ${putils} gc_probabilities --input_file ${input_csv} --out_dir . --tmp_dir ./tmp

    python ${pflow} -alpha4 1 -ag ${graph}.gz -gc ${gc_content} -out_dir . -out_file ${bins}  -score ${pbf_pred} -assembler pangenome -seed_len ${seed_len}  -seed_score ${seed_score} -min_pls_len ${min_pls_len}


    python $projectDir/bin/evaluation/transform_pbf_pred.py --input ${bins} --gfa ${graph} --output ${res} 
    """

}

process PBFPANSTAR {
    // include the penalty term for pangenomes

    input:
    tuple val(meta), path(graph), path(pred)

    output:

    tuple val(meta), path(res), emit: res

    script:
    seed_len = 1000
    seed_score = 0.5
    min_pls_len = 1000
    putils = "${binfolder}/pbfpan/code/plasbin_utils.py"
    pflow = "${binfolder}/pbfpan/code/plasbin_flow.py"

    asm = "pstar"
    input_csv = "${meta.id}.${asm}-input.csv"
    gc_content = "${meta.id}.${asm}.gc.tsv"
    bins = "${meta.id}.${asm}.pbf.bins.tsv"
    res = "${meta.id}.${asm}.pbf.pred.tab"
    asm_pred = "null.pred"
    pbf_pred = "${meta.id}.${asm}.pbf.pred"


    """
    python $projectDir/bin/mlpl.asmtopan.py --pred ${pred} --graph ${graph}  --output ${asm_pred} --pbf ${pbf_pred}

    bgzip -k ${graph}
    echo "sample,gfa,pls_score" > "${input_csv}"
    echo "${meta.id}.${asm},${graph}.gz,${pbf_pred}" >> "${input_csv}"
    
    python ${putils} gc_probabilities --input_file ${input_csv} --out_dir . --tmp_dir ./tmp

    python ${pflow} --nopenalty -alpha4 1 -ag ${graph}.gz -gc ${gc_content} -out_dir . -out_file ${bins}  -score ${pbf_pred} -assembler pangenome -seed_len ${seed_len}  -seed_score ${seed_score} -min_pls_len ${min_pls_len}


    python $projectDir/bin/evaluation/transform_pbf_pred.py --input ${bins} --gfa ${graph} --output ${res} 
    """
}

process PBF{
    input:
    tuple val(meta), path(graph), path(pred)
    val (asm)

    output:

    tuple val(meta), path(res), emit: res

    script:
    seed_len = 2000
    seed_score = 0.7
    min_pls_len = 1500

    input_csv = "${meta.id}.${asm}-input.csv"
    gc_content = "${meta.id}.${asm}.gc.tsv"
    bins = "${meta.id}.${asm}.pbf.bins.tsv"
    res = "${meta.id}.${asm}.pbf.pred.tab"
    
    putils = "${binfolder}/pbf/code/plasbin_utils.py"
    pflow = "${binfolder}/pbf/code/plasbin_flow.py"
    
    assembler = ""
    if (asm == "u") {
        assembler = "unicycler"
    }
    else if (asm == "s") {
        assembler = "skesa"
    }
    asm_pred = "null.pred"
    pbf_pred = "${meta.id}.${asm}.pbf.pred"


    """
    python $projectDir/bin/mlpl.asm.py --pred ${pred} --graph ${graph}  --output ${asm_pred} --pbf ${pbf_pred}

    bgzip -k ${graph}
    echo "sample,gfa,pls_score" > "${input_csv}"
    echo "${meta.id},${graph}.gz,${pbf_pred}" >> "${input_csv}"

    python ${putils} gc_probabilities --input_file ${input_csv} --out_dir . --tmp_dir ./tmp


    python ${pflow} -ag ${graph}.gz -gc ${gc_content} -out_dir . -out_file ${bins}  -score ${pbf_pred} -assembler ${assembler} -seed_len ${seed_len}  -seed_score ${seed_score} -min_pls_len ${min_pls_len}


    python $projectDir/bin/evaluation/transform_pbf_pred.py --input ${bins} --gfa ${graph} --output ${res} 

  """

}