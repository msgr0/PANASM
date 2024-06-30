#!/usr/bin/env nextflow

def convert_species(spec) {
  switch(spec) {
    case 'efea':
      species = 'Enterococcus faecium'
      break
    case 'kpne':
      species = 'Klebsiella pneumoniae'
      break
    case 'abau':
      species = 'Acinetobacter baumannii'
      break
    case 'ecol':
      species = 'Escherichia coli'
      break
    default:
      species = 'Other species'
  }
  return species
}

process MLPLASMIDS {
    input:
    tuple val(meta), path(mixed), path(uni), path(ske)

    output:
    tuple val(meta), path(mixedpred), emit: mixed
    tuple val(meta), path(unipred), emit: uni
    tuple val(meta), path(skepred), emit: ske

    
    script:
    mixedpred = "${meta.id}.mlplasmid.mix.pred"
    unipred = "${meta.id}.mlplasmid.uni.pred"
    skepred = "${meta.id}.mlplasmid.ske.pred"

    mlplas_threshold = '0.5'
    """
      Rscript $projectDir/bin/run_mlplasmids.R ${mixed} ${mixedpred} ${mlplas_threshold} '${meta.species}' TRUE
      Rscript $projectDir/bin/run_mlplasmids.R ${uni} ${unipred} 0.5 '${meta.species}' TRUE
      Rscript $projectDir/bin/run_mlplasmids.R ${ske} ${skepred} 0.5 '${meta.species}' TRUE
    """
}