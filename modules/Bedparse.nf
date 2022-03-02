process Gtf2Bed {

  tag "$params.annotationGTF"
  publishDir "${params.resultsDir}", mode: "copy"

  input:
  path(annotationGTF)
  
  output:
  path("*.bed"), emit: annotationBED

  script:
  """
  bedparse gtf2bed ${annotationGTF} > ${annotationGTF}.bed
  """
}