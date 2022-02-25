process Convert_GTF2BED {

  tag "$params.annotationGTF"
  publishDir "${params.genomeDir}", mode: "copy"

  input:
  path(annotationGTF)
  
  output:
  path("*.bed"), emit: annotation_bed

  script:
  """
  bedparse gtf2bed ${annotationGTF} >  ${annotationGTF}.bed
  """
}