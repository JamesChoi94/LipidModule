process Convert_GTF2BED {

  tag "$params.annotationGTF"
  publishDir "${params.genomeDir}", mode: "copy"

  input:
  path(annotationGTF)
  
  output:
  path(annotationBED), emit: annotation_bed

  script:
  """
  bedparse gtf2bed ${annotationGTF} > ${params.annotationGTF}.bed
  """
}