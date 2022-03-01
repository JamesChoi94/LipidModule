process Index_BAMs {

  tag "${alignerMethod}_${srrAccession}"
  publishDir "$params.bamsDir", mode: "copy", enabled: "$params.saveBAMs", pattern: "*.bai"

  input:
  tuple val(srrAccession), path(aligned_bams)
  val(alignerMethod)
  
  output:
  tuple val(srrAccession), path("*.bai"), emit: bam_index

  script:
  if (alignerMethod == 'STAR')
  """
  samtools index -b ${aligned_bams} ${aligned_bams}.bai
  """
}