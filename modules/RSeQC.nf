process RSeQC_bam_stat {
  
  tag "$srrAccession"
  publishDir "$params.RSeQCDir/aligned/", mode: "copy"
  label "RSeQC_bam_stat"

  input:
  tuple val(srrAccession), path(aligned_bams)

  output:
  tuple val(srrAccession), path( "*fastqc.html"), emit: raw_fastqc_report
  // path "*_.fastqc.zip", emit: raw_fastqc_zip

  script:
  """
  fastqc --threads ${task.cpus} ${fastq_reads}
  """
}