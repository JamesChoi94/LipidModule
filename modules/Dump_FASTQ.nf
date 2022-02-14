process Dump_FASTQ {
  
  tag "$srrAccession"
  publishDir "$params.readsDir", mode: "copy"
  label "Dump_FASTQ"

  input:
  val(srrAccession)

  output:
  tuple val(srrAccession), path("*.fastq"), emit: fastq_reads

  script:
  """
  fasterq-dump --progress --threads $task.cpus --split-files ${srrAccession}
  """
}

process Check_FASTQ {
  
  // tag "$srrAccession"
  // publishDir "$params.readsDir", mode: "copy"
  // label "Dump_FASTQ"

  input:
  val(srrAccession)

  output:
  val(missing_srrAccession), emit: missing_srrAccession

  script:
  """
  Rscript --vanilla bin/check_downloaded_reads.R $srrAcccession $params.readsDir
  """
}