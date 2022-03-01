process Fasterq_Dump {
  
  tag "$srrAccession"
  publishDir "$params.rawReadsDir", mode: "copy", enabled: "$params.saveRawFastq"

  input:
  val(srrAccession)

  output:
  tuple val(srrAccession), path("*.fastq"), emit: raw_reads

  script:
  """
  fasterq-dump --progress --threads $task.cpus --split-files ${srrAccession}
  """
}