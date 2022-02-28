process Dump_FASTQ {
  
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

process Subset_Testing_FASTQ {

  tag "$srrAccession"
  publishDir "$params.testReadsDir", mode: "copy", enabled: "$params.testRun"
  
  when:
  params.testRun

  input:
  tuple val(srrAccession), path(fastq_reads)

  output:
  tuple val(srrAccession), path("*_subset_{1,2}.fastq"), emit: test_reads

  script:
  """
  seqtk sample -s 100 ${fastq_reads[0]} 10000 > ${srrAccession}_subset_1.fastq
  seqtk sample -s 100 ${fastq_reads[1]} 10000 > ${srrAccession}_subset_2.fastq
  """
}
