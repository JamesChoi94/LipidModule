process Dump_FASTQ {
  
  tag "$srrAccession"
  publishDir "$params.rawReadsDir", mode: "copy", enabled: "$params.saveRawFastq"
  label "Dump_FASTQ"

  input:
  val(srrAccession)

  output:
  tuple val(srrAccession), path("*.fastq"), emit: raw_reads

  script:
  """
  fasterq-dump --progress --threads $task.cpus --split-files ${srrAccession}
  """
}

process Compress_FASTQ {

  tag "$srrAccession"
  publishDir "$params.rawReadsDir", mode: "copy", enabled: "$params.compressFastq"

  input:
  tuple val(srrAccession), path(fastq_reads)

  output:
  tuple val(srrAccession), path("*.fastq.gz"), emit: gz_raw_reads

  script:
  """
  gzip --stdout ${fastq_reads[0]} > ${srrAccession}_1.fastq.gz
  gzip --stdout ${fastq_reads[1]} > ${srrAccession}_2.fastq.gz
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
  tuple val(srrAccession), path("*.subset.fastq"), emit: test_reads

  script:
  """
  seqtk sample -s 100 ${fastq_reads[0]} 10000 > ${srrAccession}_1.subset.fastq
  seqtk sample -s 100 ${fastq_reads[1]} 10000 > ${srrAccession}_2.subset.fastq
  """
}
