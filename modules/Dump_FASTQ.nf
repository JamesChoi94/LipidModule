process Dump_FASTQ {
  
  tag "$srrAccession"
  publishDir "$params.rawReadsDir", mode: "copy"
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


process Subset_Testing_FASTQ {

  tag "$srrAccession"
  publishDir "$params.testReadsDir", mode: "copy"
  
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

// process Check_FASTQ {
  
//   // tag "$srrAccession"
//   // publishDir "$params.rawReadsDir", mode: "copy"
//   // label "Dump_FASTQ"

//   input:
//   val(srrAccession)

//   output:
//   val(missing_srrAccession), emit: missing_srrAccession

//   script:
//   """
//   Rscript --vanilla bin/check_downloaded_reads.R $srrAcccession $params.rawReadsDir
//   """
// }