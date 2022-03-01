process Sample_FASTQ {

  tag "$srrAccession"
  publishDir "$params.testReadsDir", mode: "copy", enabled: "$params.testRun"
  
  when:
  params.testRun

  input:
  tuple val(srrAccession), path(fastq_reads)

  output:
  tuple val(srrAccession), path("*_subset_{1,2}.fastq.gz"), emit: test_reads

  script:
  """
  seqtk sample -s 100 ${fastq_reads[0]} 10000 | gzip > ${srrAccession}_subset_1.fastq
  seqtk sample -s 100 ${fastq_reads[1]} 10000 | gzip > ${srrAccession}_subset_2.fastq
  """
}
