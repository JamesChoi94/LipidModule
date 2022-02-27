process GZIP_Compress {

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