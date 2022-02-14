process Raw_FastQC {
  
  tag "$fastq_reads"
  publishDir "$params.FastQCDir/raw/", mode: "copy"
  label "Raw_FastQC"

  input:
  tuple val(srrAccession), path(fastq_reads)

  output:
  tuple val(srrAccession), path( "*fastqc.html"), emit: raw_fastqc_report
  // path "*_.fastqc.zip", emit: raw_fastqc_zip

  script:
  """
  fastqc --threads ${task.cpus} ${fastq_reads}
  """
}

process Trimmed_FastQC {
  
  tag "$fastq_reads"
  publishDir "$params.FastQCDir/trimmed/", mode: "copy"
  label "Trimmed_FastQC"

  input:
  tuple val(srrAccession), path(fastq_reads)

  output:
  tuple val(srrAccession), path( "*fastqc.html"), emit: trimmed_fastqc_report

  script:
  """
  fastqc --threads ${task.cpus} ${fastq_reads}
  """
}