process Raw_FastQC {
  
  tag "$srrAccession"
  publishDir "$params.FastQCDir/raw/", mode: "copy"
  label "quality_control"

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
  
  tag "$srrAccession"
  publishDir "$params.FastQCDir/trimmed/", mode: "copy"
  label "quality_control"

  input:
  tuple val(srrAccession), path(fastq_reads)

  output:
  tuple val(srrAccession), path( "*fastqc.html"), emit: trimmed_fastqc_report

  script:
  """
  fastqc --threads ${task.cpus} ${fastq_reads}
  """
}

process Aligned_FastQC {

  tag "$srrAccession"
  publishDir "$params.FastQCDir/aligned/", mode: "copy"
  label "quality_control"

  input:
  tuple val(srrAccession), path(fastq_reads)

  output:
  tuple val(srrAccession), path( "*fastqc.html"), emit: aligned_fastqc_report

  script:
  """
  fastqc --threads ${task.cpus} ${fastq_reads}
  """
}