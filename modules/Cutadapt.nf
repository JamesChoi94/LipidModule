process Cutadapt {

  tag "$srrAccession"
  publishDir "$params.trimmedReadsDir", mode: "copy", pattern: "*.cutadapt.fastq", enabled: "$params.saveTrimmedFastq"
  publishDir "$params.CutadaptDir", mode: "copy", pattern: "*.cutadapt-report.txt"

  input:
  tuple val(srrAccession), path(fastq_reads)

  output:
  tuple val(srrAccession), path("*.cutadapt.fastq"), emit: trimmed_reads
  path("*.cutadapt-report.txt"), emit: cutadapt_report

  script:
  if ( params.adapterLibrary == 'TruSeq' )
  // TruSeq Universal Adapter R1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
  // TruSeq Universal Adapter R2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
  """
  cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -o ${srrAccession}.cutadapt.fastq \
    -p ${srrAccession}.cutadapt.fastq \
    --cores ${task.cpus} \
    ${fastq_reads[0]} ${fastq_reads[1]} \
    2> ${srrAccession}.cutadapt-report.txt
  """
}