process Cutadapt {

  tag "$srrAccession"
  publishDir "$params.trimmedReadsDir", mode: "copy", pattern: "*.cutadapt.fastq.gz", enabled: "$params.saveTrimmedFastq"
  publishDir "$params.CutadaptDir", mode: "copy", pattern: "*.cutadapt-report.txt"

  input:
  tuple val(srrAccession), path(fastq_reads)

  output:
  tuple val(srrAccession), path("*.cutadapt.fastq.gz"), emit: trimmed_reads
  path("*.cutadapt-report.txt"), emit: cutadapt_report

  script:
  if ( params.adapterLibrary == 'TruSeq' )
  // TruSeq Universal Adapter R1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
  // TruSeq Universal Adapter R2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
  """
  cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -m 5 \
    -o ${srrAccession}_1.cutadapt.fastq.gz \
    -p ${srrAccession}_2.cutadapt.fastq.gz \
    --cores ${task.cpus} \
    ${fastq_reads[0]} ${fastq_reads[1]} \
    1> ${srrAccession}.cutadapt-report.txt
  """
}