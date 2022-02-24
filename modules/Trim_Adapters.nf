process Trim_Adapters {

  tag "$srrAccession"
  publishDir "$params.trimmedReadsDir",  mode: "copy", pattern: "*.bbduk.fastq", enabled: "$params.saveTrimmedFastq"
  publishDir "$params.BBdukDir", mode: "copy", pattern: "*_bbduk-stats.txt"

  label "Trim_Adapters"

  input:
  tuple val(srrAccession), path(fastq_reads)

  output:
  tuple val(srrAccession), path( "*.bbduk.fastq"), emit: trimmed_reads
  path("*_bbduk-stats.txt"), emit: trimmed_reads_report

  script:
  """
  bbduk.sh \
    -Xmx2g \
    in=${fastq_reads[0]} \
    in2=${fastq_reads[1]} \
    out=${srrAccession}_1.bbduk.fastq \
    out2=${srrAccession}_2.bbduk.fastq \
    stats=${srrAccession}_bbduk-stats.txt \
    threads=${task.cpus} \
    ref=adapters \
    k=21 \
    hdist=1 \
    ktrim=r \
    mink=10
  """
}