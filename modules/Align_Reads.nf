process Align_Reads {

  // Run read alignment against an indexed genome. `val(genome_loaded)` input 
  // exists only to wait for Load_Genome process to complete before starting
  // Align_Reads process to ensure that genome is available in memory.

  tag "${alignerMethod}_${srrAccession}"
  publishDir "$params.bamsDir", mode: "copy", pattern: "*.[bs]am", enabled: "$params.saveBAMs"
  publishDir "$params.resultsDir/$params.alignerMethod", mode: "copy", pattern: "*.out"
  publishDir "$params.resultsDir/$params.alignerMethod/quants", mode: "copy", pattern: "*.tab", enabled: "$params.saveQuants"
  label "large_mem"

  input:
  tuple val(srrAccession), path(fastq_reads)
  val(alignerMethod)
  path(index)
  val(genome_loaded)
  
  output:
  tuple val(srrAccession), path("*.[bs]am"), emit: aligned_bams
  tuple val(srrAccession), path("*.out"), emit: aligner_outs
  tuple val(srrAccession), path("*.tab"), emit: quants
  val(true), emit: unload_genome

  script:
  if ((alignerMethod == "STAR") & genome_loaded)
  """
  STAR \
    --genomeDir ${index} \
    --runMode alignReads \
    --readFilesIn ${fastq_reads[0]} ${fastq_reads[1]} \
    --runThreadN ${task.cpus} \
    --outFileNamePrefix ${srrAccession}_ \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts
  """
}


process Index_BAMs {

  tag "${alignerMethod}_${srrAccession}"
  publishDir "$params.bamsDir", mode: "copy", enabled: "$params.saveBAMs", pattern: "*.bai"

  input:
  tuple val(srrAccession), path(aligned_bams)
  val(alignerMethod)
  
  output:
  tuple val(srrAccession), path("*.bai"), emit: bam_index

  script:
  if (alignerMethod == 'STAR')
  """
  samtools index -b ${aligned_bams} ${aligned_bams}.bai
  """
}