process Align_Reads {

  // Run read alignment against an indexed genome. `val(genome_loaded)` input 
  // exists only to wait for Load_Genome process to complete before starting
  // Align_Reads process to ensure that genome is available in memory.

  tag "${alignerMethod}_${srrAccession}"
  publishDir "$params.bamsDir", mode: "copy", enabled: "$params.saveBAMs"

  input:
  tuple val(srrAccession), path(fastq_reads)
  val(alignerMethod)
  path(index)
  val(genome_loaded)
  
  output:
  tuple val(srrAccession), path("*"), emit: aligned_bams
  val(true), emit: unload_genome

  script:
  if ((alignerMethod == "STAR") & genome_loaded))
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