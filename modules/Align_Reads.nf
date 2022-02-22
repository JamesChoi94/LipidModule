process Load_Genome {

  // Load the genome + index into memory before running read alignment. This 
  // allows the genome in-memory to be shared across parallel processes 
  // running the alignment and prevent it from being removed from memory.

  tag "${alignerMethod}"

  input:
  val(alignerMethod)
  path(index)
  
  output:
  val(true), emit: load_genome

  script:
  if(alignerMethod == "STAR")
  """
  STAR \
    --runMode alignReads \
    --genomeDir ${index} \
    --genomeLoad LoadAndExit
  """
}

process Align_Reads {

  // Run read alignment against an indexed genome. `val(genome_loaded)` input 
  // exists only to wait for Load_Genome process to complete before starting
  // Align_Reads process to ensure that genome is available in memory.

  tag "${alignerMethod}_${srrAccession}"
  publishDir path: { params.saveBAMs ? params.bamsDir: "work/dump" }, mode: "copy"

  input:
  tuple val(srrAccession), path(fastq_reads)
  val(alignerMethod)
  path(index)
  val(genome_loaded)
  
  output:
  tuple val(srrAccession), path("*"), emit: aligned_bams

  script:
  if(alignerMethod == "STAR")
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