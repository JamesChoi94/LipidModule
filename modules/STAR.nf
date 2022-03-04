process Build_Index {

  tag "${params.alignerMethod}"
  publishDir "$params.genomeDir", mode: "copy"
  label "large_mem"

  input:
  path(genome_fasta)
  path(annotation_gtf)
  val(alignerMethod)

  output:
  path(star_index), emit: index

  script:
  if(alignerMethod == "STAR")
  """
  mkdir -p star_index
   STAR --runMode genomeGenerate \
    --genomeDir star_index \
    --genomeFastaFiles ${genome_fasta} \
    --sjdbGTFfile  ${annotation_gtf} \
    --runThreadN ${task.cpus}
  """
}

process Load_Index {

  // Load the genome + index into memory before running read alignment. This 
  // allows the genome in-memory to be shared across parallel processes 
  // running the alignment and prevent it from being removed from memory.

  tag "${alignerMethod}"
  label "large_mem"

  input:
  val(alignerMethod)
  path(index)
  
  output:
  val(true), emit: genome_loaded

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
  publishDir "$params.bamsDir", mode: "copy", pattern: "*.[bs]am", enabled: "$params.saveBAMs"
  publishDir "$params.resultsDir/$params.alignerMethod/logs", mode: "copy", pattern: "*.out"
  publishDir "$params.resultsDir/$params.alignerMethod/quants", mode: "copy", pattern: "*.tab", enabled: "$params.saveQuants"
  label "large_mem"

  input:
  tuple val(srrAccession), path(fastq_reads)
  val(alignerMethod)
  path(index)
  
  output:
  tuple val(srrAccession), path("*Aligned.out.ba"), emit: aligned_bams
  tuple val(srrAccession), path("*toTranscriptome.out.bam"), emit: transcriptome_bams
  tuple val(srrAccession), path("*.out"), emit: aligner_outs
  tuple val(srrAccession), path("*.tab"), emit: quants
  val(true), emit: unload_genome

  script:
  if (alignerMethod == "STAR")
  """
  STAR \
    --genomeDir ${index} \
    --runMode alignReads \
    --readFilesIn ${fastq_reads[0]} ${fastq_reads[1]} \
    --readFilesCommand zcat \
    --runThreadN ${task.cpus} \
    --outFileNamePrefix ${srrAccession}_ \
    --outSAMtype BAM Unsorted \
    --quantMode TranscriptomeSAM GeneCounts
  """
}

process Unload_Index {

  // Unload the genome, opposite of above.

  tag "${alignerMethod}"
  label "large_mem"

  input:
  val(alignerMethod)
  path(index)
  val(unload_genome)
  
  output:
  val(true), emit: load_genome

  script:
  if(alignerMethod == "STAR" & unload_genome)
  """
  STAR \
    --runMode alignReads \
    --genomeDir ${index} \
    --genomeLoad Remove
  """
}