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


process Unload_Index {

  // Unload the genome, opposite of above.

  tag "${alignerMethod}"

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
