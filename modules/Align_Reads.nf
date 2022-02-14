process Align_Reads {

  tag "${params.aligner}_${srrAccession}"
  // publishDir "params.genomeDir", mode: "copy"

  input:
  tuple val(srrAccession), path(fastq_reads)
  
  script:
  if(params.aligner == "STAR")
  """
  STAR \
    --genomeDir ${projectDir}/${params.genomeDir}/${params.aligner}_index \
    --genomeLoad Remove
  STAR \
    --runMode alignReads \
    --genomeDir ${params.genomeDir}/${params.aligner}_index \
    --readFilesIn ${fastq_reads[0]} ${fastq_reads[1]} \
    --runThreadN ${task.cpus} \
    --genomeLoad LoadAndRemove \
    --outSAMmultNmax 1 
  """
}