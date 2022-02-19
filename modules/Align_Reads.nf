process Align_Reads {

  tag "${alignerMethod}_${srrAccession}"
  publishDir "$params.bamsDir", mode: "copy"

  input:
  tuple val(srrAccession), path(fastq_reads)
  val(alignerMethod)
  // path(index)
  
  output:
  tuple val(srrAccession), path("*"), emit: star_out

  script:
  if(alignerMethod == "STAR")
  """
  STAR \
    --runMode alignReads \
    --genomeDir ${params.genomeDir}/star_index \
    --readFilesIn ${fastq_reads[0]} ${fastq_reads[1]} \
    --runThreadN ${task.cpus} \
    --outFileNamePrefix ${srrAccession}_ \
    --genomeLoad LoadAndExit \
    --outSAMmultNmax 1 
  """
}
// STAR \
//     --genomeDir ${index} \
//     --genomeLoad Remove