process Build_Index {

  tag "${params.alignerMethod}"
  publishDir "$params.genomeDir", mode: "copy"

  input:
  path(genome_fasta)
  path(annotation_gtf)
  val(alignerMethod)

  output:
  path(star_index), emit: index

  script:
  if(${alignerMethod} == "STAR")
  """
  mkdir star_index
   STAR --runMode genomeGenerate \
    --genomeDir star_index \
    --genomeFastaFiles ${genome_fasta} \
    --sjdbGTFfile  ${annotation_gtf} \
    --runThreadN ${task.cpus}
  """
}