process Build_Index {

  tag "${params.aligner}"
  publishDir "$params.genomeDir", mode: "copy"

  input:
  path(genome_fasta)
  path(annotation_gtf)

  output:
  path(star_index), emit: index

  script:
  if(params.alignerMethod == "STAR")
  """
  mkdir star_index
   STAR --runMode genomeGenerate \
    --genomeDir star_index \
    --genomeFastaFiles ${genome_fasta} \
    --sjdbGTFfile  ${annotation_gtf} \
    --runThreadN ${task.cpus}
  """
}