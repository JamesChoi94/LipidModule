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

process Convert_GTF2BED {
  publishDir "$params.genomeDir", mode: "copy"

  input:
  path(annotation_gtf)

  output:
  path("*.bed"), emit: annotation_bed

  script:
  """
  bedparse gtf2bed $annotation_gtf > ${annotation_gtf}.bed
  """
}