process Build_Index {

  tag "${params.aligner}"
  publishDir "params.indexDir", mode: "copy"

  input:
  path(genome_fasta)
  path(annotation_gtf)

  output:
  path(star_index), emit: index

  script:
  // alignerDir = path("${params.genomeDir}/${params.aligner}_index")
  // mkdirResult = alignerDir.mkdirs()
  if(params.aligner == "STAR")
  """
  mkdir star_index
   STAR --runMode genomeGenerate \
    --genomeDir ${params.indexDir} \
    --genomeFastaFiles ${genome_fasta} \
    --sjdbGTFfile  ${annotation_gtf} \
    --runThreadN ${task.cpus}
  """
}