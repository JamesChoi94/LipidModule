process Build_Index {

  tag "${params.aligner}"
  publishDir "params.genomeDir", mode: "copy"

  output:
  path("${params.aligner}_index"), emit: index

  script:
  alignerDir = path("${params.genomeDir}/${params.aligner}_index")
  mkdirResult = alignerDir.mkdirs()
  if(params.aligner == "STAR")
  """
   STAR --runMode genomeGenerate \
    --genomeDir ${params.genomeDir}/${params.aligner}_index \
    --genomeFastaFiles ${params.genomeFasta} \
    --sjdbGTFfile  ${params.annotationGTF} \
    --runThreadN ${task.cpus}
  """
}