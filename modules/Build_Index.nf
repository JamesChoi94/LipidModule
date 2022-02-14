process Build_Index {

  tag "${params.aligner}"
  publishDir "params.genomeDir", mode: "copy"

  output:
  path("${params.aligner}_index"), emit: index

  script:
  if(params.aligner == "STAR")
  """
   STAR --runMode genomeGenerate \
    --genomeDir ${params.genomeDir} \
    --genomeFastaFiles ${params.genomeFasta} \
    --sjdbGTFfile  ${params.annotationGTF} \
    --runThreadN ${task.cpus}
  """
}
  // """
  // STAR --runMode genomeGenerate \
  //   --genomeDir ${params.genomeDir} \
  //   --genomeFastaFiles ${params.testGenomeFasta} \
  //   --sjdbGTFfile  ${params.annotationGTF} \
  //   --runThreadN ${task.cpus}
  // """
  // if(params.aligner == "hisat2")
  // """
  // echo hello
  // """