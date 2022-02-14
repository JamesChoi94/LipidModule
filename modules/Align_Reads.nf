process Align_Reads {

  tag "${params.aligner}_${srrAccession}"
  publishDir "params.genomeDir", mode: "copy"

  input:
  

  output:
  path("${params.aligner}_index"), emit: index

  script:
  if(params.aligner == "STAR")
  """
  echo ${params.genomeDir}/${aligner}_index
  """
}
  // STAR --runMode alignReads \
  //   --genomeDir ${params.genomeDir} \
  //   --genomeFastaFiles ${params.genomeFasta} \
  //   --sjdbGTFfile  ${params.annotationGTF} \
  //   --runThreadN ${task.cpus}