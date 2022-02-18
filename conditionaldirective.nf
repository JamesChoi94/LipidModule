nextflow.enable.dsl=2

params.inputFile  = "data/testdata/SRR789190_1.subset.fastq"
params.saveFile   = false
params.saveDir    = "testResult"

process TestConditionalDirective {

  publishDir path: { params.saveFile ? params.saveDir : "work/dump" }, mode: 'copy'
  
  input:
  path(input_file)

  output:
  path(output_file)

  script:
  """
  head -n 10 ${input_file} > output_file
  """
}

workflow {
  testInput = Channel.fromPath(params.inputFile)
  testInput.view()
  TestConditionalDirective(testInput)
}