nextflow.enable.dsl=2

include {Build_Index} from "./modules/Build_Index"

workflow {
  annotationGTF = Channel.fromPath(params.annotationGTF)
  genomeFasta   = Channel.fromPath(params.genomeFasta)
  alignerMethod = Channel.value(params.alignerMethod)
  index = params.existingIndex
    ? Channel.fromPath(params.existingIndex) 
    : Build_Index(genomeFasta, annotationGTF, alignerMethod)
  index.view()
}