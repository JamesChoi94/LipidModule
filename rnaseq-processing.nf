#!/usr/bin/env nextflow

nextflow.enable.dsl=2


// ##################################################################
// Import modules 
// ##################################################################

include {Query_GEO; Compile_GEO_Queries} from "./modules/Query_GEO"
include {Dump_FASTQ; Subset_Testing_FASTQ} from "./modules/Dump_FASTQ"
include {Raw_FastQC; Trimmed_FastQC} from "./modules/FastQC"
include {Trim_Adapters} from "./modules/Trim_Adapters"
include {Build_Index} from "./modules/Build_Index"
include {Align_Reads} from "./modules/Align_Reads"


// ##################################################################
// Intro message
// ##################################################################

def scriptbreak = "==========================================================================="

Date date       = new Date()
String dateNow  = date.format("yyyy-MM-dd -- ")
String timeNow  = date.format("HH:mm:ss")
def runDate     = dateNow + timeNow

println ""
println "$scriptbreak"
println "Workflow:        rnaseq_preprocess"
println "Author:          James Choi"
println "Author contact:  jsc228 at miami dot edu"
println "RunName:         $workflow.runName"
println "RunProfile:      $workflow.profile"
println "Start:           $runDate"
println "$scriptbreak"
println ""
println "Testing mode:     $params.forTesting"
// println "Using the following params:"
// println ""
// println "readsDir:              $params.readsDir"
// println "queryGEODir:           $params.queryGEODir"
// println "publishmode:           $params.publishmode"
// println "sleep:                 $params.sleep"
// println "resultsDir:            $params.resultsDir"
// println "samplesheet:           $params.samplesheet"
// println "dataDir:               $params.dataDir"
// println ""


// ##################################################################
// Main workflow run
// ##################################################################

workflow {

  alignerMethod = Channel.of('hisat2', 'STAR')
  trimmerMethod = Channel.of('bbduk.sh', 'cutadapt')
  samplesheet   = Channel.fromPath(params.samplesheet)
  genomeFasta   = Channel.fromPath(params.genomeFasta)
  annotationGTF = Channel.fromPath(params.annotationGTF)
  
  // Query GEO db for GSM to SRR mappings
  Query_GEO(samplesheet)
  geo_queries = Query_GEO.out.geo_queries
  Compile_GEO_Queries(geo_queries, samplesheet)
  samplesheet = Compile_GEO_Queries.out.samplesheet_appended
  
  // Take single srrAccession if test run
  srrAccession = Compile_GEO_Queries.out.samplesheet_appended
    .splitCsv(header:true)
    .map{ r ->
      srrAccession = r['SRR_Accession']
      return srrAccession
    }
  if ( params.forTesting ) {
    srrAccession = srrAccession.take(1)
  }

  // fasterq-dump wrapper 
  Dump_FASTQ(srrAccession)
  fastq_reads = Dump_FASTQ.out.fastq_reads
  
  // Take subset of reads if test run
  if ( params.forTesting ) {
    Subset_Testing_FASTQ(fastq_reads)
    fastq_reads = Subset_Testing_FASTQ.out.test_reads
  }

  // FastQC on raw reads
  Raw_FastQC(fastq_reads)

  // Trim adapter sequences using BBduk
  Trim_Adapters(fastq_reads)
  trimmed_reads = Trim_Adapters.out.trimmed_reads

  // FastQC on trimmed reads
  Trimmed_FastQC(trimmed_reads)

  // Build genome index. If test run, use existing index.
  if ( params.forTesting ) {
    index = Channel.fromPath(params.genomeDir + "/star_index")
  } else {
    Build_Index(genomeFasta, annotationGTF)
    index = Build_Index.out.index
  }
  
  // Align reads
  Align_Reads(trimmed_reads, index)
  aligned = Align_Reads.out.star_out.view()
  
}
