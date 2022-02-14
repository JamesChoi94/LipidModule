#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {Query_GEO; Compile_GEO_Queries} from "./modules/Query_GEO"
include {Dump_FASTQ} from "./modules/Dump_FASTQ"
include {Raw_FastQC; Trimmed_FastQC} from "./modules/FastQC"
include {Trim_Adapters} from "./modules/Trim_Adapters"
include {Build_Index} from "./modules/Build_Index"
include {Align_Reads} from "./modules/Align_Reads"

// Intro message
//-------------------------------------------------------------------------

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
println "Using the following params:"
println ""
println "readsDir:              $params.readsDir"
println "queryGEODir:           $params.queryGEODir"
println "publishmode:           $params.publishmode"
println "sleep:                 $params.sleep"
println "genome:                $params.genome"
println "resultsDir:            $params.resultsDir"
println "samplesheet:           $params.samplesheet"
println "dataDir:               $params.dataDir"
println ""


workflow {

  def readFiles = params.readsDir + "/" + "*_{1,2}.fastq"
  println "$readFiles"

  readPairs     = Channel.fromFilePairs(readFiles, checkIfExists:true)
  alignerMethod = Channel.of('hisat2', 'STAR')
  trimmerMethod = Channel.of('bbduk.sh', 'cutadapt')
  samplesheet   = Channel.fromPath(params.samplesheet)
  reads         = Channel.fromPath(params.readsDir + "/*_{1,2}.fastq")
  
  /* Query GEO db for GSM to SRR mappings */
  Query_GEO(samplesheet)
  geo_queries = Query_GEO.out.geo_queries
  Compile_GEO_Queries(geo_queries, samplesheet)
  samplesheet = Compile_GEO_Queries.out.samplesheet_appended
  
  /* Get fastq files from SRA by SRR accession */
  srrAccession = Compile_GEO_Queries.out.samplesheet_appended
  .splitCsv(header:true)
  .map{ r ->
    srrAccession = r['SRR_Accession']
    return srrAccession
  }
  .take(1)

  // Uncomment DumpFASTQ for runs ------
  // Dump_FASTQ(srrAccession)
  // fastq_reads = Dump_FASTQ.out.fastq_reads
  
  // Use this chunk for test runs ------
  test_reads = Channel
    .fromFilePairs(params.testReads)
  test_reads
    .view{r -> "key: ${r[0]} read1: ${r[1][0]} read2: ${r[1][1]}"}
  // ------------------------------------

  Raw_FastQC(test_reads)
  Trim_Adapters(test_reads)
  // trimmed_reads = Trim_Adapters.out.trimmed_reads
  Trimmed_FastQC(Trim_Adapters.out.trimmed_reads)


  // Build_Index()
  
  // test_reads.get(1).view()
  println "${params.genomeDir}/${params.aligner}_index"
  Align_Reads(test_reads)
}
