#!/usr/bin/env nextflow

nextflow.enable.dsl=2


// ##################################################################
// Import modules 
// ##################################################################

// include {Query_GEO; Compile_GEO_Queries; Extract_SRR} from "./modules/Query_GEO"
include {GZIP_Compress} from "./modules/utils.nf"
include {Dump_FASTQ; Subset_Testing_FASTQ} from "./modules/Dump_FASTQ"
include {Raw_FastQC; Trimmed_FastQC; Aligned_FastQC} from "./modules/FastQC"
include {BBDuk} from "./modules/BBTools"
include {Build_Index; Load_Index; Unload_Index} from "./modules/Build_Index"
include {Align_Reads; Index_BAMs} from "./modules/Align_Reads"
include {Convert_GTF2BED} from "./modules/Bedparse.nf"
include {BAM_Stat; Read_Distribution; Infer_Experiment; Gene_Body_Coverage} from "./modules/RSeQC.nf"
include {MultiQC} from "./modules/MultiQC.nf"


// ##################################################################
// Intro message
// ##################################################################

def scriptbreak = "==========================================================================="

Date date1          = new Date()
String dateStart    = date1.format("yyyy-MM-dd -- ")
String timeStart    = date1.format("HH:mm:ss")
def runStart        = dateStart + timeStart

println ""
println "$scriptbreak"
println "Workflow:            rnaseq_preprocess"
println "Author:              James Choi"
println "Author contact:      jsc228 at miami dot edu"
println "RunName:             $workflow.runName"
println "RunProfile:          $workflow.profile"
println "testRun:             $params.testRun"
println "Start:               $runStart"
println "$scriptbreak"
println "Using the following params:"
println ""
println "saveRawFastq                   $params.saveRawFastq"
println "saveTrimmedFastq:              $params.saveTrimmedFastq"
println "saveBAMs:                      $params.saveBAMs"
println "saveQuants:                    $params.saveQuants"
println "compressFastq:                 $params.compressFastq"
println "buildNewIndex:                 $params.buildNewIndex"
println "alignerMethod:                 $params.alignerMethod"
println "genomeDir:                     $params.genomeDir"
println "existingIndex:                 $params.existingIndex"
println "genomeFasta:                   $params.genomeFasta"
println "annotationGTF:                 $params.annotationGTF"
println "readsLength:                   $params.readsLength"
println "$scriptbreak"


// ##################################################################
// Main workflow run
// ##################################################################

workflow {

  // Set input channels ----------------------------------------

  samplesheet   = Channel.fromPath(params.samplesheet)
  alignerMethod = Channel.value(params.alignerMethod)
  genomeFasta   = Channel.fromPath(params.genomeFasta)
  annotationGTF = Channel.fromPath(params.annotationGTF)
  readsLength   = Channel.value(params.readsLength) // This will need to be modified if reads have different length
  srrAccession  = Channel.fromPath(params.srrAccessions)
    .splitText()
    .map{it -> it.trim()}
  if ( params.testRun ) {
    srrAccession = srrAccession.take(1)
  }

  // fasterq-dump wrapper -----------------------------------------

  raw_reads = params.testRun
    ? Channel.fromFilePairs(params.testReadsDir + "/*_{1,2}.subset.fastq")
    : Dump_FASTQ(srrAccession).out.raw_reads
  
  // // Take subset of reads if test run 
  // Subset_Testing_FASTQ(raw_reads)
  // raw_reads = Subset_Testing_FASTQ.out.test_reads


  // FastQC on raw reads ------------------------------------------

  Raw_FastQC(raw_reads)


  // Trim adapter sequences using BBduk ---------------------------

  BBDuk(raw_reads)
  trimmed_reads = BBDuk.out.trimmed_reads


  // FastQC on trimmed reads --------------------------------------

  Trimmed_FastQC(trimmed_reads)


  // Build genome index -------------------------------------------

  index = params.existingIndex
    ? Channel.fromPath(params.existingIndex) 
    : Build_Index(genomeFasta, annotationGTF, alignerMethod)
  
  // Create BED format annotation for RSeQC
  Convert_GTF2BED(annotationGTF)
  annotationBED = Convert_GTF2BED.out.annotationBED
  

  // Align reads ----------------------------------------------------

  Load_Index(alignerMethod, index)
  genome_loaded = Load_Index.out.genome_loaded
  Align_Reads(trimmed_reads, alignerMethod, index, genome_loaded)
  Align_Reads(trimmed_reads, alignerMethod, index)
  aligned_bams = Align_Reads.out.aligned_bams
  unload_genome = Align_Reads.out.unload_genome
  Unload_Index(alignerMethod, index, unload_genome)
  Index_BAMs(aligned_bams, alignerMethod)
  

  // Post-alignment QC ---------------------------------------------

  Aligned_FastQC(aligned_bams)
  BAM_Stat(aligned_bams)
  Read_Distribution(aligned_bams, annotationBED)
  Infer_Experiment(aligned_bams, annotationBED)
  // Gene_Body_Coverage(aligned_bams, annotationBED)


  // MultiQC: collect reports --------------------------------------
  
  combined_reports = Raw_FastQC.out.raw_fastqc_report
    .concat(
      Trimmed_FastQC.out.trimmed_fastqc_report,
      Aligned_FastQC.out.aligned_fastqc_report,
      BBDuk.out.bbduk_report,
      BAM_Stat.out.bam_stat_out,
      Read_Distribution.out.read_distribution_out,
      Infer_Experiment.out.infer_experiment_out
    )
    .collect()
  MultiQC(combined_reports)

  // Finish -------------------------------------------------------

  workflow.onComplete {
    Date date2          = new Date()
    String dateEnd      = date2.format("yyyy-MM-dd -- ")
    String timeEnd      = date2.format("HH:mm:ss")
    def runEnd          = dateEnd + timeEnd

    println ""
    println "$scriptbreak"
    println "Pipeline (rnaseq-preprocess) complete!"
    println "End: $runEnd"
    println "$scriptbreak"
    println ""
  }

}