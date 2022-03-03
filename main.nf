#!/usr/bin/env nextflow

nextflow.enable.dsl=2


// ##################################################################
// Import modules 
// ##################################################################

include {GZIP_Compress} from "./modules/utils"
include {Sample_FASTQ} from "./modules/Seqtk"
include {Gtf2Bed} from "./modules/Bedparse"
include {Raw_FastQC; Trimmed_FastQC; Aligned_FastQC} from "./modules/FastQC"
include {Cutadapt} from "./modules/Cutadapt"
include {Build_Index; Load_Index; Align_Reads; Unload_Index} from "./modules/STAR"
include {Index_BAMs} from "./modules/Samtools"
include {BAM_Stat; Read_Distribution; Infer_Experiment; Gene_Body_Coverage} from "./modules/RSeQC"
include {MultiQC} from "./modules/MultiQC"


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
println "referenceDir:                  $params.referenceDir"
println "genomeIndex:                   $params.genomeIndex"
println "genomeFasta:                   $params.genomeFasta"
println "annotationGTF:                 $params.annotationGTF"
println "annotationBED:                 $params.annotationBED"
println "readsLength:                   $params.readsLength"
println "$scriptbreak"


// ##################################################################
// Main workflow run
// ##################################################################

workflow {

  // Set input channels ----------------------------------------
  
  raw_reads = params.testRun 
    ? Channel.fromFilePairs(params.testReadsDir + "/*{1,2}.fastq.gz")
    : Channel.fromFilePairs(params.rawReadsDir + "/*{1,2}.fastq.gz")
  alignerMethod = Channel.value(params.alignerMethod)
  genomeFasta = Channel.fromPath(params.genomeFasta)
  annotationGTF = Channel.fromPath(params.annotationGTF)
  readsLength = Channel.value(params.readsLength)
  annotationBED = params.annotationBED
    ? Channel.fromPath(params.annotationBED)
    : Gtf2Bed(annotationGTF)
  index = params.genomeIndex
    ? Channel.fromPath(params.genomeIndex) 
    : Build_Index(genomeFasta, annotationGTF, alignerMethod)
  

  // FastQC on raw reads ------------------------------------------

  Raw_FastQC(raw_reads)


  // Trim adapter sequences using BBduk ---------------------------

  Cutadapt(raw_reads)
  trimmed_reads = Cutadapt.out.trimmed_reads


  // FastQC on trimmed reads --------------------------------------

  Trimmed_FastQC(trimmed_reads)


  // Align reads ----------------------------------------------------

  // Load_Index(alignerMethod, index)
  // genome_loaded = Load_Index.out.genome_loaded
  // Align_Reads(trimmed_reads, alignerMethod, index, genome_loaded)
  Align_Reads(trimmed_reads, alignerMethod, index)
  aligned_bams = Align_Reads.out.aligned_bams
  // unload_genome = Align_Reads.out.unload_genome
  // Unload_Index(alignerMethod, index, unload_genome)
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
      Cutadapt.out.cutadapt_report,
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