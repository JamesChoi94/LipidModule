process BAM_Stat {
  
  tag "$srrAccession"
  publishDir "$params.RSeQCDir/aligned/$params.alignerMethod", mode: "copy"
  label "quality_control"

  input:
  tuple val(srrAccession), path(aligned_bams)

  output:
  tuple val(srrAccession), path( "*.bam_stat.txt"), emit: bam_stat_out

  script:
  """
  bam_stat.py -i ${aligned_bams} > ${aligned_bams}.bam_stat.txt
  """
}

process Read_Distribution {
  
  tag "$srrAccession"
  publishDir "$params.RSeQCDir/aligned/$params.alignerMethod", mode: "copy"
  label "quality_control"

  input:
  tuple val(srrAccession), path(aligned_bams)
  path(annotationBED)

  output:
  tuple val(srrAccession), path( "*.read_distribution.txt"), emit: read_distribution_out

  script:
  """
  read_distribution.py \
    -i ${aligned_bams} \
    -r ${annotationBED} \
    > ${aligned_bams}.read_distribution.txt
  """
}

process Infer_Experiment {

  tag "$srrAccession"
  publishDir "$params.RSeQCDir/aligned/$params.alignerMethod", mode: "copy"
  label "quality_control"

  input:
  tuple val(srrAccession), path(aligned_bams)
  path(annotationBED)

  output:
  tuple val(srrAccession), path( "*.infer_experiment.txt"), emit: infer_experiment_out

  script:
  """
  infer_experiment.py \
    -i ${aligned_bams} \
    -r ${annotationBED} \
    > ${aligned_bams}.infer_experiment.txt
  """
}


// WORK IN PROGRESS. DO NOT IMPLEMENT. CAN TAKE MULTIPLE BAMS AS INPUT TO PRODUCE GRAPH WITH
// ALL BAMS QUANTIFIED. MAYBE BETTER TO RUN AFTER ALL BAMS PRODUCED AND TO PROVIDE SINGLE
// BAM PATH AS INPUT CHANNEL HERE.
process Gene_Body_Coverage {

  tag "$srrAccession"
  publishDir "$params.RSeQCDir/aligned/$params.alignerMethod", mode: "copy"
  label "quality_control"

  input:
  tuple val(srrAccession), path(aligned_bams)
  path(annotationBED)

  output:
  tuple val(srrAccession), path("*.gene_body_coverage.txt"), emit: gene_body_coverage

  script:
  """
  geneBody_coverage.py \
    -i ${aligned_bams} \
    -r ${annotationBED} \
    -f png \
    -o ${aligned_bams}.gene
  """
}