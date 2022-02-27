process MultiQC {

  publishDir "$params.MultiQCDir", mode: "copy"
  label "quality_control"
  
  input:
  path("logs/*")

  output:
  path("multiqc_report.html"), emit: multiqc_report
  path("multiqc_data")

  script:
  """
  multiqc .
  """
}