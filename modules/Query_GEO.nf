process Query_GEO {
  
  cache: "lenient"
  publishDir "$params.queryGEODir", mode: "copy"
  label "Query_GEO"

  input:
  path(samplesheet) 

  output:
  path "esearch-runinfo.txt", emit: geo_queries

  script:
  """
  Gsm2Srr.sh --samplesheet ${samplesheet} > esearch-runinfo.txt
  """
}

process Compile_GEO_Queries {

  cache: "lenient"
  publishDir "$params.dataDir",  mode: "copy"

  input:
  path(geo_queries)
  path(samplesheet)

  output:
  path "samplesheet_appended.csv", emit: samplesheet_appended

  script:
  """
  Rscript --vanilla ${projectDir}/bin/compile_geo_queries.R \
    ${geo_queries} ${samplesheet} > "samplesheet_appended.csv"
  """
}
