# Run esearch from Entrez Direct tools to query which SRA codes map to a given GSM accession. 
run_esearch <- function(gsm) {
  gsm2sra_runs <- paste0(gsm2sra_path, '/runs')
  dir.create(path = gsm2sra_runs)
  command_str <- paste0('esearch -db sra -query ', gsm, '| efetch -format runinfo > ', gsm, '_esearch-runinfo.txt')
  system(command = command_str)
}

# Create a summary tsv of all queried GSMs in a given 'runs' directory.
summarize_runs <- function(samples_sheet) {
  gsm2sra_runs <- paste0(gsm2sra_path, '/runs')
  queried_runs <- list.files(path = gsm2sra_runs, full.names = TRUE)
  if (length(queried_runs) == 0) {
    stop('No query results found in runs folder.')
  }
  queried_runs <- lapply(
    X = queried_runs,
    FUN = readLines
  )
  headers <- lapply(queried_runs, `[`, 1)
  identicalValue <- function(x, y) if (identical(x, y)) x else FALSE
  identical_headers <- Reduce(f = identicalValue, x = headers)
  if (isFALSE(identical_headers)) {
    stop('Headers of query result files are not identical across runs.')
  } else {
    identical_headers <- unlist(strsplit(identical_headers, split = ','))
  }
  # unlist for potential multiple SRA to single GSM mappings
  queried_runs <- unlist(sapply(queried_runs, function(a) a[2:length(a)]))
  queried_runs <- strsplit(x = queried_runs, split = ',')
  queried_runs <- do.call(what = rbind.data.frame, args = queried_runs)
  colnames(queried_runs) <- identical_headers
  queried_runs <- dplyr::left_join(queried_runs, samples_sheet, by = 'SampleName')
  write.table(x = queried_runs, file = 'gsm2sra_query_compiled.tsv', sep = '\t', quote = FALSE, row.names = FALSE)
}
