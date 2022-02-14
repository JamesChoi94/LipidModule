# Run file for querying NCBI SRA database by GSM accessions

# Import query helper functions
source(file = '/scripts/gsm2sra_functions.r')

# System Environment variables 
ENV_VAR <- Sys.getenv(names = TRUE)
gsm2sra_path <- paste0(ENV_VAR['LIPID_HOME'], '/results/gsm2sra_query')
dir.create(path = gsm2sra_path)
setwd(gsm2sra_runs)

samples_sheet <- read.csv(file = ENV_VAR['SAMPLES_SHEET'], na.strings = '')
accessions <- samples_sheet$Sample_Accession_GEO
accessions <- accessions[!is.na(accessions)]

# run
lapply(
  X = accessions,
  FUN = run_esearch
)
summarize_runs()