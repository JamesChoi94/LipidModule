# System Environment variables 
ENV_VAR <- Sys.getenv(names = TRUE)
gsm2sra_path <- paste0(ENV_VAR['LIPID_HOME'], '/results/gsm2sra_query')
dir.create(path = gsm2sra_path)
gsm2sra_runs <- paste0(gsm2sra_path, '/runs')
dir.create(path = gsm2sra_runs)
setwd(gsm2sra_runs)

samples_sheet <- read.csv(file = ENV_VAR['SAMPLES_SHEET'], na.strings = '')
accessions <- samples_sheet$Sample_Accession_GEO
accessions <- accessions[!is.na(accessions)]


# query functions
run_esearch <- function(gsm) {
    run_out <- paste(gsm, '_esearch-runinfo.txt')
    command_str <- paste0('esearch -db sra -query ', gsm, '| efetch -format runinfo > ', gsm, '_esearch-runinfo.txt')
    system(command = command_str)
}
lapply(
    X = accessions,
    FUN = run_esearch
)
stop()

dir.create('query_runs')

for (i in 1:length(accessions)) {
    run_out <- paste0(results_dir, '/', accessions[i], '_esearch-runinfo.txt')
    print(run_out)
    # run_command <- paste('esearch -db srq -query', accessions[i], '| efetch -format runinfo >', )
    # system(command = 'echo $LIPID_HOME')
}

# system(
#   command = 'esearch -db sra -query "$(awk -F , '{if (length($2) > 0 && $2 ~ /GSM/) print $2}' samples_sheet.csv)" | efetch -format runinfo > /mnt/d/MiamiProject/Lipid_module/esearch_sra_runinfo.txt'
# )