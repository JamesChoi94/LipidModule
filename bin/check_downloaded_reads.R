args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2) { 
  stop("Usage: check_SRR.R <SRR_Accession_to_check>", 
       call.=FALSE)
}
reads_check <- args[1]
readsDir <- args[2]
# reads_check <- c("SRR789190_1.fastq", "SRR789190_2.fastq")
# readsDir <- 'D:/LipidModule/data/rawdata/'
reads_done <- list.files(path = readsDir)
cat(reads_check[!reads_check %in% reads_done])
