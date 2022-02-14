args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2) { 
  stop("Usage: compile_queries.R <QueryGEO_path> <samplesheet_path>", 
       call.=FALSE)
}
query_path <- args[1]
samplesheet_path <- args[2]
queries <- readLines(con = query_path)
samplesheet <- read.csv(samplesheet_path)
tmp <- strsplit(unique(queries), ',')
queries <- do.call(rbind.data.frame, args = tmp[2:length(tmp)])
colnames(queries) <- tmp[[1]]
names(queries)[names(queries) == 'SampleName'] <- 'GSM_Accession'
names(queries)[names(queries) == 'Run'] <- 'SRR_Accession'
compile_out <- merge(x = queries,
                     y = samplesheet, 
                     by.y = 'GSM_Accession', 
                     all.x = TRUE)
first_cols <- c('GSM_Accession','SRR_Accession', colnames(samplesheet)[2:ncol(samplesheet)])
else_cols <- colnames(compile_out)[!colnames(compile_out) %in% first_cols]
col_order <- c(first_cols, else_cols)
samplesheet_out <- compile_out[col_order]
# write.csv(x = samplesheet_out,
#           file = gsub(pattern = 'samplesheet.csv',
#                       replacement = 'samplesheet_appended.csv',
#                       x = samplesheet_path),
#           row.names = FALSE)
# print(samplesheet_out)
# cat(samplesheet_out, file = gsub(pattern = 'samplesheet.csv',
#                       replacement = 'samplesheet_appended.csv',
#                       x = samplesheet_path))
write.csv(x = samplesheet_out,
          file = "",
          row.names = FALSE,
          quote = FALSE)