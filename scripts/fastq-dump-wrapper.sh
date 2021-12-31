SRR_ACCESSION=($(awk 'NR!=1 {print $1}' results/gsm2sra_query/gsm2sra_query_compiled.tsv))
echo 'Retriving .sra files for the following accessions:' ${SRR_ACCESSION[@]}
prefetch -p -O data/fastq/sra "${SRR_ACCESSION[@]}"
