SRR_ACCESSION=($(awk 'NR!=1 {print $1}' results/gsm2sra_query/gsm2sra_query_compiled.tsv))
# echo 'Retriving .sra files for the following accessions:' ${SRR_ACCESSION[@]}
for SRR in ${SRR_ACCESSION[@]:0:2}
do
  echo ${SRR}
  prefetch -p -O data/fastq/sra "${SRR}"
  # fasterq-dump
done
