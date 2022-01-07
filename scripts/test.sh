#!/bin/bash
# https://github.com/rnnh/bioinfo-notebook.git

# Help/usage text
usage="$(basename "$0") [options] -O|--outdir <path> -e|--threads <count> \n
	-s|--samplesheet <file> <SRR ID(s)> \n
\n
This script performs read trimming and downloads FASTQ reads from NCBI's SRA. It can take a single SRR \n
ID as an input, multiple SRR IDs separated by spaces, or SRR IDs listed \n
in gsm2sra_query_compiled.tsv. If -s|--samplesheet argument is provided, \n
then query_gsm2sra.R must first be run. \n
\n
Optional arguments: \n
\t      SRR ID(s)\t\t           Sequence Read Archive Run ID(s) (SRR...). \n
\t 			-s | -samplesheet\t			Whether to automatically import SRR IDs from \n
\t 			\t\t\t									a sample query sheet, e.g: \n
\t 			\t\t\t									gsm2sra_query_compiled.tsv \n
\t 			\t\t\t									This will override any SRR IDs provided. \n
\t      -O | --outdir\t\t       Output directory for fastq files. \n
\t      -h | --help\t\t         show this help text and exit \n
\t      -e | --threads\t\t      how many threads to use (dflt=2) \n
\t      --verbose\t\t           make output of script more verbose \n
"


bbduk.sh in=read1.fq in2=read2.fq out1=clean1.fq out2=clean2.fq ref=adapters \
	threads=4 stats=bbduk_stats.txt k=21 hdist=1 ktrim=r mink=10


source ENV.sh
# SRR_ACCESSION=($(awk 'NR!=1 {print $1}' ${SAMPLE_SHEET_PATH}))
