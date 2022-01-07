#!/bin/bash
# https://github.com/rnnh/bioinfo-notebook.git

# Help/usage text
usage="$(basename "$0") [options] -O|--outdir <path> -t|--threads <count> \n
	-s|--samplesheet <file> <SRR ID(s)> \n
\n
This script 

\n
Optional arguments: \n
\t      SRR ID(s)\t\t           Sequence Read Archive Run ID(s) (SRR...). \n
\t 			-s | -samplesheet\t			Whether to automatically import SRR IDs from \n
\t 			\t\t\t									a sample query sheet, e.g: \n
\t 			\t\t\t									gsm2sra_query_compiled.tsv \n
\t 			\t\t\t									This will override any SRR IDs provided. \n
\t      -h | --help\t\t         show this help text and exit \n
\t      -t | --threads\t\t      how many threads to use (dflt=2) \n
\t      --verbose\t\t           make output of script more verbose \n
"
