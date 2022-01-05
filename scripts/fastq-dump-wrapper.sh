#!/bin/bash
# https://github.com/rnnh/bioinfo-notebook.git

# Help/usage text
usage="$(basename "$0") [options] -O|--outdir <path> -e|--threads <count> \n
	-s|--samplesheet <file> <SRR ID(s)> \n
\n
This script downloads FASTQ reads from NCBI's SRA. It can take a single SRR \n
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

# Setting VERBOSE to 0
# This will be changed to "1" if --verbose is given as an argument,
# resulting in more verbose script output
VERBOSE=0

# Setting default number of PROCESSORS to use
THREADS=2

# Creating an empty variable for SRRs to be downloaded.
SRR_ACCESSION=""
SAMPLE_SHEET_PATH=""

# Set default fastq-dump directory
FASTQDUMP_OUT="data"

# Print usage instructions if script is called without any arguments
if [ "$1" = "" ] ; then
  echo -e "ERROR: please provide input files. \n"
  echo -e $usage
  exit 1
fi

# Iterating through the input arguments with a while loop. 
while [[ "$#" ]]
do
	case "$1" in
		-h|--help)
			echo -e $usage
			exit
			;;
		-t|--threads)
			THREADS=$2
			shift 2 # $1="-e" and $2=<count>
			;;
		--verbose)
			VERBOSE=1
			shift 1
			;;
		-s|--samplesheet)
			SAMPLE_SHEET_PATH=$2
			shift 2
			;;
		-O|--outdir)
			FASTQDUMP_OUT=$2
			shift 2
			;;
		--) # end argument parsing
			shift 1
			break 
			;;
		-*|--*) # unsupported flags
			echo -e "ERROR: $1 is an invalid option. \n" >&2
			echo -e $usage
			exit 1
			;;
		SRR*)
			# Preserve SRR ID(s) as positional arguments. Since SRRs are last 
			# arguments, they are concatenated to $SRR_ACCESSION
			SRR_ACCESSION=(${SRR_ACCESSION[@]} $1)
			shift
			;;
		*)
			break
			;;
	esac
done

if [[ "${SAMPLE_SHEET_PATH}" != "" ]]
then
	if [[ -f ${SAMPLE_SHEET_PATH} ]] 
	then
		SRR_ACCESSION=($(awk 'NR!=1 {print $1}' ${SAMPLE_SHEET_PATH}))
	else 
		echo Cannot find ${SAMPLE_SHEET_PATH}. Please specify another relative path.
		echo e.g. results/gsm2sra_query/gsm2sra_query_compiled.tsv
	fi
fi

# Beginning the main body of the script
# The sleep commands ("sleep 1s", "sleep 2s") slow down the script to make
# the output more readable in real-time

echo Starting fastq-dump-wrapper
echo Script started: $(date)

for SRR in ${SRR_ACCESSION[@]}
do
	echo Prefetching sra file for ${SRR}
	until prefetch -p -O ${FASTQDUMP_OUT} "${SRR}"
	do
		echo prefetch fail, retrying in 10 seconds...
		rm ${SRR}.sra.*
		sleep 10s
	done
	until fasterq-dump -p -O ${FASTQDUMP_OUT} -e ${THREADS} \
		--split-files "${SRR}"
	do
		echo fasterq-dump failed, retrying in 10 seconds...
		rm -r ${LIPID_HOME}/fasterq.tmp.*
		sleep 10s
	done
done

echo Script finished: $(date)


# Note: parallel-fastq-dump --sra-id ${SRR_ACCESSION[0]} --threads 2 --outdir data/fastq --splitfiles --gzip
# Above took ~ 20 minutes on my local machine for 26mil spots.

# SRR_ACCESSION=($(awk 'NR!=1 {print $1}' results/gsm2sra_query/gsm2sra_query_compiled.tsv))
# # echo 'Retriving .sra files for the following accessions:' ${SRR_ACCESSION[@]}
# for SRR in ${SRR_ACCESSION[@]:0:2}
# do
#   echo ${SRR}
#   prefetch -p -O data/fastq/ "${SRR}"
#   fasterq-dump -p -O data/fastq/ -e ${THREADS} --split-files "${SRR}"
# done