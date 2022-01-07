#!/bin/bash
# https://github.com/rnnh/bioinfo-notebook.git

# Help/usage text
usage="$(basename "$0") [options] -O|--outdir <path> -t|--threads <count> \n
	-s|--samplesheet <file> <SRR ID(s)> \n
\n
This script downloads FASTQ reads from NCBI's SRA and performs first-pass \n
quality control using FastQC. It can take a single SRR ID as an input, \n
multiple SRR IDs separated by spaces, or SRR IDs listed in \n
gsm2sra_query_compiled.tsv. If -s|--samplesheet argument is provided, \n
then query_gsm2sra.R must first be run. \n
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

# Adding sra-tools to PATH
export PATH=${PATH}:${LIPID_HOME}/env/sratoolkit.2.11.3-ubuntu64/bin/

# Set default output directories
# For fasterq-dump, this might be redundant since sra-tools already allows 
# specifying location of user-respository.
FASTQDUMP_OUT="${LIPID_HOME}/data/rawdata"
FASTQC_RAW_OUT="${LIPID_HOME}/results/fastqc/rawdata"
BBDUK_READS_OUT="${LIPID_HOME}/data/trimmed_reads"
BBDUK_STATS_OUT="${LIPID_HOME}/results/bbduk"
FASTQC_TRIMMED_OUT="${LIPID_HOME}/results/fastqc/trimmed_reads"

# Setting VERBOSE to 0
# This will be changed to "1" if --verbose is given as an argument,
# resulting in more verbose script output
VERBOSE=0

# Setting default number of PROCESSORS to use
THREADS=2

# Creating an empty variable for SRRs to be downloaded.
SRR_ACCESSION=""
SAMPLE_SHEET_PATH=""

# Print usage instructions if script is called without any arguments
if [[ "$1" = "" ]] ; then
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

# Verify that input sample sheet path is valid.
if [[ -f ${SAMPLE_SHEET_PATH} ]]
then
	SRR_ACCESSION=($(awk 'NR!=1 {print $1}' ${SAMPLE_SHEET_PATH}))
elif [[ ${#SRR_ACCESSION[@]} == 0 ]]
then
	echo -e "Cannot find ${SAMPLE_SHEET_PATH}. Please specify another relative file path."
	echo -e "e.g. results/gsm2sra_query/gsm2sra_query_compiled.csv\n"
fi

# %%%%%%%%%%%%% Main body of the script %%%%%%%%%%%%%

# The sleep commands ("sleep 1s", "sleep 2s") slow down the script to make
# the output more readable in real-time

echo -e "~~~~~~~~ FASTQ-DUMP-WRAPPER ~~~~~~~~\n"
echo -e "Starting fastq-dump-wrapper"
echo -e "Script started: $(date)\n"

# Create output directories if cannot be found
if [[ ! -d ${FASTQDUMP_OUT} ]]
then
	mkdir -p  ${FASTQDUMP_OUT}
fi
if [[ ! -d ${FASTQC_RAW_OUT} ]]
then
	mkdir -p ${FASTQC_RAW_OUT}
fi
if [[ ! -d ${BBDUK_READS_OUT} ]]
then
	mkdir -p  ${BBDUK_READS_OUT}
fi
if [[ ! -d ${BBDUK_STATS_OUT} ]]
then
	mkdir -p  ${BBDUK_STATS_OUT}
fi
if [[ ! -d ${FASTQC_TRIMMED_OUT} ]]
then
	mkdir -p ${FASTQC_TRIMMED_OUT}
fi

# the work
for SRR in ${SRR_ACCESSION[@]}
do
	echo -e "Prefetching sra file for ${SRR}"
	until prefetch -p -O ${FASTQDUMP_OUT} "${SRR}"
	do
		echo -e "Prefetch fail, retrying in 10 seconds..."
		rm ${SRR}.sra.*
		sleep 10s
	done
	sleep 3s
	# skip fasterq-dump if fastq files exist
	# if [[ "(ls ${FASTQDUMP_OUT}/${SRR}*.fastq)" ]]
	if ls -U ${FASTQDUMP_OUT}/${SRR}*.fastq 1> /dev/null 2>&1
	then
		echo -e "\nfastq files for ${SRR} already in output directory. Skipping fastq-dump...\n"
	else 
		until fasterq-dump -p -O ${FASTQDUMP_OUT} -e ${THREADS} \
			--split-files "${SRR}"
		do
			echo -e "fasterq-dump failed, retrying in 10 seconds..."
			rm -r ${LIPID_HOME}/fasterq.tmp.*
			sleep 10s
		done
	fi
	sleep 3s
done

# Run FastQC on raw reads
echo -e "Running FastQC on raw reads..."
FASTQC_DONE=($(ls ${FASTQC_RAW_OUT} | cut -d '_' -f 1-2 | uniq))
ALL_FASTQ=($(echo $(ls ${FASTQDUMP_OUT}/*.fastq | rev | cut -d '/' -f 1 | rev |\
	cut -d '.' -f 1)))
FASTQC_DO=($(echo ${ALL_FASTQ[@]} ${FASTQC_DONE[@]} | tr ' ' '\n' | sort | \
	uniq -u | sort -r))
FASTQC_DO=( "${FASTQC_DO[@]/#/${FASTQDUMP_OUT}/}" )
FASTQC_DO=( "${FASTQC_DO[@]/%/.fastq}" )
fastqc -o ${FASTQC_RAW_OUT} -t ${THREADS} ${FASTQC_DO[@]}
echo -e "Done running FastQC on raw reads."

# Perform adapter trimming
for SRR in ${SRR_ACCESSION[@]}
do
	READS=($(ls ${FASTQDUMP_OUT}/${SRR}*.fastq))
	echo ${READS[@]}
	if [[ "${#READS[@]}" -ge 2 ]]
	then
		bbduk.sh \
			in=${FASTQDUMP_OUT}/${SRR}_1.fastq \
			in2=${FASTQDUMP_OUT}/${SRR}_2.fastq \
			out=${BBDUK_READS_OUT}/${SRR}_1_trimmed.fastq \
			out2=${BBDUK_READS_OUT}/${SRR}_2_trimmed.fastq \
			stats=${BBDUK_STATS_OUT}/${SRR}_BBduk-stats.txt \
			threads=${THREADS} ref=adapters k=21 hdist=1 ktrim=r mink=10
	else
		bbduk.sh in=${FASTQDUMP_OUT}/${SRR}.fastq \
			out=${BBDUK_READS_OUT}/${SRR}_trimmed.fastq \
			stats=${BBDUK_STATS_OUT}/${SRR}_BBduk-stats.txt \
			threads=${THREADS} ref=adapters k=21 hdist=1 ktrim=r mink=10
	fi
done

# Run FastQC on trimmed reads
echo -e "Running FastQC on trimmed reads..."
FASTQC_DONE=($(ls ${FASTQC_TRIMMED_OUT} | cut -d '_' -f 1-2 | uniq))
ALL_FASTQ=($(echo $(ls ${FASTQDUMP_OUT}/*.fastq | rev | cut -d '/' -f 1 | rev |\
	cut -d '.' -f 1)))
FASTQC_DO=($(echo ${ALL_FASTQ[@]} ${FASTQC_DONE[@]} | tr ' ' '\n' | sort | \
	uniq -u | sort -r))
FASTQC_DO=( "${FASTQC_DO[@]/#/${FASTQDUMP_OUT}/}" )
FASTQC_DO=( "${FASTQC_DO[@]/%/.fastq}" )
fastqc -o ${FASTQC_TRIMMED_OUT} -t ${THREADS} ${FASTQC_DO[@]}
echo -e "Done running FastQC on trimmed reads."


echo Script finished: $(date)

# Note: parallel-fastq-dump --sra-id ${SRR_ACCESSION[0]} --threads 2 --outdir data/fastq --splitfiles --gzip
# Above took ~ 20 minutes on my local machine for 26mil spots.