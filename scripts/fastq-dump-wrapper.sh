# Note: parallel-fastq-dump --sra-id ${SRR_ACCESSION[0]} --threads 2 --outdir data/fastq --splitfiles --gzip
# Above took ~ 20 minutes on my local machine for 26mil spots.
SRR_ACCESSION=($(awk 'NR!=1 {print $1}' results/gsm2sra_query/gsm2sra_query_compiled.tsv))
# echo 'Retriving .sra files for the following accessions:' ${SRR_ACCESSION[@]}
for SRR in ${SRR_ACCESSION[@]:0:2}
do
  echo ${SRR}
  prefetch -p -O data/fastq/ "${SRR}"
  fasterq-dump -O data/fastq/
done




#! /bin/bash
# https://github.com/rnnh/bioinfo-notebook.git

# Help/usage text
usage="$(basename "$0") [options] -a|--annotation <annotation_file> \
	-f|--fasta <fasta_file> -L|--aLigner <aligner_program> <SRR ID(s)> \n
\n
This script downloads FASTQ reads from NCBI's SRA, aligns them to an annotated \n
genome using STAR, and generates gene count table(s) using featureCounts.\n
It can take a single SRR ID as an input, or multiple SRR IDs separated by\n
spaces.\n
\n
Required arguments: \n
\t      -a | --annotation\t     input genome annotation file \n
\t      -f | --fasta\t\t        input FASTA file for annotated genome \n
\t      -L | --aLigner\t\t      program to use for alignment. \n
\t      SRR ID(s)\t\t           Sequence Read Archive Run ID(s) (SRR...) \n
\n
Optional arguments: \n
\t      -h | --help\t\t         show this help text and exit \n
\t      -p | --processors\t	number (n) of processors to use (default: 1) \n
\t      --fastq-dump\t\t        use 'fastq-dump' instead of the 'fasterq-dump'\n
\t      --verbose\t\t           make output of script more verbose\n
\t	--removetemp\t\t	remove read and alignment files once they are\n
\t	\t\t\t  		no longer needed (minimises disk space needed) \n
\t	--log\t\t\t		redirect terminal output to log file
"

# Setting FASTQDUMP to 0
# This will be changed to "1" if --fastq-dump is given as an argument,
# resulting in fastq-dump being used instead of the default fasterq-dump
FASTQDUMP=0

# Setting VERBOSE to 0
# This will be changed to "1" if --verbose is given as an argument,
# resulting in more verbose script output
VERBOSE=0

# Setting REMOVETEMP to 0
# This will be changed to "1" if --removetemp is given as an argument,
# resulting in *.fastq, *.fastq.gz, *.sam, *.bam and *.tsv.summary, being
# removed once they are no longer needed to create a featureCounts table
# REMOVETEMP=0

# Setting LOG to 0
# This will be changed to "1" if --log is given as an argument,
# resulting in the terminal output from this script being redirected to a log
# file
LOG=0

# Setting default number of PROCESSORS to use
PROCESSORS=1

# Creating an empty variable for SRRs to be downloaded and aligned to genome
SRRs=""

# Print usage instructions if script is called without any arguments
if [ "$1" = "" ] ; then
  echo -e "ERROR: please provide input files. \n"
  echo -e $usage
  exit 1
fi

# Iterating through the input arguments with a while loop
while (( "$#" )); do
	case "$1" in
		-h|--help)
			echo -e $usage
			exit
			;;
		-a|--annotation)
			ANNOTATION=$2
			shift 2
			;;
		-f|--fasta)
			FASTA=$2
			shift 2
			;;
		-p|--processors)
			PROCESSORS=$2
			shift 2
			;;
		--fastq-dump)
			FASTQDUMP=1
			shift
			;;
		--verbose)
			VERBOSE=1
			shift
			;;
		--removetemp)
			REMOVETEMP=1
			shift
			;;
		--log)
			LOG=1
			shift
			;;
		--) # end argument parsing
			shift
			break
			;;
		-*|--*) # unsupported flags
			echo -e "ERROR: $1 is an invalid option. \n" >&2
			echo -e $usage
			exit 1
			;;
		*) # preserve SRR ID(s) as positional arguments
			SRRs="$SRRs $1"
			shift
			;;
	esac
done