#!/bin/bash

# Help/usage text
usage="$(basename "$0") [options] -O|--outdir <path> -t|--threads <count> \n
	-s|--samplesheet <file> [arguments] <SRR ID(s)> \n
\n
This script performs read mapping and quantification of reads. It performs \n
read mapping using HISAT2 and read quantification using featureCounts \n
software from the Subread package. Input can be one or more FASTQ files \n
(separated by spaces). Alternatively, it can take as input SRR IDs listed in \n
gsm2sra_query_compiled.tsv. Mapping is done against the GRCm38 mouse genome \n
assembly as made available through HISAT2 documentation: \n
http://daehwankimlab.github.io/hisat2/download/#m-musculus \n
\n
Optional arguments: \n
\t			SRR ID(s)			Sequence Read Archive Run ID(s) (SRR...). \n
\t			-s | --samplesheet\t		Whether to automatically import SRR IDs from \n
\t			\t\t\t									a sample query sheet, e.g: \n
\t			\t\t\t									gsm2sra_query_compiled.tsv \n
\t			\t\t\t									This will override any SRR IDs provided. \n
\t			\t\t\t									Matching FASTQ files are first checked in \n
\t			\t\t\t									data/trimmed_reads then data/rawdata. Exits \n
\t			\t\t\t									if not found. \n
\t			-d | --index-dir\t			Path to directory containing HISAT2 index \n
\t			\t\t\t									files. \n
\t			-x | --index-basename\t	Basename of HISAT2 index files. Includes text\n
\t			\t\t\t									up to e.g. *.1.ht2 or *.1.ht21. \n
\t			-u | --unpaired-reads\t	Use flag if reads are unpaired. Default \n
\t			\t\t\t									assumes paired reads. \n
\t			-k | --keep-bam\t				Keep intermediate BAM files. If set, they \n
\t			\t\t\t									are saved to data/aligned_bams. \n
\t			-h | --help\t						Show this help text and exit \n
\t			-t | --threads\t				How many threads to use (dflt=2) \n
\t			--verbose\t							Make output of script more verbose \n
"

# TO DO:
# STRANDEDNESS

# Default location of references
HISAT2_INDEX_DIR="${LIPID_HOME}/ref/HISAT2_index/mm10"
HISAT2_INDEX_BASENAME="genome"
GTF_PATH="${LIPID_HOME}/ref/annotations/gencode.vM10.annotation.gtf"

# Set default output directories
FASTQC_HISAT2_OUT="${LIPID_HOME}/results/fastqc/hisat2"
FASTQC__OUT="${LIPID_HOME}/results/fastqc/featureCounts"
HISAT2_LOG_OUT="${LIPID_HOME}/results/hisat2"
SAMTOOLS_FLAGSTATS_OUT="${LIPID_HOME}/results/samtools/flagstats"

SAMTOOLS_BAM_OUT="${LIPID_HOME}/data/aligned_bams"
SAMTOOLS_BAM_INDEX_OUT="${LIPID_HOME}/data/aligned_bams"

# Setting default number of PROCESSORS to use
THREADS=2

# Creating an empty variable for SRRs to be downloaded.
SRR_ACCESSION=""
SAMPLE_SHEET_PATH=""
UNPAIRED=0
KEEP_BAM=0

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
		-i|--index-dir)
			HISAT2_INDEX_DIR=$2
			shift 2
			;;
		-b|--index-basename)
			HISAT2_INDEX_BASENAME=$2
			shift 2
			;;
		-u|--unpaired-reads)
			UNPAIRED=1
			shift 1
			;;
		-k|--keep-bam)
			KEEP_BAM=1
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

# check sample sheet
# get SRRs - check trimmed - raw - else eit
# check inde
# if no directory, create and download
# if basename no match, eit
# check unpaired reads
# if paired, eit


# Verify that input sample sheet path is valid if provided.
if [[ -f ${SAMPLE_SHEET_PATH} ]]
then
	SRR_ACCESSION=($(awk 'NR!=1 {print $1}' ${SAMPLE_SHEET_PATH}))
elif [[ ${#SRR_ACCESSION[@]} == 0 ]]
then
	echo -e "Cannot find ${SAMPLE_SHEET_PATH}. Please specify another relative file path."
	echo -e "e.g. results/gsm2sra_query/gsm2sra_query_compiled.csv\n"
fi


echo $KEEP_BAM
echo $SAMPLE_SHEET_PATH



# hisat2 indices
# gencode.vM10.annotation.gtf
# sam/bam outputs
# coun

# Run
# get SRR fastqs

# hisat2 -x genome -1 ../../../data/rawdata/SRR789190_1.fastq 
# -2 ../../../data/rawdata/SRR789190_2.fastq -S ../../../data/aligned_bams/SRR789190.sam

# bwa mem genome.fa reads.fastq | samtools sort -o myfile_sorted.bam

# hisat2 (options)... | \
#   tee >(samtools flagstat - > hisat2_output.flagstat) | \
#   samtools sort -O BAM | \
#   tee hisat2_output.bam | \
#   samtools index - hisat2_output.bam.bai