#!/bin/bash
#BSUB -J LipidModule-setup
#BSUB -P lemmon
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -W 24:00
#BSUB -q general
#BSUB -n 6
#BSUB -B
#BSUB -N
#BSUB -u jsc228@miami.edu


####################################################
# Set NETHOME and SCRATCH directory variables
####################################################

export NETHOME=/nethome/jsc228/LipidModule
export SCRATCH=/scratch/projects/lemmon/jsc228/LipidModule
export SAMPLESHEET=$NETHOME/data/samplesheet.csv

####################################################
# Setup environment necessary for nextflow execution
####################################################

# All commands/dirs should be relative to NETHOME
cd ${NETHOME}
module load java/1.8.0_60 # java/1.8.0_60+ is sufficient
unset _JAVA_OPTIONS
module load R/4.1.0 # allows Rscript from command line
module load python/3.8.7

# export PATH=$PATH:/${HOME}/miniconda3/envs/LipidModule/bbtools/lib

# Make sure script files are executable 
find . -type f \( -iname "*.sh" -o -iname "*.r" \) -exec chmod +x {} \;
source setup_unix_environment.sh


####################################################
# Copy NETHOME directory to SCRATCH and cd
####################################################

rsync -rpl --exclude "*.out" --exclude "*.err" --exclude ".nextflow.log*" \
  ${NETHOME}/ ${SCRATCH}
cd ${SCRATCH}

####################################################
# Download necessary genome FASTAs and GTFs
####################################################

# Job data should be stored in SCRATCH
# FASTAs 
if [[ ! -f ${SCRATCH}/ref/GRCm39.primary_assembly.genome.fa ]]
then
  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M28/GRCm39.primary_assembly.genome.fa.gz \
    -P ${SCRATCH}/ref/
  gunzip ${SCRATCH}/ref/GRCm39.primary_assembly.genome.fa.gz
fi

# GTF
if [[ ! -f ${SCRATCH}/ref/gencode.vM28.annotation.gtf ]]
then
  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M28/gencode.vM28.annotation.gtf.gz \
    -P ${SCRATCH}/ref/
  gunzip ${SCRATCH}/ref/gencode.vM28.annotation.gtf.gz
fi

####################################################
# Setup nextflow temporary work directory
####################################################

mkdir -p ${SCRATCH}/work
export NXF_WORK=${SCRATCH}/work

####################################################
# Convert GSM to SRR Accessions from samplesheet
####################################################

GSM_ACCESSIONS=($(awk -F, 'NR!=1 {print $1}' ${SAMPLESHEET}))
for GSM in ${GSM_ACCESSIONS[@]}
do
  echo "esearchng ${GSM}..."
  esearch -db sra -query ${GSM} | efetch -format runinfo >> \
    $SCRATCH/data/esearch-runinfo.txt
  sleep 2s
done

# Also create txt file of accessions only
awk -F, 'NR%2==0 {print $1}' data/esearch-runinfo.txt > data/srrAccessions.txt

####################################################
# Download FASTQ files for SRR Accessions
####################################################

# mkdir -p $SCRATCH/data/rawdata
# SRR_ACCESSIONS=($(awk -F, 'NR%2==0 {print $1}' results/esearch-runinfo.txt))
# for SRR in ${SRR_ACCESSIONS[@]}
# do 
#   echo "fasterq-dump ${SRR}"
#   prefetch -p -O $SCRATCH/data/rawdata "${SRR}"
#   sleep 1s 
#   fasterq-dump -p -O $SCRATCH/data/rawdata --threads 6 \
# 	  --split-files "${SRR}"
#   sleep 2s
# done

####################################################
# Run main
####################################################

# nextflow run main.nf \
#   -profile lsf \
#   -resume \
#   -w ${SCRATCH}/work