#!/bin/bash
#BSUB -J LipidModule_download-data
#BSUB -P lemmon
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -W 24:00
#BSUB -q general
#BSUB -n 8
#BSUB -B
#BSUB -N
#BSUB -u jsc228@miami.edu

# Note: Job requests 8 cores for 8 GNU parallel curl jobs (see below)

####################################################
# Set directory variables and conda env
####################################################

export NETHOME=/nethome/jsc228/LipidModule
export SCRATCH=/scratch/projects/lemmon/jsc228/LipidModule

source "${HOME}/miniconda3/etc/profile.d/conda.sh"
conda activate fastq-download

cd $NETHOME
find . -type f \( -iname "*.sh" -o -iname "*.r" \) -exec chmod +x {} \;


####################################################
# Sync NETHOME and SCRATCH
####################################################

rsync -rpl --exclude "*.out" --exclude "*.err" --exclude ".nextflow.log*" \
  ${NETHOME}/ ${SCRATCH}


####################################################
# Set fastq-download conda env and env variables
####################################################

cd $SCRATCH/data/rawdata
cat $SCRATCH/bin/sra_explorer_fastq_download.sh | parallel -j 8 "{}"


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
# Script End
####################################################

conda deactivate