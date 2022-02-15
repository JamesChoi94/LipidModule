#!/bin/bash
#BSUB -J LipidModule-setup
#BSUB -P lemmon
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -W 1:00
#BSUB -q general
#BSUB -n 12
#BSUB -B
#BSUB -N
#BSUB -u jsc228@miami.edu


####################################################
# Set NETHOME and SCRATCH directory variables
####################################################

export NETHOME=/nethome/jsc228/LipidModule
export SCRATCH=/scratch/projects/lemmon/jsc228/LipidModule

####################################################
# Setup environment necessary for nextflow execution
####################################################

# All commands/dirs should be relative to NETHOME
cd ${NETHOME}
module load java # java/1.8.0_60+ is sufficient

# Make sure script files are executable
find . -type f \( -iname "*.sh" -o -iname "*.r" \) -exec chmod +x {} \;
source setup_unix_environment.sh

####################################################
# Copy NETHOME directory to SCRATCH
####################################################

rsync -rEl --exclude "*.out" --exclude "*.err" --exclude ".nextflow.log*" \
  ${NETHOME} ${SCRATCH}

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
# Run rnaseq-processing 
####################################################

nextflow run rnaseq-processing.nf \
  -profile condaLSF \
  -resume \
  -w ${SCRATCH}/work