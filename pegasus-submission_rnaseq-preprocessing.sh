#!/bin/bash
#BSUB -J LipidModule_rnaseq-preprocessing
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
# Set directory variables and environment
####################################################

export NETHOME=/nethome/jsc228/LipidModule
export SCRATCH=/scratch/projects/lemmon/jsc228/LipidModule

source "${HOME}/miniconda3/etc/profile.d/conda.sh"
conda activate rnaseq-preprocessing

cd $NETHOME
find . -type f \( -iname "*.sh" -o -iname "*.r" \) -exec chmod +x {} \;

module load java/1.8.0_60 # java/1.8.0_60+ is sufficient
unset _JAVA_OPTIONS
module load R/4.1.0 # allows Rscript from command line


####################################################
# Sync NETHOME and SCRATCH
####################################################

rsync -rpl --exclude "*.out" --exclude "*.err" --exclude ".nextflow.log*" \
  ${NETHOME}/ ${SCRATCH}


####################################################
# Setup nextflow temporary work directory
####################################################

mkdir -p ${SCRATCH}/work
export NXF_WORK=${SCRATCH}/work


####################################################
# Run main nextflows
####################################################

cd $SCRATCH
# nextflow run main.nf \
#   -profile lsf \
#   -resume \
#   -w ${SCRATCH}/work