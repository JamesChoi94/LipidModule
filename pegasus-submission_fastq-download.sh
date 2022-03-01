#!/bin/bash
#BSUB -J fastq-download
#BSUB -P lemmon
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -W 24:00
#BSUB -q general
#BSUB -n 12
#BSUB -B
#BSUB -N
#BSUB -u jsc228@miami.edu

export NETHOME=/nethome/jsc228/LipidModule
export SCRATCH=/scratch/projects/lemmon/jsc228/LipidModule

find . -type f \( -iname "*.sh" -o -iname "*.r" \) -exec chmod +x {} \;
rsync -rpl --exclude "*.out" --exclude "*.err" --exclude ".nextflow.log*" \
  ${NETHOME}/ ${SCRATCH}
cd ${SCRATCH}

cd $SCRATCH/data/rawdata
source $SCRATCH/bin/sra_explorer_fastq_download.sh