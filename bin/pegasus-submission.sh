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

export NETHOME=/nethome/jsc228/LipidModule
export LIPID_HOME=/scratch/projects/lemmon/jsc228
cd ${NETHOME}
source setup_environment.sh
cd ${NETHOME}
cd ..
scp -r LipidModule ${LIPID_HOME}
cd ${LIPID_HOME}/LipidModule

nextflow run rnaseq-processing.nf \
  -profile conda,lsf \
  -resume