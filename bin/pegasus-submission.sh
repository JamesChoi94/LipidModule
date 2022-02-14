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
export SCRATCH=/scratch/projects/lemmon/jsc228
cd ${NETHOME}
source setup_environment.sh
cd ${NETHOME}
cd ..
if [[ ! -d ${SCRATCH}/LipidModule ]]
then
  scp -r LipidModule ${LIPID_HOME}
fi
cd ${SCRATCH}/LipidModule

module load java/1.8.0_60

export NXF_WORK=${SCRATCH}/LipidModule
nextflow run rnaseq-processing.nf \
  -profile conda,lsf \
  -resume \
  -w ${SCRATCH}/LipidModule