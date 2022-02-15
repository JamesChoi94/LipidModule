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
export SCRATCH=/scratch/projects/lemmon/jsc228/LipidModule

cd ${NETHOME}
source setup_environment.sh
rsync -rEl --exclude "*.out" --exclude "*.err" --exclude ".nextflow.log*" \
  ${NETHOME} ${SCRATCH}
module load java/1.8.0_60
mkdir -p ${SCRATCH}/work
export NXF_WORK=${SCRATCH}/work
cd ${NETHOME}
nextflow run rnaseq-processing.nf \
  -profile conda,lsf \
  -resume \
  -w ${SCRATCH}/work