#!/bin/bash
#BSUB -J fastq-dump-lipidmodule
#BSUB -P lemmon
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -W 24:00
#BSUB -q general
#BSUB -n 12
#BSUB -B
#BSUB -N
#BSUB -u jsc228@miami.edu

export HOME=/scratch/projects/lemmon/jsc228
export LIPID_HOME=/scratch/projects/lemmon/jsc228/LipidModule
cd ${LIPID_HOME}
echo HOME: ${HOME}
echo LIPID_HOME: ${LIPID_HOME}
conda env create --file env/LipidModule.yml
conda activate LipidModule
scripts/fastq-dump-wrapper.sh \
  -s results/gsm2sra_query/gsm2sra_query_compiled.tsv \
  -t %n