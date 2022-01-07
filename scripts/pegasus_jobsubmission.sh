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

export LIPID_HOME=/scratch/projects/lemmon/jsc228/LipidModule/
conda activate LipidModule
scripts/fastq-dump-wrapper.sh \
  -s results/gsm2sra_query/gsm2sra_query_compiled.tsv \
  -t %n