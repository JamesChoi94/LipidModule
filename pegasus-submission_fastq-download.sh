#!/bin/bash
#BSUB -J LipidModule_fastq-download
#BSUB -P lemmon
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -W 24:00
#BSUB -q general
#BSUB -n 8
#BSUB -B
#BSUB -N
#BSUB -u jsc228@miami.edu

export NETHOME=/nethome/jsc228/LipidModule
export SCRATCH=/scratch/projects/lemmon/jsc228/LipidModule

source "${HOME}/miniconda3/etc/profile.d/conda.sh"
conda activate fastq-download

find . -type f \( -iname "*.sh" -o -iname "*.r" \) -exec chmod +x {} \;
rsync -rpl --exclude "*.out" --exclude "*.err" --exclude ".nextflow.log*" \
  ${NETHOME}/ ${SCRATCH}
cd ${SCRATCH}

cd $SCRATCH/data/rawdata
cat $SCRATCH/bin/sra_explorer_fastq_download.sh | parallel -j 8 "{}"