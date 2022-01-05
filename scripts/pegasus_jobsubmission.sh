#!/bin/bash
#BSUB -J fastq-dump-lipidmodule
#BSUB -P lemmon
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -W 1:00
#BSUB -q general
#BSUB -n 1
#BSUB -R "rusage[mem=128]"
#BSUB -B
#BSUB -N
#BSUB -u jsc228@miami.edu
#
# Run serial executable on 1 cpu of one node


export LIPID_NETHOME=/nethome/jsc228/LipidModule/
export LIPID_SCRATCH=/scratch/projects/lemmon/jsc228/LipidModule/

# 1. Setup environment.

# 2. Download SRA

# 3. 
