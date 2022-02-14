#!/bin/bash
#BSUB -J LipidModule-setup
#BSUB -P lemmon
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -W 1:00
#BSUB -q general
#BSUB -n 2
#BSUB -B
#BSUB -N
#BSUB -u jsc228@miami.edu

export NETHOME=/nethome/jsc228/LipidModule
export LIPID_HOME=/scratch/projects/lemmon/jsc228/LipidModule
cd /nethome/jsc228/LipidModule
source setup_environment.sh