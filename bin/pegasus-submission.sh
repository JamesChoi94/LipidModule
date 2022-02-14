#!/bin/bash
#BSUB -J 
#BSUB -P lemmon
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -W 1:00
#BSUB -q general
#BSUB -n 8
#BSUB -B
#BSUB -N
#BSUB -u jsc228@miami.edu

cd /nethome/jsc228
wget -qO- https://get.nextflow.io | bash
chmod +x nextflow