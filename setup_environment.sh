#!/bin/sh

# These should be modified for every new study --------------------------------

# NETHOME and LIPID_HOME should point to the `nethome` and  `scratch` file 
# systems on Pegasus. This will allow larger files (e.g. fastq's) to be saved
# in the `scratch` file system and other, smaller analysis ouputs (e.g. read 
# mappings) to be saved in `nethome` for "permanent" storage. It should be 
# assumed that anything in `scratch` can be deleted after any run.

# If NETHOME and LIPID_HOME point to same directory, the large data files (e.g.
# fastq.gz and genome index/FASTAs) will be saved in the same project 
# directory. This should only really be done for testing purposes. Do NOT 
# point both to the `nethome` file system - your account may be suspended for
# inappropriate resource usage.


export NETHOME=/nethome/jsc228/LipidModule
export LIPID_HOME=/scratch/projects/lemmon/jsc228/LipidModule
# export NETHOME=/mnt/d/LipidModule
# export LIPID_HOME=/mnt/d/LipidModule


if [ ! -d ${NETHOME} ]
then
  echo 'NETHOME cannot be found.'
  return
fi

echo "NETHOME:"
echo ${NETHOME}
echo "SCRATCH:"
echo ${LIPID_HOME}

cd {LIPID_HOME}

bash bin/setup_linux_environment.sh

export PATH=/bin/sratoolkit.2.11.3-ubuntu64/bin:$PATH
echo "PATH:"
echo $PATH

# Download necessary genome FASTAs and GTFs --------------------------

# FASTA
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M28/GRCm39.primary_assembly.genome.fa.gz \
  -P ref/
gunzip ref/GRCm39.primary_assembly.genome.fa.gz

# GTF
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M28/gencode.vM28.annotation.gtf.gz \
  -P ref/
gunzip ref/gencode.vM28.annotation.gtf.gz