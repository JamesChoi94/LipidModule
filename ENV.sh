#!/bin/sh
# export LIPID_HOME=/mnt/d/MiamiProject/Lipid_Module
export LIPID_HOME=$(pwd)
echo $LIPID_HOME
export SAMPLES_SHEET=${LIPID_HOME}/data/samples_sheet.tsv
echo $SAMPLES_SHEET
export SRA_NCBI_CONFIG=${LIPID_HOME}/env/sratoolkit-vdb-config_ec2-setup.tar.gz
echo $SRA_NCBI_CONFIG