#!/bin/sh
# export LIPID_HOME=/mnt/d/MiamiProject/Lipid_Module
export LIPID_HOME=$(pwd)
echo $LIPID_HOME
export SAMPLES_SHEET=${LIPID_HOME}/data/samples_sheet.csv
echo $SAMPLES_SHEET
export SRA_NCBI_CONFIG=${LIPID_HOME}/config/sratoolkit-vdb-config_ec2-setup.tar.gz
echo $SRA_NCBI_CONFIG

# Export executables to PATH if last environment setup was successful.
LOG_CHECK=$(find bin/ -name 'INSTALL_LOG.txt')
if [ "" != "${LOG_CHECK}" ]
then
  LAST_RUN_STATUS=$(tail -n 1 ${LIPID_HOME}/bin/INSTALL_LOG.txt | awk '{print $NF}')
  # echo "${LAST_RUN_LOG}"
  # LAST_RUN_STATUS=$(${LAST_RUN_LOG} | awk '{print $NF}')
  # echo ${LAST_RUN_STATUS}
  if [ "successful" = "${LAST_RUN_STATUS}" ]
  then
    printf "Previous source/binaries installation successful. Exporting executables to PATH\n"
    export PATH=${PATH}:\
${LIPID_HOME}/bin/sratoolkit.2.11.3-ubuntu64/bin:\
${LIPID_HOME}/bin/edirect/:\
${LIPID_HOME}/bin/STAR-2.7.9a/bin/Linux_x86_64:\
${LIPID_HOME}/bin/samtools-1.14:\
${LIPID_HOME}/bin/FastQC:\
${LIPID_HOME}/bin/hisat2-2.2.1:\
${LIPID_HOME}/bin/kallisto
  fi
fi