
# Purpose: install binaries/source packages needed for LipidModule project into clean Ubuntu 20.04 LTS (e.g. AWS ec2 instance).
# After successful installation of packages, a file called 'INSTALL_LOG' will appear in the /bin folder with info on the last time the environment setup was run.

export LIPID_HOME=$(pwd)
echo $LIPID_HOME
export SRA_NCBI_CONFIG=${LIPID_HOME}/config/sratoolkit-vdb-config_ec2-setup.tar.gz
echo $SRA_NCBI_CONFIG

# Check to make sure this environment script wasn't already run.
LOG_CHECK=$(find bin/ -name 'INSTALL_LOG.txt')
if [ "" = "${LOG_CHECK}" ]
then
  touch ${LIPID_HOME}/bin/INSTALL_LOG.txt
fi
LAST_RUN_LOG=$(tail -n 1 ${LIPID_HOME}/bin/INSTALL_LOG.txt)
LAST_RUN_STATUS=$(${LAST_RUN_LOG} | awk '{print $NF}')
if [ "successful" = "${LAST_RUN_STATUS}" ]
then
  printf "${LAST_RUN_LOG}\n"
fi

# Need java for FastQC. Can load java module available in pegasus.
module load java

##########################################################################################

### R installation ###
# Description: R is a free software environment for statistical computing and graphics. 
R_OK=$(which R)
if [ "" = "${R_OK}"]
then
  echo "R is not installed. Contact pegasus.ccs.miami.edu sysadmin."
else
  R_VERSION=$(R --version | grep "version" | head -n 1)
fi

##########################################################################################

### FastQC installation (Java application) ###
# Description: FastQC is a quality control tool for high throughput sequence data.
wget -P ${LIPID_HOME}/bin/ --no-check-certificate https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip ${LIPID_HOME}/bin/fastqc_v0.11.9.zip -d ${LIPID_HOME}/bin/
chmod 755 ${LIPID_HOME}/bin/FastQC/fastqc
FASTQC_VERSION=$(${LIPID_HOME}/bin/FastQC/fastqc --version)
if [ "" = "${FASTQC_VERSION}" ]
then
  printf "FastQC installation unsuccessful. Please install FastQC manually.\n"
  printf "Last run: $(date) Status: unsuccessful\n" >> ${LIPID_HOME}/bin/INSTALL_LOG.txt
else
  printf "FastQC successfully installed: ${FASTQC_VERSION}\n\n"
fi

##########################################################################################

### SRAtoolkit installation ### 
# Description: SRA Toolkit can be used to directly download SRA data files and reference sequences
# Specific installation instructions can be found here:
# https://github.com/ncbi/sra-tools/wiki
wget -P ${LIPID_HOME}/bin/ https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.3/sratoolkit.2.11.3-ubuntu64.tar.gz  -d ${LIPID_HOME}/bin/
tar -vxzf ${LIPID_HOME}/bin/sratoolkit.2.11.3-ubuntu64.tar.gz -C ${LIPID_HOME}/bin/
FASTQDUMP_VERSION=$(${LIPID_HOME}/bin/sratoolkit.2.11.3-ubuntu64/bin/fastq-dump --version)
# Import SRA tools config, which can be prepared on a separate machine beforehand e.g. on my local machine:
# tar -czvf sratoolkit-vdb-config_ec2-setup.tar.gz .ncbi/
tar -xvzf ${SRA_NCBI_CONFIG} -C ${HOME}
# ${LIPID_HOME}/bin/sratoolkit.2.11.3-ubuntu64/bin/vdb-config -Q yes
VDB_CONFIG_VERSION=$(${LIPID_HOME}/bin/sratoolkit.2.11.3-ubuntu64/bin/vdb-config --version)
if [ "" = "${FASTQDUMP_VERSION}" ] || [ "" = "${VDB_CONFIG_VERSION}" ]
then
  printf "SRA toolkit installation unsuccessful. Please install manually.\n"
  printf "Last run: $(date) Status: unsuccessful\n" >> ${LIPID_HOME}/bin/INSTALL_LOG.txt
else
  printf "SRA toolkit successfully installed:\nfastq-dump version: ${FASTQDUMP_VERSION}\nvdb-config version: ${VDB_CONFIG_VERSION}\n\n"
fi

##########################################################################################

### STAR installation ###
# Description: STAR is a splice-aware aligner for aligning reads to a referece genome.
wget -P ${LIPID_HOME}/bin/ https://github.com/alexdobin/STAR/archive/2.7.9a.tar.gz  -d ${LIPID_HOME}/bin/
tar -xzf ${LIPID_HOME}/bin/2.7.9a.tar.gz -C ${LIPID_HOME}/bin/
STAR_VERSION=$(${LIPID_HOME}/bin/STAR-2.7.9a/bin/Linux_x86_64/STAR --version)
if [ "" = "${STAR_VERSION}" ]
then
  printf "STAR installation unsuccessful. Please install manually.\n"
  printf "Last run: $(date) Status: unsuccessful\n" >> ${LIPID_HOME}/bin/INSTALL_LOG.txt
else
  printf "STAR successfully installed: ${STAR_VERSION}\n\n"
fi

##########################################################################################

### HISAT2 installation ###
# Description: HISAT2 is a fast and sensitive alignment program for mapping next-generation sequencing reads (both DNA and RNA) to a population of human genomes as well as to a single reference genome.
curl -s https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download > ${LIPID_HOME}/bin/hisat2-2.2.1-Linux_x86_64.zip
unzip ${LIPID_HOME}/bin/hisat2-2.2.1-Linux_x86_64.zip -d ${LIPID_HOME}/bin/
HISAT2_VERSION=$(${LIPID_HOME}/bin/hisat2-2.2.1/hisat2 --version | head -n 1)
if [ "" = "${HISAT2_VERSION}" ]
then
  printf "HISAT2 installation unsuccessful. Please install manually.\n"
  printf "Last run: $(date) Status: unsuccessful\n" >> ${LIPID_HOME}/bin/INSTALL_LOG.txt
else
  printf "HISAT2 successfully installed: ${HISAT2_VERSION}\n\n"
fi

##########################################################################################

### Kallisto installation ###
# Description: kallisto is a program for quantifying abundances of transcripts from bulk and single-cell RNA-Seq data, or more generally of target sequences using high-throughput sequencing reads.
wget -P ${LIPID_HOME}/bin/ https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_linux-v0.46.1.tar.gz  -d ${LIPID_HOME}/bin/
gunzip ${LIPID_HOME}/bin/kallisto_linux-v0.46.1.tar.gz -d ${LIPID_HOME}/bin/
tar -xvf ${LIPID_HOME}/bin/kallisto_linux-v0.46.1.tar -C ${LIPID_HOME}/bin/
KALLISTO_VERSION=$(${LIPID_HOME}/bin/kallisto/kallisto version)
if [ "" = "${KALLISTO_VERSION}" ]
then
  printf "kallisto installation unsuccessful. Please install manually.\n"
  printf "Last run: $(date) Status: unsuccessful\n" >> ${LIPID_HOME}/bin/INSTALL_LOG.txt
else
  printf "kallisto successfully installed: ${KALLISTO_VERSION}\n\n"
fi

##########################################################################################

### samtools installation ###
# Description: Samtools is a suite of programs for interacting with high-throughput sequencing data.
wget -P ${LIPID_HOME}/bin/ https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2  -d ${LIPID_HOME}/bin/
bunzip2 ${LIPID_HOME}/bin/samtools-1.14.tar.bz2
tar -xvf ${LIPID_HOME}/bin/samtools-1.14.tar -C ${LIPID_HOME}/bin/
cd ${LIPID_HOME}/bin/samtools-1.14/
./configure
make
make install
cd ${LIPID_HOME}
${LIPID_HOME}/bin/samtools-1.14/samtools --version

##########################################################################################

### Entrez Direct tools installation
# Description: Entrez Direct (EDirect) provides access to the NCBI's suite of interconnected databases (publication, sequence, structure, gene, variation, expression, etc.) from a Unix terminal window
sh -c "$(wget -q ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -P ${LIPID_HOME}/bin/)"
mv ${HOME}/edirect ${LIPID_HOME}/bin/

##########################################################################################

printf "R successfully installed: ${R_VERSION}\n\n"
printf "JRE successfully installed at: ${JAVA_PATH}\n"
printf "FastQC successfully installed: ${FASTQC_VERSION}\n"
printf "SRA toolkit successfully installed:\nfastq-dump version: ${FASTQDUMP_VERSION}\nvdb-config version: ${VDB_CONFIG_VERSION}\n"
printf "STAR successfully installed: ${STAR_VERSION}\n"
printf "HISAT2 successfully installed: ${HISAT2_VERSION}\n"
printf "kallisto successfully installed: ${KALLISTO_VERSION}\n"

printf "To export executables to PATH, please source ENV.sh\n"

# export PATH=${PATH}:\
# ${LIPID_HOME}/bin/sratoolkit.2.11.3-ubuntu64/bin:\
# ${LIPID_HOME}/bin/edirect/:${LIPID_HOME}/bin:\
# ${LIPID_HOME}/bin/STAR-2.7.9a/bin/Linux_x86_64:\
# ${LIPID_HOME}/bin/samtools-1.14:\
# ${LIPID_HOME}/bin/FastQC:\
# ${LIPID_HOME}/bin/hisat2-2.2.1:\
# ${LIPID_HOME}/bin/kallisto

printf "Last run: $(date) Status: successful\n" >> ${LIPID_HOME}/bin/INSTALL_LOG.txt