# TODO: 
# add executables to $PATH

LIPID_HOME=$(pwd)
BINARY_HOME=${LIPID_HOME}/bin

# Purpose: install binaries/source packages needed for LipidModule project into clean Ubuntu 20.04 LTS (e.g. AWS ec2 instance).
# After successful installation of packages, a file called 'INSTALL_LOG' will appear in the /bin folder with info on the last time the environment setup was run.

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
  exit
fi

##########################################################################################

### R installation ###
# Description: R is a free software environment for statistical computing and graphics. 
R_OK=$(dpkg-query -W --showformat='${Status}\n' 'r-base' | grep 'install ok installed')
if [ "" = "${R_OK}" ]
then
  REQUIRED_PKGS=('dirmngr' 'gnupg' 'apt-transport-https' 'ca-certificates' 'software-properties-common')
  # Install the dependencies necessary to add a new repository over HTTPS:
  for PKG in ${REQUIRED_PKGS[*]}
  do
    PKG_OK=$(dpkg-query -W --showformat='${Status}\n' ${PKG} | grep 'install ok installed')
    printf "Checking for ${PKG}: ${PKG_OK}\n"
    if [ "" = "${PKG_OK}" ]
    then
      printf "No ${PKG}. Setting up ${PKG}.\n"
      sudo apt-get -yes install ${PKG}
    fi
  done
  # Add the CRAN repository to your system sources' list:
  sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
  sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
  sudo apt-get -yes install r-base
fi
R_VERSION=$(R --version | grep "version" | head -n 1)
if [ "" = "${R_VERSION}" ]
then
  printf "R installation unsuccessful. Please install R manually. Exiting ec2 environment setup.\n"
  printf "Last run: $(date) Status: unsuccessful\n" >> ${LIPID_HOME}/bin/INSTALL_LOG.txt
  exit
else
  printf "R successfully installed: ${R_VERSION}\n"
fi

##########################################################################################

### FastQC installation (Java application) ###
# Description: FastQC is a quality control tool for high throughput sequence data.
JRE_OK=$(dpkg-query -W --showformat='${Status}\n' 'default-jre' | grep 'install ok installed')
if [ "" = "${JRE_OK}" ]
then
  printf "Java Runtime Environment not found in environment. Installing default-jre.\n"
  sudo apt install default-jre
fi
JAVA_PATH=$(which java)
if [ "" = "${JAVA_PATH}" ]
then
  printf "Java Runtime Environment (JRE) unsuccessful. Please install JRE manually. Exiting ec2 environment setup.\n"
  printf "Last run: $(date) Status: unsuccessful\n" >> ${LIPID_HOME}/bin/INSTALL_LOG.txt
  exit
else
  printf "JRE successfully installed at: ${JAVA_PATH}\n"
fi
wget -P ${LIPID_HOME}/bin/ https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip ${LIPID_HOME}/bin/fastqc_v0.11.0.zip
chmod 755 ${LIPID_HOME}/bin/FastQC/fastqc
FASTQC_VERSION=$(${LIPID_HOME}/bin/FastQC/fastqc --version)
if [ "" = "${FASTQC_VERSION}" ]
then
  printf "FastQC installation unsuccessful. Please install FastQC manually. Exiting ec2 environment setup.\n"
  printf "Last run: $(date) Status: unsuccessful\n" >> ${LIPID_HOME}/bin/INSTALL_LOG.txt
  exit
else
  printf "FastQC successfully installed: ${FASTQC_VERSION}\n"
fi

##########################################################################################

### SRAtoolkit installation ### 
# Description: SRA Toolkit can be used to directly download SRA data files and reference sequences
# Specific installation instructions can be found here:
# https://github.com/ncbi/sra-tools/wiki
wget -P ${LIPID_HOME}/bin/ https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.3/sratoolkit.2.11.3-ubuntu64.tar.gz
tar -vxzf ${LIPID_HOME}/bin/sratoolkit.2.11.3-ubuntu64.tar.gz
FASTQDUMP_VERSION=$(${LIPID_HOME}/bin/sratoolkit.2.11.3-ubuntu64/bin/fastq-dump --version)
# Import SRA tools config, which can be prepared on a separate machine beforehand e.g. on my local machine:
# tar -czvf sratoolkit-vdb-config_ec2-setup.tar.gz .ncbi/
cp ${LIPID_HOME}/config/sratoolkit-vdb-config_ec2-setup.tar.gz ${HOME}
tar -xvzf sratoolkit-vdb-config_ec2-setup.tar.gz
# ${LIPID_HOME}/bin/sratoolkit.2.11.3-ubuntu64/bin/vdb-config -Q yes
VDB_CONFIG_VERSION = $(${LIPID_HOME}/bin/sratoolkit.2.11.3-ubuntu64/bin/vdb-config --version)
if [ "" = "${FASTQDUMP_VERSION}" ] || [ "" = "${VDB_CONFIG_VERSION}" ]
then
  printf "SRA toolkit installation unsuccessful. Please install manually. Exiting ec2 environment setup.\n"
  printf "Last run: $(date) Status: unsuccessful\n" >> ${LIPID_HOME}/bin/INSTALL_LOG.txt
  exit
else
  printf "FastQC successfully installed: ${FASTQC_VERSION}\n"
fi
printf "SRA toolkit successfully installed:\nfastq-dump version: ${FASTQDUMP_VERSION}\nvdb-config version: ${VDB_CONFIG_VERSION}\n"

# Export to PATH for convenience and to show path to binaries.
# export PATH=$PATH:${LIPID_HOME}/sratoolkit.2.11.3-ubuntu64/bin/

##########################################################################################

### STAR installation ###
# Description: STAR is a splice-aware aligner for aligning reads to a referece genome.
wget -P ${LIPID_HOME}/bin/ https://github.com/alexdobin/STAR/archive/2.7.9a.tar.gz
tar -xzf ${LIPID_HOME}/bin/2.7.9a.tar.gz
STAR_VERSION=$(${LIPID_HOME}/bin/STAR-2.7.9a/bin/Linux_x86_64/STAR --version)
if [ "" = "${STAR_VERSION}" ]
then
  printf "STAR installation unsuccessful. Please install manually. Exiting ec2 environment setup.\n"
  printf "Last run: $(date) Status: unsuccessful\n" >> ${LIPID_HOME}/bin/INSTALL_LOG.txt
  exit
else
  printf "STAR successfully installed: ${STAR_VERSION}\n"
fi


##########################################################################################

### HISAT2 installation ###
# Description: HISAT2 is a fast and sensitive alignment program for mapping next-generation sequencing reads (both DNA and RNA) to a population of human genomes as well as to a single reference genome.
curl -s https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download > ${LIPID_HOME}/bin/hisat2-2.2.1-Linux_x86_64.zip
unzip ${LIPID_HOME}/bin/hisat2-2.2.1-Linux_x86_64.zip
HISAT2_VERSION=$(${LIPID_HOME}/bin/hisat2-2.2.1/hisat2 --version | head -n 1)
if [ "" = "${HISAT2_VERSION}" ]
then
  printf "HISAT2 installation unsuccessful. Please install manually. Exiting ec2 environment setup.\n"
  printf "Last run: $(date) Status: unsuccessful\n" >> ${LIPID_HOME}/bin/INSTALL_LOG.txt
  exit
else
  printf "HISAT2 successfully installed: ${HISAT2_VERSION}\n"
fi

##########################################################################################

### Kallisto installation ###
# Description: kallisto is a program for quantifying abundances of transcripts from bulk and single-cell RNA-Seq data, or more generally of target sequences using high-throughput sequencing reads.
wget -P ${LIPID_HOME}/bin/ https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_linux-v0.46.1.tar.gz
gunzip ${LIPID_HOME}/bin/kallisto_linux-v0.46.1.tar.gz
tar -xvf ${LIPID_HOME}/bin/kallisto_linux-v0.46.1.tar
KALLISTO_VERSION=$(${LIPID_HOME}/bin/kallisto/kallisto version)
if [ "" = "${KALLISTO_VERSION}" ]
then
  printf "kallisto installation unsuccessful. Please install manually. Exiting ec2 environment setup.\n"
  printf "Last run: $(date) Status: unsuccessful\n" >> ${LIPID_HOME}/bin/INSTALL_LOG.txt
  exit
else
  printf "kallisto successfully installed: ${KALLISTO_VERSION}\n"
fi

##########################################################################################

### samtools installation ###
# Description: Samtools is a suite of programs for interacting with high-throughput sequencing data.
wget -P ${LIPID_HOME}/bin/ https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2
bunzip2 ${LIPID_HOME}/bin/samtools-1.14.tar.bz2
tar -xvf ${LIPID_HOME}/bin/samtools-1.14.tar
cd ${LIPID_HOME}/bin/samtools-1.14/
./configure
make
make install
cd ${LIPID_HOME}
${LIPID_HOME}/bin/samtools-1.14/samtools --version

##########################################################################################

echo "Last run: $(date) Status: successful\n" >> ${LIPID_HOME}/bin/INSTALL_LOG.txt
exit