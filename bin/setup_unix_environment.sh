#!/bin/sh

####################################################
# Description
####################################################

# This script install miniconda3 and creates two conda environments.
# `fastq-download` is used to download fastq sequencing files.
# `rnaseq-preprocessing` is used to perform read QC and alignment.
# Two environments are required because of dependency conflicts
# between `sratools` in `fastq-download` and other packages in 
# `rnaseq-preprocessing`. This script should be run once (in theory).

echo "Starting setup_unix_environment.sh in $(pwd)"


####################################################
# Install Nextflow
####################################################

# Installs nextflow into a bin/ folder that can be found via $PATH

if [[ ! -f ${HOME}/bin/nextflow ]]
then
  if [[ ! -f ${HOME}/nextflow ]]
  then
    cd ${HOME}
    mkdir bin
    wget -qO- https://get.nextflow.io | bash
    chmod +x nextflow
  fi
  mv nextflow bin/
fi


####################################################
# Install Miniconda3
####################################################

echo Checking if conda is executable from PATH...
if [ "" == "$(which conda)" ]
then
  echo conda not found. Downloading Miniconda3 installation script...
  sleep 1s
  wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    -O ${HOME}/miniconda.sh

  echo Installing Miniconda3...; sleep 1s
  bash ${HOME}/miniconda.sh -b -p ${HOME}/miniconda3

  echo Miniconda3 installed, removing installation script...
  rm -f ${HOME}/miniconda.sh

  echo Setting up Miniconda3...; sleep 1s
  source "${HOME}/miniconda3/etc/profile.d/conda.sh"
  hash -r
  conda config \
    --set always_yes yes \
    --set changeps1 yes \
    --set auto_activate_base false
  conda update -q conda
  conda init

  echo Displaying information about current conda installation...; sleep 1s
  conda info -a
else
  echo conda found at: $(which conda)
  conda info -a
fi

####################################################
# Clean unused conda env packages
####################################################

echo Removing unused packages and caches using conda...
sleep 1s
conda clean --all --yes
conda config --add channels bioconda
conda config --add channels conda-forge


####################################################
# Create fastq-download conda env
####################################################

echo Installing fastq-download conda env...
sleep 1s
if [[ -d ${HOME}/miniconda3/envs/fastq-download ]]
then
  INSTALL_ENV=0
  echo fastq-download conda env already exists.; sleep 1s
  echo Checking if installed packages are current...; sleep 1s
  conda list --name fastq-download --explicit > config/tmp_env.txt
  ENV_DIFF=$(diff config/tmp_env.txt config/fastq-download.txt | wc -l)
  rm config/tmp_env.txt
  if [ "${ENV_DIFF}" -ge 1 ]
  then
    echo "Conda diff length: ${ENV_DIFF}"
    conda env update --name fastq-download --file config/fastq-download.txt --prune
  else
    echo fastq-download conda env is up to date. Exiting script.
    return 0
  fi
else
  INSTALL_ENV=1
  echo fastq-download conda env does not exist.; sleep 1s
fi
if [[ $INSTALL_ENV == 1 ]]
then
  echo Creating the fastq-download env...; sleep 1s
  conda create --name fastq-download --file config/fastq-download.txt
fi


####################################################
# Create rnaseq-preprocessing conda env
####################################################

echo Installing rnaseq-preprocessing conda env...
sleep 1s
if [[ -d ${HOME}/miniconda3/envs/rnaseq-preprocessing ]]
then
  INSTALL_ENV=0
  echo rnaseq-preprocessing conda env already exists.; sleep 1s
  echo Checking if installed packages are current...; sleep 1s
  conda list --name rnaseq-preprocessing --explicit > config/tmp_env.txt
  ENV_DIFF=$(diff config/tmp_env.txt config/rnaseq-preprocessing.txt | wc -l)
  rm config/tmp_env.txt
  if [ "${ENV_DIFF}" -ge 1 ]
  then
    echo "Conda diff length: ${ENV_DIFF}"
    conda env update --name rnaseq-preprocessing --file config/rnaseq-preprocessing.txt --prune
  else
    echo rnaseq-preprocessing conda env is up to date. Exiting script.
    return 0
  fi
else
  INSTALL_ENV=1
  echo rnaseq-preprocessing conda env does not exist.; sleep 1s
fi
if [[ $INSTALL_ENV == 1 ]]
then
  echo Creating the rnaseq-preprocessing env...; sleep 1s
  conda create --name rnaseq-preprocessing --file config/rnaseq-preprocessing.txt
fi


####################################################
# Setup SRAtools config 
####################################################

# Import SRA tools config, which can be prepared on a separate machine beforehand 
# e.g. on my local machine:
# tar -czvf sratoolkit-vdb-config_ec2-setup.tar.gz .ncbi/
echo Extracting sra-tools configuration file...
# SRA toolkit looks for configuration in $HOME/.ncbi/
tar -xvzf config/sra-tools-2.11.3-vdb-config.tar.gz -C ${HOME}


####################################################
# Script finish
####################################################

echo -e Script finished!
echo -e Please restart your Linux system for these changes to take effect.
echo conda envs can be activated using the command...
echo -e	"\t \$ conda activate rnaseq-preprocessing"
echo env can be deactivated using the command...
echo -e	"\t \$ conda deactivate"