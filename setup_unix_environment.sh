#!/bin/sh

# These should be modified for every new study --------------------------------

# NETHOME and SCRATCH should point to the `nethome` and  `scratch` file 
# systems on Pegasus. This will allow larger files (e.g. fastq's) to be saved
# in the `scratch` file system and other, smaller analysis ouputs (e.g. read 
# mappings) to be saved in `nethome` for "permanent" storage. It should be 
# assumed that anything in `scratch` can be deleted after any run.

# If NETHOME and SCRATCH point to same directory, the large data files (e.g.
# fastq.gz and genome index/FASTAs) will be saved in the same project 
# directory. This should only really be done for testing purposes. Do NOT 
# point both to the `nethome` file system - your account may be suspended for
# inappropriate resource usage.


# export NETHOME=/nethome/jsc228/LipidModule
# export SCRATCH=/scratch/projects/lemmon/jsc228/LipidModule
# export NETHOME=/mnt/d/LipidModule
# export SCRATCH=/mnt/d/LipidModule

if [ ! -d ${NETHOME} ]
then
  echo 'NETHOME cannot be found.'
  return
fi

echo "NETHOME:"
echo ${NETHOME}
echo "SCRATCH:"
echo ${SCRATCH}


# Please note!! --------------------------------------
cd ${NETHOME} 
# ----------------------------------------------------


# Download sra tools ----------------------------------------------------------------------

# Downloads sra-toolkit into bin/ directory in $NETHOME. Make these packages available 
# in the project directory bin/ folder so that they are available always.
if [[ ! -f ${NETHOME}/bin/sratoolkit.2.11.3-ubuntu64.tar.gz ]]
then
  echo sra-tools v.2.11.3 is not available through bioconda. Installing directly from NCBI...
  wget -P ${NETHOME}/bin/ https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.3/sratoolkit.2.11.3-ubuntu64.tar.gz
fi
if [[ ! -d ${NETHOME}/bin/sratoolkit.2.11.3-ubuntu64 ]]
then
  tar -vxzf ${NETHOME}/bin/sratoolkit.2.11.3-ubuntu64.tar.gz -C ${NETHOME}/bin/
fi

# Import SRA tools config, which can be prepared on a separate machine beforehand 
# e.g. on my local machine:
# tar -czvf sratoolkit-vdb-config_ec2-setup.tar.gz .ncbi/
echo Extracting sra-tools configuration file...
# SRA toolkit looks for configuration in $HOME/.ncbi/
tar -xvzf config/sra-tools-2.11.3-vdb-config.tar.gz -C ${HOME}
# But executables can be made available in $NETHOME/bin
export PATH=${NETHOME}/bin/sratoolkit.2.11.3-ubuntu64/bin:$PATH


# Install EDirect -------------------------------------------------------------------------

if [[ -d ${NETHOME}/bin/edirect ]]
then
  sh -c "$(wget -q ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O ${NETHOME}/bin)"
fi
export PATH=${PATH}:${NETHOME}/bin/edirect


# # Install BEDOPS ----------------------------------------------------------------

# if [[ ! -d ${NETHOME}/bin/bedops_linux_x86_64-v2.4.40/bin ]]
# then
#   wget https://github.com/bedops/bedops/releases/download/v2.4.40/bedops_linux_x86_64-v2.4.40.tar.bz2 \
#     -P bin
#   mkdir bin/bedops_linux_x86_64-v2.4.40
#   tar -jxvf bin/bedops_linux_x86_64-v2.4.40.tar.bz2 --directory bin/bedops_linux_x86_64-v2.4.40/
# fi
# export PATH=${PATH}:${NETHOME}/bin/bedops_linux_x86_64-v2.4.40/bin


# Install nextflow ------------------------------------------------------------------------

# Installs nextflow into a bin/ folder in home directory
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
# export PATH=${PATH}:${HOME}/bin/edirect

# # Alternatively, install nextflow into the project directory instead of making
# # the executable available to user environment.
# if [[ ! -f ${NETHOME}/bin/nextflow ]]
# then
#   if [[ ! -f ${NETHOME}/nextflow ]]
#   then
#     cd ${NETHOME}/bin
#     wget -qO- https://get.nextflow.io | bash
#     chmod +x nextflow
#   fi
# fi
# export PATH=${PATH}:${NETHOME}/bin/edirect



# Download Miniconda3 ---------------------------------------------------------------------

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


# Install LipidModule venv ----------------------------------------------------------------

echo Checking if LipidModule virtual environment is already installed...
sleep 1s
if [[ -d ${HOME}/miniconda3/envs/LipidModule ]]
# True if environment exists exists and is directory...
then
  INSTALL_ENV=0
  echo The LipidModule virtual environment already exists.; sleep 1s
  echo Checking if installed packages are current...; sleep 1s
  conda list --name LipidModule --explicit > config/tmp_env.txt
  ENV_DIFF=$(diff config/tmp_env.txt config/LipidModule.txt | wc -l)
  rm config/tmp_env.txt
  if [ "${ENV_DIFF}" -ge 1 ]
  then
    echo "Conda diff length: ${ENV_DIFF}"
    conda activate LipidModule
    conda env update --name LipidModule --file config/LipidModule.txt --prune
  else
    echo LipidModule virtual environment is up to date. Exiting script.
    return 0
  fi
# If environment does not exist...
else
  INSTALL_ENV=1
  echo LipidModule virtual environment does not exist.; sleep 1s
fi

if [[ $INSTALL_ENV == 1 ]]
then
  echo Creating the LipidModule virtual environment using conda...; sleep 1s
  conda env create -f config/LipidModule.txt
fi

# echo Removing unused packages and caches using conda...
# sleep 1s # Slows down script to make terminal output more readable
# conda clean --all --yes
# conda config --add channels conda-forge
# conda config --add channels bioconda

echo -e Script finished!

echo -e Please restart your Linux system for these changes to take effect.

echo The LipidModule environment can be activated using the command...
echo -e	"\t \$ conda activate LipidModule"
echo A conda virtual environment can be deactivated using the command...
echo -e	"\t \$ conda deactivate"