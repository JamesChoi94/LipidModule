#!/bin/bash/

# Script to setup project environment on a Linux machine. Downloads and 
# installs miniconda3, and uses conda to install the 'LipidModule' 
# virtual environment. Code adapted and modified from original script by
# Ronan Harrington (https://github.com/rnnh/bioinfo-notebook.git)

# Help/usage text
usage="$(basename "$0") \n
\n
This script downloads and installs Miniconda3, and uses conda to install \n
the 'LipidModule' virtual environment. NOTE: installation of miniconda \n
will create a new 'miniconda3' directory in your $HOME. \n
\n
Before running this script... \n
\n
\t 1. Run the following command to ensure that installed software is up \n
\t to date: \n
\t \t \$ sudo apt-get update \n
\n 
\t 2. Clone the github repository: \n
\t \t $ git clone https://github.com/JamesChoi94/LipidModule.git \n
\n
\t 3. Set LipidModule/ as the working directory and double-check path to \n
\t this directory:
\t \t $ cd LipidModule
\t \t $ pwd
\n
\n
"

# Iterate through input arguments
while (( "$#" )) # 
do
	case "$1" in
		-h|--help)
			echo -e $usage
			exit 0
			;;
	esac
done

# Download sra tools
if [[ ! -f ${NETHOME}/bin/sratoolkit.2.11.3-ubuntu64.tar.gz ]]
then
  echo sra-tools v.2.11.3 is not available through bioconda. Installing directly from NCBI...
  wget -P ${NETHOME}/bin/ https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.3/sratoolkit.2.11.3-ubuntu64.tar.gz
fi
if [[ ! -d ${NETHOME}/bin/sratoolkit.2.11.3-ubuntu64 ]]
then
  tar -vxzf ${NETHOME}/bin/sratoolkit.2.11.3-ubuntu64.tar.gz -C ${NETHOME}/bin/
fi

echo Extracting sra-tools configuration file...
# Import SRA tools config, which can be prepared on a separate machine beforehand 
# e.g. on my local machine:
# tar -czvf sratoolkit-vdb-config_ec2-setup.tar.gz .ncbi/
tar -xvzf ${NETHOME}/config/sra-tools-2.11.3-vdb-config.tar.gz -C ${HOME}
export PATH=${PATH}:${NETHOME}/bin/sratoolkit.2.11.3-ubuntu64/bin/


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
  rm -f {HOME}/miniconda.sh

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

echo Checking if LipidModule virtual environment is already installed...
sleep 1s
if [ -d ${HOME}/miniconda3/envs/LipidModule ]
# True if environment exists exists and is directory...
then
  REINSTALL_ENV=0
  echo The LipidModule virtual environment already exists.; sleep 1s
  echo Checking if installed packages are current...; sleep 1s
  conda activate LipidModule
  conda list --explicit > ${NETHOME}/tmp_env.yml
  ENV_DIFF=$(diff ${NETHOME}/tmp_env.yml \
    ${NETHOME}/config/LipidModule.yml | wc -l)
  rm ${NETHOME}/tmp_env.yml
  if [ ${ENV_DIFF} -ge 1 ]
    # True if environment specifications not identical
  then
    REINSTALL_ENV=1
  else
    echo LipidModule virtual environment is up to date. Exiting script.
    exit 0
  fi
# If environment does not exist...
else
  REINSTALL_ENV=1
  echo LipidModule virtual environment does not exist.; sleep 1s
fi

echo Creating the LipidModule virtual environment using conda...; sleep 1s
if [ REINSTALL_ENV == 1 ]
then
  conda create --name LipidModule \
    --file ${NETHOME}/config/LipidModule.yml
else
  conda env create -f ${NETHOME}/config/LipidModule.yml 
fi

echo Removing unused packages and caches using conda...
sleep 1s # Slows down script to make terminal output more readable
conda clean --all --yes

echo -e Script finished!

echo -e Please restart your Linux system for these changes to take effect.

echo The LipidModule environment can be activated using the command...
echo -e	"\t \$ conda activate LipidModule"
echo A conda virtual environment can be deactivated using the command...
echo -e	"\t \$ conda deactivate"