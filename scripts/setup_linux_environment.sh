#! /bin/bash/

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
The 'bash' command is used to run this script: \n
\t \$ bash $0 \n
\n
Optional arguments: \n
\t      -h | --help\t         show this help text and exit \n
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

# Set working directory to LipidModule home
LIPID_HOME=$(pwd)
cd $LIPID_HOME
echo $LIPID_HOME

echo Checking if conda is executable from PATH...
if [ "" == "$(which conda)" ]
then
  echo conda not found. Downloading Miniconda3 installation script...
  sleep 1s # Slows down script to make terminal output more readable
  if [ "$(uname -m)" == "x86_64" ]
  # If the Linux system is 64-bit...
  then
    # Download the script to install the 64-bit version of miniconda
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
      -O {HOME}/miniconda.sh
  # If the Linux system is not 64-bit...
  else
    # Download the script to install the 32-bit version of miniconda
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86.sh \
      -O {HOME}/miniconda.sh
  fi

  echo Installing Miniconda3...; sleep 1s
  bash miniconda.sh -b -p ${HOME}/miniconda3

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
  if [ "$(uname -m)" == "x86_64" ]
  then
    conda list --explicit > ${LIPID_HOME}/tmp_env.txt
    ENV_DIFF=$(diff ${LIPID_HOME}/tmp_env.txt \
      ${LIPID_HOME}/env/LipidModule.txt | wc -l)
    rm ${LIPID_HOME}/tmp_env.txt
    if [ ${ENV_DIFF} -ge 1 ]
    # True if environment specifications not identical
    then
      REINSTALL_ENV=1
    else
      echo LipidModule virtual environment is up to date. Exiting script.
      exit 0
    fi
  else 
    conda env export --name LipidModule > ${LIPID_HOME}/tmp_env.yml
    ENV_DIFF=$(diff ${LIPID_HOME}/tmp_env.yml \
      ${LIPID_HOME}/env/LipidModule.yml | wc -l)
    rm ${LIPID_HOME}/tmp_env.yml
    if [ ${ENV_DIFF} -ge 1 ]
    then
      REINSTALL_ENV=1
    else
      echo LipidModule virtual environment is up to date. Exiting script.
      exit 0
    fi
  fi
# If environment does not exist...
else
  REINSTALL_ENV=1
  echo LipidModule virtual environment does not exist.; sleep 1s
fi

echo Creating the LipidModule virtual environment using conda...; sleep 1s
if [ REINSTALL_ENV == 1 ]
then
  if [ "$(uname -m)" == "x86_64" ]
  then
	  conda create --name LipidModule \
      --file ${LIPID_HOME}/envs/LipidModule.txt
  else
	  conda env create -f ${LIPID_HOME}/envs/LipidModule.yml
  fi
fi

echo Removing unused packages and caches using conda...
sleep 1s # Slows down script to make terminal output more readable
conda clean --all --yes

echo sra-tools v.2.11.3 is not available through bioconda. Installing directly from NCBI...
if [ ! -f 'env/sratoolkit.2.11.3-ubuntu64.tar.gz' ]
then
  wget -P ${LIPID_HOME}/env/ https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.3/sratoolkit.2.11.3-ubuntu64.tar.gz
fi
tar -vxzf ${LIPID_HOME}/env/sratoolkit.2.11.3-ubuntu64.tar.gz -C ${LIPID_HOME}/env/

echo Extracting sra-tools configuration file...
# Import SRA tools config, which can be prepared on a separate machine beforehand 
# e.g. on my local machine:
# tar -czvf sratoolkit-vdb-config_ec2-setup.tar.gz .ncbi/
tar -xvzf ${LIPID_HOME}/env/sra-tools-2.11.3-vdb-config.tar.gz -C ${HOME}
export PATH=${PATH}:${LIPID_HOME}/env/sratoolkit.2.11.3-ubuntu64/bin/

echo -e Script finished!

echo -e Please restart your Linux system for these changes to take effect.

echo The LipidModule environment can be activated using the command...
echo -e	"\t \$ conda activate LipidModule"
echo A conda virtual environment can be deactivated using the command...
echo -e	"\t \$ conda deactivate"