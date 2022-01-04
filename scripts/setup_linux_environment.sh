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
cd ${pwd}

echo Checking if LipidModule environment is already installed...
sleep 2s # Slows down script to make terminal output more readable
if [ -d ${HOME}/miniconda/envs/LipidModule ]
# True if environment exists exists and is directory
then
  echo The LipidModule environment already exists, exiting script.
  exit 0
fi

echo Checking if conda is executable from PATH...
if [ "" != "$(which conda)" ]
then
  echo conda found at: $(which conda)
  conda info -a
else
  echo conda not found. Downloading Miniconda3 installation script...
  sleep 2s # Slows down script to make terminal output more readable
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

  echo Installing Miniconda3...
  sleep 2s # Slows down script to make terminal output more readable
  bash miniconda.sh -b -p ${HOME}/miniconda

  echo Miniconda3 installed, removing installation script...
  rm -f {HOME}/miniconda.sh

  echo Setting up Miniconda3...
  sleep 2s # Slows down script to make terminal output more readable
  source "${HOME}/miniconda/etc/profile.d/conda.sh"
  hash -r
  conda config \
    --set always_yes yes \
    --set changeps1 yes \
    --set auto_activate_base false
  conda update -q conda
  conda init

  echo Displaying information about current conda installation...
  sleep 2s # Slows down script to make terminal output more readable
  conda info -a
fi

echo Creating the LipidModule virtual environment using conda...
sleep 2s # Slows down script to make terminal output more readable
# If the Linux system is 64-bit...
if [ "$(uname -m)" == "x86_64" ];
then
	# Create the virtual environment using the explicit spec list
	conda create --name LipidModule \
		--file ~/LipidModule/envs/LipidModule.txt
# If the Linux system is not 64-bit...
else
	# Create the virtual environment using an "environment".yml file
	conda env create -f ~/LipidModule/envs/LipidModule.yml
fi

echo Removing unused packages and caches using conda...
sleep 2s # Slows down script to make terminal output more readable
conda clean --all --yes

echo -e Script finished! \n

echo -e Please restart your Linux system for these changes to take effect. \n

echo The LipidModule environment can be activated using the command...
echo -e	"\t \$ conda activate LipidModule"
echo A conda virtual environment can be deactivated using the command...
echo -e	"\t \$ conda deactivate"