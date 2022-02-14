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


# export NETHOME=/nethome/jsc228/LipidModule
# export LIPID_HOME=/scratch/projects/lemmon/jsc228/LipidModule
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

cd {NETHOME}


# Download sra tools ----------------------------------------------------------------------

if [[ ! -f bin/sratoolkit.2.11.3-ubuntu64.tar.gz ]]
then
  echo sra-tools v.2.11.3 is not available through bioconda. Installing directly from NCBI...
  wget -P bin/ https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.3/sratoolkit.2.11.3-ubuntu64.tar.gz
fi
if [[ ! -d bin/sratoolkit.2.11.3-ubuntu64 ]]
then
  tar -vxzf bin/sratoolkit.2.11.3-ubuntu64.tar.gz -C bin/
fi

echo Extracting sra-tools configuration file...
# Import SRA tools config, which can be prepared on a separate machine beforehand 
# e.g. on my local machine:
# tar -czvf sratoolkit-vdb-config_ec2-setup.tar.gz .ncbi/
tar -xvzf config/sra-tools-2.11.3-vdb-config.tar.gz -C ${HOME}
export PATH=/bin/sratoolkit.2.11.3-ubuntu64/bin:$PATH
echo "PATH:"
echo $PATH



# Download Miniconda3 ---------------------------------------------------------------------

echo Checking if conda is executable from PATH...
if [ "" == "$(which conda)" ]
then
  echo conda not found. Downloading Miniconda3 installation script...
  sleep 1s # Slows down script to make terminal output more readable
  if [[ "$(uname -m)" == "x86_64" || "$HOSTTYPE" == "X86_64" ]]
  # If the Linux system is 64-bit...
  then
    # Download the script to install the 64-bit version of miniconda
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
      -O ${HOME}/miniconda.sh
  # If the Linux system is not 64-bit...
  else
    # Download the script to install the 32-bit version of miniconda
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86.sh \
      -O ${HOME}/miniconda.sh
  fi

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


# Install LipidModule venv ----------------------------------------------------------------

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
    conda list --explicit > config/tmp_env.txt
    ENV_DIFF=$(diff config/tmp_env.txt config/LipidModule.txt | wc -l)
    rm config/tmp_env.txt
    if [ ${ENV_DIFF} -ge 1 ]
    # True if environment specifications not identical
    then
      REINSTALL_ENV=1
    else
      echo LipidModule virtual environment is up to date. Exiting script.
      exit 0
    fi
  else 
    conda env export --name LipidModule > config/tmp_env.yml
    ENV_DIFF=$(diff config/tmp_env.yml config/LipidModule.yml | wc -l)
    rm config/tmp_env.yml
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
	  conda create --name LipidModule --file config/LipidModule.txt
  else
	  conda env create -f config/LipidModule.yml
  fi
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


# Download necessary genome FASTAs and GTFs -----------------------------------------------

# FASTA
if [[ ! -f ref/GRCm39.primary_assembly.genome.fa ]]
then
  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M28/GRCm39.primary_assembly.genome.fa.gz \
    -P ref/
  gunzip ref/GRCm39.primary_assembly.genome.fa.gz
fi

# GTF
if [[ ! -f ref/gencode.vM28.annotation.gtf ]]
then
  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M28/gencode.vM28.annotation.gtf.gz \
    -P ref/
  gunzip ref/gencode.vM28.annotation.gtf.gz
fi


# Install nextflow ------------------------------------------------------------------------

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
cd ${NETHOME}

cp ${NETHOME} ${LIPID_HOME}