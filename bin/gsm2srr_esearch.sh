
####################################################
# Set NETHOME and SCRATCH directory variables
####################################################

export NETHOME=/nethome/jsc228/LipidModule
export SCRATCH=/scratch/projects/lemmon/jsc228/LipidModule
export SAMPLESHEET=$NETHOME/data/samplesheet.csv

cd ${NETHOME}

####################################################
# Convert GSM to SRR Accessions from samplesheet
####################################################

GSM_ACCESSIONS=($(awk -F, 'NR!=1 {print $1}' ${SAMPLESHEET}))
for GSM in ${GSM_ACCESSIONS[@]}
do
  echo "esearchng ${GSM}..."
  esearch -db sra -query ${GSM} | efetch -format runinfo >> \
    $SCRATCH/data/esearch-runinfo.txt
  sleep 2s
done

# Also create txt file of accessions only
awk -F, 'NR%2==0 {print $1}' data/esearch-runinfo.txt > data/srrAccessions.txt