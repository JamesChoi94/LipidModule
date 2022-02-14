#!/bin/bash
exec1>esearch-runinfo.txt
SAMPLE_SHEET_PATH=""
while [[ "$#" ]]
do
	case "$1" in
		-s|--samplesheet)
			SAMPLE_SHEET_PATH=$2
			shift 2
			;;
		--) # end argument parsing
			shift 1
			break 
			;;
		-*|--*) # unsupported flags
			echo -e "ERROR: $1 is an invalid option. \n" >&2
			usage
			exit 1
			;;
		*)
			break
			;;
	esac
done

GSM_ACCESSION=($(awk -F, 'NR>1 {print $1}' ${SAMPLE_SHEET_PATH} | sed -r 's/\s+//g' ))
# echo "Getting SRR accessions for the following GSM accessions:"
# echo ${GSM_ACCESSION[@]}
for gsm in "${GSM_ACCESSION[@]}"
do
  # echo ${gsm} > "${OUTDIR}/${gsm}_esearch-runinfo.txt"
  esearch -db sra -query ${gsm} | efetch -format runinfo 1>> \
    esearch-runinfo.txt
	sleep 2s
done