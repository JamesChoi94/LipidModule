# LipidModule

RNAseq-based detection of genes related to macrophage lipid accumulation relevant to CNS injury.

## Running the analysis

This analysis should be run on a UNIX-based machine.

1. Create a tab-delimited sample spreadsheet of study GEO accessions from which you wish to pull. See [sample spreadsheet](data/samples_sheet.csv) for an example. **Note**: Second column, which contains the GEO accession codes, is *required*. I recommend having study titles in the first column. All other columns are not required.

1. Clone the github repository. 
    ```
    git clone https://github.com/JamesChoi94/LipidModule.git
    ```
1. Setup environment and install packages.
    ```
    bash scripts/setup_linux_environment.sh
    ```
1. Setup environment variables.
    ```
    source ENV.sh
    ```   
1. Query NCBI.
    ```
    R scripts/query_gsm2sra.r
    ```
1.
