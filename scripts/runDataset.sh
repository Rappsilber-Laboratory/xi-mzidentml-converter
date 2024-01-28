#!/usr/bin/env bash

# Load environment (and make the slrum command available)
. /etc/profile.d/slurm.sh

#This job resets one document(accession based) from mongodb

##### OPTIONS
# (required) the project accession
PROJECT_ACCESSION=""

##### VARIABLES
# the name to give to the LSF job (to be extended with additional info)
JOB_NAME="xi-mzidentml-converter"
# memory limit
MEMORY_LIMIT=8G
#email notification
JOB_EMAIL="pride-report@ebi.ac.uk"

##### FUNCTIONS
printUsage() {
    echo "Description: This will parse cross-linking dataset"
    echo "$ ./scripts/xi-mzidentml-converter.sh"
    echo ""
    echo "Usage: ./xi-mzidentml-converter.sh -a|--accession [-e|--email]"
    echo "     Example: ./xi-mzidentml-converter.sh -a PXD036833 --dontdelete -w api"
    echo "     (required) accession         : the project accession of a crosslinking dataset"
    echo "     (optional) email             :  Email to send LSF notification"
}

##### PARSE the provided parameters
while [ "$1" != "" ]; do
    case $1 in
      "-a" | "--accession")
        shift
        PROJECT_ACCESSION=$1
        ;;
    esac
    shift
done

##### CHECK the provided arguments
if [ -z ${PROJECT_ACCESSION} ]; then
         echo "Need to enter a project accession"
         printUsage
         exit 1
fi

DATE=$(date +"%Y%m%d%H%M")
LOG_FILE_NAME="${JOB_NAME}/${PROJECT_ACCESSION}-${DATE}.log"
LOG_FILE_ERROR_NAME="${JOB_NAME}/${PROJECT_ACCESSION}-${DATE}-error.log"
JOB_NAME="${JOB_NAME}-${PROJECT_ACCESSION}"


##### Change directory to where the script locate
cd ${0%/*}

##### Activate conda environment
source /hps/software/users/juan/pride/anaconda3/etc/profile.d/conda.sh
conda activate xi-mzidentml-converter

#### RUN it on the cluster #####
sbatch -t 7-0 \
     --mem=${MEMORY_LIMIT} \
     --partition datamover \
     --mail-type=ALL \
     --mail-user=${JOB_EMAIL} \
     --job-name=${JOB_NAME} \
     --output ${LOG_FILE_NAME} \
     --error ${LOG_FILE_ERROR_NAME} \
     --wrap="process_dataset -p ${PROJECT_ACCESSION} --dontdelete -w api"
