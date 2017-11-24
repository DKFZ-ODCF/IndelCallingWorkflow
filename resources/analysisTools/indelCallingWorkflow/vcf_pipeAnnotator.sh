#!/bin/bash

#PBS -l walltime=0:50:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=200m

#source ${CONFIG_FILE}
#exit 0
#### Check and create output file names
#testing=FILE_NAMES
#type=create
#source ${TOOL_CHECK_FILES}
#[[ $? != 0 ]] && echo -e "\nCreation of file names for PID: ${PID} had non zero exit status, exiting pipeline\n\n" && exit 2
#
#### Remove previous downstream files
#rm ${VCF_SOMATIC} ${VCF_SOMATIC_FUNCTIONAL}

### Temp files:
filenameVCFTemp=${FILENAME_VCF}.tmp

#set -x
#set -o pipefail

############################################ Deep Annotation ####################################################

if [[ ! -f ${FILENAME_VCF} ]]
then
	echo input file ${FILENAME_VCF} does not exist
	exit 73
fi

# create a bunch of pipes using Matthias' overlapper for multiple annotation files and do a system call with this
pipe=`${PERL_BINARY} ${TOOL_CREATEPIPES} ${FILENAME_VCF} ${CONFIG_FILE} ${TOOL_ANNOTATE_VCF_FILE} ${PIPENAME}  ${TABIX_BINARY}`

if [[ "$?" != 0 ]] || [[ -z "${pipe}" ]]; then echo "problem when generating pipe: $PIPENAME. Exiting..."; exit 2; fi

eval ${pipe} > ${filenameVCFTemp}

if [[ "$?" == 0 ]]
then
	cat ${filenameVCFTemp} | ${BGZIP_BINARY} > ${FILENAME_VCF} && ${TABIX_BINARY} -f -p vcf ${FILENAME_VCF} && rm ${filenameVCFTemp} && touch ${FILENAME_CHECKPOINT}
else
    echo "There was a non-zero exit code in the ${PIPENAME} pipe; temp file ${filenameVCFTemp} not moved back"
    exit 2
fi

