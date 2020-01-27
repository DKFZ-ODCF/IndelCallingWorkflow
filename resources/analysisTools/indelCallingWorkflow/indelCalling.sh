#!/bin/bash
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/IndelCallingWorkflow).
#
#PBS -l walltime=50:00:00
#PBS -l nodes=1:ppn=12
#PBS -l mem=14g

set -o pipefail

### Create temporary filenames used in this script:
useCustomScratchDir=false
SCRATCH_DIR=${RODDY_SCRATCH}

LOG_TMP=${SCRATCH_DIR}/platyLog.tmp
####################################### Calling Script #################################################
CALL_SNP=${CALL_SNP-0}

mkdir -p ${DIR_TEMP}/indelCalling
fileECPlatypus=${DIR_TEMP}/indelCalling/ecplatypus.txt
fileECCat=${DIR_TEMP}/indelCalling/eccat.txt
fileECBGZip=${DIR_TEMP}/indelCalling/ecbgzip.txt
fileECConfAnno=${DIR_TEMP}/indelCalling/ecconfanno.txt

# Default those values to the defaults.
PLATYPUS_BUFFER_SIZE=${PLATYPUS_BUFFER_SIZE-100000}
PLATYPUS_MAX_READS=${PLATYPUS_MAX_READS-5000000}

[[ ! -r ${FILENAME_TUMOR_BAM} ]] && echo "Tumor bam is missing or not right." && exit 201
bamFiles=${FILENAME_TUMOR_BAM}
if [[ ${isControlWorkflow} == true ]]; then
    [[ ! -r ${FILENAME_CONTROL_BAM} ]] && echo "Control bam is missing or not right." && exit 200
    bamFiles="${FILENAME_CONTROL_BAM},${bamFiles}"
fi

#source ${TOOL_ANALYZE_BAM_HEADER}
#getRefGenomeAndChrPrefixFromHeader ${FILENAME_TUMOR_BAM} # Sets CHR_PREFIX and REFERENCE_GENOME

${PLATYPUS_BINARY} callVariants \
	--refFile=${REFERENCE_GENOME} \
	--output=${FILENAME_VCF_RAW}.tmp.platypus \
	--bamFiles=${bamFiles} \
	--nCPU=${CPU_COUNT} \
	--genIndels=1 \
	--genSNPs=${CALL_SNP} \
	--logFileName=${LOG_TMP} \
	--verbosity=1 \
	--bufferSize=${PLATYPUS_BUFFER_SIZE} \
	--maxReads=${PLATYPUS_MAX_READS} \
        --minFlank=0 \
	${PLATYPUS_PARAMS}

[[ $? -gt 0 ]] && echo "Error during platypus indel calling." && exit 1

lineCount=`grep -v "^#" ${FILENAME_VCF_RAW}.tmp.platypus | cut -f 12 | sort | uniq | grep -v '^$' | wc -l`
echo "Number of corrupt rows: $lineCount"

if [[ $lineCount -gt 0 ]]
then
  (grep "#" ${FILENAME_VCF_RAW}.tmp.platypus ; grep -v "^#" ${FILENAME_VCF_RAW}.tmp.platypus | awk '{if(NF == 11){print $0}}') > ${FILENAME_VCF_RAW}.tmp.platypus.11
  [[ $? -gt 0 ]] && echo "Error during platypus indel calling." && exit 3

  ${BGZIP_BINARY} -c -f ${FILENAME_VCF_RAW}.tmp.platypus.11 > ${FILENAME_VCF_RAW}.tmp && rm ${FILENAME_VCF_RAW}.tmp.platypus.11
  [[ $? -gt 0 ]] && echo "Error during platypus indel calling." && exit 4

  grep -v "^#" ${FILENAME_VCF_RAW}.tmp.platypus | awk '{if(NF > 11 ){print $0}}' > ${FILENAME_VCF_RAW}.tmp.platypus.linesCorrupt && rm ${FILENAME_VCF_RAW}.tmp.platypus
  [[ $? -gt 0 ]] && echo "Error during platypus indel calling." && exit 5

  corruptLines=`cat ${FILENAME_VCF_RAW}.tmp.platypus.linesCorrupt | wc -l` 
  echo "$corruptLine corrupt lines in the platypus indel calling raw file."

  if [[ $corruptLines -gt 10 ]]
  then 
    echo "Error: More than 10 corrupt lines in platypus indel calling." && exit 6
  fi
else 
  ${BGZIP_BINARY} -c -f ${FILENAME_VCF_RAW}.tmp.platypus > ${FILENAME_VCF_RAW}.tmp && rm ${FILENAME_VCF_RAW}.tmp.platypus

  [[ $? -gt 0 ]] && echo "Error during platypus indel calling." && exit 2

  #lineCount=`zgrep -v "^#" ${FILENAME_VCF_RAW}.tmp | cut -f 12 | sort | uniq -c | wc -l`
fi

#Sort Indel alphanumerically, needed for hg38 transfer
(zcat ${FILENAME_VCF_RAW}.tmp | grep '#' ; zcat ${FILENAME_VCF_RAW}.tmp | grep -v '#' | sort -V -k1,2) > indelSorted_pid.vcf.raw.tmp

bgzip indelSorted_pid.vcf.raw.tmp

mv ${FILENAME_VCF_RAW}.tmp indelUnsorted_pid.vcf.raw.gz #Just for control. Can be deleted in future script

mv indelSorted_pid.vcf.raw.tmp.gz ${FILENAME_VCF_RAW}.tmp

#Last line is from old code
mv ${FILENAME_VCF_RAW}.tmp ${FILENAME_VCF_RAW} && ${TABIX_BINARY} -f -p vcf ${FILENAME_VCF_RAW}
