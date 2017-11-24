#!/bin/bash

#PBS -l walltime=0:10:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=20m

source ${CONFIG_FILE}

### first check for existence of BAM files and their indexes!
#testing=FILES_TO_EVALUATE
#type=eval_job
#source ${TOOL_CHECK_FILES}
#[[ $? != 0 ]] && echo -e "\nTest for files for PID: ${PID} had non zero exit status, exiting pipeline\n\n" && exit 1
#[[ ${ok} == 0 ]] && echo -e "\nEvaluation of bam files for PID: ${PID} had non zero exit status, exiting pipeline\n\n" && exit 1

### Remove previous downstream files
#rm ${VCF_SOMATIC} ${VCF_SOMATIC_FUNCTIONAL}

########################################## Filter ###############################################

outputFilenamePrefix=${FILENAME_VCF%.vcf.gz}

if [[ "$GERMLINE_AVAILABLE" == 0 ]]; then
    FILTER_VALUES=""
    [[ ${FILTER_ExAC} == 'true' ]]         && FILTER_VALUES="${FILTER_VALUES} ${ExAC_COL} AF ${CRIT_ExAC_maxMAF}+"
    [[ ${FILTER_EVS} == 'true' ]]          && FILTER_VALUES="${FILTER_VALUES} ${EVS_COL} MAF ${CRIT_EVS_maxMAF}+"
    #[[ ${FILTER_1KGENOMES} == 'true' ]]    && FILTER_VALUES="${FILTER_VALUES} ${KGENOMES_COL} AF ${CRIT_1KGENOMES_maxMAF}+"
    #[[ ${FILTER_1KGENOMES} == 'true' ]]    && FILTER_VALUES="${FILTER_VALUES} ${KGENOMES_COL} ASN_AF ${CRIT_1KGENOMES_maxMAF}+"
    #[[ ${FILTER_1KGENOMES} == 'true' ]]    && FILTER_VALUES="${FILTER_VALUES} ${KGENOMES_COL} AMR_AF ${CRIT_1KGENOMES_maxMAF}+"
    #[[ ${FILTER_1KGENOMES} == 'true' ]]    && FILTER_VALUES="${FILTER_VALUES} ${KGENOMES_COL} AFR_AF ${CRIT_1KGENOMES_maxMAF}+"
    [[ ${FILTER_1KGENOMES} == 'true' ]]    && FILTER_VALUES="${FILTER_VALUES} ${KGENOMES_COL} EUR_AF ${CRIT_1KGENOMES_maxMAF}+"
    [[ ${FILTER_NON_CLINIC} == 'true' ]]   && FILTER_VALUES="${FILTER_VALUES} ${DBSNP_COL} CLN,COMMON nonexist,exist"
    [[ ${FILTER_LOCALCONTROL} == 'true' ]] && FILTER_VALUES="${FILTER_VALUES} ${LOCALCONTROL_COL} . ${CRIT_LOCALCONTROL_maxMAF}+"
    [[ ${FILTER_RECURRENCE} == 'true' ]]   && FILTER_VALUES="${FILTER_VALUES} ${RECURRENCE_COL} . ${CRIT_RECURRENCE}+"

    if [[ ${FILTER_VALUES} != "" ]]; then
        outputFilenamePrefix="${outputFilenamePrefix}_postFiltered"
        FILTERED_VCF="${outputFilenamePrefix}.vcf"
        ${PYPY_BINARY} -u ${TOOL_VCF_FILTER_BY_CRIT} ${FILENAME_VCF} ${FILTERED_VCF}${FILTER_VALUES}
        ${BGZIP_BINARY} -f ${FILTERED_VCF} && ${TABIX_BINARY} -f -p vcf ${FILTERED_VCF}.gz
        FILENAME_VCF=${FILTERED_VCF}.gz
    fi
fi

##### Set filenames #####
somatic_indels_vcf=${outputFilenamePrefix}_somatic_indels_conf_${MIN_CONFIDENCE_SCORE}_to_10.vcf
somatic_functional_indel_vcf=${outputFilenamePrefix}_somatic_functional_indels_conf_${MIN_CONFIDENCE_SCORE}_to_10.vcf
somatic_ncRNA_indel_vcf=${outputFilenamePrefix}_somatic_ncRNA_indels_conf_${MIN_CONFIDENCE_SCORE}_to_10.vcf
germline_functional_indel_vcf=${outputFilenamePrefix}_germline_functional_indels_conf_${MIN_CONFIDENCE_SCORE}_to_10.vcf
combined_screen_shots=indel_somatic_functional_combined.pdf

screenshot_dir=$(dirname ${FILENAME_VCF})/screenshots


set -xuv
set -o pipefail

${PERL_BINARY} ${TOOL_PLATYPUS_INDEL_EXTRACTOR} --bgzip=${BGZIP_BINARY} --tabix=${TABIX_BINARY} --infile=${FILENAME_VCF} --somout=${somatic_indels_vcf} --funcout=${somatic_functional_indel_vcf} --ncout=${somatic_ncRNA_indel_vcf} --germlineout=${germline_functional_indel_vcf} --minconf=${MIN_CONFIDENCE_SCORE} ${ADDITIONAL_FILTER_OPTS}

[[ "$?" != 0 ]] && echo "There was a non-zero exit filter step" && exit 1

[[ -d ${screenshot_dir} ]] && rm -rf ${screenshot_dir}

[[ ! -d ${screenshot_dir} ]] && mkdir ${screenshot_dir} && cd ${screenshot_dir}

[[ "$GERMLINE_AVAILABLE" == 1 ]] && ${PYTHON_BINARY} ${TOOL_SCREENSHOT} --vcf=${somatic_functional_indel_vcf} --control=${FILE_CONTROL_BAM} --tumor=${FILE_TUMOR_BAM} --ref=${REFERENCE_GENOME} --prefix=${VCF_SCREENSHOTS_PREFIX} --window=${WINDOW_SIZE} --annotations=${REPEAT_MASKER} --samtoolsbin=${SAMTOOLS_BINARY} --tabixbin=${TABIX_BINARY}
[[ "$GERMLINE_AVAILABLE" == 0 ]] && ${PYTHON_BINARY} ${TOOL_SCREENSHOT} --vcf=${somatic_functional_indel_vcf} --tumor=${FILE_TUMOR_BAM} --ref=${REFERENCE_GENOME} --prefix=${VCF_SCREENSHOTS_PREFIX} --window=${WINDOW_SIZE} --annotations=${REPEAT_MASKER} --samtoolsbin=${SAMTOOLS_BINARY} --tabixbin=${TABIX_BINARY}


pngs=(`ls *.pdf`)
sorted=$(printf "%s\n" ${pngs[@]}|sort -k1,1V)

${GHOSTSCRIPT_BINARY} -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=${combined_screen_shots} ${sorted}

#### Indel QC json file
resultBasePath=`dirname ${FILENAME_VCF}`
${PERL_BINARY} ${TOOL_PLATYPUS_INDEL_JSON} ${somatic_indels_vcf} > ${resultBasePath}/indel.json

touch ${FILENAME_CHECKPOINT}
