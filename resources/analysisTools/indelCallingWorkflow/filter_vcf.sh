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

##### Set filenames #####
somatic_indels_vcf=${FILENAME_VCF%.vcf.gz}_somatic_indels_conf_${MIN_CONFIDENCE_SCORE}_to_10.vcf
somatic_functional_indel_vcf=${FILENAME_VCF%.vcf.gz}_somatic_functional_indels_conf_${MIN_CONFIDENCE_SCORE}_to_10.vcf
somatic_ncRNA_indel_vcf=${FILENAME_VCF%.vcf.gz}_somatic_ncRNA_indels_conf_${MIN_CONFIDENCE_SCORE}_to_10.vcf
germline_functional_indel_vcf=${FILENAME_VCF%.vcf.gz}_germline_functional_indels_conf_${MIN_CONFIDENCE_SCORE}_to_10.vcf
combined_screen_shots=indel_somatic_functional_combined.pdf

screenshot_dir=$(dirname ${FILENAME_VCF})/screenshots


set -xuv
set -o pipefail

${PERL_BINARY} ${TOOL_PLATYPUS_INDEL_EXTRACTOR} --bgzip=${BGZIP_BINARY} --tabix=${TABIX_BINARY} --infile=${FILENAME_VCF} --somout=${somatic_indels_vcf} --funcout=${somatic_functional_indel_vcf} --ncout=${somatic_ncRNA_indel_vcf} --germlineout=${germline_functional_indel_vcf} --minconf=${MIN_CONFIDENCE_SCORE} ${ADDITIONAL_FILTER_OPTS}

[[ "$?" != 0 ]] && echo "There was a non-zero exit filter step" && exit 1

[[ -d ${screenshot_dir} ]] && rm -rf ${screenshot_dir}

[[ ! -d ${screenshot_dir} ]] && mkdir ${screenshot_dir} && cd ${screenshot_dir} && ${PYTHON_BINARY} ${TOOL_SCREENSHOT} --vcf=${somatic_functional_indel_vcf} --control=${FILE_CONTROL_BAM} --tumor=${FILE_TUMOR_BAM} --ref=${REFERENCE_GENOME} --prefix=${VCF_SCREENSHOTS_PREFIX} --window=${WINDOW_SIZE} --annotations=${REPEAT_MASKER} --samtoolsbin=${SAMTOOLS_BINARY} --tabixbin=${TABIX_BINARY}

pngs=(`ls *.pdf`)
sorted=$(printf "%s\n" ${pngs[@]}|sort -k1,1V)

${GHOSTSCRIPT_BINARY} -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=${combined_screen_shots} ${sorted}

touch ${FILENAME_CHECKPOINT}
