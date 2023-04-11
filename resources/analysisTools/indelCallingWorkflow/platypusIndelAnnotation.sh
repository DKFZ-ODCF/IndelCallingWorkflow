#!/bin/bash
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/IndelCallingWorkflow).
#
#PBS -l walltime=20:00:00
#PBS -l nodes=1:ppn=2
#PBS -l mem=2600m

#### first check for existence of BAM files and their indexes!
#testing=FILES_TO_EVALUATE
#type=eval_job
#source ${TOOL_CHECK_FILES}
#[[ $? != 0 ]] && echo -e "\nTest for files for PID: ${PID} had non zero exit status, exiting pipeline\n\n" && exit 1
#[[ ${ok} == 0 ]] && echo -e "\nEvaluation of bam files for PID: ${PID} had non zero exit status, exiting pipeline\n\n" && exit 1
#
#### Check and create output file names
#testing=FILENAME_NAMES
#type=create
#source ${TOOL_CHECK_FILES}
#[[ $? != 0 ]] && echo -e "\nCreation of file names for PID: ${PID} had non zero exit status, exiting pipeline\n\n" && exit 2
#
#### Test availability of previous files
#[[ ! -f ${FILENAME_VCF} ]] && echo -e "The raw file: ${VCF_RAW} was not found, exiting..." && exit 2
#
#### Remove previous downstream files
#rm ${VCF_FINAL} ${VCF_SOMATIC} ${VCF_SOMATIC_FUNCTIONAL}

### Create temporary filenames used in this script:
#TMP_FILENAME_PREFIX=${RESULTS_PER_PIDS_DIR}/${PID}/${RESULTS_SUB_DIR}/${INDELFILENAME_PREFIX}${PID}
OUTPUT_DIRECTORY=`dirname ${FILENAME_VCF_RAW}`
TMP_FILENAME_PREFIX=`dirname ${FILENAME_VCF_RAW}`/indels_${PID}
filenameVCFTemp=${FILENAME_VCF_RAW%.gz}.tmp
filenameVCFFinalUnzipped=${FILENAME_VCF_OUT%.gz}

filenameVCFPancan=${FILENAME_VCF_OUT%.vcf.gz}_pancan.vcf.gz
filenameVCFPancanTemp=${filenameVCFPancan%.gz}.tmp
filenameVCFPancanUnzipped=${filenameVCFPancan%.gz}

FOR_ANNOVAR=${TMP_FILENAME_PREFIX}.ForAnnovar.bed
ANNOVAR_SEGDUP=${TMP_FILENAME_PREFIX}.Annovar.segdup
ANNOVAR_CYTOBAND=${TMP_FILENAME_PREFIX}.Annovar.cytoband
ANNOVAR_VARIANT=${TMP_FILENAME_PREFIX}.ForAnnovar.bed.variant_function
ANNOVAR_VARIANT_EXON=${TMP_FILENAME_PREFIX}.ForAnnovar.bed.exonic_variant_function


set -x
set -o pipefail

### Annotate with polymorphisms (dbSNP, 1K genomes, ExAC, EVS, and local controls) and prepare annovar input file
cmdAnnotation="zcat ${FILENAME_VCF_RAW} | \
    ${PERL_BINARY} ${TOOL_ANNOTATE_VCF_FILE} -a - -b ${DBSNP_INDEL} --columnName=${DBSNP_COL} --reportMatchType  --bAdditionalColumn=2 --reportBFeatCoord --padding=${INDEL_ANNOTATION_PADDING} --minOverlapFraction=${INDEL_ANNOTATION_MINOVERLAPFRACTION} --maxBorderDistanceSum=${INDEL_ANNOTATION_MAXBORDERDISTANCESUM} --maxNrOfMatches=${INDEL_ANNOTATION_MAXNROFMATCHES} | \
    ${PERL_BINARY} ${TOOL_ANNOTATE_VCF_FILE} -a - -b ${KGENOME} --columnName=${KGENOMES_COL} --reportMatchType --bAdditionalColumn=2  --reportBFeatCoord --padding=${INDEL_ANNOTATION_PADDING} --minOverlapFraction=${INDEL_ANNOTATION_MINOVERLAPFRACTION} --maxBorderDistanceSum=${INDEL_ANNOTATION_MAXBORDERDISTANCESUM} --maxNrOfMatches=${INDEL_ANNOTATION_MAXNROFMATCHES} | \
    ${PERL_BINARY} ${TOOL_ANNOTATE_VCF_FILE} -a - -b ${ExAC} --columnName=${ExAC_COL} --bFileType vcf --reportLevel 4 --reportMatchType | \
    ${PERL_BINARY} ${TOOL_ANNOTATE_VCF_FILE} -a - -b ${EVS} --columnName=${EVS_COL} --bFileType vcf --reportLevel 4 --reportMatchType| \
    ${PERL_BINARY} ${TOOL_ANNOTATE_VCF_FILE} -a - -b ${GNOMAD_WES_ALL_INDEL} --columnName=${GNOMAD_WES_COL} --bFileType vcf --reportLevel 4 --reportMatchType | \
    ${PERL_BINARY} ${TOOL_ANNOTATE_VCF_FILE} -a - -b ${GNOMAD_WGS_ALL_INDEL} --columnName=${GNOMAD_WGS_COL} --bFileType vcf --reportLevel 4 --reportMatchType | \
    ${PERL_BINARY} ${TOOL_ANNOTATE_VCF_FILE} -a - -b ${LOCALCONTROL} --columnName=${LOCALCONTROL_COL} --minOverlapFraction 1 --bFileType vcf --reportLevel 4 --reportMatchType| \
    tee ${filenameVCFTemp} | perl ${TOOL_VCF_TO_ANNOVAR} ${CHR_PREFIX} ${CHR_SUFFIX} > ${FOR_ANNOVAR}.tmp"

eval ${cmdAnnotation}

[[ "$?" != 0 ]] && echo "There was a non-zero exit code in the polymorphism annotation pipe; exiting..." && exit 3

mv ${filenameVCFTemp} ${filenameVCFFinalUnzipped}
mv ${FOR_ANNOVAR}.tmp ${FOR_ANNOVAR}


###### Basic annotation with Annovar
# Gene annotation with annovar
ANNOVAR_FILE=${ANNOVAR_BUILDVER}_`eval echo $ANNOVAR_DBTYPE | cut -d " " -f 2`.txt
ANNOVAR_DBFILEPATH=${ANNOVAR_DBPATH}/${ANNOVAR_FILE}
[[ ! -f ${ANNOVAR_DBFILEPATH} ]]  && echo "Gene annotation database not found. Check ANNOVAR_DBTYPE." && exit -16

${ANNOVAR_BINARY} --buildver=${ANNOVAR_BUILDVER} ${ANNOVAR_DBTYPE} ${FOR_ANNOVAR} ${ANNOVAR_DBPATH}

# segdup annotation with annovar
${ANNOVAR_BINARY} --buildver=${ANNOVAR_BUILDVER} -regionanno -dbtype segdup --outfile=${ANNOVAR_SEGDUP} ${FOR_ANNOVAR} ${ANNOVAR_DBPATH}
av_segdup=`ls ${ANNOVAR_SEGDUP}*genomicSuperDups`

# cytoband annotation with annovar
${ANNOVAR_BINARY} --buildver=${ANNOVAR_BUILDVER} -regionanno -dbtype band --outfile=${ANNOVAR_CYTOBAND} ${FOR_ANNOVAR} ${ANNOVAR_DBPATH}
av_cytoband=`ls ${ANNOVAR_CYTOBAND}*cytoBand`

perl ${TOOL_NEW_COLS_TO_VCF} --vcfFile=${filenameVCFFinalUnzipped} \
    --newColFile="${TOOL_PROCESS_ANNOVAR} ${ANNOVAR_VARIANT} ${ANNOVAR_VARIANT_EXON} |" \
    --newColHeader=${ANNOVAR_GENEANNO_COLS} --chrPrefix=${CHR_PREFIX} --chrSuffix=${CHR_SUFFIX} --reportColumns="3,4,5,6" --bChrPosEnd="0,1,2" |
perl ${TOOL_NEW_COLS_TO_VCF} --vcfFile="-" --newColFile=${av_segdup} --newColHeader=${ANNOVAR_SEGDUP_COL} --chrPrefix=${CHR_PREFIX} --chrSuffix=${CHR_SUFFIX} \
    --reportColumns="1" --bChrPosEnd="2,7,8"  |
perl ${TOOL_NEW_COLS_TO_VCF} --vcfFile="-" --newColFile=${av_cytoband} --newColHeader=${ANNOVAR_CYTOBAND_COL} --chrPrefix=${CHR_PREFIX} --chrSuffix=${CHR_SUFFIX} \
    --reportColumns="1" --bChrPosEnd="2,7,8" > ${filenameVCFTemp}


[[ "$?" != 0 ]] && echo "There was a non-zero exit code in the Annovar annotation pipe; temp file ${filenameVCFTemp} not moved back" && exit 4

mv ${filenameVCFTemp} ${filenameVCFFinalUnzipped}

indel_reliability_pipe=`${PERL_BINARY} ${TOOL_CREATEPIPES} ${filenameVCFFinalUnzipped} ${PARAMETER_FILE} ${TOOL_ANNOTATE_VCF_FILE} INDEL_RELIABILITY  ${TABIX_BINARY}`

if [[ "$?" != 0 ]] || [[ -z "${indel_reliability_pipe}" ]]; then echo "problem when generating INDEL_RELIABILITY pipe. Exiting..."; exit 5; fi

eval ${indel_reliability_pipe} > ${filenameVCFTemp}

if [[ "$?" != 0 ]]; then
    echo "There was a non-zero exit code in the INDEL_RELIABILITY pipe; temp file ${filenameVCFTemp} not moved back" && \
    exit 5
fi

mv ${filenameVCFTemp} ${filenameVCFFinalUnzipped}


################################## Confidence annotation ############################################

### Get the real names of the columns created by Platypus
tumor_column=`${SAMTOOLS_BINARY} view -H ${FILENAME_TUMOR_BAM} | grep -m 1 SM: | ${PERL_BINARY} -ne 'chomp;$_=~m/SM:(\S+)/;print "$1\n";'`

target=/dev/null
[[ ${runOnPancan-false} == true ]] && target=${filenameVCFPancanTemp}

if [[ ${isControlWorkflow} == true ]]; then
    control_column=`${SAMTOOLS_BINARY} view -H ${FILENAME_CONTROL_BAM} \
      | grep -m 1 SM: \
      | ${PERL_BINARY} -ne 'chomp;$_=~m/SM:(\S+)/;print "$1\n";'`
    ${PYPY_BINARY} -u ${TOOL_PLATYPUS_CONFIDENCE_ANNOTATION}  \
      --infile=${filenameVCFFinalUnzipped} \
      --controlColName=${control_column} \
      --tumorColName=${tumor_column} \
      ${CONFIDENCE_OPTS_INDEL} \
      | tee ${filenameVCFTemp} \
      | cut -f 1-11 \
      > ${target}
else
    ${PYPY_BINARY} -u ${TOOL_PLATYPUS_CONFIDENCE_ANNOTATION} \
      --nocontrol \
      --infile=${filenameVCFFinalUnzipped} \
      --tumorColName=${tumor_column} \
      ${CONFIDENCE_OPTS_INDEL} \
      | tee ${filenameVCFTemp} \
      | cut -f 1-11 \
      > ${target}
fi

[[ "$?" != 0 ]] && echo "There was a non-zero exit code in the confidence annotation; temp file ${filenameVCFTemp} not moved back" && exit 6
mv ${filenameVCFTemp} ${filenameVCFFinalUnzipped}

[[ ${runOnPancan-false} == true ]] && mv ${filenameVCFPancanTemp} ${filenameVCFPancanUnzipped}

${BGZIP_BINARY} -f ${filenameVCFFinalUnzipped} && ${TABIX_BINARY} -f -p vcf ${FILENAME_VCF_OUT}
[[ ${runOnPancan-false} == true ]] && ${BGZIP_BINARY} -f ${filenameVCFPancanUnzipped} && ${TABIX_BINARY} -f -p vcf ${filenameVCFPancan}

### Remove temp files:
rm ${OUTPUT_DIRECTORY}/*Annovar*

touch ${FILENAME_CHECKPOINT}
