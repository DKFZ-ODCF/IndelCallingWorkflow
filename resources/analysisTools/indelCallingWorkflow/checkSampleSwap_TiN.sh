#!/bin/bash
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/IndelCallingWorkflow).
#
## Bam header analysis
source ${TOOL_ANALYZE_BAM_HEADER}
getRefGenomeAndChrPrefixFromHeader ${FILENAME_TUMOR_BAM} # Sets CHR_PREFIX and REFERENCE_GENOME

## Getting tumor and control header
VCF_TUMOR_HEADER_COL=`samtools view -H ${FILENAME_TUMOR_BAM} | grep -P "^@RG" | perl -ne 'chomp; @s=split(/\t/, $_) ; map{if($_=~/SM:/){$_=~/SM:(.*)/; print "$1\n"}} @s;' | sort | uniq`
VCF_NORMAL_HEADER_COL=`samtools view -H ${FILENAME_CONTROL_BAM} | grep -P "^@RG" | perl -ne 'chomp; @s=split(/\t/, $_) ; map{if($_=~/SM:/){$_=~/SM:(.*)/; print "$1\n"}} @s;' | sort | uniq`

if [[ -z $VCF_TUMOR_HEADER_COL ]] 
then
  VCF_TUMOR_HEADER_COL=`basename ${FILENAME_TUMOR_BAM} | sed 's/.bam$//'`
fi

if [[ -z $VCF_NORMAL_HEADER_COL ]]
then
  VCF_NORMAL_HEADER_COL=`basename ${FILENAME_CONTROL_BAM} | sed 's/.bam$//'`
fi

LOGFILE=$(dirname "$FILENAME_RARE_GERMLINE")/checkSampleSwap_TiN.log

${PERL_BINARY} ${TOOL_CHECK_SAMPLE_SWAP_SCRIPT} \
    --pid=${PID} \
    --raw_file=${FILENAME_VCF_RAW} \
    --annotate_vcf=${TOOL_ANNOTATE_VCF_FILE} \
    --gnomAD_genome=${GNOMAD_WGS_SNV_INDEL} \
    --gnomAD_exome=${GNOMAD_WES_SNV_INDEL} \
    --localControl_WGS=${LOCALCONTROL_PLATYPUS_WGS_SNV_INDEL} \
    --localControl_WES=${LOCALCONTROL_PLATYPUS_WES_SNV_INDEL} \
    --split_mnps_script=${TOOL_SPLIT_MNPS_SCRIPT} \
    --tumor_bam=${FILENAME_TUMOR_BAM} \
    --control_bam=${FILENAME_CONTROL_BAM} \
    --reference=${REFERENCE_GENOME} \
    --chr_prefix=${CHR_PREFIX} \
    --TiN_R_script=${TOOL_TUMOR_IN_NORMAL_PLOT_RSCRIPT} \
    --canopyFunction=${TOOL_CANOPY_CLUSTER_FUNCTION_RSCRIPT} \
    --chrLengthFile=${CHROM_SIZES_FILE} \
    --normal_header_col=${VCF_NORMAL_HEADER_COL} \
    --tumor_header_col=${VCF_TUMOR_HEADER_COL} \
    --sequenceType=${SEQUENCE_TYPE} \
    --gene_model_bed=${GENE_MODEL_BEDFILE} \
    --TiNDA_rightBorder=${TINDA_RIGHT_BORDER} \
    --TiNDA_bottomBorder=${TINDA_BOTTOM_BORDER} \
    --maf_thershold=${TINDA_MAX_MAF_CUTOFF} \
    --outfile_tindaVCF=${FILENAME_TINDA_VCF} \
    --outfile_swapJson=${FILENAME_SWAP_JSON} \
    2>&1 | tee "$LOGFILE"

### Check whether Perl run was success or not
if [[ $? -eq 0 ]]
then
    touch "$FILENAME_CHECKPOINT_SWAP"
    exit 0
elif grep -P 'Less than \d+ rare germline variants' "$LOGFILE"; then
    touch "$FILENAME_CHECKPOINT_SWAP"
    exit 0
else
    exit 1
fi
