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

${PERL_BINARY} ${TOOL_CHECK_SAMPLE_SWAP_SCRIPT} \
    --pid=${PID} \
    --raw_file=${FILENAME_VCF_RAW} \
    --annotate_vcf=${TOOL_ANNOTATE_VCF_FILE} \
    --gnomAD_commonSNV=${GNOMAD_V2_1_GENOME_SNV_INDEL} \
    --localControl_commonSNV=${LOCAL_CONTROL_2018_SNV_INDEL} \
    --bias_script=${TOOL_STRAND_BIAS_FILTER_PYTHON_FILE} \
    --tumor_bam=${FILENAME_TUMOR_BAM} \
    --control_bam=${FILENAME_CONTROL_BAM} \
    --reference=${REFERENCE_GENOME} \
    --TiN_R_script=${TOOL_TUMOR_IN_NORMAL_PLOT_RSCRIPT} \
    --canopyFunction=${TOOL_CANOPY_CLUSTER_FUNCTION_RSCRIPT} \
    --chrLengthFile=${TOOL_CANOPY_LINEAR_CHR_DATA} \
    --normal_header_col=${VCF_NORMAL_HEADER_COL} \
    --tumor_header_col=${VCF_TUMOR_HEADER_COL} \
    --sequenceType=${SEQUENCE_TYPE} \
    --exome_capture_kit_bed=${EXOME_CAPTURE_KIT_BEDFILE} \
    --bedtools_binary=${BEDTOOLS_BINARY} \
    --TiNDA_rightBorder=${TINDA_RIGHT_BORDER} \
    --TiNDA_bottomBorder=${TINDA_BOTTOM_BORDER} \
    --TiNDA_runRscript=${TINDA_RUN_RSCRIPT} \
    --maf_thershold=${TINDA_MAX_MAF_CUTOFF} 

### Check the perl run was success or not
if [[ $? == 0 ]] 
then
  touch ${FILENAME_CHECKPOINT}
  exit 0
else
  exit 1
fi
