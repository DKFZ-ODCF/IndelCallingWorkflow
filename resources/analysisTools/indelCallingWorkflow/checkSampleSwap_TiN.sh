#!/bin/bash
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/IndelCallingWorkflow).
#
## Bam header analysis
source "$TOOL_ANALYZE_BAM_HEADER"
getRefGenomeAndChrPrefixFromHeader "$FILENAME_TUMOR_BAM" # Sets CHR_PREFIX and REFERENCE_GENOME

## Getting tumor and control header
bamExtractSmTag() {
  local filename="${1:?No BAM file provided}"
  samtools view -H "$filename" \
    | grep -P "^@RG" \
    | perl -ne 'chomp; @s=split(/\t/, $_) ; map{if($_=~/SM:/){$_=~/SM:(.*)/; print "$1\n"}} @s;' \
    | sort -u
}

defaultToBamBasename() {
  local bamFile="${1:?No BAM provided}"
  local name="${2:?No argument}"
  if [[ -z "$name" ]]; then
    basename "$bamFile" .bam
  else
    echo "$name"
  fi
}

getColumnName() {
  local bamFile="${1:?No BAM provided}"
  defaultToBamBasename "$bamFile" "$(bamExtractSmTag "$bamFile")"
}

VCF_TUMOR_HEADER_COL="$(getColumnName "$FILENAME_TUMOR_BAM")"
VCF_NORMAL_HEADER_COL="$(getColumnName "$FILENAME_CONTROL_BAM")"

LOGFILE="$(dirname "$FILENAME_RARE_GERMLINE")/checkSampleSwap_TiN.log"

${PERL_BINARY} ${TOOL_CHECK_SAMPLE_SWAP_SCRIPT} \
    --pid=${PID} \
    --raw_file=${FILENAME_VCF_RAW} \
    --annotate_vcf=${TOOL_ANNOTATE_VCF_FILE} \
    --gnomAD_genome=${GNOMAD_V2_1_GENOME_SNV_INDEL} \
    --gnomAD_exome=${GNOMAD_V2_1_EXOME_SNV_INDEL} \
    --localControl_commonSNV=${LOCAL_CONTROL_PLATYPUS_SNV_INDEL} \
    --split_mnps_script=${TOOL_SPLIT_MNPS_SCRIPT} \
    --bias_script=${TOOL_STRAND_BIAS_FILTER_PYTHON_FILE} \
    --tumor_bam=${FILENAME_TUMOR_BAM} \
    --control_bam=${FILENAME_CONTROL_BAM} \
    --reference=${REFERENCE_GENOME} \
    --TiN_R_script=${TOOL_TUMOR_IN_NORMAL_PLOT_RSCRIPT} \
    --canopyFunction=${TOOL_CANOPY_CLUSTER_FUNCTION_RSCRIPT} \
    --chrLengthFile=${TOOL_CANOPY_LINEAR_CHR_DATA} \
    --normal_header_substr=${VCF_NORMAL_HEADER_COL} \
    --tumor_header_substr=${VCF_TUMOR_HEADER_COL} \
    --sequenceType=${SEQUENCE_TYPE} \
    --exome_capture_kit_bed=${EXOME_CAPTURE_KIT_BEDFILE} \
    --gene_model_bed=${GENE_MODEL_BEDFILE} \
    --TiNDA_rightBorder=${TINDA_RIGHT_BORDER} \
    --TiNDA_bottomBorder=${TINDA_BOTTOM_BORDER} \
    --maf_thershold=${TINDA_MAX_MAF_CUTOFF} \
    --outfile_rareGermline=${FILENAME_RARE_GERMLINE} \
    --outfile_somaticRescue=${FILENAME_SOMATIC_RESCUE} \
    --outfile_allSomatic=${FILENAME_ALL_SOMATIC} \
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
