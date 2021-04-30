#!/usr/bin/env bash
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/IndelCallingWorkflow).
#

# Developed for the LSF cluster. Should also work with the PBS cluster.

module load "perl/5.20.2"
module load "python/2.7.9"
module load "samtools/1.9"
module load "htslib/1.9"
module load "R/3.3.1"
module load "bedtools/2.24.0"
module load "gdc/6.3.0"

export GHOSTSCRIPT_BINARY=gs
export PYTHON_BINARY=python
export PERL_BINARY=perl
export SAMTOOLS_BINARY=samtools
export BGZIP_BINARY=bgzip
export TABIX_BINARY=tabix
export BEDTOOLS_BINARY=bedtools
export RSCRIPT_BINARY=Rscript
