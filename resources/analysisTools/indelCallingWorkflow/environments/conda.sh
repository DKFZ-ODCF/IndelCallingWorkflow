#!/usr/bin/env bash

# Copyright (c) 2017
# Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).
# Load a Conda environment.

source activate "${condaEnvironmentName:?No Conda environment name defined. Please set 'condaEnvironmentName'.}" \
    || (echo "Could not load Conda environment '$condaEnvironmentName'" && exit 100)

export GHOSTSCRIPT_BINARY=gs
export PLATYPUS_BINARY=platypus
export PERL_BINARY=perl
export PYTHON_BINARY=python
export BGZIP_BINARY=bgzip
export TABIX_BINARY=tabix
export BEDTOOLS_BINARY=bedtools
export SAMTOOLS_BINARY=samtools
export PYPY_BINARY=pypy
