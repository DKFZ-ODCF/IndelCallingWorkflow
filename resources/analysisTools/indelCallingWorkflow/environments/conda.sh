#!/usr/bin/env bash

# Copyright (c) 2017
# Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).
# Load a Conda environment.

createCleanCondaEnvironment() {
    unset LD_LIBRARY_PATH
    unset MODULE_PATH
    unset PKG_CONFIG_PATH
    unset PYTHONHOME
    unset PYTHONPATH
    unset PYTHONSTARTUP
    unset PYTHON_LIB
    unset PBS_SHARED_BIN
    unset PERL5LIB
    unset PERL_LOCAL_LIB_ROOT
    unset PERL_MB_OPT
    unset PERL_MM_OPT
    export PATH="$HOME/miniconda3/bin:/opt/torque/bin:/usr/lib64/mpi/gcc/openmpi/bin:/opt/maui/bin:/usr/local/bin:/usr/bin:/bin:/usr/bin/X11:/usr/X11R6/bin:/usr/lib/mit/bin"
}

createCleanCondaEnvironment


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

pypy-c() {
    pypy -c "$@"
}
export -f pypy-c
export PYPY_BINARY=pypy-c
