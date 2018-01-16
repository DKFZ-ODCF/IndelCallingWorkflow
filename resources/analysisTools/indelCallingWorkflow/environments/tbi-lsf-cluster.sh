#!/usr/bin/env bash

export HTSLIB_VERSION=1.4.1
export HTSLIB_VERSION_PLATYPUS=1.3.1
export PLATYPUS_VERSION=0.8.1

platypus() {
  module load "htslib/$HTSLIB_VERSION_PLATYPUS"
  module load "Platypus/$PLATYPUS_VERSION"
  Platypus.py "$@"
  module load "htslib/$HTSLIB_VERSION"
}

export -f platypus

module load "htslib/$HTSLIB_VERSION"
module load "perl/5.20.2"
module load "python/2.7.9"
module load "samtools/0.1.19"
module load "pypy/5.8.0"
module load "R/3.3.1"

export PLATYPUS_BINARY=platypus
export GHOSTSCRIPT_BINARY=gs
export PYTHON_BINARY=python
export PERL_BINARY=perl
export SAMTOOLS_BINARY=samtools
