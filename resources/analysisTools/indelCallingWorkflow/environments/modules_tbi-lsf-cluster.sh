#!/usr/bin/env bash

HTSLIB_VERSION=1.4.1
HTSLIB_VERSION_PLATYPUS=1.5.0

platypus() {
  module load "htslib/$HTSLIB_VERSION_PLATYPUS"
  module load platypus/0.8.1
  Platypus.py "$@"
  module load htslib/1.4.1
}

module load "htslib/$HTSLIB_VERSION"
module load perl/5.20.2
module load python/2.7.9
module load samtools/0.1.19

export PLATYPUS_BINARY=Platypus.py
export GHOSTSCRIPT_BINARY=gs