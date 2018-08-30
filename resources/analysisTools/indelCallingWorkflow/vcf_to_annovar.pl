#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/IndelCallingWorkflow).
#

use strict;
use warnings;
use v5.10;

my ($chr_prefix, $chr_suffix);
BEGIN {
  $chr_prefix = shift;
  $chr_suffix = shift;
  # Now make constants
  use constant;
  constant->import("CHR_PFX", $chr_prefix);
  constant->import("CHR_SFX", $chr_suffix);
}
my @fields;
my $end;
my $alt;
my ($pos_orig, $end_orig);

while (<>) {
  next if (/^\#/);
  @fields = split(/\t/);
  if (CHR_PFX) {
    $fields[0] =~ s/$chr_prefix//;
  }
  if (CHR_SFX) {
    $fields[0] =~ s/$chr_suffix//;
  }


  $alt = (split(',',$fields[4]))[0]; # in case of multiallelic variants take only first allele

  $alt = uc($alt);
  $fields[3] = uc($fields[3]);

  $pos_orig = $fields[1];
  $end_orig = $fields[1] + length($fields[3]) -1;

  if (length($fields[3]) > 1 || length($alt) > 1) { # indel: vcf reports base before event in ref and alt, annovar does not like this (otherwise detects everything as substitution)
    if (substr($fields[3],0,1) eq substr($alt,0,1)) { # should always be the case for vcf-conform indel notation, but you never know...
      $fields[3] = substr($fields[3],1); # remove first base from ref
      $alt = substr($alt,1); # remove first base from alt
      $fields[1] += 1; # adjust pos coordinate
    } else {
      warn "Found indel where 1st ref base does not equal 1st alt base:\n$_";
    }
  }

  $fields[3] = '-' if (! length($fields[3]));
  $alt = '-' if (! length($alt));

  $end = $fields[1] + length($fields[3]) - 1; # annovar has a strange understanding of the end coordinate: this works if there is '-' or one base or multiple bases in ref

  say join "\t", $fields[0],$fields[1],$end,$fields[3],$alt, $pos_orig, $end_orig;
}


### End confusion: here end is according to vcf standard (1-based, including)
