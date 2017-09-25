#!/usr/bin/env perl

###################################################################
#
# Author: Naveed Ishaque
# Date: 6th March 2017 (v1.0)
#
# Description: this script takes one input indel VCF file (indel_$pid_somatic_indels_conf_8_to_10.vcf) and prints a json file for the OTP database
#
# Output fields should include (for all): 
# - Samples details : file, pid
# - Indel stats     : numIndels, numIns, numDels, ratioInsToDels, numSize1_3, numSize4_10, numSize11plus, percentSize1_3, percentSize4_10, percentSize11plus, 
#
# CHANGELOG
#
# v1.1
# Bugix - the last entry had a comma
#
###################################################################

# Libraries
use strict;

# Hardcoded variables
my $ref="REF";
my $alt="ALT";
my $confidence="CONFIDENCE";

# Usage
my $usage="\nThis script produces a json file for indel statistics.\n\n\t$0 [indel_PID_somatic_indels_conf_8_to_10.vcf] > indel.json\n\n";

###################################################################

# PARSE INPUT

# Accept input file name
my $indel_file=shift or die "\nERROR: missing input INDEL file.\n$usage";
chomp $indel_file;

# Check if INDEL file has expected name
my $pid;
my $path;
if ($indel_file =~ m/^(.*?)indel_(.*?)_somatic_indels_conf_8_to_10.vcf$/){
  $path = $1;
  $pid = $2;
}
else{
  warn "WARNING: INDEL file did not have expected name format ($indel_file)... please script may not work as intended!\n";
}

# Read file
my $indel_file_handle;
my $numIndels=0;
my $numIns=0;
my $numDels=0;
my $numSize1_3=0;
my $numSize4_10=0;
my $numSize11plus=0;
my $numInsSize1_3=0;
my $numInsSize4_10=0;
my $numInsSize11plus=0;
my $numDelsSize1_3=0;
my $numDelsSize4_10=0;
my $numDelsSize11plus=0;
my %header_colNums;
my $seen_header=0;
my $open_statement="<$indel_file";
$open_statement="zcat $indel_file |" if $indel_file =~ m/gz$/;
open ($indel_file_handle, "$open_statement") or die "ERROR: cannot open file \"$indel_file\".\n$usage";
while (<$indel_file_handle>){
  chomp;
  if (/^\#\#/){
    # vcf header
  }
  elsif(/^\#/){
    # header line
    my @header_array = split ("\t", $_);
    @header_colNums{@header_array} = (0..$#header_array);
    $seen_header=1;
  }
  else{
    if ($seen_header==0) {
      die "ERROR: no header line for \"$indel_file\".\n$usage";
    }
    my @indel_line = split("\t", $_);
    # assume no need to check if the format of the indel file
    $numIndels++;
    my $sizeRef=length($indel_line[$header_colNums{$ref}]);
    my $sizeAlt=length($indel_line[$header_colNums{$alt}]);
    my $size;
    $size = $sizeAlt - $sizeRef if ($sizeAlt>$sizeRef);
    $size = $sizeRef - $sizeAlt if ($sizeRef>$sizeAlt);
    $numSize1_3++    if ($size <=3);
    $numSize4_10++   if ($size >=4 && $size <=10);
    $numSize11plus++ if ($size >=11);
    if ($sizeAlt>$sizeRef){
      $numIns++;
      $numInsSize1_3++    if ($size <=3);
      $numInsSize4_10++   if ($size >=4 && $size <=10);
      $numInsSize11plus++ if ($size >=11);
    }
    else{
      $numDels++;
      $numDelsSize1_3++    if ($size <=3);
      $numDelsSize4_10++   if ($size >=4 && $size <=10);
      $numDelsSize11plus++ if ($size >=11);
    }
    # indel subtypes are not supported
    # if ($function eq "exonic") ...
    # if ($classification eq "frameshift insertion") ...
  }
}
close($indel_file_handle);

###################################################################

# PRINT JSON FILE

# make basic tests or die
if ($numIndels ne ($numIns + $numDels)) {die "ERROR: number of insertions ($numIns) + deletions ($numDels) != number of indels ($numIndels) for \"$indel_file\".\n$usage"};
if ($numIndels ne ($numSize1_3 + $numSize4_10 + $numSize11plus)) {die "ERROR: number of indels of size 1-3 ($numSize1_3) + 4-10 ($numSize4_10) + >11 ($numSize11plus) != number of indels ($numIndels) for \"$indel_file\".\n$usage"};
if ($numIns ne ($numInsSize1_3 + $numInsSize4_10 + $numInsSize11plus)) {die "ERROR: number of insertions of size 1-3 ($numInsSize1_3) + 4-10 ($numInsSize4_10) + >11 ($numInsSize11plus) != number of insertions ($numIns) for \"$indel_file\".\n$usage"};
if ($numDels ne ($numDelsSize1_3 + $numDelsSize4_10 + $numDelsSize11plus)) {die "ERROR: number of deletions of size 1-3 ($numDelsSize1_3) + 4-10 ($numDelsSize4_10) + >11 ($numDelsSize11plus) != number of deletions ($numDels) for \"$indel_file\".\n$usage"};

# print json
print "{\n";
print "  \"all\": {\n";
print "    \"file\": \"$indel_file\",\n";
print "    \"numIndels\": $numIndels,\n";
print "    \"numIns\": $numIns,\n";
print "    \"numDels\": $numDels,\n";
print "    \"numSize1_3\": $numSize1_3,\n";
print "    \"numSize4_10\": $numSize4_10,\n";
print "    \"numSize11plus\": $numSize11plus,\n";
print "    \"numInsSize1_3\": $numInsSize1_3,\n";
print "    \"numInsSize4_10\": $numInsSize4_10,\n";
print "    \"numInsSize11plus\": $numInsSize11plus,\n";
print "    \"numDelsSize1_3\": $numDelsSize1_3,\n";
print "    \"numDelsSize4_10\": $numDelsSize4_10,\n";
print "    \"numDelsSize11plus\": $numDelsSize11plus,\n";
# avoid DIV 0 error
$numIndels=1 if ($numIndels==0);
print "    \"percentIns\": ".(100*$numIns/$numIndels).",\n";
print "    \"percentDels\": ".(100*$numDels/$numIndels).",\n";
print "    \"percentSize1_3\": ".(100*$numSize1_3/$numIndels).",\n";
print "    \"percentSize4_10\": ".(100*$numSize4_10/$numIndels).",\n";
print "    \"percentSize11plus\": ".(100*$numSize11plus/$numIndels).",\n";
print "    \"percentInsSize1_3\": ".(100*$numInsSize1_3/$numIndels).",\n";
print "    \"percentInsSize4_10\": ".(100*$numInsSize4_10/$numIndels).",\n";
print "    \"percentInsSize11plus\": ".(100*$numInsSize11plus/$numIndels).",\n";
print "    \"percentDelsSize1_3\": ".(100*$numDelsSize1_3/$numIndels).",\n";
print "    \"percentDelsSize4_10\": ".(100*$numDelsSize4_10/$numIndels).",\n";
print "    \"percentDelsSize11plus\": ".(100*$numDelsSize11plus/$numIndels)."\n";
print "  }\n";
print "}\n";
