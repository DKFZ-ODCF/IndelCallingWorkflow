#!/usr/bin/env perl
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/IndelCallingWorkflow).
#
#############
## Author: Nagarajan Paramasivam
## Program to 
###  1. Check for Tumor-Control sample swap from same individual
###  2. Check for Tumor in Control from sample individual (TiN)
### 
############
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use JSON::Create 'create_json';

### Input Files and parameters and paths ############################################
my ($pid, $rawFile, $ANNOTATE_VCF, $DBSNP, $tumorBAM, $controlBAM, $ref, 
  $gnomAD_genome, $gnomAD_exome, $split_mnps_script, $localControl_wgs, $localControl_wes,
  $TiN_R, $chrLengthFile, $normal_header_pattern, $tumor_header_pattern, $geneModel,
  $canopy_Function, $seqType, $captureKit, $bedtoolsBinary, $rightBorder, $chr_prefix,
  $bottomBorder, $outfile_tindaVCF, $outfile_SJ);

# Filtering setting to assign common or rare variants
my $AF_cutoff;

GetOptions ("pid=s"                      => \$pid,
            "raw_file=s"                 => \$rawFile,
            "annotate_vcf=s"             => \$ANNOTATE_VCF, 
            "gnomAD_genome=s"            => \$gnomAD_genome,
            "gnomAD_exome=s"             => \$gnomAD_exome,
            "localControl_WGS=s"         => \$localControl_wgs,
            "localControl_WES=s"         => \$localControl_wes,
            "split_mnps_script=s"        => \$split_mnps_script,            
            "tumor_bam=s"                => \$tumorBAM,
            "control_bam=s"              => \$controlBAM,
            "reference=s"                => \$ref,
            "chr_prefix=s"               => \$chr_prefix,
            "TiN_R_script=s"             => \$TiN_R,
            "canopyFunction=s"           => \$canopy_Function,
            "chrLengthFile=s"            => \$chrLengthFile,
            "normal_header_col=s"        => \$normal_header_pattern,
            "tumor_header_col=s"         => \$tumor_header_pattern,
            "sequenceType=s"             => \$seqType,
            "gene_model_bed=s"           => \$geneModel,
            "TiNDA_rightBorder:f"        => \$rightBorder,
            "TiNDA_bottomBorder:f"       => \$bottomBorder,
            "maf_thershold:f"            => \$AF_cutoff,
            "outfile_tindaVCF:s"         => \$outfile_tindaVCF,
            "outfile_swapJSON:s"         => \$outfile_SJ)
or die("Error in SwapChecker input parameters");

die("ERROR: PID is not provided\n") unless defined $pid;
die("ERROR: Raw vcf file is not provided\n") unless defined $rawFile;
die("ERROR: annotate_vcf.pl script path is missing\n") unless defined $ANNOTATE_VCF;
die("ERROR: gnomAD genome is not provided\n") unless defined $gnomAD_genome;
die("ERROR: gnomAD exome is not provided\n") unless defined $gnomAD_exome;
die("ERROR: Tumor bam is missing\n") unless defined $tumorBAM;
die("ERROR: Control bam is missing\n") unless defined $controlBAM;
die("ERROR: Genome reference file is missing\n") unless defined $ref;

# With fill path, filename in annotation
my $analysisBasePath           = dirname $rawFile;
my $snvsGT_RawFile             = $analysisBasePath."/snvs_${pid}.GTfiltered_raw.vcf"; 
my $snvsGT_gnomADFile          = $analysisBasePath."/snvs_${pid}.GTfiltered_gnomAD.vcf";
my $snvsGT_somatic             = $analysisBasePath."/snvs_${pid}.GTfiltered_gnomAD.SomaticIn.vcf";
my $snvsGT_germlineRare        = $analysisBasePath."/snvs_${pid}.GTfiltered_gnomAD.Germline.Rare.vcf";
my $snvsGT_germlineRare_txt    = $analysisBasePath."/snvs_${pid}.GTfiltered_gnomAD.Germline.Rare.txt";
my $snvsGT_germlineRare_png    = $analysisBasePath."/snvs_${pid}.GTfiltered_gnomAD.Germline.Rare.Rescue.png";
my $snvsGT_germlineRare_oFile  = $analysisBasePath."/snvs_${pid}.GTfiltered_gnomAD.Germline.Rare.Rescue.txt";
my $snvsGT_germlineRare_oVCF   = $outfile_tindaVCF;
my $jsonFile                   = $outfile_SJ; # checkSwap.json

###########################################################################################
### For JSON file

my %json = (  
  pid => $pid,
  somaticSmallVarsInTumor => 0,
  somaticSmallVarsInControl => 0,
  germlineSmallVarsInBoth => 0,  
  germlineSmallVarsInBothRare => 0,
  somaticSmallVarsInTumorCommonInGnomad => 0,
  somaticSmallVarsInTumorCommonInGnomadPer => 0,
  somaticSmallVarsInControlCommonInGnomad => 0,
  somaticSmallVarsInControlCommonInGnomadPer => 0,
  somaticSmallVarsInTumorPass => 0,
  somaticSmallVarsInTumorPassPer => 0,
  somaticSmallVarsInControlPass => 0,
  somaticSmallVarsInControlPassPer => 0,
  tindaGermlineRareAfterRescue => 0,
  tindaSomaticAfterRescue => 0,
  tindaSomaticAfterRescueMedianAlleleFreqInControl => 0
);


###########
### WES vs WGS 
# If WES filter for the seq kit
my $updated_rawFile;
# update in WES
if($seqType eq 'WES') {
   $updated_rawFile = $rawFile.".intersect.gz";
   `bedtools slop -i $geneModel -b 5 -g $chrLengthFile | \
    cut -f1-3 | awk '{if(\$3<\$2){print \$1"\t"\$3"\t"\$2}else{print \$0}}' | \
    bedtools merge -i - | \
    bedtools intersect -header -a $rawFile -b - | \
    bgzip -f > $updated_rawFile && tabix -f -p vcf $updated_rawFile`;
   #$updated_rawFile = $rawFile;
}
elsif($seqType eq 'WGS') {
   $updated_rawFile = $rawFile;
   #`zcat $rawFile | python $split_mnps_script | bgzip -f > $updated_rawFile; tabix -f -p vcf $updated_rawFile`;
}

open(my $IN, 'zcat '. $updated_rawFile.'| ')
    || die "Can't read in the '$updated_rawFile'\n";
open(JSON, ">$jsonFile")
    || die "Can't create the '$jsonFile'\n";

## Filtering for Somatic variants and germline based on platypus genotype
my ($controlCol, $tumorCol, $formatCol);

my $columnCounter;

## Creating tumor and control raw somatic snvs files 
open(GTraw, ">$snvsGT_RawFile")
    || die "Can't create the '$snvsGT_RawFile'\n";

while(!eof($IN)) {
  my $line = readline($IN)
      || die("Error reading from zcatted '$updated_rawFile': $!");
  chomp $line;
   
  if($line =~ /^#/)  {
    # Headers
    if($line =~ /^#CHROM/) {
      my @header = split(/\t/, $line);

      $columnCounter = $#header;

      for(my $i=0; $i<=$#header; $i++) {        
        $controlCol = $i, if($header[$i] =~ /$normal_header_pattern/i);
        $tumorCol = $i, if($header[$i] =~ /$tumor_header_pattern/i);
        $formatCol = $i, if($header[$i] =~ /FORMAT/);
      }

      if($controlCol =~ /^$/ || $tumorCol =~ /^$/) {
        # stop if header doesn't control or tumor in the column header
        print JSON create_json (\%json);
        print "Normal header patter provided : $normal_header_pattern\n";
        print "Tumor header patter provided : $tumor_header_pattern\n";
        die("$pid doesn't have control-tumor pair info in the column header\n");
      }
      else {
        print GTraw "$line\tControl_AF\tTumor_AF\tTumor_dpALT\tTumor_dp\tControl_dpALT\tControl_dp\tGT_Classification\n";
      }
    }
    else {
      ## Rest of the header rows
      print GTraw "$line\n";
    }
  }
  else {
    # Variants
    my @variantInfos = split(/\t/, $line);
    my $filter = $variantInfos[6];
    my $chr    = $variantInfos[0];
    my $ref    = $variantInfos[3];
    my @alt    = split(/,/, $variantInfos[4]);
    my @control = split(/:/, $variantInfos[$controlCol]);
    my @tumor   = split(/:/, $variantInfos[$tumorCol]);
    my @format  = split(/:/, $variantInfos[$formatCol]);

    my ($iGT, $iGQ, $iPL, $iNV, $iDP);
    for(my $i=0; $i<=$#format; $i++) {
      if($format[$i] eq "GT"){$iGT=$i}
      if($format[$i] eq "GQ"){$iGQ=$i}
      if($format[$i] eq "PL" || $format[$i] eq "GL"){$iPL=$i}
      if($format[$i] eq "NV"){$iNV=$i}
      if($format[$i] eq "NR"){$iDP=$i}
    }

    # Removing extra chr contigs, Indels and bad quality snvs
    # Including both indels and snvs - removed as we will have issue with bias Filter
    if($chr=~/^($chr_prefix)?(X|Y|[1-9]|1[0-9]|2[0-2])$/ && $filter =~/^(PASS|alleleBias)$/) {


      my @tumor_dp = split(/,/, $tumor[$iDP]);
      my @control_dp = split(/,/, $control[$iDP]);
      my @tumor_nv = split(/,/, $tumor[$iNV]);
      my @control_nv = split(/,/, $control[$iNV]);

      for(my $i=0;$i<=$#alt; $i++) {

        if($tumor_dp[$i] >= 5 && $control_dp[$i] >= 5) {

          $variantInfos[4] = $alt[$i] ;
          my $newLine = join("\t", @variantInfos) ;

          my $tumor_AF = $tumor_nv[$i]/$tumor_dp[$i] ;
          my $control_AF = $control_nv[$i]/$control_dp[$i] ;

          if($tumor_nv[$i] > 3 && $control_AF == 0) {
            print GTraw "$newLine\t$control_AF\t$tumor_AF\t$tumor_nv[$i]\t$tumor_dp[$i]\t$control_nv[$i]\t$control_dp[$i]\tTumor_Somatic\n";
            $json{'somaticSmallVarsInTumor'}++;
          }
          elsif($tumor_AF == 0 && $control_nv[$i] > 3) {
            print GTraw "$newLine\t$control_AF\t$tumor_AF\t$tumor_nv[$i]\t$tumor_dp[$i]\t$control_nv[$i]\t$control_dp[$i]\tControl_Somatic\n";
            $json{'somaticSmallVarsInControl'}++;
          }
          elsif($tumor_nv[$i] > 1 && $control_nv[$i] >= 1) {
            print GTraw "$newLine\t$control_AF\t$tumor_AF\t$tumor_nv[$i]\t$tumor_dp[$i]\t$control_nv[$i]\t$control_dp[$i]\tGermlineInBoth\n";
          }
        }
      }
    } 
  }
}

close GTraw;

### Cleaning up MNPs, due to the MNPs escaping in multi-SNVs lines
#`cat $snvsGT_RawFile | python $split_mnps_script | uniq > $snvsGT_RawFile.temp ; mv -f $snvsGT_RawFile.temp $snvsGT_RawFile`;
my $resolve_complex_variants = "(head -n 2000 '$snvsGT_RawFile' | grep '#' ; 
  cat '$snvsGT_RawFile' | python '$split_mnps_script' | grep -v '#' | sort -V -k1,2) | uniq > '$snvsGT_RawFile'.temp ;
  mv -f '$snvsGT_RawFile'.temp '$snvsGT_RawFile';";
 
print "\n$resolve_complex_variants\n";
my $run_resolve_complex = system($resolve_complex_variants);

if($run_resolve_complex !=0){
  die("ERROR: During cleaning up MNPs\n");
}

## Annotating with gnomAD and local control 
my $annotation_command = "cat '$snvsGT_RawFile' | \
                          perl '$ANNOTATE_VCF' -a - -b '$gnomAD_genome' --columnName='gnomAD_GENOMES' --bAdditionalColumn=2 --reportMatchType --reportLevel 1 | \
                          perl '$ANNOTATE_VCF' -a - -b '$gnomAD_exome' --columnName='gnomAD_EXOMES' --bAdditionalColumn=2 --reportMatchType --reportLevel 1 | \
                          perl '$ANNOTATE_VCF' -a - -b '$localControl_wgs' --columnName='LocalControl_WGS' --bAdditionalColumn=2  --reportMatchType --reportLevel 1 | \
                          perl '$ANNOTATE_VCF' -a - -b '$localControl_wes' --columnName='LocalControl_WES' --bAdditionalColumn=2  --reportMatchType --reportLevel 1 > '$snvsGT_gnomADFile'";

print "\n$annotation_command\n";
my $runAnnotation = system($annotation_command);

if($runAnnotation != 0 ) {
  `rm $jsonFile`;
  die("ERROR: In the allele frequency annotation step\n") ;
}

####### Germline file and rare variant filtering
open(ANN, "<$snvsGT_gnomADFile") || die "cant open the $snvsGT_gnomADFile\n";

open(GermlineRareFile, ">$snvsGT_germlineRare") || die "cant create the $snvsGT_germlineRare\n";
open(GermlineRareFileText, ">$snvsGT_germlineRare_txt") || die "cant create the $snvsGT_germlineRare_txt\n";

print GermlineRareFileText "CHR\tPOS\tREF\tALT\tControl_AF\tTumor_AF\tTumor_dpALT\tTumor_dp\tControl_dpALT\tControl_dp\tRareness\n";

open(SomaticFile, ">$snvsGT_somatic") || die "cant create the $snvsGT_somatic\n";
while(!eof(ANN)) {
  my $annLine = readline(ANN)
      || die "Error reading from '$snvsGT_gnomADFile': $!";
  chomp $annLine;

  # column counters
  my $start_col = $columnCounter+1;
  my $end_col = $columnCounter+6;
  my $gnomAD_genome_col = $columnCounter+8;
  my $gnomAD_exome_col = $columnCounter+9;
  my $localcontrol_wgs_col = $columnCounter+10;
  my $localcontrol_wes_col = $columnCounter+11;

  if($annLine =~ /^#/) {
    print GermlineRareFile "$annLine\n";
    print SomaticFile "$annLine\n";    
  }
  else {
    my @annLineSplit = split(/\t/, $annLine);

    my $germlineTextInfo = join("\t", @annLineSplit[0..1], @annLineSplit[3..4], @annLineSplit[$start_col..$end_col]);
    
    ### rare or common in gnomAD or local control
    my $AF_gnomAD_genome = 0;
    my $AF_gnomAD_exome = 0;
    my $AF_localcontrol_wgs = 0;
    my $AF_localcontrol_wes = 0;
    
    if($annLineSplit[$gnomAD_genome_col] =~/MATCH=(exact|position)/) {
      $AF_gnomAD_genome = parse_AF($annLineSplit[$gnomAD_genome_col] );      
    }
    if($annLineSplit[$gnomAD_exome_col] =~/MATCH=(exact|position)/) {
      $AF_gnomAD_exome = parse_AF($annLineSplit[$gnomAD_exome_col]);      
    }
    if($annLineSplit[$localcontrol_wgs_col] =~ /MATCH=(exact|position)/) {
      $AF_localcontrol_wgs = parse_AF($annLineSplit[$localcontrol_wgs_col]);
    }
    if($annLineSplit[$localcontrol_wes_col] =~ /MATCH=(exact|position)/) {
      $AF_localcontrol_wes = parse_AF($annLineSplit[$localcontrol_wes_col]);
    }


    my $common_rare = "RARE";
    if($AF_gnomAD_genome > $AF_cutoff || $AF_gnomAD_exome > $AF_cutoff || $AF_localcontrol_wgs > $AF_cutoff || $AF_localcontrol_wes) {
      $common_rare  = "COMMON";    
    }
    

    # Somatic rare or common 
    if($annLine=~/_Somatic/) {
      if($common_rare eq "COMMON") {
        $annLine =~ s/_Somatic/_Somatic_Common/;
        print SomaticFile "$annLine\n"; 
      }
      else {
        $annLine =~ s/_Somatic/_Somatic_Rare/;
        print SomaticFile "$annLine\n";
        ## Adding control's rare somatic variants into the germline file
        if($annLine=~/Control_Somatic/){
          print GermlineRareFile "$annLine\n";
          print GermlineRareFileText "$germlineTextInfo\tSomaticControlRare\n";
        }
      }
    }
    # Rare germline variant    
    elsif($annLine =~ /GermlineInBoth/ && $common_rare eq "RARE") {
      $json{'germlineSmallVarsInBothRare'}++;
      print GermlineRareFile "$annLine\n";
      print GermlineRareFileText "$germlineTextInfo\tRare\n";
    }
  }
}

close GermlineRareFile;
close SomaticFile;
close ANN;

### Crashing when there is less than 50 germline variants. Write whatever information has been gathered.
if($json{'germlineSmallVarsInBothRare'} < 50){
  print JSON create_json (\%json);
  close JSON;
  die "Less than 50 rare germline variants. Might be a sample swap or poor coverage on any one of the samples, please check the coverage information.\nExiting the analysis";
}

#######################################
### Finding and plotting TiN

my $runRscript_code  = join("", "Rscript '$TiN_R' -f '$snvsGT_germlineRare_txt'", 
  " --oPlot $snvsGT_germlineRare_png",
  " --oFile $snvsGT_germlineRare_oFile",
  " --pid $pid",
  " --chrLength $chrLengthFile",
  " --cFunction $canopy_Function",
  " --seqType $seqType",
  " --rightBorder $rightBorder",
  " --bottomBorder $bottomBorder",
  " --vcf $snvsGT_germlineRare",
  " --oVcf $snvsGT_germlineRare_oVCF");

print $runRscript_code, "\n";

my $runRscript = system($runRscript_code);

if($runRscript != 0) {
  `rm $jsonFile`;  
  die "Error while running $TiN_R in swapChecker, $runRscript";
}

open(TINDA_rareOutput, $snvsGT_germlineRare_oFile) || die "Can't open the $snvsGT_germlineRare_oFile: $!";
my @SomaticRescue_control_AF;

while(!eof(TINDA_rareOutput)) {
  my $line = readline(TINDA_rareOutput)
      || die "Error reading from '$snvsGT_germlineRare_oFile': $!";
  chomp $line;
  if ($line =~/Germline/) {
    $json{'tindaGermlineRareAfterRescue'}++;
  }
  elsif ($line =~ /Somatic_Rescue/) {
    $json{'tindaSomaticAfterRescue'}++;
    my @gr_ss = split(/\t/, $line);
    push(@SomaticRescue_control_AF, $gr_ss[4]);
  }  
}
close TINDA_rareOutput;

if($json{'tindaSomaticAfterRescue'} > 0) {
  $json{'tindaSomaticAfterRescueMedianAlleleFreqInControl'} = median(@SomaticRescue_control_AF);
} else {
  $json{'tindaSomaticAfterRescueMedianAlleleFreqInControl'} = 0;
}

#######################################
### Counting the somatic numbers 
open(my $SOMATIC_FH, "<$snvsGT_somatic")
    || die "Can't open the file $snvsGT_somatic\n";

while(!eof($SOMATIC_FH)) {
  my $line = readline($SOMATIC_FH)
    || die("Error reading from zcatted '$snvsGT_somatic': $!");
  chomp $line;
  
  if($line!~/^#/) {
    if($line=~/Tumor_Somatic_Common/) {
      $json{'somaticSmallVarsInTumorCommonInGnomad'}++;
    }
    elsif($line=~/Tumor_Somatic_Rare/) {
      $json{'somaticSmallVarsInTumorPass'}++;
    }

    if($line=~/Control_Somatic_Common/) {
      $json{'somaticSmallVarsInControlCommonInGnomad'}++;    
    }
    elsif($line=~/Control_Somatic_Rare/) {
      $json{'somaticSmallVarsInControlPass'}++;
    }   
  }
}
close $SOMATIC_FH;

##################
## Creating sample swap json file 

## Percentage calculations
if($json{'somaticSmallVarsInTumor'} > 0) {
  $json{'somaticSmallVarsInTumorCommonInGnomadPer'} = $json{'somaticSmallVarsInTumorCommonInGnomad'}/$json{'somaticSmallVarsInTumor'};  
  $json{'somaticSmallVarsInTumorPassPer'} = $json{'somaticSmallVarsInTumorPass'}/$json{'somaticSmallVarsInTumor'};
}
else {
  $json{'somaticSmallVarsInTumorCommonInGnomadPer'} = 0;  
  $json{'somaticSmallVarsInTumorPassPer'} = 0;
}

if($json{'somaticSmallVarsInControl'} > 0) {
  $json{'somaticSmallVarsInControlCommonInGnomasPer'} = $json{'somaticSmallVarsInControlCommonInGnomad'}/$json{'somaticSmallVarsInControl'};  
  $json{'somaticSmallVarsInControlPassPer'} = $json{'somaticSmallVarsInControlPass'}/$json{'somaticSmallVarsInControl'};
}
else {
  $json{'somaticSmallVarsInControlCommonInGnomadPer'} = 0;  
  $json{'somaticSmallVarsInControlPassPer'} = 0;

}

print JSON create_json (\%json);
close JSON;

######################################
#### Cleaning up files 
`rm $snvsGT_RawFile $snvsGT_gnomADFile`;

`rm $snvsGT_somatic $snvsGT_germlineRare $snvsGT_germlineRare_txt`;

#####################################
## Subroutine to parse AF from multiple or single match 
sub parse_AF{

  my $line = $_[0];
  my $AF = 0;
  if($line=~/\&/){
    my @lines = split(/&/, $line);
    foreach my $match(@lines) {
      my ($temp_AF) = $match =~ /;AF=(\d(\.\d+(e-\d+)?)?);/;
      if($temp_AF > $AF) {
        $AF = $temp_AF;
      }
    }
  }
  else{
    ($AF) = $line =~ /;AF=(\d(\.\d+(e-\d+)?)?);/;
  }
  return($AF);
}

## Median function
sub median {
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) {
        return $vals[int($len/2)];
    }
    else {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}