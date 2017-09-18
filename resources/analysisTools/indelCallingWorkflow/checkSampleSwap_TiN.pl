#!/usr/bin/perl 
############ 
## Author: Nagarajan Paramasivam
## Program to 
###  1. Check for Tumor-Control sample swap from same individual
###  2. Check for Tumor in Control from sample individual (TiN)
### 
############
use strict;
use File::Basename;
use Getopt::Long;
use JSON::Create 'create_json';

### Input Files and parameters and paths ############################################
my ($pid, $rawFile, $ANNOTATE_VCF, $DBSNP, $biasScript, $tumorBAM, $controlBAM, $ref, $gnomAD, $TiN_R, $localControl, $chrLengthFile, $normal_header_pattern, $tumor_header_pattern, $localControl_2, $canopy_Function, $seqType, $captureKit, $bedtoolsBinary);

GetOptions ("pid=s"                      => \$pid,
            "raw_file=s"                 => \$rawFile,
            "annotate_vcf=s"             => \$ANNOTATE_VCF, 
            "gnomAD_commonSNV=s"         => \$gnomAD,
            "localControl_commonSNV=s"   => \$localControl,
            "localControl_commonSNV_2=s" => \$localControl_2,
            "bias_script=s"              => \$biasScript,
            "tumor_bam=s"                => \$tumorBAM,
            "control_bam=s"              => \$controlBAM,
            "reference=s"                => \$ref,
            "TiN_R_script=s"             => \$TiN_R,
            "canopyFunction=s"           => \$canopy_Function,
            "chrLengthFile=s"            => \$chrLengthFile,
            "normal_header_col=s"        => \$normal_header_pattern,
            "tumor_header_col=s"         => \$tumor_header_pattern,
            "sequenceType=s"             => \$seqType,
            "exome_capture_kit_bed=s"    => \$captureKit,
            "bedtools2_24_binary=s"      => \$bedtoolsBinary)
or die("Error in SwapChecker input parameters");

die("ERROR: PID is not provided\n") unless defined $pid;
die("ERROR: Raw vcf file is not provided\n") unless defined $rawFile;
die("ERROR: annotate_vcf.pl script path is missing\n") unless defined $ANNOTATE_VCF;
die("ERROR: gnomAD common SNVs is not provided\n") unless defined $gnomAD;
die("ERROR: strand bias script path is missing\n") unless defined $ANNOTATE_VCF;
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
my $snvsGT_somaticRareBiasFile = $analysisBasePath."/snvs_${pid}.GTfiltered_gnomAD.SomaticIn.Rare.BiasFiltered.vcf";
my $jsonFile                    = $analysisBasePath."/checkSampleSwap.json"; # checkSwap.json

###########################################################################################
### For JSON file

my %json = (  
  pid => $pid,
  SomaticSmallVarsInTumor => 0,
  SomaticSmallVarsInControl => 0,
  GermlineSNVsHeterozygousInBoth => 0,  
  GermlineSNVsHeterozygousInBothRare => 0,
  SomaticSmallVarsInTumorCommonInGnomad => 0,
  SomaticSmallVarsInTumorCommonInGnomadPer => 0,
  SomaticSmallVarsInControlCommonInGnomad => 0,
  SomaticSmallVarsInControlCommonInGnomadPer => 0,
  SomaticSmallVarsInTumorInBias => 0,
  SomaticSmallVarsInTumorInBiasPer => 0,
  SomaticSmallVarsInControlInBias => 0,
  SomaticSmallVarsInControlInBiasPer => 0,
  SomaticSmallVarsInTumorPass => 0,
  SomaticSmallVarsInTumorPassPer => 0,
  SomaticSmallVarsInControlPass => 0,
  SomaticSmallVarsInControlPassPer => 0,
  TindaGermlineRareAfterResuce => 0,
  TindaSomaticAfterResuce => 0,
  TindaSomaticAfterResuceMedianAlleleFreqInControl => 0
);

###########
### WES vs WGS 
# If WES filter for the seq kit
my $updated_rawFile;
# update in WES
if($seqType eq 'WES') {
   $updated_rawFile = $rawFile.".intersect.gz";
  `$bedtoolsBinary intersect -header -a $rawFile -b $captureKit | bgzip > $updated_rawFile ; tabix -p vcf $updated_rawFile`;
}
elsif($seqType eq 'WGS') {
   $updated_rawFile = $rawFile;
}

open(my $IN, 'zcat '. $updated_rawFile.'| ') || die "Cant read in the $updated_rawFile\n";
open(JSON, ">$jsonFile") || die "Can't craete the $jsonFile\n";

## Filtering for Somatic variants and germline based on platypus genotype
my ($controlCol, $tumorCol, $formatCol);

my $columnCounter;

## Creating tumor and control raw somatic snvs files 
open(GTraw, ">$snvsGT_RawFile") || die "Can't create the $snvsGT_RawFile\n";

while(<$IN>) {
  chomp;
  my $line = $_;  
   
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
    if($chr=~/^(X|Y|[1-9]|1[0-9]|2[0-2])$/ && $filter =~/^(PASS|alleleBias)$/ && $variantInfos[4] != /,/) {


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
            $json{'SomaticSmallVarsInTumor'}++;
          }
          elsif($tumor_AF == 0 && $control_nv[$i] > 3) {
            print GTraw "$newLine\t$control_AF\t$tumor_AF\t$tumor_nv[$i]\t$tumor_dp[$i]\t$control_nv[$i]\t$control_dp[$i]\tControl_Somatic\n";
            $json{'SomaticSmallVarsInControl'}++;
          }
          elsif($tumor_AF > 0 && $control_AF > 0) {
            print GTraw "$newLine\t$control_AF\t$tumor_AF\t$tumor_nv[$i]\t$tumor_dp[$i]\t$control_nv[$i]\t$control_dp[$i]\tGermlineInBoth\n";
          }
        }
      }
    } 
  }
}

close GTraw;

## Annotating with dbSNP

my $runAnnotation = system("cat $snvsGT_RawFile | perl $ANNOTATE_VCF -a - -b '$gnomAD' --columnName='gnomAD_COMMON_SNV' --reportMatchType  --bAdditionalColumn=2 | perl $ANNOTATE_VCF -a - -b '$localControl' --columnName='LocalControl_COMMON_SNV' --reportMatchType --bAdditionalColumn=2 | perl $ANNOTATE_VCF -a - -b '$localControl_2' --columnName='LocalControl2_COMMON_SNV' --reportMatchType --bAdditionalColumn=2 > $snvsGT_gnomADFile");


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

while(<ANN>) {
  chomp;
  my $annLine = $_;
  if($annLine =~ /^#/) {
    print GermlineRareFile "$annLine\n";
    print SomaticFile "$annLine\n";    
  }
  else {
    my @annLineSplit = split(/\t/, $annLine);
    my $start_col = $columnCounter+1;
    my $end_col = $columnCounter+6;
    my $gnomAD_col = $columnCounter+8; 

    my $germlineTextInfo = join("\t", @annLineSplit[0..1], @annLineSplit[3..4], @annLineSplit[$start_col..$end_col]);

    if($annLine=~/_Somatic/) {
      if($annLine=~/MATCH/) {
        $annLine =~ s/_Somatic/_Somatic_Common/;
        print SomaticFile "$annLine\n"; 
      }
      else {
        $annLine =~ s/_Somatic/_Somatic_Rare/;
        print SomaticFile "$annLine\n";
      }
    }
    elsif($annLine=~/Germline/) {
      if($annLine =~ /GermlineInBoth/ && $annLine !~ /MATCH/ && $seqType eq 'WGS') {
        $json{'GermlineSmallVarsHeterozygousInBothRare'}++;
        print GermlineRareFile "$annLine\n";
	print GermlineRareFileText "$germlineTextInfo\tRare\n";
      }
      elsif($annLine =~ /GermlineInBoth/ && $seqType eq 'WES') {
        my $rareness;
        if($annLine !~ /MATCH/) {
          $json{'GermlineSmallVarsHeterozygousInBothRare'}++;
          $rareness = "Rare";
        }
        else {
          $rareness = "Common";
        }
        print GermlineRareFile "$annLine\n";
        print GermlineRareFileText "$germlineTextInfo\t$rareness\n";
      }  
    }
  }
}

close GermlineRareFile;
close SomaticFile;
close Ann;

#######################################
### Finding and plotting TiN

#print "Rscript-3.3.1 $TiN_R -f $snvsGT_germlineRare_txt --oPlot $snvsGT_germlineRare_png --oFile $snvsGT_germlineRare_oFile -p $pid --chrLength $chrLengthFile --cFunction $canopy_Function --SeqType $seqType\n" ;

my $runRscript = system("Rscript-3.3.1 $TiN_R -f $snvsGT_germlineRare_txt --oPlot $snvsGT_germlineRare_png --oFile $snvsGT_germlineRare_oFile -p $pid --chrLength $chrLengthFile --cFunction $canopy_Function --SeqType $seqType" ) ;

if($runRscript != 0) { 
  `rm $jsonFile`;
  die "Error while running $TiN_R in swapChecker\n";
}
 
chomp($json{'TindaGermlineRareAfterResuce'} = `cat $snvsGT_germlineRare_oFile | grep 'Germline' | wc -l`);
chomp($json{'TindaSomaticAfterResuce'}  = `cat $snvsGT_germlineRare_oFile | grep 'Somatic_Rescue' | wc -l`);

if($json{'TindaSomaticAfterResuce'} > 0) {

  $json{'TindaSomaticAfterResuceMedianAlleleFreqInControl'} = `cat $snvsGT_germlineRare_oFile | grep 'Somatic_Rescue' | cut -f5 | sort -n | awk ' { a[i++]=\$1;} END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`;

}
else
{
  $json{'TindaSomaticAfterResuceMedianAlleleFreqInControl'} = 0;
}

#######################################
## Running Bias Filters

my $runBiasScript = system("python $biasScript $snvsGT_somatic $tumorBAM $ref $snvsGT_somaticRareBiasFile --tempFolder $analysisBasePath --maxOpRatioPcr=0.34 --maxOpRatioSeq=0.34 --maxOpReadsPcrWeak=2 --maxOpReadsPcrStrong=2");

if($runBiasScript !=0) {
  die "Error while running $biasScript in swapChecker\n";
}


### Counting The Numbers 
open(SOM_RareBias, "<$snvsGT_somaticRareBiasFile") || die "Can't open the file $snvsGT_somaticRareBiasFile\n";

while(<SOM_RareBias>) {
  chomp;
  if($_!~/^#/) {
    if($_=~/Tumor_Somatic_Common/ && $_!~/bPcr|bSeq/) {
      $json{'SomaticSmallVarsInTumorCommonInGnomad'}++;
    }
    elsif($_=~/Tumor_Somatic/ && $_=~/bPcr|bSeq/) {
      $json{'SomaticSmallVarsInTumorInBias'}++;
    }
    elsif($_=~/Tumor_Somatic_Rare/) {
      $json{'SomaticSmallVarsInTumorPass'}++;
    }

    if($_=~/Control_Somatic_Common/ && $_!~/bPcr|bSeq/) {
      $json{'SomaticSmallVarsInControlCommonInGnomad'}++;    
    }
    elsif($_=~/Control_Somatic/ && $_=~/bPcr|bSeq/) {
      $json{'SomaticSmallVarsInControlInBias'}++;
    }
    elsif($_=~/Control_Somatic_Rare/) {
      $json{'SomaticSmallVarsInControlPass'}++;
    }   
  }
}

##################
## Creating sample swap json file 

## Percentage calculations
if($json{'SomaticSmallVarsInTumor'} > 0) {
  $json{'SomaticSmallVarsInTumorCommonInGnomADPer'} = $json{'SomaticSmallVarsInTumorCommonInGnomAD'}/$json{'SomaticSmallVarsInTumor'};
  $json{'SomaticSmallVarsInTumorInBiasPer'} = $json{'SomaticSmallVarsInTumorInBias'}/$json{'SomaticSmallVarsInTumor'};
  $json{'SomaticSmallVarsInTumorPassPer'} = $json{'SomaticSmallVarsInTumorPass'}/$json{'SomaticSmallVarsInTumor'};
}
else {
  $json{'SomaticSmallVarsInTumorCommonInGnomadPer'} = 0;
  $json{'SomaticSmallVarsInTumorInBiasPer'} = 0;
  $json{'SomaticSmallVarsInTumorPassPer'} = 0;
}

if($json{'SomaticSmallVarsInControl'} > 0) {
  $json{'SomaticSmallVarsInControlCommonInGnomasPer'} = $json{'SomaticSmallVarsInControlCommonInGnomad'}/$json{'SomaticSmallVarsInControl'};
  $json{'SomaticSmallVarsInControlInBiasPer'} = $json{'SomaticSmallVarsInControlInBias'}/$json{'SomaticSmallVarsInControl'};
  $json{'SomaticSmallVarsInControlPassPer'} = $json{'SomaticSmallVarsInControlPass'}/$json{'SomaticSmallVarsInControl'};
}
else {
  $json{'SomaticSmallVarsInControlCommonInGnomadPer'} = 0;
  $json{'SomaticSmallVarsInControlInBiasPer'} = 0;
  $json{'SomaticSmallVarsInControlPassPer'} = 0;

}

print JSON create_json (\%json);
close JSON;

######################################
#### Cleaning up files 
`rm $snvsGT_RawFile $snvsGT_gnomADFile`;
