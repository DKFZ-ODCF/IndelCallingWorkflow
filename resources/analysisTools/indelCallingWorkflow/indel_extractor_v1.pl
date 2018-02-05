#!usr/bin/perl

# By Ivo Buchhalter
# Hard coded script to extract indels in the snv pipeline
# Extracts: somatic indels -> indels_PID_somtatic_indels_conf_?_to_10.vcf
#           somatic coding indels -> indels_PID_somatic_coding_indels_conf_?_to_10.vcf
# optionally: germline coding indels -> indels_PID_germline_coding_indels_conf_?_to_10.vcf


use strict;
use warnings;
use Getopt::Long;

my $infile;
my $minconf=8;
my $region = 0;
my $somout;
my $funcout;
my $ncout;
my $germlineout;
my $bgzip = "bgzip";
my $tabix = "tabix";
GetOptions (	"infile=s"	 		=> \$infile,		# vcf file, can be bgzipped
				"minconf=i"			=> \$minconf,		# minimum confidence score
				"region=s"			=> \$region,		# file with the region when only a region should be extracted
				"somout=s"			=> \$somout,		# Outputfile for somatic indels
				"funcout=s"			=> \$funcout,		# Outputfile for functional somatic indels
				"ncout=s"			=> \$ncout,			# Outputfile for non coding RNAS
				"germlineout=s"		=> \$germlineout,	# germline coding indels
				"bgzip=s"			=> \$bgzip,
				"tabix=s"			=> \$tabix
) or die "Could not get the options!\n";

if($region ne "0" && (!-f $region || $infile !~ /\.gz$/ || !-f $infile.".tbi")){die "region-file: $region is not a valid file or, infile: $infile is not zipped or there is no index for the infile $infile.tbi\n";}
if(!-f $infile){die "The provided infile: $infile is not a valid file\n";}
if(-f $infile && $infile =~ /\.gz$/){open(IN, "zcat $infile |");}
else{open(IN, "<$infile") or die "Could not open the infile: $infile\n";}

if($region ne "0"){print "Region file provided: $region\n";}

open(SOM, ">$somout") or die "Could not open the file $somout\n";
open(COD, ">$funcout") or die "Could not open the file $funcout\n";
open(NCO, ">$ncout") or die "Could not open the file $ncout\n";
if (defined $germlineout)
{
	open(GER, ">$germlineout") or die "Could not open the file $germlineout\n";
}

my $head;
while(<IN>)
{
	chomp;
	$head=$_;
	last if($_ =~ /^#CHR/);
}
close IN;

my @head=split("\t", $head);
my %col;


my $i = 0;
while($i < @head)
{
	if($head[$i] eq "EXONIC_CLASSIFICATION"){$col{"EXONIC_CLASSIFICATION"} = $i; print "EXONIC_CLASSIFICATION in column ", $i+1,"\n";}
	if($head[$i] eq "ANNOVAR_FUNCTION"){$col{"ANNOVAR_FUNCTION"} = $i; print "ANNOVAR_FUNCTION in column ", $i+1,"\n";}
	if($head[$i] eq "CLASSIFICATION"){$col{"CLASSIFICATION"} = $i; print "CLASSIFICATION in column ", $i+1,"\n";}
	if($head[$i] eq "CONFIDENCE"){$col{"CONFIDENCE"} = $i; print "CONFIDENCE in column ", $i+1,"\n";}
	if($head[$i] eq "REGION_CONFIDENCE"){$col{"REGION_CONFIDENCE"} = $i; print "REGION_CONFIDENCE in column ", $i+1,"\n";}
	$i++;
}

print SOM $head, "\n";
print COD $head, "\n";
print NCO $head, "\n";
if (defined $germlineout)
{
	print GER $head, "\n";
}

if($region ne "0"){open(IN, "$tabix $infile -B $region |") or die "Could not open the file with tabix and regions\n";}
elsif($infile =~ /\.gz$/){open(IN, "zcat $infile |");}
else{open(IN, "<$infile");}

while(<IN>)
{
	chomp;
	next if($_ =~ /^#/);
	my @line = split("\t", $_);
	if (defined $germlineout && $line[$col{"REGION_CONFIDENCE"}] >= $minconf && ($line[$col{"CLASSIFICATION"}] eq "germline" || $line[$col{"CLASSIFICATION"}] eq "SNP_support_germline") && ($line[$col{"ANNOVAR_FUNCTION"}] eq "exonic" || $line[$col{"ANNOVAR_FUNCTION"}] =~ /splicing/) && $line[$col{"ANNOVAR_FUNCTION"}] !~ /ncRNA_exonic/)
	{
		print GER $_, "\n";
	}
	next if($line[$col{"CONFIDENCE"}] < $minconf);
	if($line[$col{"CLASSIFICATION"}] eq "somatic"){print SOM $_, "\n";}
	if($line[$col{"CLASSIFICATION"}] eq "somatic" && $line[$col{"ANNOVAR_FUNCTION"}] !~ /ncRNA/ && ($line[$col{"ANNOVAR_FUNCTION"}] =~ /exonic/ || $line[$col{"ANNOVAR_FUNCTION"}] =~ /splicing/)){print COD $_, "\n";}
	if($line[$col{"CLASSIFICATION"}] eq "somatic" && ($line[$col{"ANNOVAR_FUNCTION"}] =~ /ncRNA_exonic/ || $line[$col{"ANNOVAR_FUNCTION"}] =~ /ncRNA_splicing/)){print NCO $_, "\n";}
}

close IN;
close SOM;
close COD;
if (defined $germlineout)
{
	close GER;
}

