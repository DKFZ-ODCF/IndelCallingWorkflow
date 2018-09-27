#/usr/bin/env python

# confidenceAnnotation_Indels.py
# Author: Jeongbin Park (based on perl code 'platypusConfidenceAnnotation.pl')
# Purpose: Annotate 'CONFIDENCE' and 'RECLASSIFICATION' to the vcf file.
# Usage: python -u confidenceAnnotation_Indels.py --options > output.vcf
# Explanation:

import argparse
import sys
import time
from vcfparser import *  # BGZFType, LineParser, get_header_indices


def extract_info(info, keys, sep=";"):
    info_kv = {key: val for key, val in (j if len(j) == 2 else (j[0], None) for j in (i.split('=') for i in info.split(sep)))}

    if type(keys) is list:
        rtn = []
        for key in keys:
            rtn.append(info_kv.get(key, None))
        rtn = [0 if i == 'None' else i for i in rtn]
        return rtn

    if type(keys) is str:
        rtn = info_kv.get(keys, None) 
        rtn = 0 if rtn == "None" else rtn
        return rtn

def main(args):
    if not args.no_makehead:
        header = '##fileformat=VCFv4.1\n' \
                 '##fileDate=' + time.strftime("%Y%m%d") + '\n' \
                 '##pancancerversion=1.0\n' \
                 '##reference=<ID=' + args.refgenome[0] + ',Source=' + args.refgenome[1] + '>\n' \
                 '##center=' + args.center + '\n' \
                 '##workflowName=DKFZ_indel_workflow\n' \
                 '##workflowVersion=1.0.0\n'
        header += '\n'.join(args.additional_header) + '\n' if len(args.additional_header) > 0 else ""
        header += '##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Indicates if record is a somatic mutation">\n' \
                  '##INFO=<ID=GERMLINE,Number=0,Type=Flag,Description="Indicates if record is a germline mutation">\n' \
                  '##INFO=<ID=UNCLEAR,Number=0,Type=Flag,Description="Indicates if the somatic status of a mutation is unclear">\n' \
                  '##INFO=<ID=FR,Number=.,Type=Float,Description="Estimated population frequency of variant">\n' \
                  '##INFO=<ID=MMLQ,Number=1,Type=Float,Description="Median minimum base quality for bases around variant">\n' \
                  '##INFO=<ID=TCR,Number=1,Type=Integer,Description="Total reverse strand coverage at this locus">\n' \
                  '##INFO=<ID=HP,Number=1,Type=Integer,Description="Homopolymer run length around variant locus">\n' \
                  '##INFO=<ID=WE,Number=1,Type=Integer,Description="End position of calling window">\n' \
                  '##INFO=<ID=Source,Number=.,Type=String,Description="Was this variant suggested by Playtypus, Assembler, or from a VCF?">\n' \
                  '##INFO=<ID=FS,Number=.,Type=Float,Description="Fisher\'s exact test for strand bias (Phred scale)">\n' \
                  '##INFO=<ID=WS,Number=1,Type=Integer,Description="Starting position of calling window">\n' \
                  '##INFO=<ID=PP,Number=.,Type=Float,Description="Posterior probability (phred scaled) that this variant segregates">\n' \
                  '##INFO=<ID=TR,Number=.,Type=Integer,Description="Total number of reads containing this variant">\n' \
                  '##INFO=<ID=NF,Number=.,Type=Integer,Description="Total number of forward reads containing this variant">\n' \
                  '##INFO=<ID=TCF,Number=1,Type=Integer,Description="Total forward strand coverage at this locus">\n' \
                  '##INFO=<ID=NR,Number=.,Type=Integer,Description="Total number of reverse reads containing this variant">\n' \
                  '##INFO=<ID=TC,Number=1,Type=Integer,Description="Total coverage at this locus">\n' \
                  '##INFO=<ID=END,Number=.,Type=Integer,Description="End position of reference call block">\n' \
                  '##INFO=<ID=MGOF,Number=.,Type=Integer,Description="Worst goodness-of-fit value reported across all samples">\n' \
                  '##INFO=<ID=SbPval,Number=.,Type=Float,Description="Binomial P-value for strand bias test">\n' \
                  '##INFO=<ID=START,Number=.,Type=Integer,Description="Start position of reference call block">\n' \
                  '##INFO=<ID=ReadPosRankSum,Number=.,Type=Float,Description="Mann-Whitney Rank sum test for difference between in positions of variants in reads from ref and alt">\n' \
                  '##INFO=<ID=MQ,Number=.,Type=Float,Description="Root mean square of mapping qualities of reads at the variant position">\n' \
                  '##INFO=<ID=QD,Number=1,Type=Float,Description="Variant-quality/read-depth for this variant">\n' \
                  '##INFO=<ID=SC,Number=1,Type=String,Description="Genomic sequence 10 bases either side of variant position">\n' \
                  '##INFO=<ID=BRF,Number=1,Type=Float,Description="Fraction of reads around this variant that failed filters">\n' \
                  '##INFO=<ID=HapScore,Number=.,Type=Integer,Description="Haplotype score measuring the number of haplotypes the variant is segregating into in a window">\n' \
                  '##INFO=<ID=Size,Number=.,Type=Integer,Description="Size of reference call block">\n' \
                  '##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership">\n' \
                  '##INFO=<ID=1000G,Number=0,Type=Flag,Description="Indicates membership in 1000Genomes">\n' \
                  '##FILTER=<ID=GOF,Description="Variant fails goodness-of-fit test.">\n' \
                  '##FILTER=<ID=badReads,Description="Variant supported only by reads with low quality bases close to variant position, and not present on both strands.">\n' \
                  '##FILTER=<ID=alleleBias,Description="Variant frequency is lower than expected for het">\n' \
                  '##FILTER=<ID=Q20,Description="Variant quality is below 20.">\n' \
                  '##FILTER=<ID=HapScore,Description="Too many haplotypes are supported by the data in this region.">\n' \
                  '##FILTER=<ID=MQ,Description="Root-mean-square mapping quality across calling region is low.">\n' \
                  '##FILTER=<ID=strandBias,Description="Variant fails strand-bias filter">\n' \
                  '##FILTER=<ID=SC,Description="Variants fail sequence-context filter. Surrounding sequence is low-complexity">\n' \
                  '##FILTER=<ID=QD,Description="Variants fail quality/depth filter.">\n' \
                  '##FILTER=<ID=ALTC,Description="Alternative reads in control and other filter not PASS">\n' \
                  '##FILTER=<ID=VAF,Description="Variant allele frequency in tumor < 10% and other filter not PASS">\n' \
                  '##FILTER=<ID=VAFC,Description="Variant allele frequency in tumor < 5% or variant allele frequency in control > 5%">\n' \
                  '##FILTER=<ID=QUAL,Description="Quality of entry too low and/or low coverage in region">\n' \
                  '##FILTER=<ID=ALTT,Description="Less than three variant reads in tumor">\n' \
                  '##FILTER=<ID=GTQ,Description="Quality for genotypes below thresholds">\n' \
                  '##FILTER=<ID=GTQFRT,Description="Quality for genotypes below thresholds and variant allele frequency in tumor < 10%">\n' \
                  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Unphased genotypes">\n' \
                  '##FORMAT=<ID=GQ,Number=.,Type=Integer,Description="Genotype quality as phred score">\n' \
                  '##FORMAT=<ID=GOF,Number=.,Type=Float,Description="Goodness of fit value">\n' \
                  '##FORMAT=<ID=NR,Number=.,Type=Integer,Description="Number of reads covering variant location in this sample">\n' \
                  '##FORMAT=<ID=GL,Number=.,Type=Float,Description="Genotype log10-likelihoods for AA,AB and BB genotypes, where A = ref and B = variant. Only applicable for bi-allelic sites">\n' \
                  '##FORMAT=<ID=NV,Number=.,Type=Integer,Description="Number of reads containing variant in this sample">\n' \
                  '##SAMPLE=<ID=CONTROL,SampleName=control_' + args.pid + ',Individual=' + args.pid + ',Description="Control">\n' \
                  '##SAMPLE=<ID=TUMOR,SampleName=tumor_' + args.pid + ',Individual=' + args.pid + ',Description="Tumor">\n' \
                  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
        if not args.no_control:
            header += "CONTROL\t"
        header += "TUMOR\t"
        sys.stdout.write(header)

    for line in args.infile:
        if line[:2] == "##":
            continue

        if line[0] == "#":
            headers = list(line[1:].rstrip().split('\t'))
            fixed_headers = ["^QUAL$", "^INFO$", "^FILTER$" , "MAPABILITY", "HISEQDEPTH", "SIMPLE_TANDEMREPEATS",
                             "REPEAT_MASKER", "DUKE_EXCLUDED", "DAC_BLACKLIST", "SELFCHAIN", "^CONFIDENCE$",
                             "^CLASSIFICATION$", "^REGION_CONFIDENCE$", "^PENALTIES$", "^REASONS$",
                            ]
            variable_headers = { "ANNOVAR_SEGDUP_COL": "^SEGDUP$", "KGENOMES_COL": "^1K_GENOMES$", "DBSNP_COL": "^DBSNP$",
                                 "CONTROL_COL": "^" + args.controlColName + "$", "TUMOR_COL": "^" + args.tumorColName + "$"}

            if args.no_control:
                variable_headers["ExAC_COL"] = "^ExAC$"
                variable_headers["EVS_COL"] = "^EVS$"
                variable_headers["GNOMAD_EXOMES_COL"] = "^GNOMAD_EXOMES$"
                variable_headers["GNOMAD_GENOMES_COL"] = "^GNOMAD_GENOMES$"
                variable_headers["LOCALCONTROL_COL"] = "^LocalControlAF$"
            else:
                fixed_headers += [ "^INFO_control", "^ANNOTATION_control$", ]

            header_indices = get_header_indices(headers, args.configfile, fixed_headers, variable_headers)

            if args.no_control:
                if header_indices["TUMOR_COL"] == -1:
                    header_indices["TUMOR_COL"] = 9
            else:
                if header_indices["CONTROL_COL"] == -1:
                    header_indices["CONTROL_COL"] = 9

                if header_indices["TUMOR_COL"] == -1:
                    header_indices["TUMOR_COL"] = 10

            # create headers if they don't exist
            for optional_header in ["CLASSIFICATION", "CONFIDENCE", "REGION_CONFIDENCE", ]:
                if header_indices[optional_header] == -1:
                    headers.append(optional_header)

            if args.debug:
                for debug_header in ["PENALTIES", "REASONS", ]:
                    if header_indices[debug_header] == -1:
                        headers.append(debug_header)

            if args.no_makehead:
                header_str = "#" + '\t'.join(headers)
            else:
                if args.no_control:
                    header_str = '\t'.join(headers[10:])
                else:
                    header_str = '\t'.join(headers[11:])
            sys.stdout.write(header_str + '\n')

            continue

        entries = list(line.rstrip().split('\t'))

        # Ignore calls in contigs
        if entries[0][:3] == "chr":
            chrom = entries[0][3:]
        else:
            chrom = entries[0]
        if not chrom.isdigit and not chrom == "X" and not chrom == "Y":
            continue

        if not args.notonlyindel:
            # Ignore calls different from indels, should be changed later so that we also annotate replacements and maybe multiple alternatives
            if "," in entries[3] or "," in entries[4] or \
                (len(entries[3]) > 1 and len(entries[4]) > 1) or \
                (len(entries[3]) == 1 and len(entries[4]) == 1):
                continue

        help = LineParser(entries, header_indices) # This will let you access values using header name

        classification = ""
        confidence = 10 # start with maximum value and subtract something for each "bad" feature
                        # "high" = 9-10, "medium" = 6-8, low <= 5

        VAFControl = 0
        VAFTumor = 0
        filter = {}
        dbsnp_pos = None
        dbsnp_id = None
        region_conf = 10
        reasons = ""
        penalties = ""
        infos = []

        if args.no_control:
            in1KG_AF = False
            indbSNP = False

            is_commonSNP = False
            is_clinic = False
            inExAC = False
            inEVS = False
            inGnomAD_WES = False
            inGnomAD_WGS = False
            inLocalControl = False

        # 1) external information of if these SNPs have already been found (incl. false positives from 1000 genomes!)
        # dbSNP
        if help["DBSNP_COL_VALID"] and "MATCH=exact" in help["DBSNP_COL"]:
            indbSNP = True
            infos.append("DB")
            dbsnp_pos = entries[1]
            dbsnp_id = extract_info(help["DBSNP_COL"].split("&")[0], "ID")
            # precious!
            #INFO=<ID=CLN,Number=0,Type=Flag,Description="SNP is Clinical(LSDB,OMIM,TPA,Diagnostic)">
            #INFO=<ID=PM,Number=0,Type=Flag,Description="SNP is Precious(Clinical,Pubmed Cited)">
            if ";CLN;" in help["DBSNP_COL"]:
                is_clinic = True
            if "COMMON=1" in help["DBSNP_COL"]:
                is_commonSNP = True

        # 1000 genomes
        if help["KGENOMES_COL_VALID"] and "MATCH=exact" in help["KGENOMES_COL"]:
            if args.no_control:
                af = extract_info(help["KGENOMES_COL"].split("&")[0], "EUR_AF")
                if af is not None and any(af > 0.01 for af in map(float, af.split(','))) > 0.01:
                    in1KG_AF = True
            infos.append("1000G")

        if args.no_control:
            if help["ExAC_COL_VALID"] and any(af > 0.001 for af in map(float, extract_info(help["ExAC_COL"], "AF").split(','))):
                inExAC = True
                infos.append("ExAC")
            if help["EVS_COL_VALID"] and any(af > 1.0 for af in map(float, extract_info(help["EVS_COL"], "MAF").split(','))):
                inEVS = True
                infos.append("EVS")
            if help["GNOMAD_EXOMES_COL_VALID"] and any(af > 0.001 for af in map(float, extract_info(help["GNOMAD_EXOMES_COL"], "AF").split(','))):
                inGnomAD_WES = True
                infos.append("gnomAD_Exomes")
            if help["GNOMAD_GENOMES_COL_VALID"] and any(af > 0.001 for af in map(float, extract_info(help["GNOMAD_GENOMES_COL"], "AF").split(','))):
                inGnomAD_WGS = True
                infos.append("gnomAD_Genomes")
            if help["LOCALCONTROL_COL_VALID"] and any(af > 0.02 for af in map(float, extract_info(help["LOCALCONTROL_COL"], "AF").split(','))):
                inLocalControl = True
                infos.append("LOCALCONTROL")

        qual = help["QUAL"]
        ### variants with more than one alternative are still skipped e.g. chr12	19317131	.	GTT	GT,G	...
        if (args.no_control or help["CONTROL_COL"].split(":")[0] == "0/0") and any(word in help["TUMOR_COL"].split(":")[0] for word in ("1/0", "0/1", "1/1")):
            classification = "somatic" # All calls with this genotype are called somatic
            if not args.no_control:
                controlGT, controlGP, C_GOF, controlGQ, controlDP, controlDP_V = help["CONTROL_COL"].strip().split(":")
                C_GOF, controlGQ, controlDP, controlDP_V = int(C_GOF), int(controlGQ), int(controlDP), int(controlDP_V)
            tumorGT, tumorGP, T_GOF, tumorGQ, tumorDP, tumorDP_V = help["TUMOR_COL"].strip().split(":")
            T_GOF, tumorGQ, tumorDP, tumorDP_V = int(T_GOF), int(tumorGQ), int(tumorDP), int(tumorDP_V)

            if (args.no_control or controlDP > 0) and tumorDP > 0:
                if not args.no_control:
                    VAFControl= (controlDP_V/float(controlDP))*100.0
                VAFTumor  = (tumorDP_V/float(tumorDP))*100.0

            if not "PASS" in help["FILTER"]:
                if "alleleBias" in help["FILTER"]:
                    confidence -= 2
                    penalties += "alleleBias_-2_"
                    filter["alleleBias"] = 1
                    region_conf -= 2
                    reasons += "alleleBias(-2)"
                if "badReads" in help["FILTER"]:
                    confidence -= 3
                    penalties += "badReads_-3_"
                    filter["badReads"] = 1
                    region_conf -= 3
                    reasons+="badReads(-3)"
                if "MQ" in help["FILTER"]:
                    confidence -= 1
                    penalties += "MQ_-1_"
                    filter["MQ"] = 1
                    region_conf -= 1
                    reasons += "MQ(-1)"
                if "SC" in help["FILTER"]:
                    confidence -= 1
                    penalties += "SC_-1_"
                    filter["SC"] = 1
                    region_conf -= 1
                    reasons += "SC(-1)"
                if "GOF" in help["FILTER"]:
                    confidence -= 1
                    penalties += "GOF_-1_"
                    filter["GOF"] = 1
                    region_conf -= 1
                    reasons += "GOF(-1)"
                if "QD" in help["FILTER"]:
                    confidence -= 1
                    penalties += "QD_-1_"
                    filter["QD"] = 1
                    region_conf -= 1
                    reasons += "QD(-1)"
                if "strandBias" in help["FILTER"]:
                    confidence -= 2
                    penalties += "strandBias_-2_"
                    filter["strandBias"] = 1
                    region_conf -= 2
                    reasons += "strandBias(-2)"
                if not args.no_control and controlDP_V > 0:
                    confidence -= 1
                    penalties += "alt_reads_in_control_-1_"
                    filter["ALTC"] = 1
                if VAFTumor < 10:
                    confidence -= 1
                    penalties += "VAF<10_-1_"
                    filter["VAF"] = 1

            # Minimum base quality, read depth and genotype quality
            if qual > 40 and (args.no_control or controlDP >=10) and tumorDP >= 10 and \
                    (args.no_control or controlGQ >= 20) and tumorGQ >= 20: # All quality filters are OK
                # Do nothing
                pass
            elif qual > 20 and (args.no_control or controlDP >= 5) and tumorDP >= 5 and \
                    (args.no_control or controlGQ >= 20) and tumorGQ >= 20:
                confidence -= 2
                penalties += "medium_quality_values_-2_"
                filter["QUAL"] = 1
            else: # Bad quality values (is equal to medium bad values)(-2)
                confidence -= 2
                penalties += "bad_quality_values_-2_"
                filter["QUAL"] = 1

            if tumorDP_V < 3: # Less than three variant reads in tumor (-2)
                confidence -= 2
                penalties += "<3_reads_in_tumor_-2_"
                filter["ALTT"] = 1

            if not args.no_control and controlDP_V > 0:
                if VAFControl < 5 and VAFTumor > 5:	# VAF in control below 0.05 and bigger 0.05 in tumor (-1)
                    confidence -= 1
                    penalties += "alt_reads_in_control(VAF<0.05)_-1_"
                    filter["ALTC"] = 1
                else:	# Almost always bad (-3)
                    confidence -= 3
                    penalties += "alt_reads_in_control(VAF>=0.05_or_tumor_VAF<=0.05)_-3_"
                    filter["VAFC"] = 1

            if not args.no_control:
                ssControlGP = list(map(float, controlGP.split(",")))
            sstumorGP = list(map(float, tumorGP.split(",")))

            # Genotype probability
            if tumorGT == "1/0" or tumorGT == "0/1":
                if not ((args.no_control or (ssControlGP[1] < args.hetcontr and ssControlGP[2] < args.homcontr)) and sstumorGP[0] < args.homreftum and sstumorGP[2] < args.tumaltgen):
                    # At least one genotype quality is bad (-2)
                    confidence -= 2
                    penalties += "bad_genotype_quality_-2_"
                    filter["GTQ"] = 1
                    if not args.no_control and controlDP_V > 0:	# Bad genotype and alternative reads in control (-1)
                        confidence -= 1
                        penalties += "alt_reads_in_control_-1_"
                        filter["ALTC"] = 1
                    if VAFTumor < 10:	# Bad genotype and VAF below 0.1 (-1)
                        confidence -= 1
                        penalties += "VAF<10_-1_"
                        filter["GTQFRT"] = 1
            elif tumorGT == ("1/1"):
                if not ((args.no_control or (ssControlGP[1] < args.hetcontr and ssControlGP[2] < args.homcontr)) and sstumorGP[0] < args.homreftum and sstumorGP[1] < args.tumaltgen):
                    # At least one genotype quality is bad (-2)
                    confidence -= 2
                    penalties += "bad_genotype_quality_-2_"
                    filter["GTQ"] = 1
                    if not args.no_control and controlDP_V > 0:	# Bad genotype and alternative reads in control (-1)
                        confidence -= 1
                        penalties += "alt_reads_in_control_-1_"
                        filter["ALTC"] = 1
                    if VAFTumor < 10:	# Bad genotype and VAF below 0.1 (-1)
                        confidence -= 1
                        penalties += "VAF<10_-1_"
                        filter["GTQFRT"] = 1
        else:
            if not args.no_control and any(help["CONTROL_COL"][:3] == word for word in ("1/0", "0/1", "1/1", )):
                classification = "germline"
            else:
                classification = "unclear"
            confidence = 1

        if args.no_control and (in1KG_AF or (indbSNP and is_commonSNP and not is_clinic) or inExAC or inEVS or inGnomAD_WES or inGnomAD_WGS or inLocalControl):
            classification = "SNP_support_germline"

        if confidence < 1:	# Set confidence to 1 if it is below one
            confidence = 1
        if confidence > 10:	# Set confidence to 10 if above (will not happen at the moment as we never give a bonus)
            confidence = 10

        ###### Stuff from Barbara
        # more filters to assess how good the region is, i.e. if indel overlaps with strange regions
        # the blacklists have few entries; the HiSeqDepth has more "reads attracting" regions,
        # often coincide with tandem repeats and CEN/TEL, not always with low mapability
        # Duke excluded and ENCODE DAC blacklist, only consider if not already annotated as suspicious repeat
        if help["DUKE_EXCLUDED_VALID"] or help["DAC_BLACKLIST_VALID"] or help["HISEQDEPTH_VALID"]:
            region_conf -= 3 # really bad region, usually centromeric repeats
            reasons += "Blacklist(-3)"

        if help["ANNOVAR_SEGDUP_COL_VALID"] or help["SELFCHAIN_VALID"]:
            region_conf -= 1
            reasons += "SelfchainAndOrSegdup(-1)"

        if any(word in help["REPEAT_MASKER"] for word in ["Simple_repeat", "Low_", "Satellite", ]) or help["SIMPLE_TANDEMREPEATS_VALID"]:
            region_conf -= 2
            reasons += "Repeat(-2)"
        # other repeat elements to penalize at least a bit
        elif help["REPEAT_MASKER_VALID"]:
            region_conf -= 1
            reasons += "Other_repeat(-1)"

        # Mapability is 1 for unique regions, 0.5 for regions appearing twice, 0.33... 3times, ...
        # Everything with really high number of occurences is artefacts
        # does not always correlate with the above regions
        # is overestimating badness bc. of _single_ end read simulations
        if help["MAPABILITY"] == ".":
            # in very rare cases (CEN), there is no mapability => ".", which is not numeric but interpreted as 0
            region_conf -= 5
            reasons += "Not_mappable(-5)"
        else:
            reduce = 0
            # can have several entries for indels, e.g. 0.5&0.25 - take worst (lowest) or best (highest)?!
            mapability = min(map(float, help["MAPABILITY"].split("&"))) # just chose one - here: min
            if mapability < 0.5:
                # 0.5 does not seem to be that bad: region appears another time in
                # the genome and we have paired end data!
                region_conf -= 1
                reduce += 1

                #if mapability < 0.4: # 3-4 times appearing region is worse but still not too bad
                #    region_conf -= 1
                #    reduce += 1

                if mapability < 0.25: # > 4 times appearing region
                    region_conf -= 1
                    reduce += 1

                if mapability < 0.1: # > 5 times is bad
                    region_conf -= 2
                    reduce += 2

                if mapability < 0.05: # these regions are clearly very bad (Lego stacks)
                    region_conf -= 3
                    reduce += 3

                reasons += "Low_mappability(%s=>-%d)"%(help["MAPABILITY"], reduce)

        if classification != "somatic" and not "PASS" in help["FILTER"]:
            # such an indel is probably also "unclear" and not "germline"
            # all filters were already punished before but for germline we really want to be strict
            region_conf -= 3
            reasons += "notPASS(-3)"

        if region_conf < 1:
            region_conf = 1

        idx_classification = header_indices["CLASSIFICATION"]
        idx_confidence = header_indices["CONFIDENCE"]
        idx_region_confidence = header_indices["REGION_CONFIDENCE"]

        if args.debug:
            idx_penalties = header_indices["PENALTIES"]
            idx_reasons = header_indices["REASONS"]

        if idx_classification == -1:
            entries.append(classification)
        else:
            entries[idx_classification] = classification

        if idx_confidence == -1:
            entries.append(str(confidence))
        else:
            entries[idx_confidence] = str(confidence)

        if idx_region_confidence == -1:
            entries.append(str(region_conf))
        else:
            entries[idx_region_confidence] = str(region_conf)

        if args.debug:
            if idx_penalties == -1:
                entries.append(penalties)
            else:
                entries[idx_penalties] = penalties
            if idx_reasons == -1:
                entries.append(reasons)
            else:
                entries[idx_reasons] = reasons

        if classification == "somatic":
            entries[header_indices["INFO"]] = 'SOMATIC;' + entries[header_indices["INFO"]]
            if confidence >= 8:
                entries[header_indices["FILTER"]] = "PASS"
            else:
                filter_list = []
                filteroptions = ["GOF","badReads","alleleBias","MQ","strandBias","SC","QD","ALTC","VAF","VAFC","QUAL","ALTT","GTQ","GTQFRT",]
                for filteroption in filteroptions:
                    if filter.get(filteroption, 0) == 1:
                        filter_list.append(filteroption)
                entries[header_indices["FILTER"]] = ';'.join(filter_list)
        elif classification == "germline" or (args.no_control and classification == "SNP_support_germline"):
            entries[header_indices["INFO"]] = 'GERMLINE;' + entries[header_indices["INFO"]]
        else:
            entries[header_indices["INFO"]] = 'UNCLEAR;' + entries[header_indices["INFO"]]

        entries[header_indices["QUAL"]] = "."
        if dbsnp_id is not None and dbsnp_pos is not None:
            entries[2] = dbsnp_id + "_" + dbsnp_pos

        print '\t'.join(entries)

    args.infile.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate 'CONFIDENCE' and 'CLASSIFICATION' to the vcf file.")
    parser.add_argument("-b", "--configfile", dest="configfile", nargs="?", type=argparse.FileType('r'), default=None,
                        help='Specify the config file which contains header names.')
    parser.add_argument("-w", "--nocontrol", dest="no_control", action="store_true", default=False,
                        help='Set this flag if input vcf file does not contain control information.')
    parser.add_argument("-m", "--nomakehead", dest="no_makehead", action="store_true", default=False,
                        help="Set this flag if you do NOT want to produce a pancancer conform head.")
    parser.add_argument("-i", "--infile", dest="infile", nargs="?",
                        type=BGZFType('r'), default=sys.stdin,
                        help="Specify the path of vcf file. If not specified, then stdin is used instead.")
    parser.add_argument("-c", "--controlColName", dest="controlColName", default="CONTROL_COL", help="Column header name of Control.")
    parser.add_argument("-t", "--tumorColName", dest="tumorColName", default="TUMOR_COL", help="Column header name of Tumor.")
    parser.add_argument("-d", "--debug", dest="debug", action="store_true", default=False,
                        help="Print justification for confidence score.")
    parser.add_argument("-a", "--anno", dest="anno", type=int, nargs="?", default=0,
                        help="Set this flag if you want to print the original annotation.")
    parser.add_argument("-o", "--notonlyindel", dest="notonlyindel", action="store_true", default=False,
                        help="Set this flag if you want to also annotate not only real indels.")
    parser.add_argument("-H", "--addhead", dest="additional_header", nargs="+", default=[],
                        help="String with additional header line infer multiple times for multiple additional lines.")
    parser.add_argument("-p", "--pid", dest="pid", nargs="?", help="Patient ID (default NA).", default="NA")
    parser.add_argument("-g", "--refgenome", dest="refgenome", nargs=2,
                        default=["hs37d5", "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/" \
                                           "phase2_reference_assembly_sequence/hs37d5.fa.gz", ],
                        help="reference genome used for calling ID, path (default hs37d5, ftp://ftp.1000genomes.ebi" \
                             ".ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz)")
    parser.add_argument("-z", "--center", dest="center", nargs="?", default="DKFZ",
                        help="Center (unclear if the center where the data was produced or the center of the " \
                             "calling pipeline; default DKFZ).")
    parser.add_argument("--hetcontr", dest="hetcontr", type=float, default=-4.60517,
                        help="Score that a 0/0 call in the control is actually 0/1 or 1/0 (the more negative, the less likely).")
    parser.add_argument("--homcontr", dest="homcontr", type=float, default=-4.60517,
                        help="Score that a 0/0 call in the control is actually 1/1 (the more negative, the less likely).")
    parser.add_argument("--homreftum", dest="homreftum", type=float, default=-4.60517,
                        help="Score that a 0/1 or 1/0 or 1/1 in tumor is actually 0/0 (the more negative, the less likely).")
    parser.add_argument("--tumaltgen", dest="tumaltgen", type=float, default=0,
                        help="Score that a 0/1 or 1/0 call in tumor is actually 1/1 or that a 1/1 call in tumor is " \
                             "actually 1/0 or 0/1 (the more negative, the less likely).")
args = parser.parse_args()
main(args)
