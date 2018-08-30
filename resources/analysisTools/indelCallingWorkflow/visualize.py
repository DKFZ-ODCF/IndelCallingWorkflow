#!/usr/bin/env python
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/IndelCallingWorkflow).
#


# This script generates a png file for each entry in a VCF file. The file
# displays the reads in a bam file around the VCF entry (+-WINDOW_SIZE).
# Reallocation of memry works now

import sys
import subprocess
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
from matplotlib import gridspec


argument_parser = argparse.ArgumentParser(
  description='This script generates a png file for each entry in a vcf file, a bed file or a manually specified region.' )
argument_parser.add_argument('--control', metavar='FILE', type=str,
    default=None, help='input bam file of the control')
argument_parser.add_argument('--tumor', metavar='FILE', type=str, 
    required=True, help='input bam file of the tumor')
argument_parser.add_argument('--ref', metavar='FILE', type=str,
    required=True, help='input reference genome file (fastq format)')
argument_parser.add_argument('--vcf', metavar='FILE', type=str,
    default=None, help='input vcf file ( as an alternative use --bed )')
argument_parser.add_argument('--bed', metavar='FILE', type=str,
    default=None, help='input bed file ( as an alternative use --vcf )')
argument_parser.add_argument('--annotations', metavar='FILE', type=str,
    default=None, help='annotation track indexed with tabix')
argument_parser.add_argument('--prefix', metavar='PREFIX', type=str,
    default="./", help='target directory and file name prefix for generated output files')
argument_parser.add_argument('--window', metavar='N', type=int,
    default=100, help='the output file for position X will show the region [X-window,X+window]')
argument_parser.add_argument('--samtoolsbin', metavar='N', type=str,
    default="samtools", help='the path to the samtools binary, default is \'samtools\'')
argument_parser.add_argument('--tabixbin', metavar='N', type=str,
    default="tabix", help='the path to the tabix binary, default is \'tabix\'')
argument_parser.add_argument('region', nargs='?', type=str,
    default=None, help='syntax either \'chr:start-end\' or \'chr:center\', use --vcf or --bed for more convenience')
parsed_arguments = argument_parser.parse_args()


basepair_colors = { 'A':"#009600", 'C':"#3030fe", 'G':"#d17105", 'T':"#ff0000", 'N':"#00ffff" }



class ReferenceBuffer(object):

	def __init__( self, filename, chromosome ):

		self.samtools_call = [ parsed_arguments.samtoolsbin, "faidx", filename ]
		self.chromosome    = chromosome;
		self.offset        = 0
		self.sequence      = ""

	def __getitem__( self, pos ):

		if len(self.sequence) > pos - self.offset >= 0:

			return self.sequence[ pos - self.offset ]

		else:

			region = "%s:%i-%i" % ( self.chromosome, pos - 1000, pos + 1000 )
			call   = self.samtools_call + [ region ]
			output = subprocess.check_output( call )

			self.offset   = max( 0, pos - 1000 )
			self.sequence = "".join( output.split('\n')[1:] ).upper()

			return self.sequence[ pos - self.offset ]



def get_annotations( region ):

	if parsed_arguments.annotations:

		call   = [ parsed_arguments.tabixbin, parsed_arguments.annotations, region ]
		output = subprocess.check_output( call )

		if not output:

			if region[:3] == "chr":

				call[2] = call[2][3:]

			else:

				call[2] = "chr" + call[2]

		output = subprocess.check_output( call )

		return [ [line.split('\t')[3],int(line.split('\t')[1]),int(line.split('\t')[2])] for line in output.split('\n') if len(line.split('\t')) >= 3 ]

	else:

		return None



def parse_cigar( cigar, pos ):

	cigar_struct = [["", None]]

	for char in cigar:

		if "0" <= char <= "9":

			if type(cigar_struct[-1][0]) != str:

				cigar_struct.append([ "", None ])

			cigar_struct[-1][0] = cigar_struct[-1][0] + char

		else:

			cigar_struct[-1][0] = int( cigar_struct[-1][0] )
			cigar_struct[-1][1] = char

	if cigar_struct[-1][1] == None:

		del cigar_struct[-1]

	abs_pos = pos
	rel_pos = 0

	if cigar_struct[0][1] == 'S':
		abs_pos -= cigar_struct[0][0]

	cigar_struct2 = []

	for n,t in cigar_struct:

		if t in ['M','S']:

			for m in range(n):
				cigar_struct2.append( (t,abs_pos+m, rel_pos+m) )

			abs_pos += n
			rel_pos += n

		elif t == 'D':

			for m in range(n):
				cigar_struct2.append( (t,abs_pos+m, rel_pos) )

			abs_pos += n

		elif t == 'I':

			cigar_struct2.append( (t,abs_pos, rel_pos+m) )

			for m in range(1,n):
				cigar_struct2.append( ('i',abs_pos, rel_pos+m) )

			rel_pos += n

	return cigar_struct2



def plot_histogram( cigars, ax ):

	points_deleted   = []
	points_clipped   = []
	points_unclipped = []

	for cigar in cigars:

		for (c_type, abs_pos, rel_pos) in cigar:

			if c_type in ['M']:
				points_deleted.append( abs_pos )
				points_clipped.append( abs_pos )
				points_unclipped.append( abs_pos )
			elif c_type in 'D':
				points_clipped.append( abs_pos )
				points_unclipped.append( abs_pos )
			elif c_type in ['S']:
				points_unclipped.append( abs_pos )

	points_deleted = sorted( points_deleted )
	points_clipped = sorted( points_clipped )
	points_unclipped = sorted( points_unclipped )

	min_unclipped = min( points_unclipped )
	max_unclipped = max( points_unclipped )

	histogram_deleted = [0] * ( 1 + max_unclipped - min_unclipped )
	histogram_clipped = [0] * ( 1 + max_unclipped - min_unclipped )
	histogram_unclipped = [0] * ( 1 + max_unclipped - min_unclipped )

	for p in points_deleted:

		histogram_deleted[ p - min_unclipped ] += 1

	for p in points_clipped:

		histogram_clipped[ p - min_unclipped ] += 1

	for p in points_unclipped:

		histogram_unclipped[ p - min_unclipped ] += 1


	original_x         = range( min_unclipped, max_unclipped + 1 )
	subdivided_x       = [ a+b for a in original_x for b in [-0.33,0.33]][1:-2]

	smooth_unclipped = [ x for x in histogram_unclipped for y in [0,1] ][1:-2]
	smooth_clipped   = [ x for x in histogram_clipped for y in [0,1] ][1:-2]
	smooth_deleted   = [ x for x in histogram_deleted for y in [0,1] ][1:-2]

	ax.fill_between( subdivided_x, 0, smooth_unclipped, color="lightblue" )
	ax.fill_between( subdivided_x, 0, smooth_clipped, color="darkblue" )
	ax.fill_between( subdivided_x, 0, smooth_deleted, color="blue" )
	ax.set_ylim( ymin=0 )



def plot_cigars( cigars, sequences, reverses, ax, reference_function ):

	patches = []

	right_limits = [0] * len(cigars)

	for cigar,sequence,reverse in zip(cigars,sequences,reverses):

		for j in range(len(cigars)):
			if right_limits[j] <= cigar[0][1]:
				line = j
				break

		clipped_cigar = [ (t,a,r) for (t,a,r) in cigar if t != 'S' ]

		right_limits[line] = cigar[-1][1] + 10

		ax.barh( -line, clipped_cigar[-1][1]-clipped_cigar[0][1]+1, height=1, left=clipped_cigar[0][1]-0.5, color="#c8c8c8", linewidth=0 )

		for (c_type,abs_pos,rel_pos) in cigar:
			if c_type == 'D':
				ax.barh( -line, 1, height=1, left=abs_pos-0.5, color="#505050", linewidth=0 )
			elif c_type == 'M' and sequence[rel_pos] != reference_function[abs_pos]:
				ax.barh( -line, 1, height=1, left=abs_pos-0.5, color=basepair_colors[sequence[rel_pos]], linewidth=0, alpha=0.5 )
			elif c_type == 'S' and sequence[rel_pos] != reference_function[abs_pos]:
				ax.barh( -line, 1, height=1, left=abs_pos-0.5, color=basepair_colors[sequence[rel_pos]], linewidth=0, alpha=0.1 )

		for (c_type,abs_pos,rel_pos) in cigar:
			if c_type == 'I':
				ax.barh( -line, 0.6, height=1, left=abs_pos-0.8, fill=False, linewidth=1 )

		ax.barh( -line, cigar[-1][1]-cigar[0][1]+1, height=1, left=cigar[0][1]-0.5, fill=False, linewidth=0 ) # linewidth changed to 0

		if reverse:
			patches.append( matplotlib.patches.Polygon( [(cigar[0][1]-0.5,-line),(cigar[0][1]-0.5,1-line),(cigar[0][1]-1.5,0.5-line)] ) )
			patches.append( matplotlib.patches.Polygon( [(cigar[-1][1]+0.5,1-line),(cigar[-1][1]+1,1-line),(cigar[-1][1]+0.5,0.5-line)] ) )
			patches.append( matplotlib.patches.Polygon( [(cigar[-1][1]+0.5,-line),(cigar[-1][1]+1,-line),(cigar[-1][1]+0.5,0.5-line)] ) )
		else:
			patches.append( matplotlib.patches.Polygon( [(cigar[-1][1]+0.5,-line),(cigar[-1][1]+0.5,1-line),(cigar[-1][1]+1.5,0.5-line)] ) )
			patches.append( matplotlib.patches.Polygon( [(cigar[0][1]-0.5,1-line),(cigar[0][1]-1,1-line),(cigar[0][1]-0.5,0.5-line)] ) )
			patches.append( matplotlib.patches.Polygon( [(cigar[0][1]-0.5,-line),(cigar[0][1]-1,-line),(cigar[0][1]-0.5,0.5-line)] ) )

	collection = matplotlib.collections.PatchCollection( patches, linewidths=0.5, edgecolors="black", facecolors="yellow" )
	ax.add_collection( collection )

	ax.set_ylim( ymin=( 1 - len([ l for l in right_limits if l > 0 ]) ), ymax=1 )
	ax.set_yticks([])



def  plot_region( region_chrom, region_center, region_left, region_right, plot_title ):

	region_string = "%s:%i-%i" % ( region_chrom, region_left, region_right )
	annotations = get_annotations( region_string )
	#	print( "region %s annotations: " % region_string, annotations )

	bams = [ parsed_arguments.tumor ]
	if parsed_arguments.control is not None:
		fig = plot.figure(figsize=(19.2, 10.8))
		bams.insert(0, parsed_arguments.control)
		rows = 7
		cols = 3
		h_ratios = [1,5,15,2,5,15]
		ax_indices = [1,2,4,5]
		if annotations:
			ax_indices += [7]
	else:
		fig = plot.figure(figsize=(19.2, 5.4))
		rows = 4
		cols = 3
		h_ratios = [1,5,15]
		ax_indices = [1,2]
		if annotations:
			ax_indices += [4]

	if annotations:
		rows += 2
		h_ratios += [2,1]

	h_ratios += [1]

	grid = gridspec.GridSpec( rows, cols, height_ratios=h_ratios, hspace=0,
							width_ratios=[1,42,1], wspace=0,
							left=0, right=1, bottom=0, top=1 )
	
	ax = [ plot.subplot( grid[i,1] ) for i in ax_indices ]

	reference_buffer = ReferenceBuffer( parsed_arguments.ref, region_chrom )
	visible_basepairs = [ reference_buffer[i] for i in range( region_left, region_right + 1 ) ]

	for idx, bam in enumerate(bams):
		samtools_call   = ( parsed_arguments.samtoolsbin, "view", "-q", "20", "-F", "1024", bam, region_string )
		samtools_output = subprocess.check_output( samtools_call )
		samtools_output = [ line.split('\t') for line in samtools_output.split('\n') ]
		samtools_reads = [ line for line in samtools_output if len(line) > 5 ]

		plot_histogram( [ parse_cigar( read[5], int(read[3]) ) for read in samtools_reads if read[5] != "*" ], ax[idx*2] )
		plot_cigars( [ parse_cigar( read[5], int(read[3]) ) for read in samtools_reads if read[5] != "*"  ],
	    	         [ read[9] for read in samtools_reads if read[5] != "*"  ],
	        	     [ bool(int(read[1])&0x10) for read in samtools_reads if read[5] != "*"  ],
	            	 ax[idx*2+1], reference_buffer )

		if len(bams) == 2:
			if idx == 0:
				ax[idx*2].set_title( "%s - control" % plot_title )
			else:
				ax[idx*2].set_title( "%s - tumor" % plot_title)
		else:
			ax[idx*2].set_title( "%s - tumor" % plot_title )
		ax[idx*2].set_xticks([])
		ax[idx*2].yaxis.set_tick_params( labelleft=True, labelright=True )
		ax[idx*2].ticklabel_format( style='plain', axis='x', useOffset=False )

		ax[idx*2+1].xaxis.set_tick_params( width=0 )
		ax[idx*2+1].set_xticks([ i for i in range( region_left, region_right + 1 ) ])
		ax[idx*2+1].xaxis.set_ticklabels( visible_basepairs )

		for tick in ax[idx*2+1].get_xticklabels():
			tick.set_color( basepair_colors[tick._text] )

		for axis in ax[idx*2:(idx+1)*2]:
			axis.axvline( region_center - 0.5, color="black", linewidth=0.5 )
			axis.axvline( region_center + 0.5, color="black", linewidth=0.5 )
			for x in range( 10, parsed_arguments.window, 10 ):
				axis.axvline( region_center - x - 0.5, color="black", linewidth=0.25 )
				axis.axvline( region_center + x + 0.5, color="black", linewidth=0.25 )

	for axis in ax:
		axis.set_xlim( xmin=region_left-0.5, xmax=region_right+0.5 )

	if annotations:
		ax[-1].set_title( "annotations from " + parsed_arguments.annotations.split('/')[-1] )
		ax[-1].set_xticks([])
		ax[-1].set_yticks([])
		ax[-1].axis('off')

		for ann in annotations:

			ann[1]=max( ann[1], region_left )
			ann[2]=min( ann[2], region_right )

			ax[-1].barh( 0, ann[2]-ann[1]+1, height=1, left=ann[1]-0.5, color="#c8c8c8" )
			ax[-1].text( float( ann[1] + ann[2] ) / 2.0, 0.5, ann[0], ha='center', va='center' )



if parsed_arguments.region:

	region_chrom = parsed_arguments.region.split(':')[0]

	if len(parsed_arguments.region.split('-')) == 2:

		region_left   = int(parsed_arguments.region.split(':')[1].split('-')[0])
		region_right  = int(parsed_arguments.region.split(':')[1].split('-')[1])
		region_center = ( region_left + region_right ) // 2

	else:

		region_center = int(parsed_arguments.region.split(':')[1])
		region_left   = region_center - parsed_arguments.window
		region_right  = region_center + parsed_arguments.window

	plot_title = "%s:%s" % ( region_chrom, region_center )

	print region_chrom, region_left, region_center, region_right

	plot_region( region_chrom, region_center, region_left, region_right, plot_title )

	plot.savefig( "%s%s_%i.pdf" % ( parsed_arguments.prefix, region_chrom, region_center ) )
	plot.clf()
	plot.cla()
	plot.close()



if parsed_arguments.vcf:

	vcf_columns = {}

	for line in open(parsed_arguments.vcf, 'r' ):
		if line[:1] == "#":

			for i,col in enumerate( line[1:].rstrip('\n').split('\t') ):

				vcf_columns[col] = i

		elif line[:1] != "#" and "CHROM" in vcf_columns and "POS" in vcf_columns:

			region_chrom  = line.split('\t')[ vcf_columns["CHROM"] ]
			region_center = int( line.split('\t')[ vcf_columns["POS"] ] )
			region_left   = region_center - parsed_arguments.window
			region_right  = region_center + parsed_arguments.window

			plot_title = "%s:%s" % ( region_chrom, region_center )
			if "GENE" in vcf_columns and line.split('\t')[ vcf_columns["GENE"] ] != '.':
				plot_title += " variant name %s" % line.split('\t')[ vcf_columns["GENE"] ].rstrip('\n')		# .rstrip('\n') added to avoid new line in title
			if "EXONIC_CLASSIFICATION" in vcf_columns and line.split('\t')[ vcf_columns["EXONIC_CLASSIFICATION"] ] != '.':
				plot_title += " with exonic classification %s" % line.split('\t')[ vcf_columns["EXONIC_CLASSIFICATION"] ].rstrip('\n')		# .rstrip('\n') added to avoid new line in title

			plot_region( region_chrom, region_center, region_left, region_right, plot_title )

			plot.savefig( "%s%s_%i.pdf" % ( parsed_arguments.prefix, region_chrom, region_center ) )
			plot.clf()
			plot.cla()
			plot.close()



if parsed_arguments.bed:

	for line in open(parsed_arguments.bed, 'r' ):

		if line[:1] != "#":

			region_chrom  = line.split('\t')[0]
			region_left   = int(line.split('\t')[1])
			region_right  = int(line.split('\t')[2])
			region_center = ( region_left + region_right ) // 2

			plot_title = "%s:%s" % ( region_chrom, region_center )

			plot_region( region_chrom, region_center, region_left, region_right, plot_title )

			plot.savefig( "%s%s_%i.pdf" % ( parsed_arguments.prefix, region_chrom, region_center ) )
			plot.clf()
			plot.cla()
			plot.close()

