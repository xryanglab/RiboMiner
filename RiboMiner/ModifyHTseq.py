#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
Usage: ModifyHTseq -i bamFile -g gtfFile -o countsFile -t exon -m union -q 10 --minLen 25 --maxLen 35 --exclude-first 45 --exclude-last 15 --id-type gene_id
Notes: This script only used for strand specific library.
'''


from __future__ import division
import numpy
import sys
import pysam
import os
import itertools
import HTSeq
from optparse import OptionParser
from .__init__ import __version__
# __version__="0.1"


def create_parse_for_htseq():
	'''argument parser.'''
	usage="usage: python %prog [options]" + '\n' + __doc__ + "\n"
	parser=OptionParser(usage=usage,version=__version__)
	parser.add_option("-i","--input", action="store",type="string",default=None,dest="bamFile",
			help="Input file in bam format. default=%default")
	parser.add_option("-g","--gtfFile",action="store",type="string",dest="gtfFile",
			help="Standard GTF file of a specific species.")
	parser.add_option("-o","--outputFile",action="store",type="string",dest="outputFile",
			help="File name of output files.")
	parser.add_option("-t","--type",action="store",type="string",dest="Type",default='exon',
			help="Feature type (3rd column in GFF file) to be used. [exon or CDS]")
	parser.add_option("-m","--mode",action="store",type="string",default='union',dest="mode",
			help="mode to handle reads overlapping more than one feature, the same as htseq-count [union,intersection-strict,intersection-nonempty]. default=%default")
	parser.add_option("-q","--min-quality",action="store",type="int",dest="minQuality",default=10,
			help="The minimum quality of base to be required! default=%default")
	parser.add_option('--minLen',action="store",type="int",dest="minLen",default=25,
			help="The minimum length of reads to be considered. default=%default(nt)")
	parser.add_option('--maxLen',action="store",type="int",dest="maxLen",default=35,
			help="The max length of reads to be considered. default=%default(nt)")
	parser.add_option('--exclude-first',action="store",type="int",dest="excludeFirst",default=45,
			help="The number of nucleotides need to be excluded from start codon. default=%default(nt)=15(codon)")
	parser.add_option('--exclude-last',action="store",type="int",dest="excludeLast",default=15,
			help="The number of nucleotides need to be excluded from stop codon. default=%default(nt)=5(codon)")
	parser.add_option('--id-type',action="store",type="string",dest="id_type",default="gene_id",
			help="define the id type users input. the default is gene id, the same as '-i' in htseq-count. default=%default")
	return parser

def modifHTSeq(bam_filename, gff_filename, out_file, overlap_mode, feature_type, id_attribute, minaqual, exclude_start_distance, exclude_stop_distance, min_len, max_len):
	#feature GenomicArrayOfSets
	features = HTSeq.GenomicArrayOfSets( "auto", stranded=True)
	counts = {}
	start_codon_sites = {}
	stop_codon_sites = {}
	#GTF
	gff = HTSeq.GFF_Reader( gff_filename, end_included=True)
	i = 0
	for f in gff:
		if f.type == feature_type:
			if id_attribute in f.attr : #the same to the f.attr.keys()
				feature_id = f.attr[ id_attribute ] # f.attr will return the 9-th colum of the input gtf file as {}
			else:
				feature_id = f.attr[ 'gene_id' ] #in the gtf file of Rat, there are some CDS/exon dont have gene_name ,but every items have gene_id
			features[ f.iv ] += feature_id #label the chrmosome with gene_name, if dont have gene_name,replaced by gene_id
			#counts[ f.attr[ id_attribute ] ] = 0 #only counts reads for genes with id_attribute, so cant repaced by counts[ feature_id ] = 0
			counts[ feature_id ] = 0
		### if there are multiple TIS, use the most 5' end start codon and the most 3' end stop codon
		if f.type == "start_codon":
			if id_attribute in f.attr :
				gname=f.attr[ id_attribute ]
			if gname not in start_codon_sites :
				start_codon_sites[gname] = f.iv.start_d
			else :
				if f.iv.strand == "+":
					start_codon_sites[gname] = min(f.iv.start_d, start_codon_sites[gname])
				else :
					start_codon_sites[gname] = max(f.iv.start_d, start_codon_sites[gname])
		#
		if f.type == "stop_codon":
			if id_attribute in f.attr :
				gname=f.attr[ id_attribute ]
			if gname not in stop_codon_sites :
				stop_codon_sites[gname] = f.iv.end_d
			else :
				if f.iv.strand == "+":
					stop_codon_sites[gname] = max(f.iv.end_d, stop_codon_sites[gname])
				else :
					stop_codon_sites[gname] = min(f.iv.end_d, stop_codon_sites[gname])
		i += 1
		if i % 100000 == 0 :
			sys.stderr.write( "%d GFF lines processed.\n" % i )
	#bam
	read_seq = HTSeq.BAM_Reader( bam_filename )
	#counts
	empty = 0
	ambiguous = 0
	notaligned = 0
	lowqual = 0
	nonunique = 0
	i = 0
	for r in read_seq:
		if i > 0 and i % 100000 == 0 :
			sys.stderr.write( "%d SAM alignment record processed.\n" % i )
		i += 1
		if not r.aligned:
			notaligned += 1
			continue
		if r.optional_field( "NH" ) > 1:
			nonunique += 1
			continue
		if r.aQual < minaqual:
			lowqual += 1
			continue
		###
		if len(r.read.seq) < min_len or len(r.read.seq) > max_len:
			continue
		iv_seq = ( co.ref_iv for co in r.cigar if co.type == "M" and co.size > 0 )
		if overlap_mode == "union":
			fs = set()
			for iv in iv_seq:
				for iv2, fs2 in features[ iv ].steps():
					fs = fs.union( fs2 )
		elif overlap_mode == "intersection-strict" or overlap_mode == "intersection-nonempty":
			fs = None
			for iv in iv_seq:
				for iv2, fs2 in features[ iv ].steps():
					if len(fs2) > 0 or overlap_mode == "intersection-strict":
						if fs is None:
							fs = fs2.copy()
						else:
							fs = fs.intersection( fs2 )
		else:
			sys.exit( "Illegal overlap mode." )
		if fs is None or len( fs ) == 0:
			empty += 1
		elif len( fs ) > 1:
			ambiguous += 1
		else:
			try : #some genes may dont have start or stop codon
				if abs(start_codon_sites[ list(fs)[0] ] - r.iv.start_d) < exclude_start_distance :
					continue
				elif abs(r.iv.end_d - stop_codon_sites[ list(fs)[0] ]) < exclude_stop_distance:
					continue
				else :
					counts[ list(fs)[0] ] += 1
			except :
				counts[ list(fs)[0] ] += 1
	#output
	with open(out_file,"w") as fout:
		fout.write("%s\t%s\n" %(id_attribute.strip(),"count"))
		for fn in sorted( counts.keys() ):
			fout.write("%s\t%s\n" % (fn,counts[fn]))
		fout.write("__no_feature\t%d\n" % empty)
		fout.write("__ambiguous\t%d\n" % ambiguous)
		fout.write("__too_low_aQual\t%d\n" % lowqual)
		fout.write("__not_aligned\t%d\n" % notaligned)
		fout.write("__alignment_not_unique\t%d\n" % nonunique)


def main():
    parser=create_parse_for_htseq()
    (options,args)=parser.parse_args()
    (bamFile,gtfFile,outFile,modeType,featureType,idAttribute,minQual,excludeStartDistance,excludeStopDistance,minLen,maxLen)=(options.bamFile,options.gtfFile,options.outputFile,
    options.mode,options.Type,options.id_type,options.minQuality,options.excludeFirst,options.excludeLast,options.minLen,options.maxLen)
    print("Start read counting...",file=sys.stderr)
    modifHTSeq(bamFile,gtfFile,outFile,modeType,featureType,idAttribute,minQual,excludeStartDistance,excludeStopDistance,minLen,maxLen)
    print("Finish read counting!",file=sys.stderr)

if __name__=="__main__":
    main()
