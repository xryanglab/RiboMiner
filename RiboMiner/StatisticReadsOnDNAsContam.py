#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
@Author: Li Fajin
@Date: 2019-08-22 16:58:19
@LastEditors: Li Fajin
@LastEditTime: 2019-09-17 15:48:56
@Description: This script is used for statistic reads mapped to DNA  based on bam files mapped to transcriptome.
Becaused the reads mapped to DNA may be contaminations. Code part from Xiao Zhengzao.
'''
import sys
import HTSeq
from collections import Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from optparse import OptionParser
from .__init__ import __version__

def create_parse_for_DNA_mapped_reads():
	'''argument parser'''
	usage="usage: python %prog [options]"
	parser=OptionParser(usage=usage,version=__version__)
	parser.add_option('-i','--input',action="store",type="string",dest="bamFile",help="Input files mapped to transcriptome with bam format. [required]")
	parser.add_option("-g","--gtfFile",action="store",type="string",dest="gtfFile",help="geome annotation file with gtf format.[required]")
	parser.add_option('-o',"--otput_prefix",action="store",type="string",dest="output_prefix",help="Prefix of output files.[required]")
	return parser

def get_gene_features(gtfFile,id_type,feature_type):
	'''
	get exon features and gene interval features
	'''
	features=HTSeq.GenomicArrayOfSets( "auto", stranded = "yes" )
	geneFeatures=HTSeq.GenomicArrayOfSets( "auto", stranded = "yes" )
	geneRange={}
	gtf=HTSeq.GFF_Reader(gtfFile)
	i=0
	for line in gtf:
		if line.type == feature_type:
			feature_id=line.attr[id_type]
			features[line.iv]+=feature_id
			if feature_id not in geneRange:
				geneRange[feature_id]=[line.iv.chrom,0,0,line.iv.strand]
			if geneRange[feature_id][1] != 0:
				geneRange[feature_id][1]=min(geneRange[feature_id][1],line.iv.start)
			else:
				geneRange[feature_id][1]=line.iv.start
			geneRange[feature_id][2]=max(geneRange[feature_id][2],line.iv.end)
		i+=1
		if i % 100000 == 0:
			print("%d GFF lines processed.\n" % i,file=sys.stderr)
	for g,v in geneRange.items():
		chrom,start,end,strand=v
		tmp_iv=HTSeq.GenomicInterval(chrom,start,end,strand)
		geneFeatures[tmp_iv]+=g
	return features,geneFeatures

def statistic_mapped_reads(bamFile,gtfFile,id_type,feature_type,output_prefix):
	'''
	statistic the reads number mapped to RNA, DNA or introns
	'''
	uniqueGene=0
	uniqueNotGene=0
	uniqueIntron=0
	uniqueAmbiguous=0
	RNA_len=Counter()
	DNA_len=Counter()
	Intron_len=Counter()
	bam=HTSeq.BAM_Reader(bamFile)
	features,geneFeatures=get_gene_features(gtfFile,id_type,feature_type)
	i=0
	for r in bam:
		i+=1
		if i%1000000==0:
			print("%d BAM alignments records processed.\n" % i,file=sys.stderr)
		if not r.aligned:
			continue
		if r.optional_field("NH")==1:
			r.read.seq=r.read_as_aligned.seq
			iv_seq=(co.ref_iv for co in r.cigar if co.type=="M" and co.size >0)
			feature_sets=set()
			intron_featuresets=set()
			for iv in iv_seq:
				for iv2,fs in features[iv].steps():
					feature_sets=feature_sets.union(fs)
				for iv3,fs in geneFeatures[iv].steps():
					intron_featuresets=intron_featuresets.union(fs)
			if feature_sets is None or len(feature_sets) == 0:
				if intron_featuresets is None or len(intron_featuresets)==0:
					uniqueNotGene+=1
					DNA_len[len(r.read.seq)]+=1
				else:
					uniqueIntron+=1
					Intron_len[len(r.read.seq)]+=1
			elif len(feature_sets) ==1:
				uniqueGene+=1
				RNA_len[len(r.read.seq)]+=1
			else:
				uniqueAmbiguous+=1
				RNA_len[len(r.read.seq)]+=1
	with open(output_prefix+"_reads_distribution.txt",'w') as f:
		f.write("unique mapped reads of exon: %i\n" % uniqueGene)
		f.write("unique mapped reads of intergenic region: %i\n" % uniqueNotGene)
		f.write("unique mapped reads of intron: %i\n" % uniqueIntron)
		f.write("unique mapped ambiguous reads of RNA: %i\n" % uniqueAmbiguous)
	return DNA_len,RNA_len,Intron_len

def plot_reads_distribution(lengthDict,output_prefix,text_font={"size":20,"family":"Arial","weight":"bold"}):
	'''
	plot the reads length distribution
	'''
	labels=lengthDict.keys()
	values=lengthDict.values()
	plt.rc('font',weight='bold')
	fig=plt.figure(figsize=(5,4))
	ax=fig.add_subplot(111)
	plt.bar(labels,values,color="b",width=0.5,alpha=0.9)
	ax.spines["top"].set_linewidth(2)
	ax.spines["right"].set_linewidth(2)
	ax.spines["bottom"].set_linewidth(2)
	ax.spines["left"].set_linewidth(2)
	ax.set_xlabel("Length of reads",fontdict=text_font)
	ax.set_ylabel("Read Counts",fontdict=text_font)
	ax.tick_params(which="both",width=2,labelsize=10)
	plt.tight_layout()
	plt.savefig(output_prefix+".pdf")
	plt.close()


def main():
	parser=create_parse_for_DNA_mapped_reads()
	(options,args)=parser.parse_args()
	if not options.bamFile or not options.gtfFile or not options.output_prefix:
		raise IOError("Please reset your parameters!")
	id_type="gene_id"
	feature_type="exon"
	DNA_len,RNA_len,Intron_len=statistic_mapped_reads(options.bamFile,options.gtfFile,id_type,feature_type,options.output_prefix)
	print("Start the step of plot...",file=sys.stderr)
	plot_reads_distribution(DNA_len,options.output_prefix+"_DNA",text_font={"size":15,"family":"Arial","weight":"bold"})
	plot_reads_distribution(RNA_len,options.output_prefix+"_RNA",text_font={"size":15,"family":"Arial","weight":"bold"})
	plot_reads_distribution(Intron_len,options.output_prefix+"_Intron",text_font={"size":15,"family":"Arial","weight":"bold"})
	print("Finish the step Reads statistic!",file=sys.stderr)

if __name__=="__main__":
	main()
