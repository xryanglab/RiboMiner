#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
@Author: Li Fajin
@Date: 2019-08-22 11:42:42
@LastEditors: Li Fajin
@LastEditTime: 2019-09-04 20:53:15
@Description: This script is used for statistic the length distribution of sequence reads based on a fastq file.
'''



import sys
import numpy as np
from itertools import groupby
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from optparse import OptionParser
from .__init__ import __version__


def create_parse_for_plot_reads_length():
	'''argument parser'''
	usage="usage: python %prog [options]"
	parser=OptionParser(usage=usage,version=__version__)
	parser.add_option('-i','--input',action="store",type="string",dest="fastqFile",help="Sequence fastq file.[required]")
	parser.add_option('-o',"--otput_prefix",action="store",type="string",dest="output_prefix",help="Prefix of output files.[required]")
	return parser

def fq2seqDict(fqFile):
	'''
	This function is used to get a dict of transcript sequence based on fastq file
	'''
	fastaDict={}
	f=open(fqFile,'r')
	faiter=(x[1] for x in groupby(f,lambda line: line[0]=="@")) ## groupby returns a tuple (key, group)
	for header in faiter:
		read_name=header.__next__().strip("@").split(" ")[0]
		seq='\t'.join(s.strip() for s in faiter.__next__()).strip().split("\t")[0]
		fastaDict[read_name]=seq
	return fastaDict

def get_read_length(fqFile,output_prefix):
	'''
	This function is used for getting length of all reads from a sequence fastq file.
	'''
	lengths_list=[]
	lengths_dict=defaultdict(int)
	zeros=0
	read_sequences=fq2seqDict(fqFile)
	with open(output_prefix+"_reads_length.txt",'w') as f:
		f.write("%s\t%s\n" %("read_name","read_length"))
		for read in read_sequences.keys():
			read_seq=read_sequences[read]
			read_length=len(read_seq)
			if read_length == 0:
				zeros+=1
			lengths_dict[read_length]+=1
			lengths_list.append(read_length)
			f.write("%s\t%s\n" % (str(read),str(read_length)))
	lengths_list=np.array(lengths_list)
	print("The total reads number is: " + str(len(lengths_list)),file=sys.stderr)
	print("The zero length number is: " + str(zeros),file=sys.stderr)
	print("The mean of reads length is: " + str(lengths_list.mean()),file=sys.stderr)
	print("The sd of reads length is: " + str(lengths_list.std()),file=sys.stderr)
	return lengths_list,lengths_dict

def plot_reads_length(lengths_dict,output_prefix,text_font={"size":20,"family":"Arial","weight":"bold"}):
	'''
	plot length distribution
	'''
	labels=lengths_dict.keys()
	values=lengths_dict.values()
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
	# ax.set_xticks(np.arange(0,50,5))
	# ax.set_xticklabels((np.arange(0,50,5)))
	ax.tick_params(which="both",width=2,labelsize=10)
	plt.tight_layout()
	plt.savefig(output_prefix+"_reads_length.pdf")
	plt.close()

def main():
	parser=create_parse_for_plot_reads_length()
	(options,args)=parser.parse_args()
	if not options.fastqFile or not options.output_prefix:
		raise IOError("Please reset your parameters!")
	print("Start the step of length statistics...",file=sys.stderr)
	lengths_list,lengths_dict=get_read_length(options.fastqFile,options.output_prefix)
	plot_reads_length(lengths_dict,options.output_prefix,text_font={"size":15,"family":"Arial","weight":"bold"})
	print("Finish the step of length statistics!",file=sys.stderr)

if __name__=="__main__":
	main()