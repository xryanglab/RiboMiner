#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
@Author: Li Fajin
@Date: 2019-08-21 20:55:35
@LastEditors: Li Fajin
@LastEditTime: 2019-08-30 17:01:33
@Description: This script is used for statistic GC contents of specific transcripts or GC contents on different reading frames of given transcripts.

Notice:
1) input any DNA or RNA sequences, and set --mode noraml. This will output GC contents of sequences you input.
2) input cds sequences and set --mode frames. This will output GC contents on different reading frames for those cds sequences.
'''

import sys
from itertools import groupby
from optparse import OptionParser
from .__init__ import __version__

def create_parser_for_GC_content():
	'''argument parser'''
	usage="usage: python %prog [options]"
	parser=OptionParser(usage=usage,version=__version__)
	parser.add_option("-i","--input", action="store",type="string",dest="sequences",
			help="Input file(s) in fasta format.")
	parser.add_option("-o","--otput_prefix",action="store",type="string",dest="output_prefix",
			help="Prefix of output files.[required]")
	parser.add_option("--mode",action="store",type="string",dest="mode",default="normal",
			help="The type of GC content you want to statistic. Either the normal type or GC content from each reading frame. [normal or frames]. defaul=%default")
	return parser

def fastaIter(transcriptFile):
	'''
	This function is used to get a dict of transcript sequence
	'''
	fastaDict={}
	f=open(transcriptFile,'r')
	faiter=(x[1] for x in groupby(f,lambda line: line.strip()[0]==">")) ## groupby returns a tuple (key, group)
	for header in faiter:
		geneName=header.__next__().strip(">").split(" ")[0]
		seq=''.join(s.strip() for s in faiter.__next__())
		fastaDict[geneName]=seq
	return fastaDict

def write_GC_content(sequences,output_prefix):
    sequences_dict=fastaIter(sequences)
    with open(output_prefix+"_GC_content.txt","w") as f:
        f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %('transcripts','A','T/U','C','G','GC','total','GC%'))
        for trans in sequences_dict.keys():
            seq=sequences_dict[trans]
            A,T_or_U,C,G,GC,total,GC_pct=GC_content_for_single_sequence(seq)
            f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(trans,str(A),str(T_or_U),str(C),str(G),str(GC),str(total),str(GC_pct)))
def GC_content_for_single_sequence(seq):
    A=seq.count('A')
    G=seq.count('G')
    C=seq.count('C')
    T_or_U=len(seq)-A-G-C
    GC=G+C
    total=len(seq)
    if total==0:
        raise IOError("There is empty sequence in your input file.")
    GC_pct=GC/total
    return A,T_or_U,C,G,GC,total,GC_pct

def GC_content_for_different_frame(sequences,output_prefix):
    ''' Notice: the input file must be cds sequences with fasta format.'''
    sequences_dict=fastaIter(sequences)
    i=0
    with open(output_prefix+"_GC_content_frames.txt","w") as f:
        f.write("%s\t%s\t%s\t%s\t%s\n" %('transcripts','frame0','frame1','frame2','frameSum'))
        for trans in sequences_dict.keys():
            seq=sequences_dict[trans]
            if len(seq) % 3!=0:
                seq = seq[:-(len(seq) % 3)]
                i+=1
            frame0=[seq[i] for i in range(0,len(seq),3)]
            frame1=[seq[i] for i in range(1,len(seq),3)]
            frame2=[seq[i] for i in range(2,len(seq),3)]
            frame0_GC_pct=GC_content_for_single_sequence(frame0)[-1]
            frame1_GC_pct=GC_content_for_single_sequence(frame1)[-1]
            frame2_GC_pct=GC_content_for_single_sequence(frame2)[-1]
            frameSum=GC_content_for_single_sequence(seq)[-1]
            f.write("%s\t%s\t%s\t%s\t%s\n" %(trans,frame0_GC_pct,frame1_GC_pct,frame2_GC_pct,frameSum))

def main():
    parser=create_parser_for_GC_content()
    (options,args)=parser.parse_args()
    if not options.sequences or not options.output_prefix:
        raise IOError("Please your input sequences and prefix of your output files.")
    if options.mode == 'normal':
        print("Start statistic GC content...",file=sys.stderr)
        write_GC_content(options.sequences,options.output_prefix)
        print("Finish the step of GC content!",file=sys.stderr)
    elif options.mode == 'frames':
        print("Start statistic GC content...",file=sys.stderr)
        GC_content_for_different_frame(options.sequences,options.output_prefix)
        print("Finish the step of GC content!",file=sys.stderr)
    else:
        raise IOError("Please reset your --mode parameter [normal/frames]")

if __name__=="__main__":
    main()


