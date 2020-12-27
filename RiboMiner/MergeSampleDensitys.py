#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
usage: python MergeSampleDensity.py inputFile outputFile
1) inputFile: input files separated by comma. e.g. test1_dataframe.txt,test2_dataframe.txt
2) outputFile: an output file after merging.
'''
import os
import sys
import numpy as np
import pandas as pd
from optparse import OptionParser
from .__init__ import __version__

def create_parser_for_merge_density():
	'''argument parser'''
	usage="usage: python %prog [options]"
	parser=OptionParser(usage=usage,version=__version__)
	parser.add_option("-i","--input",action="store",type="string",dest="densityFiles",help="Density files in txt format separated by comma. e.g. test1_dataframe.txt,test2_dataframe.txt")
	parser.add_option("-o","--output",action="store",type="string",dest="outputFile",help="Output filename.[required]")
	return parser

def MergeSampleData(inputFile,outputFile):
    inputFiles=inputFile.strip().split(',')
    data=[pd.read_csv(File,sep="\t",header=0) for File in inputFiles]
    Final_data=pd.concat(data,axis=0)
    Final_data.to_csv(outputFile,sep="\t",index=0)

def main():
    parser=create_parser_for_merge_density()
    (options,args)=parser.parse_args()
    print("Start merging...",file=sys.stderr)
    MergeSampleData(options.densityFiles,options.outputFile)
    print("Finish the files merging!",file=sys.stderr)
if __name__=="__main__":
    main()