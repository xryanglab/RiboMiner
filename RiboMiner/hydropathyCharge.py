#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
@Author: Li Fajin
@Date: 2019-08-18 22:24:43
@LastEditors: Li Fajin
@LastEditTime: 2019-08-30 17:03:21
@Description: The script is used for hydropathy or charge calculaiton for a specific region
input:
1) cds sequences with fasta format you are interested in. Different fasta files would be separated by comma.
2) hydropathy index file download from AAindex. Three columns: AA(single character), amino acids (three characters), and hydropathy, just lile this:
AA	amio_acids	hydropathy
A	Ala 	1.8
R	Arg 	-4.5
N	Asn 	-3.5
D	Asp 	-3.5
charge index file is just like:
AA	amio_acids	charges
A	Ala 	0
R	Arg 	1
N	Asn 	0
D	Asp 	-1
output:
1) hydropathy or charge values at each position along transcripts. One fasta file, one hydropathy or charge value file.
2) a meta results in a specific regions.
'''

from .FunctionDefinition import *
import re
from itertools import groupby




def get_hydropathy_or_charge_vector(index_dict,AA_seq):
	vector=np.zeros(len(AA_seq),dtype="float")
	for i in np.arange(len(AA_seq)):
		vector[i]+=index_dict[AA_seq[i]]
	return vector


def hydropathy_or_charge(index_dict,transcriptFile,in_extendRegionLengthParma,in_regionLengthParma,table):
	trans_seq_dict=fastaIter(transcriptFile)
	in_selectTrans=trans_seq_dict.keys()
	passTransSet=set()
	valuePerCodon={}
	startDensity=np.zeros(int(in_regionLengthParma+in_extendRegionLengthParma+1),dtype="float64")
	stopDensity =np.zeros(int(in_regionLengthParma+in_extendRegionLengthParma+1),dtype="float64")
	startNormedWindowsList=[]
	startWindowsExistPosList=[]
	stopNormedWindowsList=[]
	stopWindowsExistPosList=[]
	for trans in in_selectTrans:
		cds_seq=trans_seq_dict[trans]
		AA_seq=translation(cds_seq,table=table,cds=False)
		vector=get_hydropathy_or_charge_vector(index_dict,AA_seq)
		###get trans density vector
		(tmpStartWin,tmpStartPos)=getWindowsVector(in_extendRegionLengthParma,in_regionLengthParma,vector,0) #start codon coor is 0 (0-based)
		(tmpStopWin, tmpStopPos) =getWindowsVector(in_regionLengthParma,in_extendRegionLengthParma,vector,(len(vector)-1)) #stop codon coor is len-1 (0-based)
		startNormedWindowsList.append(tmpStartWin)
		startWindowsExistPosList.append(tmpStartPos)
		stopNormedWindowsList.append(tmpStopWin)
		stopWindowsExistPosList.append(tmpStopPos)
		passTransSet.add(trans)
		valuePerCodon[trans]=vector
	###norm per position by mean or median
	startNormedWindowsList=np.array(startNormedWindowsList)
	startWindowsExistPosList=np.array(startWindowsExistPosList)
	stopNormedWindowsList=np.array(stopNormedWindowsList)
	stopWindowsExistPosList=np.array(stopWindowsExistPosList)
	#
	for terms in range(in_extendRegionLengthParma+in_regionLengthParma+1):
		startDensity[terms]=np.mean(startNormedWindowsList[np.where(startWindowsExistPosList[:,terms]==1),terms])
		stopDensity[terms] =np.mean(stopNormedWindowsList[np.where(stopWindowsExistPosList[:,terms]==1),terms])
	return startDensity,stopDensity,valuePerCodon


def write_hydropathy_or_charge_dataframe(inFastaAttr,outFile):
	data=[]
	for fasta in inFastaAttr:
		k=pd.DataFrame([fasta.fastaLegend]*len(fasta.startDensity))
		start=pd.DataFrame(fasta.startDensity)
		stop=pd.DataFrame(fasta.stopDensity)
		values=pd.merge(start,stop,how="left",left_index=True,right_index=True)
		values=pd.merge(k,values,how="left",left_index=True,right_index=True)
		data.append(values)
	temp=data[0]
	if len(data) < 1:
		raise EOFError("Empty file, there is nothing in the file.")
	if len(data) == 1:
		temp.columns=["sample","start_codon","stop_codon"]
		temp.to_csv(outFile,sep="\t",index=0)
	else:
		for i in np.arange(1,len(data)):
			temp=np.vstack((temp,data[i]))
		temp=pd.DataFrame(temp,columns=["sample","start_codon","stop_codon"])
		temp.to_csv(outFile,sep="\t",index=0)

def write_hydropathy_or_charge__per_codon(inFastaAttr,outFile):
	for fasta in inFastaAttr:
		with open(outFile+"_"+fasta.fastaLegend+"_values_at_each_position.txt",'w') as f:
			f.write("%s\t%s\n" %("transcripts","values"))
			for trans,perCodon in fasta.valuePerCodon.items():
				f.write("%s\t" %(trans))
				for pos in range(len(perCodon)):
					f.write("%s\t" %(str(perCodon[pos])))
				f.write("\n")

def main():
	"""main program"""
	parsed=create_parser_for_hydropathy_or_charge()
	(options,args)=parsed.parse_args()
	(transcriptFiles, output_prefix, upstream_codon, downstream_codon,index,trans_file_legend,table) = (options.transcriptFiles.strip().split(','),
	options.output_prefix,options.upstream_codon,options.downstream_codon,options.index,options.trans_file_legend.strip().split(','),options.geneticCode)
	if not index:
		raise IOError("Please input the hydropathy or charge index file")
	print("your input : "+ str(len(transcriptFiles))+" transcript files",file=sys.stderr)
	## handle bam file attr
	fasta_attr=[]
	for ii,jj in zip(transcriptFiles,trans_file_legend):
		fasta=fasta_attrbution(ii,jj)
		fasta_attr.append(fasta)
	data=pd.read_csv(index,sep="\t")
	index_dict={i:j for i,j in zip(data.iloc[:,0],data.iloc[:,2])}
	print("Start calculation...",file=sys.stderr)
	for fasta in fasta_attr:
		(fasta.startDensity,fasta.stopDensity,fasta.valuePerCodon) = hydropathy_or_charge(index_dict,fasta.fastaName,upstream_codon,downstream_codon,table)
	write_hydropathy_or_charge_dataframe(fasta_attr,output_prefix+"_values_dataframe.txt")
	write_hydropathy_or_charge__per_codon(fasta_attr,output_prefix)
	print("Finish the calculation!",file=sys.stderr)




if __name__ == "__main__":
		main()













