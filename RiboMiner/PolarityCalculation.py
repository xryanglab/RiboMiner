#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
@Author: Li Fajin
@Date: 2019-08-14 17:46:02
@LastEditors: Li Fajin
@LastEditTime: 2019-08-30 17:35:12
@Description: file content
'''


from .FunctionDefinition import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.colors as cols
import seaborn as sns


def calculate_polarity(in_bamFile,in_selectTrans,in_transLengthDict,in_startCodonCoorDict,in_stopCodonCoorDict,in_readLengths,in_readOffset,min_cds_counts):
	""" calculate the polarity of every gene (the longest transcript as for the final gene)
		1) calculate the density of ribosome at all positions of a gene , d
		2) calculate the normalized distance from the center of a gene at position i, w
		3) calculate the polarity at all positions of a gene ,pi and sum up the pi as the final polarity of the gene

		pi = di*wi/sum(di)
		wi = (2*i-(len(gene)+1))/(len(gene)-1)
	"""
	pysamFile=pysam.AlignmentFile(in_bamFile,"rb")
	pysamFile_trans=pysamFile.references
	in_selectTrans=set(pysamFile_trans).intersection(in_selectTrans)
	GeneNum = len(in_selectTrans)
	polarityScore=[]
	polarityScorePerTrans=[]
	passTransSet={}
	filter=0
	for trans in in_selectTrans:
		GeneNum-+1
		leftCoor =int(in_startCodonCoorDict[trans])-1 #the first base of start codon 0-base
		rightCoor=int(in_stopCodonCoorDict[trans])-3 #the first base of stop codon 0-base
		(read_counts,read_counts_frameSum,trans_reads,cds_reads)=get_trans_frame_counts(pysamFile, trans, in_readLengths, in_readOffset, in_transLengthDict[trans], leftCoor, rightCoor)
		## get read density in nt unit
		read_counts_frame=read_counts[(leftCoor+15):(rightCoor-15)] ##  exclude the first 5 and last f codons
		tmpDistanceNormed=np.arange(1,len(read_counts_frame)+1)
		DistanceNormed=(2*tmpDistanceNormed-(len(read_counts_frame)+1))/(len(read_counts_frame)-1)
		sumValue=np.sum(read_counts_frame)
		if sumValue < min_cds_counts :
			filter+=1
			continue
		if sumValue == 0:
			tmpPolarityScorePerTrans=read_counts_frame
		else:
			tmpPolarityScorePerTrans=np.multiply(read_counts_frame,DistanceNormed)/np.sum(read_counts_frame)
		polarityScorePerTrans.append(tmpPolarityScorePerTrans)
		passTransSet[trans]=np.sum(tmpPolarityScorePerTrans)
	## calculate polarity of each gene
	for i in range(len(polarityScorePerTrans)):
		polarityScore.append(np.sum(polarityScorePerTrans[i]))
	pysamFile.close()
	print("Metaplots Transcript Number for bam file"+in_bamFile+" is :"+str(len(passTransSet)),file=sys.stderr)
	print("The number of genes filtered by read counts equals to "+ str(min_cds_counts)+" is: " + str(filter),file=sys.stderr)
	print("There are "+str(len(polarityScore))+" genes used for plot polarity value",file=sys.stderr)
	return polarityScore,passTransSet

def DrawPolarity(data,inOutPrefix):
	"""plot polarity scores"""
	text_font={"size":20,"family":"Arial","weight":"bold"}
	legend_font={"size":10,"family":"Arial","weight":"bold"}
	samples=np.unique(data.columns)
	plt.rc('font',weight='bold')
	if len(samples) <=8:
		colors=["b","orangered","green","c","m","y","k","w"]
	else:
		colors=sns.color_palette('husl',len(samples))
	fig=plt.figure(figsize=(5,4))
	ax=fig.add_subplot(111)
	## try to use for loop to re-write the plot function
	for i in np.arange(len(samples)):
		lst=data.iloc[:,i].values
		lst=lst[~np.isnan(lst)]
		sns.distplot(lst,hist=False,rug=False,label=samples[i],color=colors[i])
	ax.set_xlabel("Polarity score",fontdict=text_font)
	ax.set_ylabel("Relative gene numbers",fontdict=text_font)
	ax.spines["top"].set_visible(False)
	ax.spines["right"].set_visible(False)
	ax.spines["bottom"].set_linewidth(2)
	ax.spines["left"].set_linewidth(2)
	ax.tick_params(which="both",width=2,length=2)
	plt.legend(loc="best",prop=legend_font)
	plt.tight_layout()
	plt.savefig(inOutPrefix+"_polarity.pdf")
	plt.close()

def write_bam_file_polarity_dataframe(inBamAttr,outFile):
	data=[]
	data_index=[]
	for bms in inBamAttr:
		d=bms.passTransSet
		i=bms.bamLegend
		data.append(d)
		data_index.append(i)
	data=pd.DataFrame(data,index=data_index)
	data=data.T
	data.to_csv(outFile,sep="\t")
	return data

def parse_args_for_polarity_calculation():
	parsed=create_parser_for_polarity_calculation()
	(options,args)=parsed.parse_args()
	if options.bamListFile and (options.bam_files or options.read_length or options.read_offset or options.bam_file_legend):
		raise IOError("'-f' parameter and '-i -r -s -t' are mutually exclusive.")
	if options.bamListFile:
		bamFiles,readLengths,Offsets,bamLegends=parse_bamListFile(options.bamListFile)
	elif options.bam_files:
		bamFiles,readLengths,Offsets,bamLegends=options.bam_files.split(","),options.read_length.split("_"),options.read_offset.split("_"),options.bam_file_legend.split(",")
	else:
		raise IOError("Please check you input files!")
	print("your input : "+ str(len(bamFiles))+" bam files",file=sys.stderr)
	bam_attr=[]
	for ii,jj,mm,nn in zip(bamFiles,readLengths,Offsets,bamLegends):
		bam=bam_file_attr(ii,jj,mm,nn)
		bam_attr.append(bam)
	## calculate density for each bam files
	selectTrans,transLengthDict,startCodonCoorDict,stopCodonCoorDict,transID2geneID,transID2geneName,cdsLengthDict=reload_transcripts_information(options.coorFile)
	geneID2transID={v:k for k,v in transID2geneID.items()}
	geneName2transID={v:k for k,v in transID2geneName.items()}
	if options.in_selectTrans:
		select_trans=pd.read_csv(options.in_selectTrans,sep="\t")
		select_trans=set(select_trans.iloc[:,0].values)
		if options.id_type == 'transcript_id':
			select_trans=select_trans.intersection(selectTrans)
			print("There are " + str(len(select_trans)) + " transcripts from "+options.in_selectTrans+" used for following analysis.",file=sys.stderr)
		elif options.id_type == 'gene_id':
			tmp=[geneID2transID[gene_id] for gene_id in select_trans]
			select_trans=set(tmp)
			select_trans=select_trans.intersection(selectTrans)
			print("There are " + str(len(select_trans))+" gene id could be transformed into transcript id and used for following analysis.",file=sys.stderr)
		elif options.id_type == 'gene_name' or options.id_type=='gene_symbol':
			tmp=[geneName2transID[gene_name] for gene_name in select_trans]
			select_trans=set(tmp)
			select_trans=select_trans.intersection(selectTrans)
			print("There are " + str(len(select_trans))+" gene symbol could be transformed into transcript id and used for following analysis.",file=sys.stderr)
		else:
			raise IOError("Please input a approproate id_type parameters.[transcript_id/gene_id/gene_name/]")
	else:
		select_trans=selectTrans
	print("There are "+str(len(select_trans))+" transcripts with both start and stop codon will be used for following analysis.",file=sys.stderr)

	for bamfs in bam_attr:
		(bamfs.polarity,bamfs.passTransSet)=calculate_polarity(bamfs.bamName,select_trans,transLengthDict,startCodonCoorDict,stopCodonCoorDict,bamfs.bamLen,bamfs.bamOffset,options.min_cds_counts)
	print("Finish the step of calculate the polarity of genes",file=sys.stderr)
	data_polarity=write_bam_file_polarity_dataframe(bam_attr,options.output_prefix+"_polarity_dataframe.txt")
	print("Finish the step of write_bam_file_polarity_dataframe",file=sys.stderr)
	if options.plot.upper() in ['YES','Y']:
		DrawPolarity(data_polarity,options.output_prefix)
		print("Finish the step of plotting!",file=sys.stderr)
	elif options.plot.upper() in ['NO','N','NONE']:
		pass
	else:
		raise IOError("Please input a correct --plot parameter! [yes/no]")

def main():
	"""main program"""
	parse_args_for_polarity_calculation()

if __name__ == "__main__":
		main()

