#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
@Author: Li Fajin
@Date: 2019-08-18 21:19:13
@LastEditors: Li Fajin
@LastEditTime: 2019-08-30 17:10:01
@Description: This script is used for calculating RPFdist vaules of each transcript
RPFdist=(read counts in 5UTR)/(read counts in CDS region) or
RPFdist=(density in 5UTR)/(density in CDS region)
'''


from .FunctionDefinition import *


def RPFdist(in_bamFile,in_selectTrans,in_transLengthDict,in_startCodonCoorDict,in_stopCodonCoorDict,inCDS_lengthFilterParma,inCDS_countsFilterParma,in_excludeLengthParma,in_excludeCodonCountsParma,in_readLengths,in_readOffset,Mode):

	pysamFile=pysam.AlignmentFile(in_bamFile,"rb")
	pysamFile_trans=pysamFile.references
	in_selectTrans=set(pysamFile_trans).intersection(in_selectTrans)
	all_counts=0
	filter_1=0
	filter_2=0
	filter_3=0
	filter_4=0
	filter_5=0
	RPFdist_dict={}
	passTransSet=set()
	for trans in in_selectTrans:
		leftCoor =int(in_startCodonCoorDict[trans])-1
		rightCoor=int(in_stopCodonCoorDict[trans])-3
		(trans_counts,read_counts_frameSum,total_reads,cds_reads)=get_trans_frame_counts(pysamFile, trans, in_readLengths, in_readOffset, in_transLengthDict[trans], leftCoor, rightCoor)
		all_counts+=total_reads ## total_reads for transcript level

	for trans in in_selectTrans:
		leftCoor =int(in_startCodonCoorDict[trans])-1 #the first base of start codon 0-base
		rightCoor=int(in_stopCodonCoorDict[trans])-3 #the first base of stop codon 0-base
		# cdsLength should obey : (1) >= minFilterLength codon unit (2) CDSlength=3n
		if ((rightCoor-leftCoor) < inCDS_lengthFilterParma*3):
			filter_1+=1
			continue
		if (rightCoor-leftCoor)%3 !=0:
			filter_2+=1
			continue
		# get transNum and trans 5Pos offset counts vector
		(trans_counts,read_counts_frameSum,total_reads,cds_reads)=get_trans_frame_counts(pysamFile, trans, in_readLengths, in_readOffset, in_transLengthDict[trans], leftCoor, rightCoor)
		if total_reads==0 or in_transLengthDict[trans]==0:
			filter_3+=1
			continue
		cds_reads_normed=10**9*(cds_reads/(all_counts*len(read_counts_frameSum)))
		# cds num
		if Mode=="RPKM":
			if cds_reads_normed < inCDS_countsFilterParma:
				filter_4+=1
				continue
		elif Mode=="counts":
			if cds_reads < inCDS_countsFilterParma:
				filter_4+=1
				continue
		else:
			raise KeyError("There is no such modes, please check it.['RPKM' or 'counts']")

		# normalize the counts for per transcript
		normValue=np.mean(read_counts_frameSum[int(in_excludeLengthParma):])
		sumValue=np.sum(read_counts_frameSum[int(in_excludeLengthParma):])
		sumValue_normed=10**9*(sumValue/(all_counts*len(read_counts_frameSum)))
		if Mode == 'RPKM':
			if normValue == 0 or sumValue_normed < in_excludeCodonCountsParma:
				filter_5+=1
				continue
		elif Mode == 'counts':
			if normValue == 0 or sumValue < in_excludeCodonCountsParma:
				filter_5+=1
				continue
		else:
			raise KeyError("There is no such modes, please check it.['RPKM' or 'counts']")
		Five_UTR_counts=trans_counts[:leftCoor]
		CDS_counts=trans_counts[leftCoor:rightCoor]
		trans_counts_normed=10**9*(trans_counts/(all_counts*len(trans_counts)))
		Five_UTR_counts_normed=trans_counts_normed[:leftCoor]
		CDS_counts_normed=trans_counts_normed[leftCoor:rightCoor]
		if Mode == 'RPKM':
			RPFdist_dict[trans]=np.sum(Five_UTR_counts_normed)/np.sum(CDS_counts_normed)
			passTransSet.add(trans)
		elif Mode == 'counts':
			RPFdist_dict[trans]=np.sum(Five_UTR_counts)/np.sum(CDS_counts)
			passTransSet.add(trans)
		else:
			raise KeyError("There is no such modes, please check it.['RPKM' or 'counts']")
	pysamFile.close()
	print("Lenght filter(-l)---Transcripts number filtered by criterion one is : "+str(filter_1),file=sys.stderr)
	print("Lenght filter (3n)---Transcripts number filtered by criterion two is : "+str(filter_2),file=sys.stderr)
	print("Total counts filter---Transcripts number filtered by criterion three is : "+str(filter_3),file=sys.stderr)
	print("CDS density filter(RPKM-n or counts-n)---Transcripts number filtered by criterion four is : "+str(filter_4),file=sys.stderr)
	print("CDS density filter(normed-m)---Transcripts number filtered by criterion five is : "+str(filter_5),file=sys.stderr)
	print("Metaplots Transcript Number for bam file"+in_bamFile+" is :"+str(len(passTransSet)),file=sys.stderr)
	return RPFdist_dict

def write_bam_file_RPFdist_dataframe(inBamAttr,outFile):
	data=[]
	data_index=[]
	for bms in inBamAttr:
		d=bms.RPFdist_dict
		i=bms.bamLegend
		data.append(d)
		data_index.append(i)
	data=pd.DataFrame(data,index=data_index)
	data=data.T
	data.to_csv(outFile,sep="\t")


def main():
	parsed=create_parser_for_RPFdist()
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

	for bamfs in bam_attr:
		(bamfs.RPFdist_dict) = RPFdist(bamfs.bamName,select_trans,transLengthDict,startCodonCoorDict,stopCodonCoorDict,options.min_cds_codon,
		options.min_cds_counts,options.norm_exclude_codon,options.min_norm_region_counts,bamfs.bamLen,bamfs.bamOffset,options.mode)
		print("Finish the step of RPFdist calculation",file=sys.stderr)
	## write density
	write_bam_file_RPFdist_dataframe(bam_attr,options.output_prefix+"_RPFdist.txt")
	print("Finish the step of write_bam_file_RPFdist_dataframe",file=sys.stderr)

if __name__ == "__main__":
    main()