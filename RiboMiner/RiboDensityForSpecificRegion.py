#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
@Author: Li Fajin
@Date: 2019-08-23 15:14:11
@LastEditors: Li Fajin
@LastEditTime: 2019-08-30 17:09:01
@Description:

This script is used for getting ribosome density at specific region. For example, if there are ribosomes enriched on codon 25 to codon 75, we could
use this script to calculate the mean read density at 25-75 codons based on which we could split genes into different gene sets [up regulated or down regulated or unblocked].
The input is almost the same as MetageneAnalysis.py without filtering, so the output is the mean density of all longest transcripts.
'''


from .FunctionDefinition import *


def RibosomeDensity_for_specific_region(in_bamFile,in_selectTrans,in_transLengthDict,in_startCodonCoorDict,in_stopCodonCoorDict,in_readLengths,in_readOffset,left_position,right_position,mode,unit):
		pysamFile=pysam.AlignmentFile(in_bamFile,'rb')
		pysamFile_trans=pysamFile.references
		in_selectTrans=set(pysamFile_trans).intersection(in_selectTrans)
		local_mean_density={}
		local_density={}
		all_counts=0
		for trans in in_selectTrans:
				leftCoor =int(in_startCodonCoorDict[trans])-1
				rightCoor=int(in_stopCodonCoorDict[trans])-3
				(trans_counts,read_counts_frameSum,total_reads,cds_reads)=get_trans_frame_counts(pysamFile, trans, in_readLengths, in_readOffset, in_transLengthDict[trans], leftCoor, rightCoor)
				all_counts+=total_reads
		for trans in in_selectTrans:
				leftCoor =int(in_startCodonCoorDict[trans])-1
				rightCoor=int(in_stopCodonCoorDict[trans])-3
				(trans_counts,read_counts_frameSum,total_reads,cds_reads)=get_trans_frame_counts(pysamFile, trans, in_readLengths, in_readOffset, in_transLengthDict[trans], leftCoor, rightCoor)
				local_cds_counts=np.zeros(int(right_position-left_position+1),dtype="float64")
				local_cds_counts_normed=np.zeros(int(right_position-left_position+1),dtype="float64")
				## codon level
				if unit == 'codon':
						read_counts_frameSum_normed=10**6*(read_counts_frameSum/all_counts) ## RPM
						temp_local_cds_counts=read_counts_frameSum[int(left_position-1):int(right_position)] ## could offer two parameters later on
						temp_local_cds_counts_normed=read_counts_frameSum_normed[int(left_position-1):int(right_position)]
						local_cds_counts[np.arange(len(temp_local_cds_counts))]=temp_local_cds_counts
						local_cds_counts_normed[np.arange(len(temp_local_cds_counts_normed))]=temp_local_cds_counts_normed
						if mode == 'RPKM':
								if np.sum(local_cds_counts_normed) == 0:
										tmp_local_cds_counts_normed=0
										local_mean_density[trans]=tmp_local_cds_counts_normed
								else:
										local_mean_density[trans]=np.mean(local_cds_counts_normed)
								local_density[trans]=local_cds_counts_normed
						if mode == 'counts':
								if np.sum(local_cds_counts) == 0:
										tmp_local_cds_counts=0
										local_mean_density[trans]=tmp_local_cds_counts
								else:
										local_mean_density[trans]=np.mean(local_cds_counts)
								local_density[trans]=local_cds_counts
				if unit == 'nt':
						# trans_counts_normed=10**9*(trans_counts/all_counts*len(trans_counts))
						trans_counts_normed=10**6*(trans_counts/all_counts)
						temp_local_cds_counts=trans_counts[int(left_position-1):int(right_position)] ## could offer two parameters later on
						temp_local_cds_counts_normed=trans_counts_normed[int(left_position-1):int(right_position)]
						local_cds_counts[np.arange(len(temp_local_cds_counts))]=temp_local_cds_counts
						local_cds_counts_normed[np.arange(len(temp_local_cds_counts_normed))]=temp_local_cds_counts_normed
						if mode == 'RPKM':
								if np.sum(local_cds_counts_normed) == 0:
										tmp_local_cds_counts_normed=0
										local_mean_density[trans]=tmp_local_cds_counts_normed
								else:
										local_mean_density[trans]=np.mean(local_cds_counts_normed)
								local_density[trans]=local_cds_counts_normed
						if mode == 'counts':
								if np.sum(local_cds_counts) == 0:
										tmp_local_cds_counts=0
										local_mean_density[trans]=tmp_local_cds_counts
								else:
										local_mean_density[trans]=np.mean(local_cds_counts)
								local_density[trans]=local_cds_counts

		return local_mean_density,local_density

def write_local_codon_units_density(inBamAttr,outFile,left_position,right_position):
	for bms in inBamAttr:
		data=pd.DataFrame(bms.local_density)
		data=data.T
		data.columns=["codon_"+str(i+1) for i in np.arange(int(right_position)+1-int(left_position))]
		data.to_csv(outFile+"_"+bms.bamLegend+"_local_density.txt",sep='\t')


def write_bam_file_local_mean_cds_counts_dataframe(inBamAttr,outFile):
	data=[]
	data_index=[]
	for bms in inBamAttr:
		d=bms.local_mean_density
		i=bms.bamLegend
		data.append(d)
		data_index.append(i)
	data=pd.DataFrame(data,index=data_index)
	data=data.T
	data.to_csv(outFile,sep="\t")

def parse_args_for_specific_region_metagene():
	parsed=creat_parser_for_specific_region()
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
			(bamfs.local_mean_density,bamfs.local_density)=RibosomeDensity_for_specific_region(bamfs.bamName,select_trans,transLengthDict,startCodonCoorDict,stopCodonCoorDict,
			bamfs.bamLen,bamfs.bamOffset,options.left_position,options.right_position,options.mode,options.unit)
	print("Finish the step of ribosomeDensityNormPerTrans",file=sys.stderr)
	## write density
	write_bam_file_local_mean_cds_counts_dataframe(bam_attr,options.output_prefix+"_local_mean_density.txt")
	write_local_codon_units_density(bam_attr,options.output_prefix,options.left_position,options.right_position)
	print("Finish the step of write_bam_file_density",file=sys.stderr)

def main():
		'''main funciton'''
		parse_args_for_specific_region_metagene()

if __name__=="__main__":
		main()
