#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
@Author: Li Fajin
@Date: 2019-08-16 15:09:23
@LastEditors: Li Fajin
@LastEditTime: 2019-08-30 17:08:36
@Description: The script is used for calculating ribosome density at each position along transcript and as well as "coverage" of each transcript.
'''


from .FunctionDefinition import *



def ribosomeDensityAtEachPosition(in_bamFile,in_bamLegend,in_selectTrans,in_transLengthDict,in_startCodonCoorDict,in_stopCodonCoorDict,in_readLengths,in_readOffset,Unit,output_prefix):
	pysamFile=pysam.AlignmentFile(in_bamFile,"rb")
	pysamFile_trans=pysamFile.references
	in_selectTrans=set(pysamFile_trans).intersection(in_selectTrans)
	## codon unit
	if Unit=="codon":
		with open(output_prefix+"_"+in_bamLegend.strip()+"_cds_codon_density.txt",'w') as f1,open(output_prefix+"_"+in_bamLegend.strip()+"_cds_codon_coverage.txt",'w') as f2:
			f1.write("%s\t%s\n" %("transcript",in_bamLegend))
			f2.write("%s\t%s\n" %("transcript",in_bamLegend))
			for trans in in_selectTrans:
				leftCoor =int(in_startCodonCoorDict[trans])-1 #the first base of start codon 0-base
				rightCoor=int(in_stopCodonCoorDict[trans])-3 #the first base of stop codon 0-base
				(trans_counts,read_counts_frameSum,total_reads,cds_reads)=get_trans_frame_counts(pysamFile, trans, in_readLengths, in_readOffset, in_transLengthDict[trans], leftCoor, rightCoor)
				cds_codon_density=np.array(read_counts_frameSum)
				cds_coverage=np.sum(cds_codon_density!=0)/len(read_counts_frameSum)
				f2.write("%s\t%s\n" %(trans,str(cds_coverage)))
				f1.write("%s\t" %(trans))
				for i in range(len(read_counts_frameSum)):
					# print(trans,str(i+1),str(read_counts_frameSum[i]))
					# f1.write("%s\t%s\t%s\n" %(trans,str(i+1),str(read_counts_frameSum[i])))
					f1.write("%s\t" %(str(read_counts_frameSum[i])))
				f1.write("\n")
	elif Unit=="nt":
		with open(output_prefix+"_"+in_bamLegend.strip()+"_cds_nt_density.txt",'w') as f1, open(output_prefix+"_"+in_bamLegend.strip()+"_cds_nt_coverage.txt",'w') as f2:
			f1.write("%s\t%s\n" %("transcript",in_bamLegend))
			f2.write("%s\t%s\n" %("transcript",in_bamLegend))
			for trans in in_selectTrans:
				leftCoor =int(in_startCodonCoorDict[trans])-1 #the first base of start codon 0-base
				rightCoor=int(in_stopCodonCoorDict[trans])-3 #the first base of stop codon 0-base
				(trans_counts,read_counts_frameSum,total_reads,cds_reads)=get_trans_frame_counts(pysamFile, trans, in_readLengths, in_readOffset, in_transLengthDict[trans], leftCoor, rightCoor)
				trans_counts_cds=np.array(trans_counts[leftCoor:rightCoor+3])
				cds_coverage=sum(trans_counts_cds!=0)/len(trans_counts_cds)
				f2.write("%s\t%s\n" %(trans,str(cds_coverage)))
				f1.write("%s\t" %(trans))
				for i in range(len(trans_counts_cds)):
					# print(trans,str(i+1),str(trans_counts[i]))
					# f1.write("%s\t%s\t%s\n" %(trans,str(i+1),str(trans_counts_cds[i])))
					f1.write("%s\t" %(str(trans_counts_cds[i])))
				f1.write("\n")
	else:
		raise IOError("Please enter a proper -U parameter![codon or nt]")

def parse_args_for_riboDensity_atEachPosition():
	parsed=create_parser_for_riboDensity_atEachPosition()
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
	select_trans=sorted(list(select_trans))
	for bamfs in bam_attr:
		ribosomeDensityAtEachPosition(bamfs.bamName,bamfs.bamLegend,select_trans,transLengthDict,startCodonCoorDict,stopCodonCoorDict,bamfs.bamLen,bamfs.bamOffset,options.unit,options.output_prefix)
	print("Finish the step of ribosomeDensityAtEachPosition",file=sys.stderr)

def main():
	"""main program"""
	## parse the options and args and output
	parse_args_for_riboDensity_atEachPosition()


if __name__ == "__main__":
	main()


