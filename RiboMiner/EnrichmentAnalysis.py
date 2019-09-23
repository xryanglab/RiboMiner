#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
@Author: Li Fajin
@Date: 2019-08-27 21:13:08
@LastEditors: Li Fajin
@LastEditTime: 2019-08-30 16:59:46
@Description: This script is used for calculating enrichment ratio and outputing files used for plot.

'''
from .FunctionDefinition import *
from functools import reduce
from operator import mul
from scipy import stats
from collections import defaultdict


def get_density_dict(densityFile):
	density_dict={}
	with open(densityFile,'r') as f:
		i=0
		for line in f:
			i+=1
			if i==1:
				continue
			tmp=line.strip().split("\t")
			trans_id=tmp[0]
			density_dict[trans_id]=[float(i) for i in tmp[1:]]
	return density_dict

def get_all_counts(density_dict):
	all_counts=0
	for trans,counts in density_dict.items():
		all_counts+=sum(counts)
	return all_counts
def CI_for_t_distribution(data,confidence=0.95):
	'''
	scripts source:
	https://www.jianshu.com/p/6cfce4cc2f7f
	'''
	mean=np.mean(data)
	std=np.std(data,ddof=1)
	size=len(data)
	alpha=1-confidence
	t_score=stats.t.isf(alpha/2,df=(size-1))
	ME=t_score*std/np.sqrt(size)
	lower_CI=mean-ME
	upper_CI=mean+ME
	return lower_CI,upper_CI

def enrichment_ratio(ctrl_dict,treat_dict,in_selectTrans,inCDS_lengthFilterParma,inCDS_countsFilterParma,in_regionLengthParma,in_extendRegionLengthParma,in_excludeLengthParma,in_excludeCodonCountsParma,mode,unit,confidence):
	ctrl_all_counts=get_all_counts(ctrl_dict)
	treat_all_counts=get_all_counts(treat_dict)
	passTransSet=set()
	ratio_dict={}
	startratio=np.zeros(int(in_regionLengthParma+in_extendRegionLengthParma+1),dtype="float64")
	stopratio =np.zeros(int(in_regionLengthParma+in_extendRegionLengthParma+1),dtype="float64")
	start_lower_CI=np.zeros(int(in_regionLengthParma+in_extendRegionLengthParma+1),dtype="float64")
	start_upper_CI=np.zeros(int(in_regionLengthParma+in_extendRegionLengthParma+1),dtype="float64")
	stop_lower_CI=np.zeros(int(in_regionLengthParma+in_extendRegionLengthParma+1),dtype="float64")
	stop_upper_CI=np.zeros(int(in_regionLengthParma+in_extendRegionLengthParma+1),dtype="float64")
	startNormedWindowsList=[]
	startWindowsExistPosList=[]
	stopNormedWindowsList=[]
	stopWindowsExistPosList=[]
	filter_1=0
	filter_2=0
	filter_3=0
	filter_4=0
	filter_5=0
	for trans in in_selectTrans:
		ctrl_trans_counts=np.array(ctrl_dict[trans])
		treat_trans_counts=np.array(treat_dict[trans])
		cds_length=len(ctrl_trans_counts)
		if unit == 'codon':
			if cds_length < inCDS_lengthFilterParma:
				filter_1+=1
				continue
			if (cds_length*3) % 3 != 0:
				filter_2+=1
				continue
		elif unit == 'nt':
			if cds_length < inCDS_lengthFilterParma:
				filter_1+=1
				continue
			if cds_length % 3 != 0:
				filter_2+=1
				continue
		else:
			raise KeyError("There is no such unit, please check it.['codon' or 'nt']")
		ctrl_trans_counts_rpm=10**6*(ctrl_trans_counts/ctrl_all_counts)
		treat_trans_counts_rpm=10**6*(treat_trans_counts/treat_all_counts)
		ctrl_trans_counts_sum=sum(ctrl_trans_counts)
		treat_trans_counts_sum=sum(treat_trans_counts)
		ctrl_trans_counts_sum_norm=10**9*(ctrl_trans_counts_sum/(ctrl_all_counts*cds_length))
		treat_trans_counts_sum_norm=10**9*(treat_trans_counts_sum/(treat_all_counts*cds_length))
		if ctrl_trans_counts_sum == 0 or treat_trans_counts_sum ==0:
			filter_3+=1
			continue
		if mode == 'RPKM':
			if ctrl_trans_counts_sum_norm < inCDS_countsFilterParma or treat_trans_counts_sum_norm < inCDS_countsFilterParma:
				filter_4+=1
				continue
		elif mode == 'counts':
			if ctrl_trans_counts_sum < inCDS_countsFilterParma or treat_trans_counts_sum < inCDS_countsFilterParma:
				filter_4+=1
				continue
		else:
			raise KeyError("There is no such modes, please check it.['RPKM' or 'counts']")
		ctrl_normValue=np.mean(ctrl_trans_counts_rpm[in_excludeLengthParma:])
		treat_normValue=np.mean(treat_trans_counts_rpm[in_excludeLengthParma:])
		ctrl_sumValue=sum(ctrl_trans_counts[in_excludeLengthParma:])
		treat_sumValue=sum(treat_trans_counts[in_excludeLengthParma:])
		ctrl_sumValue_norm=10**9*(ctrl_sumValue/(ctrl_all_counts*cds_length))
		treat_sumValue_norm=10**9*(treat_sumValue/(treat_all_counts*cds_length))
		if mode == 'RPKM':
			if ctrl_normValue==0 or treat_normValue==0 or ctrl_sumValue_norm <in_excludeCodonCountsParma or treat_sumValue_norm <in_excludeCodonCountsParma:
				filter_5+=1
				continue
		elif mode == 'counts':
			if ctrl_normValue==0 or treat_normValue==0 or ctrl_sumValue <in_excludeCodonCountsParma or treat_sumValue <in_excludeCodonCountsParma:
				filter_5+=1
				continue
		else:
			raise KeyError("There is no such modes, please check it.['RPKM' or 'counts']")
		ratio=np.zeros(cds_length,dtype="float64")
		for i in np.arange(len(ratio)):
			ctrl=ctrl_trans_counts_rpm[i]+1 ## ctrl=ctrl+1
			treat=treat_trans_counts_rpm[i]+1 ## treat=treat+1
			# if ctrl ==0 or treat ==0:
			# 	ctrl+=1
			# 	treat+=1
			ratio[i]+=treat/ctrl
		####get trans density vector
		(tmpStartWin,tmpStartPos)=getWindowsVector(in_extendRegionLengthParma,in_regionLengthParma,ratio,0) #start codon coor is 0 (0-based)
		(tmpStopWin, tmpStopPos) =getWindowsVector(in_regionLengthParma,in_extendRegionLengthParma,ratio,(len(ratio)-1)) #stop codon coor is len-1 (0-based)
		startNormedWindowsList.append(tmpStartWin)
		startWindowsExistPosList.append(tmpStartPos)
		stopNormedWindowsList.append(tmpStopWin)
		stopWindowsExistPosList.append(tmpStopPos)
		passTransSet.add(trans)
		ratio_dict[trans]=ratio

	startNormedWindowsList=np.array(startNormedWindowsList)
	startWindowsExistPosList=np.array(startWindowsExistPosList)
	stopNormedWindowsList=np.array(stopNormedWindowsList)
	stopWindowsExistPosList=np.array(stopWindowsExistPosList)
	for terms in range(in_extendRegionLengthParma+in_regionLengthParma+1):
		startratio[terms]=np.mean(startNormedWindowsList[np.where(startWindowsExistPosList[:,terms]==1),terms])
		stopratio[terms] =np.mean(stopNormedWindowsList[np.where(stopWindowsExistPosList[:,terms]==1),terms])
		if len(in_selectTrans) > 1:
			start_lower_CI[terms],start_upper_CI[terms]=CI_for_t_distribution(startNormedWindowsList[np.where(startWindowsExistPosList[:,terms]==1),terms][0],confidence=confidence)
			stop_lower_CI[terms],stop_upper_CI[terms]=CI_for_t_distribution(stopNormedWindowsList[np.where(stopWindowsExistPosList[:,terms]==1),terms][0],confidence=confidence)
		else:
			pass
	print("filter_1: Lenght filter(-l)---Transcripts number filtered by criterion one is : "+str(filter_1),file=sys.stderr)
	print("filter_2: Lenght filter (3n)---Transcripts number filtered by criterion two is : "+str(filter_2),file=sys.stderr)
	print("filter_3: Trans counts filter (total counts=0)---Transcripts number filtered by criterion three is : "+str(filter_3),file=sys.stderr)
	print("filter_4: CDS counts filter(RPKM-n)---Transcripts number filtered by criterion four is : "+str(filter_4),file=sys.stderr)
	print("filter_5: RPKM(-m) ---Transcripts number filtered by criterion five is : "+str(filter_5),file=sys.stderr)
	print("The number of genes used for following analysis is: " + str(len(passTransSet)),file=sys.stderr)
	if len(in_selectTrans) > 1:
		return (startratio,stopratio,passTransSet,ratio_dict,start_lower_CI,start_upper_CI,stop_lower_CI,stop_upper_CI)
	else:
		return (startratio,stopratio,passTransSet,ratio_dict)


def write_enrichment_transcripts(passTransSet,output_prefix):
	data=[]
	data_index=[]
	trans=passTransSet
	sample=output_prefix
	data.append(trans)
	data_index.append(sample)
	data=pd.DataFrame(data,index=data_index)
	data=data.T
	data.to_csv(output_prefix+"_enrichment_transcript_id.txt",sep="\t",index=0)

def write_mean_density_dataframe(startMeanDensityDict,stopMeanDensityDict,startLowerCI,startUpperCI,stopLowerCI,stopUpperCI,outFile):
	data=[]
	for sample in startMeanDensityDict:
		k=pd.DataFrame([sample]*len(startMeanDensityDict[sample]))
		startDensity=pd.DataFrame(startMeanDensityDict[sample])
		stopDensity=pd.DataFrame(stopMeanDensityDict[sample])
		startLCI=pd.DataFrame(startLowerCI[sample])
		startUCI=pd.DataFrame(startUpperCI[sample])
		stopLCI=pd.DataFrame(stopLowerCI[sample])
		stopUCI=pd.DataFrame(stopUpperCI[sample])
		density=pd.merge(startDensity,stopDensity,how="left",left_index=True,right_index=True)
		startCI=pd.merge(startLCI,startUCI,how="left",left_index=True,right_index=True)
		stopCI=pd.merge(stopLCI,stopUCI,how="left",left_index=True,right_index=True)
		CI=pd.merge(startCI,stopCI,how="left",left_index=True,right_index=True)
		density=pd.merge(density,CI,how="left",left_index=True,right_index=True)
		density=pd.merge(k,density,how="left",left_index=True,right_index=True)
		data.append(density)
	temp=data[0]
	if len(data) < 1:
		raise EOFError("Empty file, there is nothing in the file.")
	if len(data) == 1:
		temp.columns=['sample','start_density','stop_density','start_lower_CI','start_upper_CI','stop_lower_CI','stop_upper_CI']
		temp.to_csv(outFile,sep="\t",index=0)
	else:
		for i in np.arange(1,len(data)):
			temp=np.vstack((temp,data[i]))
		temp=pd.DataFrame(temp,columns=["sample","start_density","stop_density",'start_lower_CI','start_upper_CI','stop_lower_CI','stop_upper_CI'])
		temp.to_csv(outFile,sep="\t",index=0)
def write_enrichment_dataframe(startratio,stopratio,output_prefix):
		data=[]
		k=pd.DataFrame([output_prefix]*len(startratio))
		start=pd.DataFrame(startratio)
		stop=pd.DataFrame(stopratio)
		density=pd.merge(start,stop,how="left",left_index=True,right_index=True)
		density=pd.merge(k,density,how="left",left_index=True,right_index=True)
		data.append(density)
		temp=data[0]
		if len(data) < 1:
				raise EOFError("Empty file, there is nothing in the file.")
		if len(data) == 1:
				temp.columns=['sample','start_density','stop_density']
				temp.to_csv(output_prefix+"_enrichment_dataframe.txt",sep="\t",index=0)
		else:
				for i in np.arange(1,len(data)):
						temp=np.vstack((temp,data[i]))
				temp=pd.DataFrame(temp,columns=["sample","start_density","stop_density"])
				temp.to_csv(output_prefix+"_enrichment_dataframe.txt",sep="\t",index=0)

def write_ratio_dict(ratio_dict,output_prefix):
	with open(output_prefix+"_codon_ratio.txt",'w') as f:
		f.write("%s\t%s\n" %("transcript",output_prefix))
		for trans,ratio in ratio_dict.items():
			f.write("%s\t" % trans)
			for i in range(len(ratio)):
				f.write("%s\t" % str(ratio[i]))
			f.write("\n")


def main():
	parser=create_parser_for_enrichment_analysis()
	(options,args)=parser.parse_args()
	(ctrlDensity,treatDensity,output_prefix, unit,mode,min_cds_codon, min_cds_counts,
	min_norm_region_counts, upstream_codon, downstream_codon, norm_exclude_codon,confidence) = (options.ctrlDensity,options.treatDensity,options.output_prefix,
	options.unit,options.mode,options.min_cds_codon,options.min_cds_counts,options.min_norm_region_counts,options.upstream_codon,
	options.downstream_codon,options.norm_exclude_codon,options.confidence)
	if not ctrlDensity or not treatDensity:
		raise IOError("Please set your --ctrl and --treat parameters!")
	ctrl_dict=get_density_dict(ctrlDensity)
	treat_dict=get_density_dict(treatDensity)
	selectTrans=set(ctrl_dict.keys()).intersection(set(treat_dict))
	## trans id transformation
	transID2geneID,transID2geneName=reload_transcripts_information(options.coorFile)[4:6]
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
	if len(select_trans) > 1:
		startMeanDensityDict={}
		stopMeanDensityDict={}
		startLowerCI={}
		startUpperCI={}
		stopLowerCI={}
		stopUpperCI={}
		(startratio,stopratio,passTransSet,ratio_dict,start_lower_CI,start_upper_CI,stop_lower_CI,stop_upper_CI) = enrichment_ratio(ctrl_dict,treat_dict,select_trans,
		min_cds_codon,min_cds_counts,downstream_codon,upstream_codon,norm_exclude_codon,min_norm_region_counts,mode,unit,confidence)
		startMeanDensityDict[output_prefix]=startratio
		stopMeanDensityDict[output_prefix]=stopratio
		startLowerCI[output_prefix]=start_lower_CI
		startUpperCI[output_prefix]=start_upper_CI
		stopLowerCI[output_prefix]=stop_lower_CI
		stopUpperCI[output_prefix]=stop_upper_CI
		print("Finish the step of enrichment ratio calculation",file=sys.stderr)
		## write density
		write_mean_density_dataframe(startMeanDensityDict,stopMeanDensityDict,startLowerCI,startUpperCI,stopLowerCI,stopUpperCI,output_prefix+"_enrichment_dataframe.txt")
	else:
		(startratio,stopratio,passTransSet,ratio_dict) =enrichment_ratio(ctrl_dict,treat_dict,select_trans,
		min_cds_codon,min_cds_counts,downstream_codon,upstream_codon,norm_exclude_codon,min_norm_region_counts,mode,unit,confidence)
		print("Finish the step of enrichment ratio calculation",file=sys.stderr)
		## write density
		write_enrichment_dataframe(startratio,stopratio,output_prefix)
	write_enrichment_transcripts(passTransSet,output_prefix)
	write_ratio_dict(ratio_dict,output_prefix)
	print("Finish the step of enrichment analysis!",file=sys.stderr)

if __name__=="__main__":
	main()