#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
@Description:
	Perform metagene plot around start and stop codon in codon units or nt units
	Note:
	1) Only input sorted and indexed BAM file(s). SAM format is not supported.
	2) Acoording to one gene may has serveral isoforms , we select the longest isoform.
	3) The selected isoform whose CDS length is not multiple to 3 duing to programmed ribosome frameshifting were excluded also.
'''

from __future__ import division
from .FunctionDefinition import *
import itertools
from scipy import stats
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

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

def ribosomeDensityNormPerTrans(in_bamFile,in_selectTrans,in_transLengthDict,in_startCodonCoorDict,in_stopCodonCoorDict,inCDS_lengthFilterParma,inCDS_countsFilterParma,in_regionLengthParma,in_extendRegionLengthParma,in_excludeLengthParma,in_excludeCodonCountsParma,in_readLengths,in_readOffset,Mode,Unit,Type,confidence,norm):
	pysamFile=pysam.AlignmentFile(in_bamFile,"rb")
	pysamFile_trans=pysamFile.references
	in_selectTrans=set(pysamFile_trans).intersection(in_selectTrans)
	passTransSet=set()
	startDensity=np.zeros(int(in_regionLengthParma+in_extendRegionLengthParma+1),dtype="float64")
	stopDensity =np.zeros(int(in_regionLengthParma+in_extendRegionLengthParma+1),dtype="float64")
	start_lower_CI=np.zeros(int(in_regionLengthParma+in_extendRegionLengthParma+1),dtype="float64")
	start_upper_CI=np.zeros(int(in_regionLengthParma+in_extendRegionLengthParma+1),dtype="float64")
	stop_lower_CI=np.zeros(int(in_regionLengthParma+in_extendRegionLengthParma+1),dtype="float64")
	stop_upper_CI=np.zeros(int(in_regionLengthParma+in_extendRegionLengthParma+1),dtype="float64")
	startNormedWindowsList=[]
	startWindowsExistPosList=[]
	stopNormedWindowsList=[]
	stopWindowsExistPosList=[]
	i=len(in_selectTrans)
	filter_1=0
	filter_2=0
	filter_3=0
	filter_4=0
	filter_5=0
	all_counts=0
	startNormedWindowsDict={}
	for trans in in_startCodonCoorDict.keys():
		leftCoor =int(in_startCodonCoorDict[trans])-1
		rightCoor=int(in_stopCodonCoorDict[trans])-3
		(trans_counts,read_counts_frameSum,total_reads,cds_reads)=get_trans_frame_counts(pysamFile, trans, in_readLengths, in_readOffset, in_transLengthDict[trans], leftCoor, rightCoor)
		all_counts+=total_reads ## total_reads for transcript level
	## codon unit
	if Unit=="codon":
		for trans in in_selectTrans:
			i-=1
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
			read_counts_frameSum_normed=10**6*(read_counts_frameSum/(all_counts))
			# read_counts_frameSum_normed=10**9*(read_counts_frameSum/(all_counts*len(read_counts_frameSum)))
			trans_counts_normed=10**6*(trans_counts/(all_counts))
			# trans_counts_normed=10**9*(trans_counts/(all_counts*len(trans_counts)))
			sumValue=np.sum(read_counts_frameSum[int(in_excludeLengthParma):])
			sumValue_normed=10**9*(sumValue/(all_counts*len(read_counts_frameSum)))
			NormValue=np.mean(read_counts_frameSum_normed[int(in_excludeLengthParma):])
			if Mode == 'RPKM':
				if NormValue == 0 or sumValue_normed < in_excludeCodonCountsParma:
					filter_5+=1
					continue
			elif Mode == 'counts':
				if NormValue == 0 or sumValue < in_excludeCodonCountsParma:
					filter_5+=1
					continue
			else:
				raise KeyError("There is no such modes, please check it.['RPKM' or 'counts']")
			#get trans density vector
			if Type.upper() == 'CDS':
				(tmpStartWin,tmpStartPos)=getWindowsVector(in_extendRegionLengthParma,in_regionLengthParma,read_counts_frameSum_normed,0) #start codon coor is 0 (0-based)
				(tmpStopWin, tmpStopPos) =getWindowsVector(in_regionLengthParma,in_extendRegionLengthParma,read_counts_frameSum_normed,(len(read_counts_frameSum_normed)-1)) #stop codon coor is len-1 (0-based)
			elif Type.upper() == 'UTR':
				# trans_counts_codon_level=[np.sum(trans_counts[i:i+3]) for i in range(0,len(trans_counts),3)]
				trans_counts_codon_level_normed=[np.sum(trans_counts_normed[i:i+3]) for i in range(0,len(trans_counts_normed),3)]
				(tmpStartWin,tmpStartPos)=getWindowsVector(in_extendRegionLengthParma,in_regionLengthParma,trans_counts_codon_level_normed,int(leftCoor/3)) #start codon coor is 0 (0-based)
				(tmpStopWin, tmpStopPos) =getWindowsVector(in_regionLengthParma,in_extendRegionLengthParma,trans_counts_codon_level_normed,int(rightCoor/3)) #stop codon coor is len-1 (0-based)
			else:
				raise IOError("Please input an correct --type parameters!")
			if norm in ['yes','y','Y','YES','Yes']:
				startNormedWindowsList.append(tmpStartWin/NormValue)
				startWindowsExistPosList.append(tmpStartPos)
				stopNormedWindowsList.append(tmpStopWin/NormValue)
				stopWindowsExistPosList.append(tmpStopPos)
				startNormedWindowsDict[trans]=tmpStartWin/NormValue
				passTransSet.add(trans)
			elif norm in ['no','n','N','NO','No','not','Not','NOT']:
				startNormedWindowsList.append(tmpStartWin)
				startWindowsExistPosList.append(tmpStartPos)
				stopNormedWindowsList.append(tmpStopWin)
				stopWindowsExistPosList.append(tmpStopPos)
				passTransSet.add(trans)
				startNormedWindowsDict[trans]=tmpStartWin
			else:
				raise IOError("Please enter a proper --norm parameter [yes or not]")
	## nt unit
	if Unit=="nt":
		for trans in in_selectTrans:
			i-=1
			leftCoor =int(in_startCodonCoorDict[trans])-1 #the first base of start codon 0-base
			rightCoor=int(in_stopCodonCoorDict[trans])-3 #the first base of stop codon 0-base
			# cdsLength should obey : (1) >= minFilterLength nt unit (2) CDSlength=n (n could be divided by 3)
			if ((rightCoor-leftCoor) < inCDS_lengthFilterParma):
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
			trans_counts_normed=10**6*(trans_counts/(all_counts))
			# trans_counts_normed=10**9*(trans_counts/(all_counts*len(trans_counts)))
			NormValue=np.mean(trans_counts_normed[leftCoor:rightCoor][int(in_excludeLengthParma):])
			sumValue=np.sum(trans_counts[leftCoor:rightCoor][int(in_excludeLengthParma):]) # for nt unit
			sumValue_normed=10**9*(sumValue/(all_counts*len(trans_counts[leftCoor:rightCoor])))
			if Mode == 'RPKM':
				if NormValue == 0 or sumValue_normed < in_excludeCodonCountsParma:
					filter_5+=1
					continue
			elif Mode == 'counts':
				if NormValue == 0 or sumValue < in_excludeCodonCountsParma:
					filter_5+=1
					continue
			else:
				raise KeyError("There is no such modes, please check it.['RPKM' or 'counts']")
			if Type.upper() == 'CDS':
				(tmpStartWin,tmpStartPos)=getWindowsVector(in_extendRegionLengthParma,in_regionLengthParma,trans_counts_normed[leftCoor:rightCoor],0) #for nt unit
				(tmpStopWin, tmpStopPos) =getWindowsVector(in_regionLengthParma,in_extendRegionLengthParma,trans_counts_normed[leftCoor:rightCoor],(len(trans_counts_normed[leftCoor:rightCoor])-1)) #for nt unit
			elif Type.upper() == 'UTR':
				(tmpStartWin,tmpStartPos)=getWindowsVector(in_extendRegionLengthParma,in_regionLengthParma,trans_counts_normed,leftCoor) #for nt unit
				(tmpStopWin, tmpStopPos) =getWindowsVector(in_regionLengthParma,in_extendRegionLengthParma,trans_counts_normed,rightCoor) #for nt unit
			else:
				raise IOError("Please input an correct --type parameters!")
			if norm in ['yes','y','Y','YES','Yes']:
				startNormedWindowsList.append(tmpStartWin/NormValue)
				startWindowsExistPosList.append(tmpStartPos)
				stopNormedWindowsList.append(tmpStopWin/NormValue)
				stopWindowsExistPosList.append(tmpStopPos)
				startNormedWindowsDict[trans]=tmpStartWin/NormValue
				passTransSet.add(trans)
			elif norm in ['no','n','N','NO','No','not','Not','NOT']:
				startNormedWindowsList.append(tmpStartWin)
				startWindowsExistPosList.append(tmpStartPos)
				stopNormedWindowsList.append(tmpStopWin)
				stopWindowsExistPosList.append(tmpStopPos)
				passTransSet.add(trans)
				startNormedWindowsDict[trans]=tmpStartWin
			else:
				raise IOError("Please enter a proper --norm parameter [yes or not]")
	###norm per position by mean or median
	startNormedWindowsList=np.array(startNormedWindowsList)
	startWindowsExistPosList=np.array(startWindowsExistPosList)
	stopNormedWindowsList=np.array(stopNormedWindowsList)
	stopWindowsExistPosList=np.array(stopWindowsExistPosList)
	#
	for terms in range(in_extendRegionLengthParma+in_regionLengthParma+1):
		startDensity[terms]=np.mean(startNormedWindowsList[np.where(startWindowsExistPosList[:,terms]==1),terms])
		stopDensity[terms] =np.mean(stopNormedWindowsList[np.where(stopWindowsExistPosList[:,terms]==1),terms])
		start_lower_CI[terms],start_upper_CI[terms]=CI_for_t_distribution(startNormedWindowsList[np.where(startWindowsExistPosList[:,terms]==1),terms][0],confidence=confidence)
		stop_lower_CI[terms],stop_upper_CI[terms]=CI_for_t_distribution(stopNormedWindowsList[np.where(stopWindowsExistPosList[:,terms]==1),terms][0],confidence=confidence)
	pysamFile.close()
	print("Lenght filter(-l)---Transcripts number filtered by criterion one is : "+str(filter_1),file=sys.stderr)
	print("Lenght filter (3n)---Transcripts number filtered by criterion two is : "+str(filter_2),file=sys.stderr)
	print("Total counts filter---Transcripts number filtered by criterion three is : "+str(filter_3),file=sys.stderr)
	print("CDS density filter(RPKM-n or counts-n)---Transcripts number filtered by criterion four is : "+str(filter_4),file=sys.stderr)
	print("CDS density filter(normed-m)---Transcripts number filtered by criterion five is : "+str(filter_5),file=sys.stderr)
	print("Metaplots Transcript Number for bam file"+in_bamFile+" is :"+str(len(passTransSet)),file=sys.stderr)
	return (startDensity,stopDensity,passTransSet,startNormedWindowsDict,start_lower_CI,start_upper_CI,stop_lower_CI,stop_upper_CI)

def metagenePlot(meanDensityDataFrame,in_regionLengthParma,in_extendRegionLengthParma,inOutPrefix,unit):
	plt.rc('font',weight='bold')
	samples=np.unique(meanDensityDataFrame.iloc[:,0].values)
	text_font={"size":30,"family":"Arial","weight":"bold"}
	legend_font={"size":30,"family":"Arial","weight":"bold"}
	winLen=in_regionLengthParma+in_extendRegionLengthParma+1
	if len(samples) <=8:
		colors=["b","orangered","green","c","m","y","k","w"]
	else:
		colors=sns.color_palette('husl',len(samples))
	## start codon
	with PdfPages(inOutPrefix + "_metagenePlot.pdf") as pdf:
		fig1=plt.figure(figsize=(16,8))
		ax1=fig1.add_subplot(111)
		for i in np.arange(len(samples)):
			plt.plot(np.arange(0,winLen),meanDensityDataFrame.iloc[np.where(meanDensityDataFrame.iloc[:,0]==samples[i])].iloc[:,1],color=colors[i],label=samples[i],linewidth=1)
			ax1.set_xticks(np.arange(0,winLen,50))
			ax1.set_xticklabels((np.arange(0,winLen,50)-in_extendRegionLengthParma))
			if unit == 'codon':
				ax1.set_xlabel("Distance from start codon"+' (codon)',fontdict=text_font)
			elif unit == 'nt':
				ax1.set_xlabel("Distance from start codon "+' (nt)',fontdict=text_font)
			else:
				raise KeyError("KeyError:"+str(unit))
			ax1.set_ylabel("Relative footprint density (AU)",fontdict=text_font)
			ax1.set_title("Ribosome footprint density profiles",fontdict=text_font)
			ax1.spines["top"].set_visible(False)
			ax1.spines["right"].set_visible(False)
			ax1.spines["bottom"].set_linewidth(2)
			ax1.spines["left"].set_linewidth(2)
			ax1.tick_params(which="both",width=2,labelsize=20)
			ax1.axhline(1,color="r",dashes=(2,3),clip_on=False,alpha=0.5)
			plt.legend(loc="best",prop=legend_font)
		## stop codon
		fig2=plt.figure(figsize=(16,8))
		ax2=fig2.add_subplot(111)
		for i in np.arange(len(samples)):
			##  stop codon
			plt.plot(np.arange(0,winLen),meanDensityDataFrame.iloc[np.where(meanDensityDataFrame.iloc[:,0]==samples[i])].iloc[:,2],color=colors[i],label=samples[i],linewidth=1)
			ax2.set_xticks(np.arange(0,winLen,50))
			ax2.set_xticklabels((np.arange(0,winLen,50)-in_regionLengthParma))
			if unit == 'codon':
				ax2.set_xlabel("Distance from stop codon"+' (codon)',fontdict=text_font)
			elif unit == 'nt':
				ax2.set_xlabel("Distance from stop codon "+' (nt)',fontdict=text_font)
			else:
				raise KeyError("KeyError:"+str(unit))
			ax2.set_ylabel("Relative footprint density (AU)",fontdict=text_font)
			ax2.set_title("Ribosome footprint density profiles",fontdict=text_font)
			ax2.spines["top"].set_visible(False)
			ax2.spines["right"].set_visible(False)
			ax2.spines["bottom"].set_linewidth(2)
			ax2.spines["left"].set_linewidth(2)
			ax2.tick_params(which="both",width=2,labelsize=20)
			ax2.axhline(1,color="r",dashes=(2,3),clip_on=False,alpha=0.5)
			plt.legend(loc="best",prop=legend_font)
		pdf.savefig(fig1)
		pdf.savefig(fig2)
		plt.close()


def write_passed_transcripts(inBamAttr,outFile):
	data=[]
	data_index=[]
	for bms in inBamAttr:
		trans=bms.passTransSet
		sample=bms.bamLegend
		data.append(trans)
		data_index.append(sample)
	data=pd.DataFrame(data,index=data_index)
	data=data.T
	data.to_csv(outFile,sep="\t",index=0)

def write_codon_units_density(inBamAttr,upstream,downstream,num,outFile,Type):
	try:
		for bms in inBamAttr:
			with open(outFile+"_"+bms.bamLegend+"_codon_density.txt","w") as fout:
				fout.write("%s\t%s\n" % (bms.bamLegend,"start_codon"))
				fout.write("%s" % ("transcript_id"))
				if Type.upper()=="CDS":
					for ss in range(int(num)):
						fout.write("\t%s" % ("codon_"+str(ss+1)))
				elif Type.upper()=="UTR":
					for ss in range(int(num)*2+1):
						fout.write("\t%s" %("codon_"+str(ss-int(num))))
				else:
					raise IOError("Please enter a correct Type [CDS or UTR]")

				fout.write("\n")
				for ts in bms.startNormedWindowsDict:
					fout.write("%s" % ( ts ) ) #the first num codons
					if Type.upper()=="CDS":
						for codons in range(int(num)):
							fout.write("\t%f" % ( bms.startNormedWindowsDict[ts][codons] ) )
					elif Type.upper()=="UTR":
						for nt in range(int(upstream)+int(downstream)+1):
							if nt < (int(upstream)-int(num)) or (nt > (int(upstream)+int(downstream)+1-int(num))):
								pass
							else:
								fout.write("\t%f" %(bms.startNormedWindowsDict[ts][nt]))
					fout.write("\n")
	except IOError:
			print("IOError! Please check your parameters!")

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
		return temp
	else:
		for i in np.arange(1,len(data)):
			temp=np.vstack((temp,data[i]))
		temp=pd.DataFrame(temp,columns=["sample","start_density","stop_density",'start_lower_CI','start_upper_CI','stop_lower_CI','stop_upper_CI'])
		temp.to_csv(outFile,sep="\t",index=0)
		return temp

def parse_args_for_CDS_metagene_analysis():
	parsed=create_parser_for_metagene_analysis()
	(options,args)=parsed.parse_args()
	if options.bamListFile and (options.bam_files or options.read_length or options.read_offset or options.bam_file_legend):
		raise IOError("'-f' parameter and '-i -r -s -t' are mutually exclusive.")
	if options.bamListFile:
		bamFiles,readLengths,Offsets,bamLegends=parse_bamListFile(options.bamListFile)
	elif options.bam_files:
		bamFiles,readLengths,Offsets,bamLegends=options.bam_files.split(","),options.read_length.split("_"),options.read_offset.split("_"),options.bam_file_legend.split(",")
	else:
		raise IOError("Please check you input files!")
	if options.norm_codon_density_num > options.min_cds_codon:
		print("Notice! The -y parameter should not be bigger than -l parameter!")
	if options.norm_codon_density_num > abs(int(options.downstream_codon)+int(options.upstream_codon)):
		raise IOError("Output codon density num bigger than the selective region length [-y > -u + -d]. Reset -y or -u and -d  parameters.")
	if options.type.upper()=="UTR":
		if options.norm_codon_density_num> options.downstream_codon or options.norm_codon_density_num>options.upstream_codon:
			raise IOError("Output codon density num bigger than the selective region length [-y > -u or -y > -d]. Reset -y or -u and -d  parameters.")
	else :
		if options.norm_codon_density_num>options.norm_codon_density_num> options.downstream_codon:
			raise IOError("Output codon density num bigger than the selective region length [-y > -d]. Reset -y or -u and -d  parameters.")

	print("your input : "+ str(len(bamFiles))+" bam files",file=sys.stderr)
	bam_attr=[]
	for ii,jj,mm,nn in zip(bamFiles,readLengths,Offsets,bamLegends):
		bam=bam_file_attr(ii,jj,mm,nn)
		bam_attr.append(bam)
	## calculate density for each bam files
	selectTrans,transLengthDict,startCodonCoorDict,stopCodonCoorDict,transID2geneID,transID2geneName,cdsLengthDict,transID2ChromDict=reload_transcripts_information(options.coorFile)
	geneID2transID={v:k for k,v in transID2geneID.items()}
	geneName2transID={v:k for k,v in transID2geneName.items()}
	if options.in_selectTrans:
		select_trans=pd.read_csv(options.in_selectTrans,sep="\t")
		select_trans=set(select_trans.iloc[:,0].values)
		if options.id_type == 'transcript_id':
			select_trans=select_trans.intersection(selectTrans)
			print("There are " + str(len(select_trans)) + " transcripts from "+options.in_selectTrans+" used for following analysis.",file=sys.stderr)
		elif options.id_type == 'gene_id':
			tmp=[geneID2transID[gene_id] for gene_id in select_trans if gene_id in geneID2transID]
			select_trans=set(tmp)
			select_trans=select_trans.intersection(selectTrans)
			print("There are " + str(len(select_trans))+" gene id could be transformed into transcript id and used for following analysis.",file=sys.stderr)
		elif options.id_type == 'gene_name' or options.id_type=='gene_symbol':
			tmp=[geneName2transID[gene_name] for gene_name in select_trans if gene_name in geneName2transID]
			select_trans=set(tmp)
			select_trans=select_trans.intersection(selectTrans)
			print("There are " + str(len(select_trans))+" gene symbol could be transformed into transcript id and used for following analysis.",file=sys.stderr)
		else:
			raise IOError("Please input a approproate id_type parameters.[transcript_id/gene_id/gene_name/]")
	else:
		select_trans=selectTrans
	startMeanDensityDict={}
	stopMeanDensityDict={}
	startLowerCI={}
	startUpperCI={}
	stopLowerCI={}
	stopUpperCI={}
	for bamfs in bam_attr:
		(bamfs.start_density,bamfs.stop_density,bamfs.passTransSet,bamfs.startNormedWindowsDict,bamfs.start_lower_CI,bamfs.start_upper_CI,bamfs.stop_lower_CI,bamfs.stop_upper_CI) = ribosomeDensityNormPerTrans(bamfs.bamName,select_trans,transLengthDict,
		startCodonCoorDict,stopCodonCoorDict,options.min_cds_codon,options.min_cds_counts,options.downstream_codon,options.upstream_codon,options.norm_exclude_codon,
		options.min_norm_region_counts,bamfs.bamLen,bamfs.bamOffset,options.mode,options.unit,options.type,options.confidence,options.norm)
		print("Finish the step of ribosomeDensityNormPerTrans",file=sys.stderr)
		startMeanDensityDict[bamfs.bamLegend]=bamfs.start_density
		stopMeanDensityDict[bamfs.bamLegend]=bamfs.stop_density
		startLowerCI[bamfs.bamLegend]=bamfs.start_lower_CI
		startUpperCI[bamfs.bamLegend]=bamfs.start_upper_CI
		stopLowerCI[bamfs.bamLegend]=bamfs.stop_lower_CI
		stopUpperCI[bamfs.bamLegend]=bamfs.stop_upper_CI
	## write density
	mean_density_dataframe=write_mean_density_dataframe(startMeanDensityDict,stopMeanDensityDict,startLowerCI,startUpperCI,stopLowerCI,stopUpperCI,options.output_prefix+"_dataframe.txt")
	write_passed_transcripts(bam_attr,options.output_prefix+"_transcript_id.txt")
	write_codon_units_density(bam_attr,options.upstream_codon,options.downstream_codon,options.norm_codon_density_num, options.output_prefix,options.type)
	print("Finish the step of MetageneAnalysis!",file=sys.stderr)
	if options.plot.upper() in ['YES','Y']:
		metagenePlot(mean_density_dataframe,options.downstream_codon,options.upstream_codon,options.output_prefix,options.unit)
		print("Finish the step of metagenePlot!",file=sys.stderr)
	elif options.plot.upper() in ['NO','N','NONE']:
		pass
	else:
		raise IOError("Please input a correct --plot parameter! [yes/no]")


def main():
	"""main program"""
	## parse the options and args and output
	parse_args_for_CDS_metagene_analysis()


if __name__ == "__main__":
		main()


