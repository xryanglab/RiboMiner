#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
@Author: Li Fajin
@Date: 2019-08-29 09:25:05
@LastEditors: Li Fajin
@LastEditTime: 2019-09-01 19:01:08
@Description: This script is used for plot of enrichment ratio for a single transcript
usage: python EnrichmentAnalysisForSingleTrans.py -i <all_ratio.txt> -o <output_prefix> -c <coorFile> -s transcript_id --id-type transcript_id --unit codon [-S| --ymin|--ymax...]
'''

from .FunctionDefinition import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def plot_ratio_for_single_trans(ratio,trans,inOutPrefix,ymin,ymax,unit,axvline,axhline,color='b',text_font={"size":20,"family":"Arial","weight":"bold"},legend_font={"size":20,"family":"Arial","weight":"bold"}):
	'''plot the ratio dsitribution'''
	plt.rc("font",weight="bold")
	winLen=len(ratio)
	with PdfPages(inOutPrefix + "_enrichment_ratio.pdf") as pdf:
		fig=plt.figure(figsize=(16,8))
		ax=fig.add_subplot(111)
		plt.plot(np.arange(0,winLen),ratio,color=color,label=trans,linewidth=1)
		if axvline:
			ax.axvline(axvline,color="r",dashes=[1,2],clip_on=False,linewidth=2)
		else:
			pass

		ax.set_xlabel("Distance from start codon" + "("+unit+")",fontdict=text_font)
		ax.set_ylabel("Mean enrichment (A.U)",fontdict=text_font)
		if axhline:
			ax.axhline(1,color="r",dashes=[2,3],clip_on=False,linewidth=2)
		else:
			pass
		ax.spines["top"].set_visible(False)
		ax.spines["right"].set_visible(False)
		ax.spines["bottom"].set_linewidth(2)
		ax.spines["left"].set_linewidth(2)
		ax.tick_params(which="both",width=2,labelsize=20)
		if not ymin and not ymax:
			pass
		elif not ymin and ymax:
			ax.set_ylim(0,ymax)
		elif ymin and not ymax:
			raise IOError("Please offer the ymax parameter as well!")
		elif ymin and ymax:
			ax.set_ylim(ymin,ymax)
		else:
			raise IOError("Please enter correct ymin and ymax parameters!")
		plt.legend(loc="best",prop=legend_font)
		plt.tight_layout()
		pdf.savefig(fig)
		plt.close()

def plot_ratio_for_all_trans(ratio_dict,transList,inOutPrefix,ymin,ymax,unit,axvline,axhline,start,window,step,slideWindow,color='b',text_font={"size":20,"family":"Arial","weight":"bold"},legend_font={"size":20,"family":"Arial","weight":"bold"}):
	'''plot the ratio dsitribution'''
	plt.rc("font",weight="bold")
	with PdfPages(inOutPrefix + "_enrichment_ratio.pdf") as pdf:
		for trans in transList:
			ratio=ratio_dict[trans]
			if slideWindow:
				ratio=slide_window_average(ratio,start,window,step)
			else:
				pass
			winLen=len(ratio)
			fig=plt.figure(figsize=(16,8))
			ax=fig.add_subplot(111)
			plt.plot(np.arange(0,winLen),ratio,color=color,label=trans,linewidth=1)
			if axvline:
				ax.axvline(axvline,color="r",dashes=[1,2],clip_on=False,linewidth=2)
			else:
				pass

			ax.set_xlabel("Distance from start codon" + "("+unit+")",fontdict=text_font)
			ax.set_ylabel("Mean enrichment (A.U)",fontdict=text_font)
			if axhline:
				ax.axhline(1,color="r",dashes=[2,3],clip_on=False,linewidth=2)
			else:
				pass
			ax.spines["top"].set_visible(False)
			ax.spines["right"].set_visible(False)
			ax.spines["bottom"].set_linewidth(2)
			ax.spines["left"].set_linewidth(2)
			ax.tick_params(which="both",width=2,labelsize=20)
			if not ymin and not ymax:
				pass
			elif not ymin and ymax:
				ax.set_ylim(0,ymax)
			elif ymin and not ymax:
				raise IOError("Please offer the ymax parameter as well!")
			elif ymin and ymax:
				ax.set_ylim(ymin,ymax)
			else:
				raise IOError("Please enter correct ymin and ymax parameters!")
			plt.legend(loc="best",prop=legend_font)
			plt.tight_layout()
			pdf.savefig(fig)
			plt.close()

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


def ID_transformation(transcript,coorFile,Type='singleTrans',id_type="transcript_id"):
	transID2geneID,transID2geneName=reload_transcripts_information(coorFile)[4:6]
	geneID2transID={v:k for k,v in transID2geneID.items()}
	geneName2transID={v:k for k,v in transID2geneName.items()}
	if Type=='transList':
		if id_type == 'transcript_id':
			select_trans=set(transcript)
		elif id_type == 'gene_id':
			tmp=[geneID2transID[gene_id] for gene_id in transcript]
			select_trans=set(tmp)
			print("There are " + str(len(select_trans))+" gene id could be transformed into transcript id and used for following analysis.",file=sys.stderr)
		elif id_type == 'gene_name' or id_type=='gene_symbol':
			tmp=[geneName2transID[gene_name] for gene_name in transcript]
			select_trans=set(tmp)
			print("There are " + str(len(select_trans))+" gene symbol could be transformed into transcript id and used for following analysis.",file=sys.stderr)
		else:
			raise IOError("Please input a approproate id_type parameters.[transcript_id/gene_id/gene_name/]")
	elif Type == 'singleTrans':
		if id_type == 'transcript_id':
			select_trans=transcript
		elif id_type == 'gene_id':
			select_trans=geneID2transID[transcript]
		elif id_type == 'gene_name' or id_type=='gene_symbol':
			select_trans=geneName2transID[transcript]
		else:
			raise IOError("Please input a approproate id_type parameters.[transcript_id/gene_id/gene_name/]")
	else:
		pass
	return select_trans

def slide_window_average(ratio,start,window,step):
	''' Used for calculating mean density with a slide window'''

	winLen=len(ratio)
	tmp_data=np.zeros(winLen)
	tmp_data[0:int(start)]+=ratio[0:int(start)]
	tmp_data[-int(start):]+=ratio[-int(start):]
	for j in np.arange(start,winLen-start,step):
		tmp_data[j]+=np.mean(ratio[(j-int((window-1)/2)):(j+int((window-1)/2))])
	return tmp_data

def main():
	parser=create_parser_for_single_ratio_plot()
	(options,args)=parser.parse_args()
	if not options.ratioFile or not options.output_prefix or not options.coorFile:
		raise IOError("Please reset your parameters!")
	if options.window%2 == 0:
		raise IOError("Please reset your --window parameter. It must be a odd number.")
	if (options.start_position-1) < (options.window-1)/2:
		raise IOError("Please reset your --step and --window parameters. The (window-1)/2 must be less than start-1")
	ratio_dict=get_density_dict(options.ratioFile)
	text_font={"size":40,"family":"Arial","weight":"bold"}
	legend_font={"size":30,"family":"Arial","weight":"bold"}
	if options.singleTrans and not options.in_selectTrans:
		select_trans=ID_transformation(options.singleTrans,options.coorFile,Type='singleTrans',id_type=options.id_type)
		if select_trans not in ratio_dict.keys():
			raise IOError("Please reset your -s parameter! "+options.singleTrans+" not in " + options.ratioFile+"!")
		ratio=ratio_dict[select_trans]
		if not options.slideWindow:
			pass
		else:
			ratio=slide_window_average(ratio,options.start_position,options.window,options.step)
		plot_ratio_for_single_trans(ratio,select_trans,options.output_prefix,options.ymin,options.ymax,options.unit,options.axvline,options.axhline,color='b',text_font=text_font,legend_font=legend_font)
		print("Finish the step of ratio plot!",file=sys.stderr)
	elif options.in_selectTrans and not options.singleTrans:
		transList=pd.read_csv(options.in_selectTrans,sep="\t")
		select_trans=ID_transformation(transList,options.coorFile,Type='transList',id_type=options.id_type)
		select_trans=select_trans.intersection(set(ratio_dict.keys()))
		plot_ratio_for_all_trans(ratio_dict,select_trans,options.output_prefix,options.ymin,options.ymax,options.unit,options.axvline,options.axhline,options.start_position,options.window,options.step,options.slideWindow,color='b',text_font=text_font,legend_font=legend_font)
		print("Finish the step of ratio plot!",file=sys.stderr)
	elif options.in_selectTrans and options.singleTrans:
		raise IOError("The -s and -S are mutually exclusive!")
	else:
		plot_ratio_for_all_trans(ratio_dict,ratio_dict.keys(),options.output_prefix,options.ymin,options.ymax,options.unit,options.axvline,options.axhline,options.start_position,options.window,options.step,options.slideWindow,color='b',text_font=text_font,legend_font=legend_font)
		print("Finish the step of ratio plot!",file=sys.stderr)

if __name__=="__main__":
	main()