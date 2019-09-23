#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
@Author: Li Fajin
@Date: 2019-08-28 14:02:54
@LastEditors: Li Fajin
@LastEditTime: 2019-08-30 16:42:02
@Description: This script is used for plot of enrichment ratio.
1) the input file must be python DataFrame format, and has required three columns:
		1. columns one : sample name
		2. columns two : enrichment ratio values from the start codon
		3. columns three: enrichment ratio values from the stop codon
2) the output file could be pdf/png/jpg format
'''

import numpy as np
import pandas as pd
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from optparse import OptionParser
from functools import reduce
from itertools import chain
from collections import defaultdict
from .__init__ import __version__

def create_parser():
	'''argument parser'''
	usage="usage: python %prog [options]" +"\n" + __doc__+"\n"
	parser=OptionParser(usage=usage,version=__version__)
	parser.add_option("-i","--input",action="store",type="string",dest="density_file",help="Input file in txt format.And the files has three columns; column 1: sample;columns 2: start_density; column 3: stop_density")
	parser.add_option("-d","--downstream_codon",action="store",type="int",default=500, dest="downstream_codon", help="Downstream codon corresponding to start codon (codon unit). While corresponding to stop codon, it is the upstream codon.")
	parser.add_option("-u","--upstream_codon",action="store",type="int",default=0,dest="upstream_codon",help="Upstream codon corresponding to start codon (codon unit). While corresponding to stop codon, it is the downstream codon.")
	parser.add_option("-o","--output",action="store",type="string",dest="output_prefix",help="Prefix of output files.[required]")
	parser.add_option("-f","--format",action="store",type="string",dest="output_format",default='pdf',help="Output file format,'pdf','png' or 'jpg'. default=%default")
	parser.add_option("--ymin",action="store",type="float",dest="ymin",default=None,help="The max of ylim. default=%default")
	parser.add_option("--ymax",action="store",type="float",dest="ymax",default=None,help="The max of ylim. default=%default")
	parser.add_option("--unit",action="store",type="string",dest="unit",default='codon',help="Unit for density calculation.[codon or nt]")
	parser.add_option("--mode",action="store",type="string",dest="mode",default='all',help="Mode for plot, you can chose plot all samples or single sample. [all or single]")
	parser.add_option("--axvline",action="store",type="float",dest="axvline",default=None,help="Position to plot a vertical line in x axis. default=%default")
	parser.add_option("--axhline",action="store",type="float",dest="axhline",default=None,help="Position to plot a vertical line in y axis. default=%default")
	parser.add_option("--label",action="store",type='string',dest='label_list',default=None,help="The legend labels used for plot. default=%default")
	parser.add_option("--slide-window",action="store",type="string",dest="slideWindow",default=None,help="Using slide window to average the density.Input a	true strings such as yes, y or 1. %default=default")
	parser.add_option("--start",action="store",type="int",dest="start_position",default=5,help="The start position need to be averaged.default=%default")
	parser.add_option("--window",action="store",type="int",dest="window",default=7,help="The length of silde window. ddefault=%default")
	parser.add_option("--step",action="store",type='int',dest="step",default=1,help="The step length of slide window. default=%default")
	parser.add_option("--CI",action="store",type='float',dest="CI",default=None,help="plot the confidence intervals or not. If yes, plot the CI region(95% CI default the same as metageneAnalysis.py). else, no. default=%default")


	return parser

def slide_window_average(data,samples,in_regionLengthParma,in_extendRegionLengthParma,inOutPrefix,start,window,step):
	''' Used for calculating mean density with a slide window'''
	data_average=defaultdict(dict)
	winLen=in_regionLengthParma+in_extendRegionLengthParma+1
	columns=data.columns
	flag=0
	for column in columns:
		data_average[column]=[]
		for i in np.arange(len(samples)):
			if flag == 0:
				data_average[column].extend([samples[i]]*winLen)
			else:
				tmp_data=np.zeros(winLen)
				tmp_data[0:int(start)]+=data.iloc[np.where(data.iloc[:,0]==samples[i])].loc[:,column][0:int(start)]
				tmp_data[-int(start):]+=data.iloc[np.where(data.iloc[:,0]==samples[i])].loc[:,column][-int(start):]
				for j in np.arange(start,winLen-start,step):
					tmp_data[j]+=np.mean(data.iloc[np.where(data.iloc[:,0]==samples[i])].loc[:,column][(j-int((window-1)/2)):(j+int((window-1)/2))])
				data_average[column].extend(tmp_data)
		flag+=1
	data_average=pd.DataFrame(data_average)
	data_average.to_csv(inOutPrefix+"_average_denisty.txt",sep="\t",index=0)
	return data_average

def plot_density_for_all_samples(data,samples,Type,in_regionLengthParma,in_extendRegionLengthParma,inOutPrefix,inOutFomat,ymin,ymax,unit,axvline,axhline,confidence,text_font={"size":20,"family":"Arial","weight":"bold"},legend_font={"size":20,"family":"Arial","weight":"bold"}):
	'''plot the density dsitribution'''
	plt.rc("font",weight="bold")
	fig=plt.figure(figsize=(16,8))
	ax=fig.add_subplot(111)
	winLen=in_regionLengthParma+in_extendRegionLengthParma+1
	if len(samples) <=8:
		colors=["b","orangered","green","c","m","y","k","w"]
	else:
		colors=colors=sns.color_palette('husl',len(samples))

	if unit == 'codon':
		for i in np.arange(len(samples)):
			if Type=="start codon":
				plt.plot(np.arange(0,winLen),data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,1],color=colors[i],label=samples[i],linewidth=1)
				if confidence:
					ax.fill_between(np.arange(0,winLen),data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,3],data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,4],color=colors[i],alpha=0.2)
				else:
					pass
				if axvline:
					ax.axvline(axvline,color="r",dashes=[1,2],clip_on=False,linewidth=2)
				else:
					pass

			else:
				plt.plot(np.arange(-winLen,0),data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,2],color=colors[i],label=samples[i],linewidth=1)
				if confidence:
					ax.fill_between(np.arange(-winLen,0),data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,5],data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,6],color=colors[i],alpha=0.2)
				else:
					pass

		ax.set_xlabel("Distance from "+Type + " (codon)",fontdict=text_font)
		ax.set_ylabel("Mean enrichment (A.U)",fontdict=text_font)
		if axhline:
			ax.axhline(axhline,color="r",dashes=[2,3],clip_on=False,linewidth=2)
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
		plt.savefig(inOutPrefix+"_"+Type+"."+inOutFomat,format=inOutFomat)
		# plt.show()
		plt.close()
	elif unit =='nt':
		for i in np.arange(len(samples)):
			if Type=="start codon":
				plt.plot(np.arange(0,winLen),data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,1],color=colors[i],label=samples[i],linewidth=1)
				if confidence:
					ax.fill_between(np.arange(0,winLen),data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,3],data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,4],color=colors[i],alpha=0.2)
				else:
					pass
				if axvline:
					ax.axvline(axvline,color="r",dashes=[1,2],clip_on=False,linewidth=2)
				else:
					pass

			else:
				plt.plot(np.arange(-winLen,0),data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,2],color=colors[i],label=samples[i],linewidth=1)
				if confidence:
					ax.fill_between(np.arange(-winLen,0),data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,5],data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,6],color=colors[i],alpha=0.2)
				else:
					pass
		ax.set_xlabel("Distance from "+Type + " (nt)",fontdict=text_font)
		ax.set_ylabel("Mean enrichment (A.U)",fontdict=text_font)
		if axhline:
			ax.axhline(axhline,color="r",dashes=[2,3],clip_on=False,linewidth=2)
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
		plt.savefig(inOutPrefix+"_"+Type+"."+inOutFomat,format=inOutFomat)
		# plt.show()
		plt.close()
	else:
		raise IOError("Please reset your --unit parameter [codon or nt]")

def plot_density_for_each_sample(data,samples,Type,in_regionLengthParma,in_extendRegionLengthParma,inOutPrefix,inOutFomat,ymin,ymax,unit,axvline,axhline,confidence,text_font={"size":20,"family":"Arial","weight":"bold"},legend_font={"size":20,"family":"Arial","weight":"bold"}):
	'''plot the density dsitribution'''
	plt.rc("font",weight="bold")
	winLen=in_regionLengthParma+in_extendRegionLengthParma+1
	if len(samples) <=8:
		colors=["b","orangered","green","c","m","y","k","w"]
	else:
		colors=colors=sns.color_palette('husl',len(samples))
	sample_dict=dict([i for i in enumerate(samples)])
	sample_dict={j:i for i,j in sample_dict.items()}
	if unit == 'codon':
		for i in sample_dict.keys():
			fig=plt.figure(figsize=(16,8))
			ax=fig.add_subplot(111)
			if Type=="start codon":
				plt.plot(np.arange(0,winLen),data.iloc[np.where(data.iloc[:,0]==i)].iloc[:,1],color=colors[sample_dict[i]],label=i,linewidth=1)
				if confidence:
					ax.fill_between(np.arange(0,winLen),data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,3],data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,4],color=colors[i],alpha=0.2)
				else:
					pass
				if axvline:
					ax.axvline(axvline,color="r",dashes=[1,2],clip_on=False,linewidth=2)
				else:
					pass

			else:
				plt.plot(np.arange(-winLen,0),data.iloc[np.where(data.iloc[:,0]==i)].iloc[:,2],color=colors[sample_dict[i]],label=i,linewidth=1)
				if confidence:
					ax.fill_between(np.arange(-winLen,0),data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,5],data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,6],color=colors[i],alpha=0.2)
				else:
					pass
			ax.set_xlabel("Distance from "+Type + " (codon)",fontdict=text_font)
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
			plt.savefig(inOutPrefix+"_"+i+"_"+Type+"."+inOutFomat,format=inOutFomat)
			# plt.show()
			plt.close()

	elif unit =='nt':
		for i in sample_dict.keys():
			fig=plt.figure(figsize=(16,8))
			ax=fig.add_subplot(111)
			if Type=="start codon":
				plt.plot(np.arange(0,winLen),data.iloc[np.where(data.iloc[:,0]==i)].iloc[:,1],color=colors[sample_dict[i]],label=i,linewidth=1)
				if confidence:
					ax.fill_between(np.arange(0,winLen),data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,3],data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,4],color=colors[i],alpha=0.2)
				else:
					pass
				if axvline:
					ax.axvline(axvline,color="r",dashes=[1,2],clip_on=False,linewidth=2)
				else:
					pass

			else:
				plt.plot(np.arange(-winLen,0),data.iloc[np.where(data.iloc[:,0]==i)].iloc[:,2],color=colors[sample_dict[i]],label=i,linewidth=1)
				if confidence:
					ax.fill_between(np.arange(-winLen,0),data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,5],data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,6],color=colors[i],alpha=0.2)
				else:
					pass
			ax.set_xlabel("Distance from "+Type + " (nt)",fontdict=text_font)
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
			plt.savefig(i+"_"+Type+"."+inOutFomat,format=inOutFomat)
			# plt.show()
			plt.close()
	else:
		raise IOError("Please reset your --unit parameter [codon or nt]")

def change_motif_names(data,type_list):
		samples=data.iloc[:,0]
		if len(type_list) < 1:
			pass
		if len(type_list)==1:
			new_sample=pd.DataFrame([type_list[0]]*len(samples),columns=[samples.name])
		if len(type_list) > 1:
			new_sample=list(reduce(chain,[[i]*(len(samples)//len(type_list)) for i in type_list]))
			new_sample=pd.DataFrame(new_sample,columns=[samples.name])
		data.iloc[:,0]=new_sample
		return data

def main():
	parsed=create_parser()
	(options,args)=parsed.parse_args()
	(data,in_regionLengthParma,in_extendRegionLengthParma,output_prefix,output_format,ymin,ymax,label_list,unit,axvline,axhline,mode,start,window,step,slideWindow,confidence,)=(options.density_file,options.downstream_codon,options.upstream_codon,
	options.output_prefix,options.output_format,options.ymin,options.ymax,options.label_list,options.unit,options.axvline,options.axhline,options.mode,options.start_position,options.window,options.step,options.slideWindow,options.CI)
	if window%2 == 0:
		raise IOError("Please reset your --window parameter. It must be a odd number.")
	if (start-1) < (window-1)/2:
		raise IOError("Please reset your --step and --window parameters. The (window-1)/2 must be less than start-1")
	print("your input file is: "+str(data),file=sys.stderr)
	data=pd.read_csv(data,sep="\t")
	samples=np.unique(data.iloc[:,0])
	if label_list:
		label_list=label_list.strip().split(',')
		data=change_motif_names(data,label_list)
		samples=np.unique(data.iloc[:,0])
	else:
		pass
	## plot density
	text_font={"size":40,"family":"Arial","weight":"bold"}
	legend_font={"size":30,"family":"Arial","weight":"bold"}
	if not slideWindow:
		pass
	else:
		data=slide_window_average(data,samples,in_regionLengthParma,in_extendRegionLengthParma,output_prefix,start,window,step)
	if mode == 'all':
		plot_density_for_all_samples(data,samples,"start codon",in_regionLengthParma,in_extendRegionLengthParma,output_prefix,output_format,ymin,ymax,unit,axvline,axhline,confidence,text_font=text_font,legend_font=legend_font)
		plot_density_for_all_samples(data,samples,"stop codon",in_regionLengthParma,in_extendRegionLengthParma,output_prefix,output_format,ymin,ymax,unit,axvline,axhline,confidence,text_font=text_font,legend_font=legend_font)
	elif mode == 'single':
		plot_density_for_each_sample(data,samples,"start codon",in_regionLengthParma,in_extendRegionLengthParma,output_prefix,output_format,ymin,ymax,unit,axvline,axhline,confidence,text_font=text_font,legend_font=legend_font)
		plot_density_for_each_sample(data,samples,"stop codon",in_regionLengthParma,in_extendRegionLengthParma,output_prefix,output_format,ymin,ymax,unit,axvline,axhline,confidence,text_font=text_font,legend_font=legend_font)
	else:
		raise IOError("Please reset your --mode parameter.[all or single]")

	print("finished plot the ribosome footprint density",file=sys.stderr)

if __name__ =="__main__":
	main()

