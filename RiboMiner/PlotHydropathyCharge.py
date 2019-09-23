#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
@Author: Li Fajin
@Date: 2019-08-19 09:49:32
@LastEditors: Li Fajin
@LastEditTime: 2019-09-01 21:51:34
@Description: This script is used for hydropathy or charge plot
'''

import numpy as np
import pandas as pd
import sys
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from optparse import OptionParser
from functools import reduce
from .__init__ import __version__

def create_parser_for_hydropathy_plot():
	'''argument parser'''
	usage="usage: python %prog [options]" +"\n"
	parser=OptionParser(usage=usage,version=__version__)
	parser.add_option("-i","--input",action="store",type="string",dest="density_file",help="Input file in txt format.And the files has three columns; column 1: sample;columns 2: start_density; column 3: stop_density")
	parser.add_option("-d","--downstream_codon",action="store",type="int",default=500, dest="downstream_codon", help="Downstream codon corresponding to start codon (codon unit). While corresponding to stop codon, it is the upstream codon.")
	parser.add_option("-u","--upstream_codon",action="store",type="int",default=0,dest="upstream_codon",help="Upstream codon corresponding to start codon (codon unit). While corresponding to stop codon, it is the downstream codon.")
	parser.add_option("-o","--output",action="store",type="string",dest="output_prefix",help="Prefix of output files.[required]")
	parser.add_option("-f","--format",action="store",type="string",dest="output_format",default='pdf',help="Output file format,'pdf','png' or 'jpg'. default=%default")
	parser.add_option("-l",'--left-position',action="store",type="int",dest='left_position',default=None,help="The left position of the region that you are interested in. default=%default")
	parser.add_option("-r",'--right-position',action='store',type='int',dest='right_position',default=None,help="The right position of the region that you are interested in. default=%default")
	parser.add_option("--ymin",action="store",type="float",dest="ymin",default=None,help="The max of ylim. default=%default")
	parser.add_option("--ymax",action="store",type="float",dest="ymax",default=None,help="The max of ylim. default=%default")
	parser.add_option("--mode",action="store",type="string",dest="mode",default='all',help="Control the mode for plot.[all or single]. default=%default")
	parser.add_option("--axvline",action="store",type="float",dest="axvline",default=None,help="Position to plot vetical line")
	parser.add_option("--ylab",action="store",type="string",dest="ylab",default="Average Index",help="labels for y axis.'Average hydrophobicity' or 'Average charges'")
	parser.add_option("--start",action="store",type="int",dest="start_position",default=5,help="The start position need to be averaged.default=%default")
	parser.add_option("--window",action="store",type="int",dest="window",default=7,help="The length of silde window. ddefault=%default")
	parser.add_option("--step",action="store",type='int',dest="step",default=1,help="The step length of slide window. default=%default")

	return parser

def plot_all_density(data,samples,type,in_regionLengthParma,in_extendRegionLengthParma,inOutPrefix,inOutFomat,ymin,ymax,axvline,ylab=None,text_font={"size":20,"family":"Arial","weight":"bold"},legend_font={"size":20,"family":"Arial","weight":"bold"}):
	'''plot the density dsitribution'''
	plt.rc("font",weight="bold")
	fig=plt.figure(figsize=(16,8))
	ax=fig.add_subplot(111)
	winLen=in_regionLengthParma+in_extendRegionLengthParma+1
	if len(samples) <=8:
		colors=["b","orangered","green","c","m","y","k","w"]
	else:
		colors=colors=sns.color_palette('husl',len(samples))

	for i in np.arange(len(samples)):
		if type=="start codon":
			plt.plot(np.arange(0,winLen),data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,1],color=colors[i],label=samples[i],linewidth=1)
			ax.set_xticks(np.arange(0,winLen,50))
			ax.set_xticklabels((np.arange(0,winLen,50)-in_extendRegionLengthParma))
			if axvline:
				ax.axvline(axvline,color='red',dashes=[3,2])
			else:
				pass

		else:
			plt.plot(np.arange(0,winLen),data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,2],color=colors[i],label=samples[i],linewidth=1)
			ax.set_xticks(np.arange(0,winLen,50))
			ax.set_xticklabels((np.arange(0,winLen,50)-in_regionLengthParma))

	ax.set_xlabel("Distance from "+type + ' (codon)',fontdict=text_font)
	ax.set_ylabel(ylab,fontdict=text_font)
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
	# plt.legend(loc="best",prop={'size':10,'weight':'bold'})
	plt.tight_layout()
	plt.savefig(inOutPrefix+"_"+type+"."+inOutFomat,format=inOutFomat)
	plt.close()

def plot_density_for_each_sample(data,samples,type,in_regionLengthParma,in_extendRegionLengthParma,inOutPrefix,inOutFomat,ymin,ymax,axvline,ylab=None,text_font={"size":20,"family":"Arial","weight":"bold"},legend_font={"size":20,"family":"Arial","weight":"bold"}):
	'''plot the hydropathy'''
	plt.rc('font',weight='bold')
	winLen=in_regionLengthParma+in_extendRegionLengthParma+1
	colors=["b","orangered","green","c","m","y","k","w"]
	for i in np.arange(len(samples)):
		if type=="start codon":
			fig=plt.figure(figsize=(16,8))
			ax=fig.add_subplot(111)
			plt.plot(np.arange(0,winLen),data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,1],color=colors[i],label=samples[i],linewidth=1)
			ax.set_xticks(np.arange(0,winLen,50))
			ax.set_xticklabels((np.arange(0,winLen,50)-in_extendRegionLengthParma))
			ax.set_xlabel("Distance from "+type+" (codon)",fontdict=text_font)
			ax.set_ylabel(ylab,fontdict=text_font)
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
			plt.savefig(inOutPrefix+"_"+str(samples[i])+"_start_codon"+"."+inOutFomat,format=inOutFomat)
			plt.close()

		else:
			fig=plt.figure(figsize=(16,8))
			ax=fig.add_subplot(111)
			plt.plot(np.arange(0,winLen),data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,2],color=colors[i],label=samples[i],linewidth=1)
			ax.set_xticks(np.arange(0,winLen,50))
			ax.set_xticklabels((np.arange(0,winLen,50)-in_regionLengthParma))
			ax.set_xlabel("Distance from "+type+" (codon)",fontdict=text_font)
			ax.set_ylabel(ylab,fontdict=text_font)
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
			plt.savefig(inOutPrefix+"_"+str(samples[i])+"_stop_codon"+"."+inOutFomat,format=inOutFomat)
			plt.close()
def plot_density_for_specific_region(data,samples,type,left_position,right_position,inOutPrefix,inOutFomat,ymin,ymax,axvline,mode,ylab=None,text_font={"size":20,"family":"Arial","weight":"bold"},legend_font={"size":20,"family":"Arial","weight":"bold"}):
	'''plot the density dsitribution for a specific region'''
	plt.rc("font",weight="bold")
	winLen=right_position-left_position+1
	if len(samples) <=8:
		colors=["b","orangered","green","c","m","y","k","w"]
	else:
		colors=colors=sns.color_palette('husl',len(samples))
	if mode == "all":
		fig=plt.figure(figsize=(16,8))
		ax=fig.add_subplot(111)
		for i in np.arange(len(samples)):
			if type=="start codon":
				plt.plot(np.arange(0,winLen),data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,1][(left_position-1):right_position],color=colors[i],label=samples[i],linewidth=1)
				ax.set_xticks(np.arange(0,winLen+5,5))
				ax.set_xticklabels(np.arange(left_position-1,right_position+5,5))
				if axvline:
					ax.axvline(axvline,color='red',dashes=[3,2])
				else:
					pass

			else:
				plt.plot(np.arange(0,winLen),data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,2][(left_position-1):right_position],color=colors[i],label=samples[i],linewidth=1)
				ax.set_xticks(np.arange(0,winLen,5))
				ax.set_xticklabels((np.arange(0,winLen,5)-right_position))

		ax.set_xlabel("Distance from "+type + ' (codon)',fontdict=text_font)
		ax.set_ylabel(ylab,fontdict=text_font)
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
		# plt.legend(loc="best",prop={'size':10,'weight':'bold'})
		plt.tight_layout()
		plt.savefig(inOutPrefix+"__specific_region_"+type+"."+inOutFomat,format=inOutFomat)
		plt.close()
	elif mode == "single":
		for i in np.arange(len(samples)):
			if type=="start codon":
				fig=plt.figure(figsize=(16,8))
				ax=fig.add_subplot(111)
				plt.plot(np.arange(0,winLen),data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,1][left_position-1:right_position],color=colors[i],label=samples[i],linewidth=1)
				ax.set_xticks(np.arange(0,winLen+5,5))
				ax.set_xticklabels(np.arange(left_position-1,right_position+5,5))
				if axvline:
					ax.axvline(axvline,color='red',dashes=[3,2])
				else:
					pass
				ax.set_xlabel("Distance from "+type+" (codon)",fontdict=text_font)
				ax.set_ylabel(ylab,fontdict=text_font)
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
				plt.savefig(inOutPrefix+"_"+str(samples[i])+"_start_codon"+"."+inOutFomat,format=inOutFomat)
				plt.close()

			else:
				fig=plt.figure(figsize=(16,8))
				ax=fig.add_subplot(111)
				plt.plot(np.arange(0,winLen),data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,2][left_position-1:right_position],color=colors[i],label=samples[i],linewidth=1)
				ax.set_xticks(np.arange(0,winLen+5,5))
				ax.set_xticklabels(np.arange(left_position-1,right_position+5,5))
				if axvline:
					ax.axvline(axvline,color='red',dashes=[3,2])
				else:
					pass
				ax.set_xlabel("Distance from "+type+" (codon)",fontdict=text_font)
				ax.set_ylabel(ylab,fontdict=text_font)
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
				plt.savefig(inOutPrefix+"_"+str(samples[i])+"_stop_codon"+"."+inOutFomat,format=inOutFomat)
				plt.close()
	else:
		raise IOError("Please enter a correct --mode parameter![all or single]")

def slide_window_average(data,samples,in_regionLengthParma,in_extendRegionLengthParma,inOutPrefix,start,window,step):
	start_average=[]
	stop_average=[]
	label=[]
	winLen=in_regionLengthParma+in_extendRegionLengthParma+1
	for i in np.arange(len(samples)):
			tmp1_data=np.zeros(winLen)
			tmp1_data[0:int(start)]+=data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,1][0:int(start)]
			tmp1_data[-int(start):]+=data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,1][-int(start):]
			for j in np.arange(start,winLen-start,step):
				tmp1_data[j]+=np.mean(data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,1][(j-int((window-1)/2)):(j+int((window-1)/2))])
			start_average.extend(tmp1_data)
			label.extend([samples[i]]*winLen)

			tmp2_data=np.zeros(winLen)
			tmp2_data[0:int(start)]+=data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,2][0:int(start)]
			tmp2_data[-int(start):]+=data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,2][-int(start):]
			for j in np.arange(start,winLen-start,step):
				tmp2_data[j]+=np.mean(data.iloc[np.where(data.iloc[:,0]==samples[i])].iloc[:,2][(j-int((window-1)/2)):(j+int((window-1)/2))])
			stop_average.extend(tmp2_data)

	data_average=pd.DataFrame([label,start_average,stop_average],index=['sample','start_density','stop_density'])
	data_average=data_average.T
	data_average.to_csv(inOutPrefix+"_average_frame.txt",sep="\t",index=0)
	return data_average
def main():
	parsed=create_parser_for_hydropathy_plot()
	(options,args)=parsed.parse_args()
	(inputFile,in_regionLengthParma,in_extendRegionLengthParma,output_prefix,output_format,ymin,ymax,left_position,right_position,axvline,mode,ylab,start,window,step)=(options.density_file,options.downstream_codon,options.upstream_codon,
	options.output_prefix,options.output_format,options.ymin,options.ymax,options.left_position,options.right_position,options.axvline,options.mode,options.ylab,options.start_position,options.window,options.step)
	if window%2 == 0:
		raise IOError("Please reset your --window parameter. It must be a odd number.")
	if (start-1) < (window-1)/2:
		raise IOError("Please reset your --step and --window parameters. The (window-1)/2 must be less than start-1")
	data=pd.read_csv(inputFile,sep="\t")
	samples=np.unique(data.iloc[:,0])
	text_font={"size":30,"family":"Arial","weight":"bold"}
	legend_font={"size":30,"family":"Arial","weight":"bold"}
	data_average=slide_window_average(data,samples,in_regionLengthParma,in_extendRegionLengthParma,output_prefix,start,window,step)
	if mode == 'all':
		if (not left_position and right_position) or (left_position and not right_position):
			raise IOError("Do you want to plot in a specific region? Input both the left and right position please!")
		elif ((not left_position) and (not right_position)):
			plot_all_density(data,samples,"start codon",in_regionLengthParma,in_extendRegionLengthParma,output_prefix,output_format,ymin,ymax,axvline,ylab=ylab,text_font=text_font,legend_font=legend_font)
			plot_all_density(data,samples,"stop codon",in_regionLengthParma,in_extendRegionLengthParma,output_prefix,output_format,ymin,ymax,axvline,ylab=ylab,text_font=text_font,legend_font=legend_font)
			plot_all_density(data_average,samples,"start codon",in_regionLengthParma,in_extendRegionLengthParma,output_prefix+"_average",output_format,ymin,ymax,axvline,ylab=ylab,text_font=text_font,legend_font=legend_font)
			plot_all_density(data_average,samples,"stop codon",in_regionLengthParma,in_extendRegionLengthParma,output_prefix+"_average",output_format,ymin,ymax,axvline,ylab=ylab,text_font=text_font,legend_font=legend_font)
			print("finished plot the hydropathy or charge density",file=sys.stderr)
		else:
			if left_position>=right_position:
				raise IOError("The left position must less than the right position. Please reset you -l and -r parameters!")
			if right_position>abs((in_extendRegionLengthParma-in_regionLengthParma)):
				raise IOError("The right position is out of region. Please reset your -r parameer!")
			plot_density_for_specific_region(data,samples,"start codon",left_position,right_position,output_prefix,output_format,ymin,ymax,axvline,mode,ylab=ylab,text_font=text_font,legend_font=legend_font)
			plot_density_for_specific_region(data,samples,"stop codon",left_position,right_position,output_prefix,output_format,ymin,ymax,axvline,mode,ylab=ylab,text_font=text_font,legend_font=legend_font)
			plot_density_for_specific_region(data_average,samples,"start codon",left_position,right_position,output_prefix+"_average",output_format,ymin,ymax,axvline,mode,ylab=ylab,text_font=text_font,legend_font=legend_font)
			plot_density_for_specific_region(data_average,samples,"stop codon",left_position,right_position,output_prefix+"_average",output_format,ymin,ymax,axvline,mode,ylab=ylab,text_font=text_font,legend_font=legend_font)
			print("finished plot the hydropathy or charge density",file=sys.stderr)
	elif mode == 'single':
		if (not left_position and right_position) or (left_position and not right_position):
			raise IOError("Do you want to plot in a specific region? Input both the left and right position please!")
		elif ((not left_position) and (not right_position)):
			plot_density_for_each_sample(data,samples,"start codon",in_regionLengthParma,in_extendRegionLengthParma,output_prefix,output_format,ymin,ymax,axvline,ylab=ylab,text_font=text_font,legend_font=legend_font)
			plot_density_for_each_sample(data,samples,"stop codon",in_regionLengthParma,in_extendRegionLengthParma,output_prefix,output_format,ymin,ymax,axvline,ylab=ylab,text_font=text_font,legend_font=legend_font)
			print("finished plot the hydropathy or charge density",file=sys.stderr)
		else:
			if left_position>=right_position:
				raise IOError("The left position must less than the right position. Please reset you -l and -r parameters!")
			if right_position>abs((in_extendRegionLengthParma-in_regionLengthParma)):
				raise IOError("The right position is out of region. Please reset your -r parameer!")
			plot_density_for_specific_region(data,samples,"start codon",left_position,right_position,output_prefix,output_format,ymin,ymax,axvline,mode,ylab=ylab,text_font=text_font,legend_font=legend_font)
			plot_density_for_specific_region(data,samples,"stop codon",left_position,right_position,output_prefix,output_format,ymin,ymax,axvline,mode,ylab=ylab,text_font=text_font,legend_font=legend_font)
			print("finished plot the hydropathy or charge density",file=sys.stderr)
	else:
		raise IOError("Please enter a correct --mode parameter![all or single]")

if __name__ =="__main__":
	main()

