#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
@Author: Li Fajin
@Date: 2019-08-19 15:56:32
@LastEditors: Li Fajin
@LastEditTime: 2019-08-30 17:14:20
@Description: The script is used for plot the density around a specific di- or tri-AA motif.
'''

import sys
import numpy as np
import pandas as pd
from functools import reduce
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from optparse import OptionParser
from .__init__ import __version__


def create_parse_for_AAmotifs():
	usage="usage: python %prog [options]" + "\n"
	parser=OptionParser(usage=usage,version=__version__)
	parser.add_option('-i','--input',action='store',type='string',dest='density_file',
		help='Input file in dataframe format (pd.DataFrame).')
	parser.add_option('-o','--output',action='store',type='string',dest='output_prefix',
		help='Output files')
	parser.add_option('-g','--group',action="store",type="string",dest="group_name",
		help="Group name of each group separated by comma. e.g. 'si-control,si-eIF3e'")
	parser.add_option('-r','--replicate',action="store",type="string",dest="replicate_name",
		help="Replicate name of each group separated by comma. e.g. 'si_3e_1_80S,si_3e_2_80S__si_cttl_1_80S,si_ctrl_2_80S'")
	parser.add_option("--ymax",action="store",type="float",dest="ymax",default=None,
		help="The max of ylim. default=%default")
	parser.add_option("--ymin",action="store",type="float",dest="ymin",default=None,
		help="The min of ylim. default=%default")
	parser.add_option('--mode',action="store",type="string",dest="mode",default="mean",
		help="Control the mode for plot. if '--mode single', return the plot of each samples;else return the mean plot. default=%default.")
	return parser

def DrawMotifDensity_for_replicates_of_different_groups(data,motifs,groups,replicates,output,ymin,ymax,distance=101):
	''' plot the motif density of dfiierent replicates from different groups'''
	## prepare data for plot
	label_dict={}
	data_dict={}
	for g,r in zip(groups,replicates):
		label_dict[g]=r.strip().split(',')
	## separate the data into different groups
	for g in groups:
		columns=['motif']+label_dict[g]
		data_dict[g]=data.loc[:,columns]
	plt.rc('font',weight='bold')
	# colors=['blue','orangered']
	if len(groups) <=8:
		colors=["b","orangered","green","c","m","y","k","w"]
	else:
		colors=colors=sns.color_palette('husl',len(groups))
	text_font={"size":15,"family":"Arial","weight":"bold"}
	with PdfPages(output + "_density_on_triAAMotifs.pdf") as pdf:
		x=np.arange(distance,dtype="int64")
		for motif in motifs:
			fig=plt.figure(figsize=(5,4))
			ax=fig.add_subplot(111)
			for g in groups:
				for i in np.arange(len(label_dict[g])):
					plt.plot(x,data_dict[g][data_dict[g].iloc[:,0]==motif].loc[:,label_dict[g][i]],linewidth=1,color=colors[i],label=label_dict[g][i],alpha=0.8)
			ax.set_xlabel('Distance of P/E site to ' + str(motif)+' (nt)',fontdict=text_font)
			ax.set_ylabel('Relative footprint density(A.U)',fontdict=text_font)
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.spines['left'].set_linewidth(2)
			ax.spines['bottom'].set_linewidth(2)
			ax.tick_params(which="both",width=2)
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
			plt.xticks(np.arange(0,distance,25),(np.arange(0,distance,25)-int((distance-1)/2)))
			plt.legend(loc='best',prop=text_font)
			plt.tight_layout()
			pdf.savefig(fig)
			plt.close()



def DrawMotifDensity_for_mean_denisty(data_mean,motifs,output,ymin,ymax,distance=101):
	''' plot the motif density of dfiierent replicates from different groups'''
	plt.rc('font',weight='bold')
	groups_name=data_mean.columns[1:]
	if len(groups_name) <=8:
		colors=["b","orangered","green","c","m","y","k","w"]
	else:
		colors=colors=sns.color_palette('husl',len(groups_name))
	text_font={"size":15,"family":"Arial","weight":"bold"}
	with PdfPages(output + "_mean_density_on_triAAMotifs.pdf") as pdf:
		x=np.arange(distance,dtype="int64")
		for motif in motifs:
			fig=plt.figure(figsize=(5,4))
			ax=fig.add_subplot(111)
			for g in np.arange(len(groups_name)):
				plt.plot(x,data_mean[data_mean.iloc[:,0]==motif].loc[:,groups_name[g]],linewidth=1,color=colors[g],label=groups_name[g],alpha=0.8)
			ax.set_xlabel('Distance of P/E site to ' + str(motif)+' (nt)',fontdict=text_font)
			ax.set_ylabel('Relative footprint density(A.U)',fontdict=text_font)
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.spines['left'].set_linewidth(2)
			ax.spines['bottom'].set_linewidth(2)
			ax.tick_params(which="both",width=2)
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
			plt.xticks(np.arange(0,distance,25),(np.arange(0,distance,25)-int((distance-1)/2)))
			plt.legend(loc='best',prop=text_font)
			plt.tight_layout()
			pdf.savefig(fig)
			plt.close()


def calculate_mean_density_for_triplete_AA_motif(data,groups,replicates,output_prefix):
	'''calculate the mean values among different replicates. the input data have N+1 columns, N is the number of samples
	1) col 1: motifs
	2) col2-colN: different samples
	'''
	label_dict={}
	data_dict={}
	data_mean_dict=defaultdict(dict)
	motifs=data.iloc[:,0]
	for g,r in zip(groups,replicates):
		label_dict[g]=r.strip().split(',')
	## separate the data into different groups
	for g in groups:
		data_dict[g]=data.loc[:,label_dict[g]]

	for g in groups:
		for i in np.arange(data.shape[0]):
			density=np.mean(data_dict[g].loc[i,label_dict[g]])
			data_mean_dict[g][i]=density
	## transform the dict to a python dataframe
	for g in groups:
		data_mean_dict[g]=pd.DataFrame(data_mean_dict[g],index=[g]).T
	## concatenate different data frame
	data_mean=pd.concat([v for v in data_mean_dict.values()],axis=1)
	data_mean=pd.concat((motifs,data_mean),axis=1)
	## write the mean density file
	data_mean.to_csv(output_prefix+"_mean_motifDensity_dataframe.txt",sep="\t",index=0)
	return data_mean


def parse_args_for_plot_density_of_triplete_AA_motif():
	parser=create_parse_for_AAmotifs()
	(options,args)=parser.parse_args()
	(data,output,groups,replicates,ymin,ymax,mode)=(options.density_file,options.output_prefix,options.group_name.strip().split(','),
	options.replicate_name.strip().split('__'),options.ymin,options.ymax,options.mode)
	##  input data
	data=pd.read_csv(data,sep='\t')
	motifs=np.unique(data.iloc[:,0].values)
	if mode == 'single':
		print("Start the step of plotting...",file=sys.stderr)
		DrawMotifDensity_for_replicates_of_different_groups(data,motifs,groups,replicates,output,ymin,ymax,distance=101)
		print('finishing the plot step!',file=sys.stderr)
	elif mode == 'mean':
		## calculate the mean value of different replicates
		print("Start the step of mean value calculation...",file=sys.stderr)
		data_mean=calculate_mean_density_for_triplete_AA_motif(data,groups,replicates,output)
		print('finishing the mean value calculation!',file=sys.stderr)
		print("Start the step of plotting...",file=sys.stderr)
		DrawMotifDensity_for_mean_denisty(data_mean,motifs,output,ymin,ymax,distance=101)
		print('finishing the plot step!',file=sys.stderr)
	else:
		raise IOError("Please reset your --mode parameter.[single or mean]")

def main():
	parse_args_for_plot_density_of_triplete_AA_motif()


if __name__=="__main__":
	main()
