#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
@Author: Li Fajin
@Date: 2019-08-15 21:29:12
@LastEditors: Li Fajin
@LastEditTime: 2019-09-18 20:38:23
@Description: This script is used for checking periodicity of ribosome profiling data, but without P-site identification.
And the part code are adapted from RiboCode our lab developed before. [Xiao, et al. NAR.2018]
usage: python Periodicity -i bam -a RiboCode_annote -c longest.trans.info.txt -o outprefix -L 25 -R 35 --id-type transcript-id
'''


from .FunctionDefinition import *
from functools import reduce
import pickle
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from optparse import OptionParser


def load_transcripts_pickle(pickle_file):
	if os.path.exists(pickle_file):
		sys.stdout.write("\tLoading transcripts.pickle ...\n")
		with open(pickle_file,"rb") as fin:
			gene_dict, transcript_dict = pickle.load(fin)
	else:
		raise IOError("\tError, %s file not found\n" % pickle_file)

	return gene_dict,transcript_dict

def periodicity(in_bamFile,in_selectTrans,transcript_dict,left_length,right_length):
	pysamFile=pysam.AlignmentFile(in_bamFile,"rb")
	pysamFile_trans=pysamFile.references
	in_selectTrans=set(pysamFile_trans).intersection(in_selectTrans)
	start_density=defaultdict(lambda:np.zeros(101,dtype="int64"))
	stop_density=defaultdict(lambda:np.zeros(101,dtype="int64"))
	total_reads=0
	specific_counts=defaultdict(int)
	i=0
	for trans in in_selectTrans:
		i+=1
		for record in pysamFile.fetch(trans):
			if record.flag == 16 or record.flag == 272:
				continue
			total_reads += 1
			specific_counts[record.query_length]+=1
			if trans not in transcript_dict.keys():
				continue
			if left_length<=record.query_length<=right_length:
				distance_to_start=record.reference_start-transcript_dict[trans].startcodon.start
				dsitance_to_stop=record.reference_start-transcript_dict[trans].stopcodon.end
				if abs(distance_to_start)<=50:
					start_density[record.query_length][50+distance_to_start]+=1
				if abs(dsitance_to_stop) <= 50:
					stop_density[record.query_length][50+dsitance_to_stop]+=1

	specific_counts['sum']=total_reads
	print("There are " + str(i)+" used for statistics.",file=sys.stderr)
	return start_density,stop_density,total_reads,specific_counts

def plot_periodicity_start_codon(start_density,specific_counts,output_prefix):
	start_density=pd.DataFrame(start_density)
	start_density=start_density.reindex(columns=sorted(start_density.columns))
	start_density['sum']=start_density.apply(sum,axis=1)
	start_density.loc['sum']=np.array([specific_counts[key] for key in start_density.columns])
	start_density.to_csv(output_prefix+"_from_start_codon.txt",index=0,sep="\t")
	text_font={"size":15,"family":"Arial","weight":"bold"}
	plt.rc('font',weight='bold')
	with PdfPages(output_prefix + "_start_periodicity.pdf") as pdf:
		x = np.arange(-50,51,dtype=int)
		colors = sns.color_palette('husl',3)*34
		for L in start_density.columns:
			xticks=[-40,-20,0,20,40]
			perct = '{:.2%}'.format(specific_counts[L] / specific_counts['sum'])
			fig,ax1 = plt.subplots(figsize=(5.6,2.8))
			y1=start_density.loc[:,L]
			ax1.vlines(x,ymin=np.zeros(101),ymax=y1,colors=colors[:-1],linewidth=2)
			ax1.tick_params(axis='x',which="both",top=False,direction='out')
			ax1.set_xticks(xticks)
			ax1.set_xlim((-50,50))
			ax1.spines["top"].set_visible(False)
			ax1.spines["right"].set_visible(False)
			ax1.spines["bottom"].set_linewidth(2)
			ax1.spines["left"].set_linewidth(2)
			ax1.set_xlabel("Distance from start codon (nt)",fontdict=text_font)
			ax1.set_ylabel("Alignments",fontdict=text_font)
			if not L == 'sum':
				ax1.set_title("{} nt reads,proportion:{}".format(L,perct),fontdict=text_font)
			else:
				ax1.set_title("Total reads, proporation:{}".format(perct),fontdict=text_font)
			fig.tight_layout()
			pdf.savefig(fig)
			plt.close()
	return None
def plot_periodicity_stop_codon(stop_density,specific_counts,output_prefix):
	stop_density=pd.DataFrame(stop_density)
	stop_density=stop_density.reindex(columns=sorted(stop_density.columns))
	stop_density['sum']=stop_density.apply(sum,axis=1)
	stop_density.loc['sum']=np.array([specific_counts[key] for key in stop_density.columns])
	stop_density.to_csv(output_prefix+"_from_stop_codon.txt",index=0,sep="\t")
	text_font={"size":15,"family":"Arial","weight":"bold"}
	plt.rc('font',weight='bold')
	with PdfPages(output_prefix + "_stop_periodicity.pdf") as pdf:
		x = np.arange(-50,51,dtype=int)
		colors = sns.color_palette('husl',3)*34
		for L in stop_density.columns:
			xticks=[-40,-20,0,20,40]
			perct = '{:.2%}'.format(specific_counts[L] / specific_counts['sum'])
			fig,ax1 = plt.subplots(figsize=(5.6,2.8))
			y1=stop_density.loc[:,L]
			ax1.vlines(x,ymin=np.zeros(101),ymax=y1,colors=colors[:-1],linewidth=2)
			ax1.tick_params(axis='x',which="both",top=False,direction='out')
			ax1.set_xticks(xticks)
			ax1.set_xlim((-50,50))
			ax1.spines["top"].set_visible(False)
			ax1.spines["right"].set_visible(False)
			ax1.spines["bottom"].set_linewidth(2)
			ax1.spines["left"].set_linewidth(2)
			ax1.set_xlabel("Distance from stop codon (nt)",fontdict=text_font)
			ax1.set_ylabel("Alignments",fontdict=text_font)
			if not L == 'sum':
				ax1.set_title("{} nt reads,proportion:{}".format(L,perct),fontdict=text_font)
			else:
				ax1.set_title("Total reads, proporation:{}".format(perct),fontdict=text_font)
			fig.tight_layout()
			pdf.savefig(fig)
			plt.close()
	return None

def main():
	parsed=create_parser_for_periodicity()
	(options,args)=parsed.parse_args()
	## calculate density for each bam files
	if not options.annot_dir:
		raise IOError("Please run RiboCode::prepare_transcripts first")
	else:
		gene_dict,transcript_dict = load_transcripts_pickle(os.path.join(options.annot_dir,"transcripts.pickle"))
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
	print("Start calculating read density...",file=sys.stderr)
	start_density,stop_density,total_reads,specific_counts=periodicity(options.bamFile,select_trans,transcript_dict,options.left_length,options.right_length)
	plot_periodicity_start_codon(start_density,specific_counts,options.output_prefix)
	plot_periodicity_stop_codon(stop_density,specific_counts,options.output_prefix)
	print("Finish the step of plot periodicity",file=sys.stderr)


if __name__=="__main__":
	main()