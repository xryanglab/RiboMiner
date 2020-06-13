#########################################################################
# File Name: fileTransfer.sh
# Author: lifj
# mail: lfj17@mails.tsinghua.edu.cn
# Created Time: Tue 09 Jun 2020 10:08:02 PM CST
#########################################################################
#!/bin/bash
#BSUB -J fileTransfer.sh
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -n 4
#BSUB -q TEST-A

mkdir -p ../results/output

cp ../results/periodicity/eIF5A_SRR5008135*.pdf ../results/output
cp ../results/length/eIF5A_SRR5008135_reads_length.pdf ../results/output
cp ../results/RPF_statisitcs/* ../results/output
cp ../results/RiboDensityFromDiffFrames/eIF5A_si-Ctrl-2_reading_frames.pdf ../results/output
cp ../results/MetageneAnalysis/eIF5A_whole_mean_metaplot.pdf ../results/MetageneAnalysis/eIF5A_mean_polarity.pdf ../results/output
cp ../results/MetageneAnalysis/eIF5A_CDS_normed_mean_start*.pdf ../results/MetageneAnalysis/eIF5A_CDS_normed_mean_stop*.pdf ../results/output
cp ../results/MetageneAnalysis/eIF5A_UTR_normed_mean_start*.pdf ../results/MetageneAnalysis/eIF5A_UTR_normed_mean_stop*.pdf ../results/output
cp ../results/FeatureAnalysis/up_density_on_each_kind_of*.pdf ../results/output
cp ../results/FeatureAnalysis/up*.pdf ../results/output
cp ../results/FeatureAnalysis/eIF5A*_average_*.pdf ../results/output

echo "files transfered finished!"
