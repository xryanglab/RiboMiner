#########################################################################
# File Name: MA.sh
# Author: lifj
# mail: lfj17@mails.tsinghua.edu.cn
# Created Time: Tue 09 Jun 2020 02:47:46 PM CST
#########################################################################
#!/bin/bash
#BSUB -J MA.sh
#BSUB -o MA.out
#BSUB -e MA.err
#BSUB -n 1
#BSUB -q TEST-A


results=../results
mkdir -p $results/MetageneAnalysis

MetageneAnalysisForTheWholeRegions -f ../data/attributes_P.txt -c ../results/Ref/longest.transcripts.info.txt -o $results/MetageneAnalysis/test -b 15,90,15 -l 100 -n 10 -m 1 -e 5 

PlotMetageneAnalysisForTheWholeRegions -i $results/MetageneAnalysis/test_scaled_density_dataframe.txt -o $results/MetageneAnalysis/test_whole -g Ctrl,ASNase -r Ctrl-1__ASNase-1 -b 15,90,15 --mode all

MetageneAnalysis -f ../data/attributes_P.txt -c ../results/Ref/longest.transcripts.info.txt -o $results/MetageneAnalysis/test_CDS_normed -U codon -M RPKM -u 0 -d 500 -l 100 -n 10 -m 1 -e 5 --norm yes -y 100 --CI 0.95 --type CDS
MetageneAnalysis -f ../data/attributes_P.txt -c ../results/Ref/longest.transcripts.info.txt -o $results/MetageneAnalysis/test_UTR_normed -U nt -M RPKM -u 50 -d 50 -l 50 -n 10 -m 1 -e 5 --norm yes -y 50 --CI 0.95 --type UTR

PolarityCalculation -f ../data/attributes_P.txt -c ../results/Ref/longest.transcripts.info.txt -o $results/MetageneAnalysis/test -n 64
PlotMetageneAnalysis -i $results/MetageneAnalysis/test_CDS_normed_dataframe.txt -o $results/MetageneAnalysis/test_CDS_normed -g Ctrl,ASNase -r Ctrl-1__ASNase-1 -U codon -u 0 -d 500 --mode mean --CI 0.95 --axhline 1

PlotMetageneAnalysis -i $results/MetageneAnalysis/test_UTR_normed_dataframe.txt -o $results/MetageneAnalysis/test_UTR_normed -g Ctrl,ASNase -r Ctrl-1__ASNase-1 -U nt -u 50 -d 50 --mode mean --CI 0.95 --axhline 1

PlotPolarity -i  $results/MetageneAnalysis/test_polarity_dataframe.txt -o $results/MetageneAnalysis/test -g Ctrl,ASNase -r Ctrl-1__ASNase-1 -y 5 --mode all

RiboDensityForSpecificRegion -f ../data/attributes_P.txt -c ../results/Ref/longest.transcripts.info.txt -o $results/MetageneAnalysis/test_all -U codon -M RPKM -L 1 -R 100

Rscript DivideGeneSets.R

MetageneAnalysis -f ../data/attributes_P.txt -c ../results/Ref/longest.transcripts.info.txt -o $results/MetageneAnalysis/test_up_CDS_normed -U codon -M RPKM -u 0 -d 500 -l 100 -n 5 -m 1 -e 5 --norm yes -y 100 --CI 0.95 --type CDS -S $results/MetageneAnalysis/up_trans.txt
MetageneAnalysis -f ../data/attributes_P.txt -c ../results/Ref/longest.transcripts.info.txt -o $results/MetageneAnalysis/test_down_CDS_normed -U codon -M RPKM -u 0 -d 500 -l 100 -n 5 -m 1 -e 5 --norm yes -y 100 --CI 0.95 --type CDS -S $results/MetageneAnalysis/down_trans.txt
MetageneAnalysis -f ../data/attributes_P.txt -c ../results/Ref/longest.transcripts.info.txt -o $results/MetageneAnalysis/test_unblocked_CDS_normed -U codon -M RPKM -u 0 -d 500 -l 100 -n 5 -m 1 -e 5 --norm yes -y 100 --CI 0.95 --type CDS -S $results/MetageneAnalysis/unblocked_trans.txt

Rscript GetCommonGeneSets.R
