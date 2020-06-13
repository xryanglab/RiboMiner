#########################################################################
# File Name: EA.sh
# Author: lifj
# mail: lfj17@mails.tsinghua.edu.cn
# Created Time: Tue 09 Jun 2020 02:48:15 PM CST
#########################################################################
#!/bin/bash
#BSUB -J EA.sh
#BSUB -o EA.out
#BSUB -e EA.err
#BSUB -n 1
#BSUB -q TEST-A


results=../results

mkdir -p $results/EnrichmentAnalysis

RiboDensityAtEachPosition -c ../results/Ref/longest.transcripts.info.extended.txt -f ../data/attributes.txt -o $results/EnrichmentAnalysis/eIF5A  -U codon

EnrichmentAnalysis --ctrl $results/EnrichmentAnalysis/eIF5A_si-Ctrl-2_cds_codon_density.txt --treat $results/EnrichmentAnalysis/eIF5A_si-eIF5A-2_cds_codon_density.txt -c ../results/Ref/longest.transcripts.info.extended.txt -o $results/EnrichmentAnalysis/eIF5A -U codon -M RPKM -l 150 -n 10 -m 1 -e 30 --CI 0.95 -u 0 -d 500





