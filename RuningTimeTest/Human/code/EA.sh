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

RiboDensityAtEachPosition -c ../results/Ref/longest.transcripts.info.txt -f ../data/attributes_P.txt -o $results/EnrichmentAnalysis/test  -U codon

#enrichmentMeanDensity -i $results/EnrichmentAnalysis/test_Ctrl-1_cds_codon_density.txt,$results/EnrichmentAnalysis/test_Ctrl-2_cds_codon_density.txt -o $results/EnrichmentAnalysis/Ctrl
#enrichmentMeanDensity -i $results/EnrichmentAnalysis/test_ASNase-1_cds_codon_density.txt,$results/EnrichmentAnalysis/test_ASNase-2_cds_codon_density.txt -o $results/EnrichmentAnalysis/ASNase


#EnrichmentAnalysis --ctrl $results/EnrichmentAnalysis/Ctrl_mean_density.txt --treat $results/EnrichmentAnalysis/ASNase_mean_density.txt -c ../results/Ref/longest.transcripts.info.txt -o $results/EnrichmentAnalysis/test -U codon -M RPKM -l 150 -n 10 -m 1 -e 30 --CI 0.95 -u 0 -d 500

EnrichmentAnalysis --ctrl $results/EnrichmentAnalysis/test_Ctrl-1_cds_codon_density.txt --treat $results/EnrichmentAnalysis/test_ASNase-1_cds_codon_density.txt -c ../results/Ref/longest.transcripts.info.txt -o $results/EnrichmentAnalysis/test -U codon -M RPKM -l 150 -n 10 -m 1 -e 30 --CI 0.95 -u 0 -d 500





