#########################################################################
# File Name: FA.sh
# Author: lifj
# mail: lfj17@mails.tsinghua.edu.cn
# Created Time: Tue 09 Jun 2020 02:47:59 PM CST
#########################################################################
#!/bin/bash
#BSUB -J FA.sh
#BSUB -o FA.out
#BSUB -e FA.err
#BSUB -n 1
#BSUB -q TEST-A



results=../results

mkdir -p $results/FeatureAnalysis

RiboDensityAtEachKindAAOrCodon -f ../data/attributes.txt -c ../results/Ref/longest.transcripts.info.txt -o ../results/FeatureAnalysis/up -M RPKM -S ../results/up_common_trans.txt -l 100 -n 10 --table 1 -F ../results/Ref/longest_cds_sequences.fa
RiboDensityAtEachKindAAOrCodon -f ../data/attributes.txt -c ../results/Ref/longest.transcripts.info.txt -o ../results/FeatureAnalysis/down -M RPKM -S ../results/down_common_trans.txt -l 100 -n 10 --table 1 -F ../results/Ref/longest_cds_sequences.fa
RiboDensityAtEachKindAAOrCodon -f ../data/attributes.txt -c ../results/Ref/longest.transcripts.info.txt -o ../results/FeatureAnalysis/unblocked -M RPKM -S ../results/unblocked_common_trans.txt -l 100 -n 10 --table 1 -F ../results/Ref/longest_cds_sequences.fa
PlotRiboDensityAtEachKindAAOrCodon -i ../results/FeatureAnalysis/up_all_codon_density.txt -o ../results/FeatureAnalysis/up -g Ctrl,Severe-HS-2hrs -r Ctrl-2__Severe-HS-2hrs-2  --level codon
PlotRiboDensityAtEachKindAAOrCodon -i ../results/FeatureAnalysis/up_all_codon_density.txt -o ../results/FeatureAnalysis/up -g Ctrl,Severe-HS-2hrs -r Ctrl-2__Severe-HS-2hrs-2  --level AA
PlotRiboDensityAtEachKindAAOrCodon -i ../results/FeatureAnalysis/down_all_codon_density.txt -o ../results/FeatureAnalysis/down -g Ctrl,Severe-HS-2hrs -r Ctrl-2__Severe-HS-2hrs-2  --level codon
PlotRiboDensityAtEachKindAAOrCodon -i ../results/FeatureAnalysis/down_all_codon_density.txt -o ../results/FeatureAnalysis/down -g Ctrl,Severe-HS-2hrs -r Ctrl-2__Severe-HS-2hrs-2  --level AA

PausingScore -f ../data/attributes.txt -c ../results/Ref/longest.transcripts.info.txt -o ../results/FeatureAnalysis/up -M RPKM -S ../results/up_common_trans.txt -l 100 -n 10 --table 1 -F ../results/Ref/longest_cds_sequences.fa
ProcessPausingScore -i ../results/FeatureAnalysis/up_Ctrl-2_pausing_score.txt,../results/FeatureAnalysis/up_Severe-HS-2hrs-2_pausing_score.txt -o ../results/FeatureAnalysis/up -g Ctrl,Severe-HS-2hrs -r Ctrl-2__Severe-HS-2hrs-2  --mode raw --ratio_filter 1.5 --pausing_score_filter 0.5

conda activate python27
./seq2logo-2.1/Seq2Logo.py -f ../results/FeatureAnalysis/up_pwm.txt  -u probability -I 5 -o ../results/FeatureAnalysis/up --format PDF
conda deactivate

RiboDensityAroundTripleteAAMotifs -f ../data/attributes.txt -c ../results/Ref/longest.transcripts.info.txt -o ../results/FeatureAnalysis/up_DD -M RPKM -S ../results/up_common_trans.txt -l 100 -n 10 --table 1 -F ../results/Ref/longest_cds_sequences.fa --type2 DDD --type1 DD
RiboDensityAroundTripleteAAMotifs -f ../data/attributes.txt -c ../results/Ref/longest.transcripts.info.txt -o ../results/FeatureAnalysis/up_PP -M RPKM -S ../results/up_common_trans.txt -l 100 -n 10 --table 1 -F ../results/Ref/longest_cds_sequences.fa --type2 PPP --type1 PP
RiboDensityAroundTripleteAAMotifs -f ../data/attributes.txt -c ../results/Ref/longest.transcripts.info.txt -o ../results/FeatureAnalysis/up_KK -M RPKM -S ../results/up_common_trans.txt -l 100 -n 10 --table 1 -F ../results/Ref/longest_cds_sequences.fa --type2 KKK --type1 KK

PlotRiboDensityAroundTriAAMotifs -i ../results/FeatureAnalysis/up_PP_motifDensity_dataframe.txt -o ../results/FeatureAnalysis/up_PPP -g Ctrl,Severe-HS-2hrs -r Ctrl-2__Severe-HS-2hrs-2  --mode mean --ymax 0.2
PlotRiboDensityAroundTriAAMotifs -i ../results/FeatureAnalysis/up_DD_motifDensity_dataframe.txt -o ../results/FeatureAnalysis/up_DDD -g Ctrl,Severe-HS-2hrs -r Ctrl-2__Severe-HS-2hrs-2  --mode mean --ymax 0.2  
PlotRiboDensityAroundTriAAMotifs -i ../results/FeatureAnalysis/up_KK_motifDensity_dataframe.txt -o ../results/FeatureAnalysis/up_KKK -g Ctrl,Severe-HS-2hrs -r Ctrl-2__Severe-HS-2hrs-2   --mode mean --ymax 0.2

bash GetCdsSeqForGeneSets.sh
tAI -i ../results/Ref/up_cds_sequences.fa,../results/Ref/down_cds_sequences.fa,../results/Ref/unblocked_cds_sequences.fa -N ../data/tRNA_GCNs_mouse.txt -o ../results/FeatureAnalysis/test -u 0 -d 500 -t up-regulated-genes,down-regulated-genes,unblocked-genes
cAI -i ../results/Ref/up_cds_sequences.fa,../results/Ref/down_cds_sequences.fa,../results/Ref/unblocked_cds_sequences.fa -o  ../results/FeatureAnalysis/test -u 0 -d 500 -t up-regulated-genes,down-regulated-genes,unblocked-genes --reference ../results/Ref/longest_cds_sequences.fa

hydropathyCharge -i ../results/Ref/up_cds_sequences.fa,../results/Ref/down_cds_sequences.fa,../results/Ref/unblocked_cds_sequences.fa -t up-regulated-genes,down-regulated-genes,unblocked-genes -o ../results/FeatureAnalysis/test_hydropathy -u 0 -d 500 --index ../data/hydropathy.txt
hydropathyCharge -i ../results/Ref/up_cds_sequences.fa,../results/Ref/down_cds_sequences.fa,../results/Ref/unblocked_cds_sequences.fa -t up-regulated-genes,down-regulated-genes,unblocked-genes -o ../results/FeatureAnalysis/test_Charge -u 0 -d 500 --index ../data/AA_charge.txt

PlotHydropathyCharge -i ../results/FeatureAnalysis/test_Charge_values_dataframe.txt -o ../results/FeatureAnalysis/test_charge -u 0 -d 500 --mode all --ylab "Average Charges"
PlotHydropathyCharge -i ../results/FeatureAnalysis/test_hydropathy_values_dataframe.txt -o ../results/FeatureAnalysis/test_hydro -u 0 -d 500 --mode all --ylab "Average Hydrophobicity"
cAIPlot -i ../results/FeatureAnalysis/test_local_cAI_dataframe.txt -o ../results/FeatureAnalysis/test_cAI -u 0 -d 500 --mode all --start 5 --window 7 --step 1
tAIPlot -i ../results/FeatureAnalysis/test_tAI_dataframe.txt -o ../results/FeatureAnalysis/test_tAI -u 0 -d 500 --mode all --start 5 --window 7 --step 1
