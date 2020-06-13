#########################################################################
# File Name: QC.sh
# Author: lifj
# mail: lfj17@mails.tsinghua.edu.cn
# Created Time: Tue 09 Jun 2020 02:47:28 PM CST
#########################################################################
#!/bin/bash
#BSUB -J QC.sh
#BSUB -o QC.out
#BSUB -e QC.err
#BSUB -n 1
#BSUB -q TEST-A



## periodicity
results=../results
mkdir -p $results/periodicity
cd $results/periodicity
i=SRR3297803
Periodicity -a ../../results/Ref/RiboCode_annot -i ../../data/00.STAR/${i}_STAR/$i.Aligned.toTranscriptome.out.sorted.bam -o test_${i} -c ../../results/Ref/longest.transcripts.info.txt -L 25 -R 35
cd ../../code
echo "periodicity finished!"

## length
mkdir -p $results/length
cd $results/length
i=SRR3297803
LengthDistribution -i ../../data/00.STAR/${i}_STAR/${i}.Aligned.sortedByCoord.out.bam -o test_$i -f bam
cd ../../code
echo "length finished!"

## RPF statistics
mkdir -p $results/RPF_statisitcs
cd $results/RPF_statisitcs
i=SRR3297803
StatisticReadsOnDNAsContam -i ../../data/00.STAR/${i}_STAR/${i}.Aligned.sortedByCoord.out.bam -g ../../results/Ref/Drosophila_melanogaster.BDGP6.22.98.gtf -o test_${i}
cd ../../code
echo "RPF statistics finished!"

## RiboDensityOfDiffFrames
mkdir -p $results/RiboDensityFromDiffFrames
RiboDensityOfDiffFrames -f ../data/attributes.txt -c ../results/Ref/longest.transcripts.info.txt -o $results/RiboDensityFromDiffFrames/test --plot yes
echo "RiboDensityFromDiffFrames.sh finished!"


#bash periodicity.sh
#echo "periodicity.sh finished!"
#bash length.sh
#echo "length.sh finished!"
#bash RPF_statistics.sh
#echo "RPF_statisitcs.sh finished!"
#bash RiboDensityFromDiffFrames.sh
#echo "RiboDensityFromDiffFrames.sh finished!"
