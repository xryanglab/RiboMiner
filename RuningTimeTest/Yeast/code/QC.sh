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
i=SRR5008135
Periodicity -a ../../results/Ref/RiboCode_annot -i ../../data/00.STAR/${i}_STAR/$i.Aligned.toTranscriptome.out.sorted.bam -o eIF5A_${i} -c ../../results/Ref/longest.transcripts.info.extended.txt -L 25 -R 35
cd ../../code
echo "periodicity finished!"

## length
mkdir -p $results/length
cd $results/length
i=SRR5008135
LengthDistribution -i ../../data/noncontam_${i}.fastq -o eIF5A_$i -f fastq
cd ../../code
echo "length finished!"

## RPF statistics
mkdir -p $results/RPF_statisitcs
cd $results/RPF_statisitcs
i=SRR5008135
StatisticReadsOnDNAsContam -i ../../data/00.STAR/${i}_STAR/${i}.Aligned.sortedByCoord.out.bam -g ../../results/Ref/Saccharomyces_extend.gtf -o eIF5A_${i}
cd ../../code
echo "RPF statistics finished!"

## RiboDensityOfDiffFrames
mkdir -p $results/RiboDensityFromDiffFrames
RiboDensityOfDiffFrames -f ../data/attributes.txt -c ../results/Ref/longest.transcripts.info.extended.txt -o $results/RiboDensityFromDiffFrames/eIF5A --plot yes
echo "RiboDensityFromDiffFrames.sh finished!"


#bash periodicity.sh
#echo "periodicity.sh finished!"
#bash length.sh
#echo "length.sh finished!"
#bash RPF_statistics.sh
#echo "RPF_statisitcs.sh finished!"
#bash RiboDensityFromDiffFrames.sh
#echo "RiboDensityFromDiffFrames.sh finished!"
