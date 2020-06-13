#########################################################################
# File Name: prepare_references.sh
# Author: lifj
# mail: lfj17@mails.tsinghua.edu.cn
# Created Time: Thu 30 Apr 2020 03:22:02 PM CST
#########################################################################
#!/bin/bash

workdir=.
genome_fasta=../data/Ref/Saccharomyces_genome.fa
gtf=../data/Ref/Saccharomyces.gtf

mkdir -p ../results/Ref
cp -r ../data/Ref/* ../results/Ref
#cp $genome_fasta ../results/Ref
#cp $gtf ../results/Ref
## extend gtf
#python ExtendGTF.py ../results/Ref/Saccharomyces.gtf ../results/Ref/Saccharomyces_extend.gtf

## run prepare_transcripts
#prepare_transcripts -g ../results/Ref/Saccharomyces_extend.gtf -f ../results/Ref/Saccharomyces_genome.fa -o ../results/Ref/RiboCode_annot

## modify transcript_cds.txt
#python ModifyTransCDS.py ../results/Ref/RiboCode_annot/transcripts_cds.txt ../results/Ref/RiboCode_annot/transcripts_cds_extended.txt

## output the longest transcript info
OutputTranscriptInfo -c ../results/Ref/RiboCode_annot/transcripts_cds_extended.txt -g ../results/Ref/Saccharomyces_extend.gtf -f ../results/Ref/RiboCode_annot/transcripts_sequence.fa -o ../results/Ref/longest.transcripts.info.extended.txt -O ../results/Ref/all.transcripts.info.extended.txt

## prepare fasta sequences
GetProteinCodingSequence -i ../results/Ref/RiboCode_annot/transcripts_sequence.fa  -c ../results/Ref/longest.transcripts.info.extended.txt -o ../results/Ref/longest --mode whole --table 1 

