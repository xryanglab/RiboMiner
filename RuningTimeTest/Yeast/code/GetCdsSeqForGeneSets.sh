#########################################################################
# File Name: GetCdsSeqForGeneSets.sh
# Author: lifj
# mail: lfj17@mails.tsinghua.edu.cn
# Created Time: Fri 01 May 2020 06:26:21 PM CST
#########################################################################
#!/bin/bash


GetProteinCodingSequence -i ../results/Ref/RiboCode_annot/transcripts_sequence.fa  -c ../results/Ref/longest.transcripts.info.extended.txt -o ../results/Ref/up  --mode whole --table 1  -S ../results/up_common_trans.txt
GetProteinCodingSequence -i ../results/Ref/RiboCode_annot/transcripts_sequence.fa  -c ../results/Ref/longest.transcripts.info.extended.txt -o ../results/Ref/down  --mode whole --table 1  -S ../results/down_common_trans.txt
GetProteinCodingSequence -i ../results/Ref/RiboCode_annot/transcripts_sequence.fa  -c ../results/Ref/longest.transcripts.info.extended.txt -o ../results/Ref/unblocked  --mode whole --table 1  -S ../results/unblocked_common_trans.txt



