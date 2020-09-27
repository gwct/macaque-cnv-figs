#!/bin/bash
############################################################
# For macaque CNV gene overlaps, 09.2020
# Uses bedtools to get overlaps between genic regions and
# CNVs.
############################################################

echo "01 MACAQUE GENES"
bedtools intersect -a bed/macaque-cnvs-filtered-noalu.bed -b bed/macaque-genes.bed -wb | bedtools overlap -i stdin -cols 2,3,6,7 > bed/macaque-cnvs-to-genes.bed

echo "02 MACAQUE GENES 10KB UPSTREAM"
bedtools intersect -a bed/macaque-cnvs-filtered-noalu.bed -b bed/macaque-genes-10kb-up.bed -wb | bedtools overlap -i stdin -cols 2,3,6,7 > bed/macaque-cnvs-to-genes-10kb-up.bed

echo "03 MACAQUE GENES 10KB DOWNSTREAM"
bedtools intersect -a bed/macaque-cnvs-filtered-noalu.bed -b bed/macaque-genes-10kb-down.bed -wb | bedtools overlap -i stdin -cols 2,3,6,7 > bed/macaque-cnvs-to-genes-10kb-down.bed

echo "04 MACAQUE TRANSCRIPTS"
bedtools intersect -a bed/macaque-cnvs-filtered-noalu.bed -b bed/macaque-transcripts.bed -wb | bedtools overlap -i stdin -cols 2,3,6,7 > bed/macaque-cnvs-to-transcripts.bed

echo "05 MACAQUE EXONS"
bedtools intersect -a bed/macaque-cnvs-filtered-noalu.bed -b bed/macaque-exons.bed -wb | bedtools overlap -i stdin -cols 2,3,6,7 > bed/macaque-cnvs-to-exons.bed


echo "06 HUMAN GENES"
bedtools intersect -a bed/human-cnvs-filtered-noalu.bed -b bed/human-genes.bed -wb | bedtools overlap -i stdin -cols 2,3,6,7 > bed/human-cnvs-to-genes.bed

echo "07 HUMAN GENES GENES 10KB UPSTREAM"
bedtools intersect -a bed/human-cnvs-filtered-noalu.bed -b bed/human-genes-10kb-up.bed -wb | bedtools overlap -i stdin -cols 2,3,6,7 > bed/human-cnvs-to-genes-10kb-up.bed

echo "08 HUMAN GENES 10KB DOWNSTREAM"
bedtools intersect -a bed/human-cnvs-filtered-noalu.bed -b bed/human-genes-10kb-down.bed -wb | bedtools overlap -i stdin -cols 2,3,6,7 > bed/human-cnvs-to-genes-10kb-down.bed

echo "09 HUMAN TRANSCRIPTS"
bedtools intersect -a bed/human-cnvs-filtered-noalu.bed -b bed/human-transcripts.bed -wb | bedtools overlap -i stdin -cols 2,3,6,7 > bed/human-cnvs-to-transcripts.bed

echo "10 HUMAN EXONS"
bedtools intersect -a bed/human-cnvs-filtered-noalu.bed -b bed/human-exons.bed -wb | bedtools overlap -i stdin -cols 2,3,6,7 > bed/human-cnvs-to-exons.bed