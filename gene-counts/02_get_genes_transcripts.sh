#!/bin/bash
############################################################
# For macaque CNV gene overlaps, 09.2020
# Uses awk to convert certain entries in GTF files to BED format.
############################################################


echo "01 MACAQUE GENES: bed/macaque-genes.bed"
awk '{if($3=="gene"){print "chr"$1"\t"$4"\t"$5"\t"$10}}' gtf/macaque-chromes.gtf | sort | sed 's/"//g; s/;//g' > bed/macaque-genes.bed
# Get macaque genes in BED format.

echo "02 MACAQUE TRANSCRIPTS: bed/macaque-transcripts.bed"
awk '{if($3=="transcript"){print "chr"$1"\t"$4"\t"$5"\t"$14}}' gtf/macaque-chromes.gtf | sort | sed 's/"//g; s/;//g' > bed/macaque-transcripts.bed
# Get macaque transcripts in BED format.

echo "03 HUMAN GENES: bed/humab-genes.bed"
awk '{if($3=="gene"){print "chr"$1"\t"$4"\t"$5"\t"$10}}' gtf/human-chromes.gtf | sort | sed 's/"//g; s/;//g' > bed/human-genes.bed
# Get human genes in BED format.

echo "04 HUMAN TRANSCRIPTS: bed/human-transcripts.bed"
awk '{if($3=="transcript"){print "chr"$1"\t"$4"\t"$5"\t"$14}}' gtf/human-chromes.gtf | sort | sed 's/"//g; s/;//g' > bed/human-transcripts.bed
# Get human transcripts in BED format.