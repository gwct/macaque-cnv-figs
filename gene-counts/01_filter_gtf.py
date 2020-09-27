#!/usr/bin/python3
############################################################
# For macaque CNV gene overlaps, 09.2020
# This script simply takes the raw GTF file and filters out
# all annotations not on one of the assembled chromosomes.
############################################################
import sys

spec = "macaque";
if spec not in ["human", "macaque"]:
    sys.exit("ERROR: spec must be human or macaque.")

if spec == "macaque":
    gtfin = "gtf/Macaca_mulatta.Mmul_8.0.1.97.gtf";
    gtfout = "gtf/macaque-chromes.gtf";
if spec == "human":
    gtfin = "gtf/Homo_sapiens.GRCh38.84.gtf";
    gtfout = "gtf/human-chromes.gtf"; 

chromes = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"];

with open(gtfout, "w") as out:
    for line in open(gtfin):
        if line[0] == "#":
            continue;
        # Skips the header lines

        c = line.strip("\t")[0];
        if c in chromes:
            out.write(line);
        # Writes the line to output only if it is in the chromes list.

