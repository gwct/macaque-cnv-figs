#!/usr/bin/python3
############################################################
# For macaque CNV gene overlaps, 09.2020
# This script gets the exons from the GTF files.
# Necessary because the field that contains the exon id
# is inconsistent.
############################################################
import sys, re

spec = "human";
if spec not in ["human", "macaque"]:
    sys.exit("ERROR: spec must be human or macaque.")

inp="gtf/" + spec + "-chromes.gtf"
out="bed/" + spec + "-exons.bed"

with open(out, "w") as o:
    for line in open(inp):
        if line[0] == "#":
            continue;
        ll = line.strip().split("\t")
        if ll[2] == "exon":
            if spec == "human":
                outline = ["chr" + ll[0], ll[3], ll[4], re.findall("ENSE[\d]+", line)[0]]
            elif spec == "macaque":
                outline = ["chr" + ll[0], ll[3], ll[4], re.findall("ENSMMUE[\d]+", line)[0]]
            # Uses regex to find the exon id depending on the Ensembl species identifier

            o.write("\t".join(outline) + "\n")