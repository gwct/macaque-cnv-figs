#!/usr/bin/python3
############################################################
# For macaque CNV gene overlaps, 09.2020
# This script gets the number of CNVs per chromosome.
############################################################
import sys

inp="macaque-cnvs-filtered.csv"
out="macaque-cnv-chrome-counts.csv"
ind="../gene-counts/indices/Macaca_mulatta.Mmul_8.0.1.97.fai";
# The index files contain chromosome lengths.

chromes = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","X","Y"];
sizes = {};
for line in open(ind):
    line = line.strip().split("\t");
    if line[0] not in chromes:
        continue;
    c = "chr" + line[0];
    s = line[1];
    sizes[c] = s
# Read the index files to get the chromosome lengths to correct for any regions
# that may exceed them.

counts = { "chr" + c : { 'del' : 0, 'dup' : 0, 'len' : sizes["chr"+c], 'lab' : c } for c in chromes };
first = True;
for line in open(inp):
    if first:
        first = False;
        continue;
    line = line.strip().split(",");
    cnv_type = line[1].replace("<","").replace(">","").lower();
    c = line[2];
    counts[c][cnv_type] += 1;

with open(out, "w") as outfile:
    headers = "Chromosome,Length,Num dels,Num dups, Num CNVs,label";
    outfile.write(headers + "\n");
    for c in counts:
        outline = [c, counts[c]['len'], str(counts[c]['del']), str(counts[c]['dup']), str(counts[c]['del'] + counts[c]['dup']), str(counts[c]['lab'])];
        outfile.write(",".join(outline) + "\n");