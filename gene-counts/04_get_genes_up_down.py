
#!/usr/bin/python3
############################################################
# For macaque CNV gene overlaps, 09.2020
# This script gets the regions 10kb up and downstream from
# genes. Corrects for chromosome length.
############################################################
import sys

spec = "human";
if spec not in ["human", "macaque"]:
    sys.exit("ERROR: spec must be human or macaque.")

inp="bed/" + spec + "-genes.bed"
out_up="bed/" + spec + "-genes-10kb-up.bed"
out_dn="bed/" + spec + "-genes-10kb-down.bed"

if spec == "macaque":
    ind="indices/Macaca_mulatta.Mmul_8.0.1.97.fai";
elif spec == "human":
    ind="indices/Homo_sapiens.GRCh38.84.fai";
# The index files contain chromosome lengths.

sizes = {};
for line in open(ind):
    line = line.strip().split("\t");
    c = "chr" + line[0];
    s = int(line[1]);
    sizes[c] = s
# Read the index files to get the chromosome lengths to correct for any regions
# that may exceed them.

with open(out_up, "w") as ou, open(out_dn, "w") as od:
    for line in open(inp):
        line = line.strip().split("\t");
        c, start, end = line[0], int(line[1]), int(line[2]);

        up_start = start - 10000;
        if up_start < 1:
            up_start = 1;
        # For upstream regions, if the start position is negative (i.e. the gene is within
        # 10kb of the start of the chromosome) set the start position to 1.

        down_end = end + 10000;
        if down_end > sizes[c]:
            down_end = sizes[c];
        # For downstream regions, if the end position exceeds the length of the chromosome,
        # set the end position to the length of the chromosome.

        up_line = [c, str(up_start), str(start), line[3]];
        down_line = [c, str(end), str(down_end), line[3]];

        ou.write("\t".join(up_line) + "\n");
        od.write("\t".join(down_line) + "\n");
