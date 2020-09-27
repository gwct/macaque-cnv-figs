#!/usr/bin/python3
############################################################
# For macaque revisions, 04.2020
# Gets proteins from the CAFE analysis and checks if they
# are overlapped by CNVs
############################################################

import sys, core, re
from collections import defaultdict

############################################################

if "-h" in sys.argv:
    sys.exit("Usage: python gene_count.py <species: macaque or human (required)> <max CNV length to annotate (optional)>");

if len(sys.argv) < 2:
    sys.exit(" * ERROR: Species must be provided: macaque or human");
species = sys.argv[1];
if sys.argv[1] not in ["macaque", "human"]:
    sys.exit(" * ERROR: Species must be provided: macaque or human");

if species == "human":
    regstr = "";
    cafestr = "Hsapi";
elif species == "macaque":
    regstr = "MMU"
    cafestr = "Mmula";
# Species string checking.

logfilename = "get-cafe-genes-" + species + ".log";
#sys.exit(logfilename);

with open(logfilename, "w") as logfile:
    core.runTime("# Macaque CNV and CAFE gene counter\n# Species:   " + species, logfile);

    cnvfile = "../data/sv-calls/" + species + "-cnvs-filtered.csv";
    dumpfile = "../data/dump.out.mlemur-blast.txt.I30";
    if species == "macaque":
        gtffile = "../data/ens-gene-tables/Macaca_mulatta.Mmul_8.0.1.97.gtf";   
    elif species == "human":
        gtffile = "../data/ens-gene-tables/Homo_sapiens.GRCh38.84.gtf";
    overlap_file = "../data/sv-calls/" + species + "-gene-counts-features.csv";
    outfilename = "../data/sv-calls/" + species + "-cafe-genes" + ".csv";
    # File names.

    core.PWS("# Log file:    " + logfilename, logfile);
    core.PWS("# GTF file:    " + gtffile, logfile);
    core.PWS("# Overlap file:" + overlap_file, logfile);
    core.PWS("# Dump file:   " + dumpfile, logfile);
    core.PWS("# Output file: " + outfilename, logfile);
    # I/O options and info.
    core.PWS("# ----------------", logfile);

    core.PWS("# " + core.getDateTime() + " Reading features from GTF...", logfile);
    ttp = defaultdict(str);
    for line in open(gtffile):
        if line[0] == "#":
            continue;
        line = line.strip().split("\t");
        feature_type, chrome, start, end, feature_info = line[2], line[0], int(line[3]), int(line[4]), line[8];
        if feature_type == "CDS":
            tid = re.findall('ENS' + regstr + 'T[\d]+', feature_info)[0];
            pid = re.findall('ENS' + regstr + 'P[\d]+', feature_info)[0];

            if pid not in ttp[tid]:
                ttp[tid] = pid;
    # Read features from GTF file in to memory
    core.PWS("# Transcripts/proteins read:    " + str(len(ttp)), logfile);
    core.PWS("# ----------------", logfile);

    core.PWS("# " + core.getDateTime() + " Reading CAFE proteins...", logfile);
    cafe_proteins = [];
    for line in open(dumpfile):
        line = line.strip().split("\t");
        for pid in line:
            if pid.startswith(cafestr):
                pid = pid[:pid.index(".")].replace(cafestr + "_", "");
                cafe_proteins.append(pid);
    core.PWS("# Proteins read:    " + str(len(cafe_proteins)), logfile);
    core.PWS("# ----------------", logfile);    

    core.PWS("# " + core.getDateTime() + " Writing transcripts in CAFE dataset overlapped by CNVs...", logfile);
    transcripts_written = 0;
    with open(outfilename, "w") as outfile:
        headers = "Transcript,Protein,Num del,Num dup";
        outfile.write(headers + "\n");
        cnv_overlaps = {};
        first = True;
        for line in open(overlap_file):
            if first:
                first = False;
                continue;

            line = line.split(",");
            if line[1] != "transcript":
                continue;

            tid = line[0];
            pid = ttp[tid];
            if pid not in cafe_proteins:
                continue;

            outline = [tid, pid, line[10], line[12]];
            outfile.write(",".join(outline) + "\n");
            transcripts_written += 1;

    core.PWS("# Transcripts written: " + str(transcripts_written), logfile);
    core.PWS("# ----------------", logfile);