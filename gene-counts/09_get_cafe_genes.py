#!/usr/bin/python3
############################################################
# For macaque revisions, 09.2020
# Gets proteins from the CAFE analysis and checks if they
# are overlapped by CNVs.
############################################################

import sys, re
from collections import defaultdict
sys.path.append("../lib/");
import mqcore as MQ

############################################################

if "-h" in sys.argv:
    sys.exit("Usage: python gene_count.py <species: macaque or human (required)>");

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
    MQ.runTime("# Macaque CNV and CAFE gene counter\n# Species:   " + species, logfile);

    dumpfile = "../cafe-data/dump.out.mlemur-blast.txt.I30";
    overlap_file = "bed/" + species + "-cnvs-filtered-noalu.csv";
    gtf_file = "gtf/" + species + "-chromes.gtf";
    overlap_file = "bed/" + species + "-cnvs-to-transcripts.bed";
    outfilename = "../cafe-data/" + species + "-cafe-genes" + ".csv";
    # File names.

    MQ.PWS("# Log file:    " + logfilename, logfile);
    MQ.PWS("# GTF file:    " + gtf_file, logfile);
    MQ.PWS("# Overlap file:" + overlap_file, logfile);
    MQ.PWS("# Dump file:   " + dumpfile, logfile);
    MQ.PWS("# Output file: " + outfilename, logfile);
    MQ.PWS("# ----------------", logfile);
    # I/O options and info.

    MQ.PWS("# " + MQ.getDateTime() + " Reading features from GTF...", logfile);
    ttp = defaultdict(str);
    overlaps = {};
    for line in open(gtf_file):
        if line[0] == "#":
            continue;
        line = line.strip().split("\t");
        feature_type, chrome, start, end, feature_info = line[2], line[0], int(line[3]), int(line[4]), line[8];
        if feature_type == "CDS":
            tid = re.findall('ENS' + regstr + 'T[\d]+', feature_info)[0];
            pid = re.findall('ENS' + regstr + 'P[\d]+', feature_info)[0];

            if pid not in ttp[tid]:
                ttp[tid] = pid;
                overlaps[tid] = { 'dup' : 0, 'del' : 0 };
    MQ.PWS("# Transcripts/proteins read:    " + str(len(ttp)), logfile);
    MQ.PWS("# ----------------", logfile);
    # Read features from GTF file in to memory. Necessary to associate transcript IDs with protein IDs.

    MQ.PWS("# " + MQ.getDateTime() + " Reading CAFE proteins...", logfile);
    cafe_proteins = [];
    for line in open(dumpfile):
        line = line.strip().split("\t");
        for pid in line:
            if pid.startswith(cafestr):
                pid = pid[:pid.index(".")].replace(cafestr + "_", "");
                cafe_proteins.append(pid);
    MQ.PWS("# Proteins read:    " + str(len(cafe_proteins)), logfile);
    MQ.PWS("# ----------------", logfile);
    # Read the proteins from the CAFE dataset to restrict counting overlaps to them.

    MQ.PWS("# " + MQ.getDateTime() + " Getting overlaps...", logfile);
    transcripts_skipped = 0;
    for line in open(overlap_file):
        line = line.strip().split("\t");
        
        tid = line[7];
        if tid not in overlaps:
            transcripts_skipped += 1;
            continue;

        pid = ttp[tid];

        if "<DEL>" in line[3]:
            overlaps[tid]['del'] += 1;
        elif "<DUP>" in line[3]:
            overlaps[tid]['dup'] += 1;
    MQ.PWS("# Transcripts skipped (no associated protein): " + str(transcripts_skipped), logfile);
    MQ.PWS("# ----------------", logfile);
    # Count the number of overlaps per transcript for transcripts that are included in the CAFE analysis.

    MQ.PWS("# " + MQ.getDateTime() + " Writing transcripts in CAFE dataset overlapped by CNVs...", logfile);
    transcripts_written = 0;
    with open(outfilename, "w") as outfile:
        headers = "Transcript,Protein,Num del,Num dup";
        outfile.write(headers + "\n");
        for tid in overlaps:
            pid = ttp[tid];
            if pid not in cafe_proteins:
                continue;
            # Skip if not in the set of CAFE proteins.

            outline = [tid, pid, str(overlaps[tid]['del']), str(overlaps[tid]['dup'])];
            outfile.write(",".join(outline) + "\n");
            transcripts_written += 1;            
    MQ.PWS("# Transcripts written: " + str(transcripts_written), logfile);
    MQ.PWS("# ----------------", logfile);
    # Write the counts to a file.