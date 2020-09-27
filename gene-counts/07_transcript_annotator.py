#!/usr/bin/python3
############################################################
# For macaque revisions, 09.2020
# Given a file with transcript IDs that overlap CNVs (from
# 06_gene_count_bed.py) this script retrieves all GO terms
# for all transcripts in the input file.
# Also writes all GO terms not associated with input transcripts
# to another file as the background file for a Fisher's test.
############################################################

import sys, gzip
sys.path.append("../lib/");
import mqcore as MQ

############################################################

if len(sys.argv) < 3:
    sys.exit(" * ERROR: Species must be provided: macaque or human\n *        Mode must be provided:    dup or del");
species = sys.argv[1];
if species not in ["macaque", "human"]:
    sys.exit(" * ERROR: Species must be provided: macaque or human");
mode = sys.argv[2];
if mode not in ["dup", "del"]:
    sys.exit(" * ERROR: Mode must be provided: dup or del");
# Species string checking.

MQ.runTime("# Macaque transcript ID annotator");

gofile = "go/" + species + "-go-terms-uniq.txt";
transcripts_file = "go/" + species + "-noalu-transcripts-" + mode + ".txt";
# The input files: A GO term database from Ensembl (with identical lines removed with bash's sort | uniq commands)
# and a file containing transcript IDs that overlap CNVs (from 06_gene_count_bed.py)

queryoutfile = "go/" + species + "-cnv-" + mode + "-go-query.tab";
bgoutfile = "go/" + species + "-cnv-" + mode + "-go-bg.tab";

MQ.PWS("# Go file          : " + gofile);
MQ.PWS("# Transcripts file : " + transcripts_file);
MQ.PWS("# Mode             : " + mode);
MQ.PWS("# Query out        : " + queryoutfile);
MQ.PWS("# BG out           : " + bgoutfile);
MQ.PWS("# ----------------");

MQ.PWS( "# " + MQ.getDateTime() + " Getting annotation info...");
transcript_go = {};
go_accs = {};
first = True;
for line in gzip.open(gofile):
    if first:
        first = False;
        continue;
    
    line = line.decode().strip().split("\t");
    #print(line);

    chrome = line[4];
    if chrome == "MT":
        continue;

    tid = line[1];

    if tid not in transcript_go:
        transcript_go[tid] = [];

    if tid not in go_accs:
        go_accs[tid] = [];

    if len(line) > 7:
        go_acc, go_name, go_def, go_dom = "", "", "", "";
        go_acc = line[7];
        if go_acc not in go_accs[tid]:
            go_accs[tid].append(go_acc);
            if len(line) > 8:
                go_name, go_def = line[8], line[9] 
                if len(line) > 10:
                    go_dom = line[10];
            # go = "|".join([go_acc,go_name,go_def,go_dom]);
            go = [go_acc,go_name,go_def,go_dom]
            transcript_go[tid].append(go);
# Get GO terms for transcripts

lines_written, tids_not_found = 0,0;
MQ.PWS( "# " + MQ.getDateTime() + " Annotating transcripts...");
all_overlap_tids = [];
with open(queryoutfile, "w") as queryout:
    first = True;
    cnvs = {};
    i = 0;
    for line in open(transcripts_file):
        tid = line.strip();
        if tid in transcript_go:
            for go in transcript_go[tid]:
                outline = [tid, go[0], go[1], go[2], go[3]];
                queryout.write("\t".join(outline) + "\n");
        else:
            tids_not_found += 1;

MQ.PWS( "# Transcripts not found: " + str(tids_not_found));
MQ.PWS( "# " + MQ.getDateTime() + " Writing background annotations...");
with open(bgoutfile, "w") as bgout:
    for tid in transcript_go:
        for go in transcript_go[tid]:
            if tid in all_overlap_tids:
                continue;
            outline = [tid, go[0], go[1], go[2], go[3]];
            bgout.write("\t".join(outline) + "\n");

