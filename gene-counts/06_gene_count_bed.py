#!/usr/bin/python3
############################################################
# For macaque revisions, 09.2020
# Counts overlaps between CNVs and genic regions based on
# bedtools intersect files.
############################################################

import sys, datetime
from collections import defaultdict
sys.path.append("../lib/");
import mqcore as MQ

############################################################
# Functions
def countOverlaps(cnvs, bedfile, ol_type):
# Given a list of CNVs and a bedfile with overlaps, this function counts the number of times a CNV overlaps a region and vice versa
# and reports relevant counts.
    cnvs_genes = { cnv : { 'len' : int(cnv.split(":")[3]), 'part' : [], 'full' : [] } for cnv in cnvs };
    for line in open(bedfile):
        line = line.strip().split("\t");
        cnvs_start, cnv_end, cnv, gene_start, gene_end, gene, overlap = int(line[1]), int(line[2]), line[3], int(line[5]), int(line[6]), line[7], int(line[8]);
        if cnv not in cnvs:
            continue;
        gene_len = gene_end - gene_start;    

        if overlap == gene_len:
            cnvs_genes[cnv]['full'].append(gene);
        else:
            cnvs_genes[cnv]['part'].append(gene);
        # CNV gene overlap
    # Get the genes overlapped for each CNV from the bed file.

    num_full_overlapped, num_full_overlapped_cond, num_part_overlapped, num_part_overlapped_cond, num_overlapped, num_overlapped_cond = [], [], [], [], [], [];
    # Lists of COUNTS of overlaps for each type.

    del_part_genes,del_full_genes, dup_part_genes, dup_full_genes = [], [], [], [];
    # Lists of genes that have been overlapped.

    del_full, del_part, del_total, dup_full, dup_part, dup_total = 0,0,0,0,0,0;
    # Total counts.

    cnv_at_least_one = 0;
    del_at_least_one, del_at_least_one_full, del_at_least_one_part = 0,0,0;
    dup_at_least_one, dup_at_least_one_full, dup_at_least_one_part = 0,0,0;
    # Conditional counts

    gene_counts = {};
    del_genes = [];
    dup_genes = [];
    # To count the number of CNVs that overlap each gene region.

    for cnv in cnvs_genes:
        if "<DEL>" in cnv:
            cnv_type = "del";
        elif "<DUP>" in cnv:
            cnv_type = "dup";
        else:
            sys.exit("Invalid CNV type: " + cnv);
        # Get the CNV type for later classifcation.

        full_overlaps = len(cnvs_genes[cnv]['full']);
        part_overlaps = len(cnvs_genes[cnv]['part']);
        # Get the number of full and partial gene overlaps.

        cur_genes = cnvs_genes[cnv]['full'] + cnvs_genes[cnv]['part'];
        if cur_genes != []:
            for gene in cur_genes:
                if gene not in gene_counts:
                    gene_counts[gene] = { 'del' : [], 'dup' : [], 'count' : 0 };

                gene_counts[gene]['count'] += 1;

                if cnv_type == "del":
                    gene_counts[gene]['del'].append(cnv);
                
                if cnv_type == "dup":
                    gene_counts[gene]['dup'].append(cnv);
        # This block counts the number of overlaps by genic region (i.e. how many CNVs overlap this gene, rather than 
        # how many genes does this CNV overlap).

        num_full_overlapped.append(full_overlaps);
        num_part_overlapped.append(part_overlaps);
        num_overlapped.append(full_overlaps + part_overlaps);
        # Append the number of current overlaps to the appropriate list

        if cnv_type == "del":
            del_full += full_overlaps;
            del_part += part_overlaps;
            del_total += full_overlaps + part_overlaps;
            # Increment counts of each type of overlap.

            del_part_genes += cnvs_genes[cnv]['part'];
            del_full_genes += cnvs_genes[cnv]['full'];
            # Append the current deleted genes to the list of total genes that have been deleted.
        # For deletions

        if cnv_type == "dup":
            dup_full += full_overlaps;
            dup_part += part_overlaps;
            dup_total += full_overlaps + part_overlaps;
            # Increment counts of each type of overlap.

            dup_part_genes += cnvs_genes[cnv]['part'];
            dup_full_genes += cnvs_genes[cnv]['full'];
            # Append the current duplicated genes to the list of total genes that have been duplicated.
        # For duplications

        if full_overlaps + part_overlaps != 0:
        # If at least one gene has been overlapped
            cnv_at_least_one += 1;
            # Increment the conditional count.

            num_full_overlapped_cond.append(full_overlaps);
            num_part_overlapped_cond.append(part_overlaps);
            num_overlapped_cond.append(full_overlaps + part_overlaps);
            # Append the number of current overlaps to the appropriate list

            if cnv_type == "del":
                if full_overlaps != 0:
                    del_at_least_one_full += 1;
                # If there is at least one full deletion, count it.
                if part_overlaps != 0:
                    del_at_least_one_part += 1;
                # If there is at least one partial deletion, count it.
                del_at_least_one += 1;
                # Count that this is a deletion with at least one overlap.
            # For deletions

            if cnv_type == "dup":
                if full_overlaps != 0:
                    dup_at_least_one_full += 1;
                # If there is at least one full deletion, count it.
                if part_overlaps != 0:
                    dup_at_least_one_part += 1;
                # If there is at least one partial deletion, count it.
                dup_at_least_one += 1;
                # Count that this is a deletion with at least one overlap.
            # For duplications

    MQ.PWS("\n" + ol_type, logfile);
    MQ.PWS("AT LEAST ONE GENIC REGION: " + str(cnv_at_least_one), logfile);
    MQ.PWS("NUM REGIONS OVERLAPPED AT LEAST ONCE: " + str(len(gene_counts)), logfile);

    MQ.PWS("\n\tDel\tDup", logfile);
    MQ.PWS("Full\t" + str(len(set(del_full_genes))) + "\t" + str(len(set(dup_full_genes))), logfile);
    MQ.PWS("Part\t" + str(len(set(del_part_genes))) + "\t" + str(len(set(dup_part_genes))), logfile);
    # Table 1

    MQ.PWS("\n\tDel\tDup", logfile);
    MQ.PWS("Full\t" + str(del_at_least_one_full) + "\t" + str(dup_at_least_one_full), logfile);
    MQ.PWS("Part\t" + str(del_at_least_one_part) + "\t" + str(dup_at_least_one_part), logfile);
    # Table 2

    MQ.PWS("\nAll CNVs avg. overlap", logfile);
    MQ.PWS("\tDel\tDup", logfile);
    MQ.PWS("Full\t" + str(round(float(del_full) / float(num_cnvs), 3)) + "\t" + str(round(float(dup_full) / float(num_cnvs), 3)), logfile);
    MQ.PWS("Part\t" + str(round(float(del_part) / float(num_cnvs), 3)) + "\t" + str(round(float(dup_part) / float(num_cnvs), 3)), logfile);
    MQ.PWS("All\t" + str(round(float(del_full + del_part) / float(num_cnvs), 3)) + "\t" + str(round(float(dup_full + dup_part) / float(num_cnvs), 3)), logfile);
    # Table 3

    MQ.PWS("\nOverlapping CNVs avg. overlap", logfile);
    MQ.PWS("\tDel\tDup", logfile);
    MQ.PWS("Full\t" + str(round(float(del_full) / float(del_at_least_one_full), 3)) + "\t" + str(round(float(dup_full) / float(dup_at_least_one_full), 3)), logfile);
    MQ.PWS("Part\t" + str(round(float(del_part) / float(del_at_least_one_part), 3)) + "\t" + str(round(float(dup_part) / float(dup_at_least_one_part), 3)), logfile);
    MQ.PWS("All\t" + str(round(float(del_full + del_part) / float(del_at_least_one), 3)) + "\t" + str(round(float(dup_full + dup_part) / float(dup_at_least_one), 3)), logfile);

    MQ.PWS("\nTOTAL AVG. OVERLAP CONDITIONAL: " + str(float(del_full + del_part) / float(cnv_at_least_one)), logfile);
    # Table 3

    MQ.PWS("\nNUMBER OF TIMES EACH GENE IS OVERLAPPED:");
    hist = { i : 0 for i in range(1,11) };
    for gene in gene_counts:
        if gene_counts[gene]['count'] <= 10:
            hist[gene_counts[gene]['count']] += 1;
    for i in hist:
        print(i, "\t", hist[i]);

    if ol_type == "GENES":
        MQ.PWS("\nEnsembl ID\tDels\tDups")
        for gene in gene_counts:
            if gene_counts[gene]['count'] > 5:
                MQ.PWS(gene + "\t" + str(len(gene_counts[gene]['del'])) + "\t" + str(len(gene_counts[gene]['dup'])), logfile);
    # Table 4

    if ol_type == "TRANSCRIPTS" and species == "macaque" and outfile_suffix == "":
        del_genes = list(set(del_part_genes + del_full_genes));
        dup_genes = list(set(dup_part_genes + dup_full_genes));
        with open("go/" + species + outfile_suffix + "-noalu-transcripts-del.txt", "w") as out:
            for g in del_genes:
                out.write(g + "\n");
        with open("go/" + species + outfile_suffix + "-noalu-transcripts-dup.txt", "w") as out:
            for g in dup_genes:
                out.write(g + "\n");
    # Write the transcripts to annotate with GO terms next.

############################################################

if "-h" in sys.argv:
    sys.exit("Usage: python gene_count.py <species: macaque or human (required)> <max CNV length to annotate (optional)>");

if len(sys.argv) < 2:
    sys.exit(" * ERROR: Species must be provided: macaque or human");
species = sys.argv[1];
if sys.argv[1] not in ["macaque", "human"]:
    sys.exit(" * ERROR: Species must be provided: macaque or human");
if len(sys.argv) > 2:
    max_cnv_len = sys.argv[2];
    outfile_suffix = "-" + sys.argv[2];
    try:
        max_cnv_len = int(max_cnv_len);
    except:
        sys.exit(" * ERROR: Max CNV length must be an integer.");
else:
    max_cnv_len = 999999999999999;
    outfile_suffix = "";
# Max CNV length checking.

# longest_iso = False;
# Option to run only on longest transcripts.

cnvs_file = "bed/" + species + "-cnvs-filtered-noalu.bed";
gene_file = "bed/" + species + "-cnvs-to-genes.bed";
gene_up_file = "bed/" + species + "-cnvs-to-genes-10kb-up.bed";
gene_down_file = "bed/" + species + "-cnvs-to-genes-10kb-down.bed";
transcript_file = "bed/" + species + "-cnvs-to-transcripts.bed";
exon_file = "bed/" + species + "-cnvs-to-exons.bed";
# File names.

logfilename = "gene-counts-" + species + outfile_suffix + "-noalu.log";
#sys.exit(logfilename);

with open(logfilename, "w") as logfile:
# Log file.

    MQ.runTime("# Macaque CNV gene counter\n# Species:   " + species, logfile);
    MQ.PWS("# Log file:        " + logfilename, logfile);
    MQ.PWS("# Max CNV length:  " + str(max_cnv_len), logfile);
    MQ.PWS("# CNVs file:       " + cnvs_file, logfile);
    MQ.PWS("# Gene file:       " + gene_file, logfile);
    MQ.PWS("# 10kb up file:    " + gene_up_file, logfile);
    MQ.PWS("# 10kb down file:  " + gene_down_file, logfile);
    MQ.PWS("# Transcript file: " + transcript_file, logfile);
    MQ.PWS("# Exon file:       " + exon_file, logfile);
    MQ.PWS("# ----------------", logfile);
    # I/O options and info.

    feature_types = ["gene", "gene-10kb-up", "gene-10kb-down", "transcript", "exon"];
    cnv_types = ["del", "dup"];
    overlap_types = ["full", "partial"];
    # Categories for features, CNVs, and overlaps.

    MQ.PWS("# " + MQ.getDateTime() + " Reading CNVs...", logfile);
    cnvs = [];
    for line in open(cnvs_file):
        cnv = line.strip().split("\t")[3];
        cnv_len = int(cnv.split(":")[3]);
        if cnv_len < max_cnv_len:
            cnvs.append(cnv);
    num_cnvs = len(cnvs);
    MQ.PWS("# CNVs read: " + str(num_cnvs), logfile);

    MQ.PWS("\n# " + MQ.getDateTime() + " Counting gene overlaps...", logfile);
    countOverlaps(cnvs, gene_file, "GENES");

    MQ.PWS("\n# " + MQ.getDateTime() + " Counting 10kb upstream gene overlaps...", logfile);
    countOverlaps(cnvs, gene_up_file, "10KB UP");

    MQ.PWS("\n# " + MQ.getDateTime() + " Counting 10kb downstream gene overlaps...", logfile);
    countOverlaps(cnvs, gene_down_file, "10KB DOWN");

    MQ.PWS("\n# " + MQ.getDateTime() + " Counting transcript overlaps...", logfile);
    countOverlaps(cnvs, transcript_file, "TRANSCRIPTS");

    MQ.PWS("\n# " + MQ.getDateTime() + " Counting exon overlaps...", logfile);
    countOverlaps(cnvs, exon_file, "EXONS");

