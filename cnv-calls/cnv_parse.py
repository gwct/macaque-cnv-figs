#!/usr/bin/python
############################################################
# For macaque, 04.19
# Takes a VCF with SV calls and a CSV file with 
# pedigree/sample info and retrieves information in a format 
# more suitable for R
############################################################

import sys, os, argparse, numpy as np
sys.path.append("../lib/");
import mqcore as MQ
import cnvlib

############################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Script to parse an CNV VCF to a CSV file.");
    parser.add_argument("-i", dest="input", help="A VCF file with called SVs.", default="vcf/mq-duphold-pruned.vcf.gz");
    parser.add_argument("-s", dest="samples", help="A CSV with info on the input samples.", default="../sample-info/macaque-ped-info.csv");
    parser.add_argument("-f", dest="prefix", help="A prefix string to add to all output file names.", default=False);
    args = parser.parse_args();
    # Input options.

    if not os.path.isfile(args.input):
        sys.exit(" ** ERROR: input (-i) file not found!");
    # Check the input file
    if not os.path.isfile(args.samples):
        sys.exit(" ** ERROR: sample (-s) file not found!");
    # Check the input file

    print(" --> Prepping outputs...");
    if not args.prefix:
        prefix = "";
    else:
        prefix = args.prefix + "-";
    outfilename = prefix + "cnvs.csv";
    # Prep the output file.

    with open(outfilename, "w") as outfile:
        MQ.runTime("# CNV VCF parsing", outfile);
        MQ.PWS("# VCF File:      " + args.input, outfile);
        MQ.PWS("# Samples File:  " + args.samples, outfile);
        MQ.PWS("# Output prefix: " + args.prefix, outfile);
        MQ.PWS("# Output file:   " + outfilename, outfile);
        MQ.PWS("# ----------", outfile);

        MQ.PWS("# " + MQ.getDateTime() + " Reading sample info...", outfile);
        samples = cnvlib.readSamples(args.samples);
        MQ.PWS("# " + MQ.getDateTime() + " Samples read: " + str(len(samples)), outfile);
        # Read the sample info.

        MQ.PWS("# " + MQ.getDateTime() + " Reading VCF...", outfile);
        vcf, vcf_headers, vcf_format, human_flag = cnvlib.readVCF(args.input);
        MQ.PWS("# " + MQ.getDateTime() + " Variants read: " + str(len(vcf)), outfile);
        # Read the VCF file.

        MQ.PWS("# " + MQ.getDateTime() + " Fixing VCF headers...", outfile);
        vcf_headers = [ h.replace("NHP-", "") for h in vcf_headers ];
        if "39239A" in vcf_headers:
            vcf_headers[vcf_headers.index("39239A")] = "39239";
        # Header replacement stuff for the macaque data.
                        
        MQ.PWS("# " + MQ.getDateTime() + " Checking CNV transmission...", outfile);
        vcf_passed = cnvlib.childCheck(vcf, vcf_headers, samples);
        # Get the CNVs passed on to F2s.

        MQ.PWS("# " + MQ.getDateTime() + " Parsing variants and writing to output...", outfile);
        cnvlib.outCount(vcf, vcf_passed, vcf_headers, vcf_format, outfile, samples, human_flag);
        # Make the counts for every call and output to a CSV file.

    ind_counts = { i : [0,0,0,0] for i in vcf_headers };
    gt_ind = vcf_format.index("GT");
    totals = [0,0,0,0];
    sv_types = ["<DEL>","<DUP>","<INV>"];
    for key in vcf:
        cur_type = key.split(":")[2];
        #print cur_type;
        if cur_type == "<DUP:TANDEM>":
            cur_type = "<DUP>";
        if cur_type not in sv_types:
            continue;
        sv_ind = sv_types.index(cur_type);

        for s in range(len(vcf_headers)):
            ind = vcf_headers[s];
            gt = vcf[key][s].split(":")[gt_ind];
            if "1" in gt:
                totals[sv_ind] += 1;
                totals[3] += 1;
                ind_counts[ind][sv_ind] += 1;
                ind_counts[ind][3] += 1;

    print("-" * 60); 
    pad = 10;
    print(MQ.spacedOut("Ind",24) + MQ.spacedOut("<DEL>",pad) + MQ.spacedOut("<DUP>",pad) + MQ.spacedOut("<INV>",pad) + MQ.spacedOut("Total",pad));
    for ind in vcf_headers:
        print (MQ.spacedOut(ind,24) + MQ.spacedOut(str(ind_counts[ind][0]),pad) + MQ.spacedOut(str(ind_counts[ind][1]),pad) + MQ.spacedOut(str(ind_counts[ind][2]),pad) + MQ.spacedOut(str(ind_counts[ind][3]),pad));
    print("-" * 60);
    print(MQ.spacedOut("Total:",24) + MQ.spacedOut(str(totals[0]),pad) + MQ.spacedOut(str(totals[1]),pad) + MQ.spacedOut(str(totals[2]),pad) + MQ.spacedOut(str(totals[3]),pad));
    print("-" * 60);    
    print("\nDone!");
    # This block makes a brief count for output to screen.
############################################################