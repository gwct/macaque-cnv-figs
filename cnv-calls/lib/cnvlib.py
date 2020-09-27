############################################################
# sv support modules -- used by sv_parse.py
# 04.19
############################################################

import sys, os, gzip
from collections import defaultdict

############################################################
def readSamples(samplefile):
# Read CSV file with pedigree info.
    samples = {};
    first = True;
    for line in open(samplefile):
        line = line.strip().split(",");
        if first:
            header_inds = { i : line[i] for i in range(len(line)) };
            first = False;
            continue;

        ind = line[1];
        samples[ind] = {};
        for col in range(len(line)):
            if line[col] == ind:
                continue;
    
            samples[ind][header_inds[col]] = line[col];

    return samples;
############################################################
def getFileReader(i_name):
# Check if a file is gzipped, and if so set gzip as the file reader. Otherwise, read as a normal text file.
	try:
		gzip_check = gzip.open(i_name).read(1);
		reader = gzip.open;
	except:
		reader = open;
	return reader;
############################################################
def readVCF(vcf_input):
# Reads the input VCF files
    vcf = defaultdict(list);
    human_flag = False;
    reader = getFileReader(vcf_input);

    for line in reader(vcf_input):
        if reader == gzip.open:
            line = line.decode();

        line = line.strip().split("\t");
        if line[0].startswith("##"):
            if line[0].startswith("##INFO=<ID=GENES,"):
                human_flag = True;
            continue;
        elif line[0].startswith("#"):
            vcf_headers = line[line.index("FORMAT")+1:];
            continue;
        info = line[7].split(";");
        svlen = info[2][info[2].index("=")+1:].replace("-","");
        svtype = line[4];
        if svtype == "<DUP:TANDEM>":
            svtype = "<DUP>";
        qual = line[5];
        if svtype not in ["<DEL>", "<DUP>", "<INV>"]:
            continue;

        if human_flag:
            for i in info:
                if i.startswith("ALGORITHM"):
                    algo = i.split("=")[1];
                    if "_" in algo:
                        algo = "Both";               
                if i.startswith("GENOTYPER"):
                    genotyper = i.split("=")[1];
                    if "_" in genotyper:
                        genotyper = "Both";
                if i.startswith("PARENT_AF"):
                    paf = i.split("=")[1];
            key = ":".join([line[0], line[1], svtype, svlen, qual, algo, genotyper, paf]);
        # This parses the extra info in the Brandler VCF.

        else:
            key = ":".join([line[0], line[1], svtype, svlen, qual]);
        # Info for the macaque VCF.

        vcf_format = line[8].split(":");
        if vcf_format[-1] == "DHSP":
            vcf_format = vcf_format[:-1];
        samples = line[9:];
        vcf[key] = samples;

    return vcf, vcf_headers, vcf_format, human_flag;
############################################################
def getDenovo(vcf_line, vcf_headers, samples, vcf_format, k, debug):
    denovo = [];
    ref_event = "N";
    gt_ind = vcf_format.index("GT");
    presence = [];
    for col in range(len(vcf_line)):
        ind = vcf_headers[col];
        gt = vcf_line[col].split(":")[gt_ind];
        if "1" in gt:
            presence.append(ind);

    sv_af = float(len(presence)) / float(len(samples));

    if sv_af >= 0.95:
        ref_event = "Y";

    if presence == []:
        denovo = "N";
        return denovo, 0.0, ref_event;

    if len(presence) == 1 and (samples[presence[0]]['Focal'] == "T" or samples[presence[0]]['F2'] == "T"):
        denovo.append(presence[0]);

    else:
        for ind in samples:
            relatives = [];
            if samples[ind]['Focal'] == "T" or samples[ind]['F2'] == "T":
                relatives.append(ind);
                childs = samples[ind]['Child ID'].split(";");
                sibs = samples[ind]['Sibling ID'].split(";");
                relatives += childs
                relatives += sibs;
                relatives = [ r for r in relatives if r != "NA" ];

            # if debug and ind in ['F0122-05|REACH000312', 'F0122-01|REACH000145', 'F0122-06|REACH000313']:
            #     print ind;
            #     print presence;
            #     print relatives;
            #     print "-"*10;
            
            if all(p in relatives for p in presence):
                #denovo = ind;
                denovo.append(ind);

    if len(denovo) > 1:
        indcheck = denovo[0];
        sibcheck = denovo[1:];
        sibs = samples[indcheck]['Sibling ID'].split(";");
        if not all(p in sibcheck for p in sibs):
            print("More than one denovo found in non-siblings O_o")
            print(vcf_line)
            print(k)
            print(denovo);
            print(indcheck);
            print(sibcheck);
            print(samples[indcheck]['Sibling ID'])
            print(sibs);
            print(presence);
            sys.exit();

    if denovo == []:
        denovo = "N"

    # if len(presence) == 1 and samples[presence[0]]['Focal'] == "T" and samples[presence[0]]['F2'] == "F":
    #     denovo = presence[0];
    # elif len(presence) == 2:
    #     if samples[presence[0]]['Focal'] == "T" and samples[presence[0]]['Child ID'] == presence[1]:
    #         denovo = presence[0];
    #     elif samples[presence[1]]['Focal'] == "T" and samples[presence[1]]['Child ID'] == presence[0]:
    #         denovo = presence[1];

    return denovo, str(sv_af), ref_event;
############################################################
def childCheck(vcf, vcf_headers, samples):
## Checks if a mutation given by vcf_key for a given trio has been passed on to the F1s children.

    passed = defaultdict(list);
    for key in vcf:
        for s in range(len(vcf_headers)):
            ind = vcf_headers[s];
            ind_ind = vcf_headers.index(ind);
            if samples[ind]['Focal'] == "F":
                continue;
            cur_ped = samples[ind]['Family ID'];
            cur_f2 = samples[ind]['Child ID'];
            if cur_f2 == "NA":
                continue;
            cur_f2_ind = vcf_headers.index(cur_f2);
            f1_gt = vcf[key][ind_ind].split(":")[0];
            f2_gt = vcf[key][cur_f2_ind].split(":")[0];
            if "1" in f1_gt and "1" in f2_gt:
                passed[key].append(cur_ped);
    return passed;

############################################################
def cnvParse(sv_item):
    sv_key, vcf, vcf_passed, vcf_headers, samples, vcf_format, gt_ind, debug, human_flag = sv_item;

    outlines = [];

    var = sv_key.split(":");   
    # The current SV is encoded in the key. Split to get the separate pieces of info.

    denovo, sv_af, ref_event = getDenovo(vcf[sv_key], vcf_headers, samples, vcf_format, sv_key, debug);
    # Check whether the current variant is de novo.

    for s in range(len(vcf_headers)):
        ind = vcf_headers[s];
        # Go through every individual in the vcf header. 

        cur_trio = samples[ind]['Family ID'];
        # The current trio is encoded in the sample data.

        outline = ",".join([sv_key,cur_trio,ind,samples[ind]['Sex'],var[0],var[2],var[1],var[3],var[4]]) + ",";
        # Basic info about the SV.

        gt_info = vcf[sv_key][s].split(":");
        # Getting the genotype info from the current sample field.

        if len(gt_info) == len(vcf_format) + 1:
            gt_info = gt_info[:-1];
        # This is because some of the entries in the macaque data didn't have DHSP calculated by duphold for some reason...

        if len(gt_info) != len(vcf_format):
            print(" * NOT ENOUGH FORMAT COLUMNS:");
            print(sv_key);
            print(len(vcf_format), "\t".join(vcf_format));
            print(len(gt_info), "\t".join(gt_info));
            sys.exit(1);
        # A catch in case something is off for the given sample.

        if gt_info[gt_ind] in ["0/0","./.", "0"]:
            continue;
        # We're only interested in samples with the SV.

        for field_col in range(len(gt_info)):
            if field_col == gt_ind:
                if gt_info[field_col] in ["0/1","1/0"]:
                    gt_info[field_col] = "Het";
                else:
                    gt_info[field_col]  = "HomAlt";
                # If the current field is the genotype field, re-encode as this string.

            if gt_info[field_col] == ".":
                gt_info[field_col] = "NA";
            # If the current field doesn't have an entry, re-encode as NA.

            gt_info[field_col] = gt_info[field_col].replace(",",";");
            # Replace any commas with semicolons in the current field since we're saving to a csv.

            outline += gt_info[field_col] + ",";
        # Go through every field in the current sample and add it to the output.
        #print(outline);

        trans = "NA";
        f1 = "N";
        if samples[ind]['Focal'] == "T":
            f1 = "Y";
            if cur_trio in vcf_passed[sv_key]:
                trans = "Y";
            else:
                trans = "N";
        # This checks if the current sample is an F1 and if the variant has been transmitted to the F2.
        
        if denovo != "N" and ind in denovo:
            ind_denovo = "Y";
            denovo = "N";
        else:
            ind_denovo = "N";

        trio_f = samples[ind]['Father ID'];
        trio_m = samples[ind]['Mother ID'];
        
        if samples[ind]['Relationship'] in ["Father", "Mother"]:
            trio_f1 = samples[ind]['Child ID'];
        else:
            trio_f1 = ind;

        pat_age = samples[ind]['Father age'];
        mat_age = samples[ind]['Mother age'];
        # Saving redundant trio information so we don't have to look up other rows later.

        outline += ",".join([f1,trans,ind_denovo,trio_f,trio_m,trio_f1,pat_age,mat_age,sv_af,ref_event]);
        #print(outline);
        if human_flag:
            outline += "," + ",".join([var[5], var[6], var[7]]);
        # Add everything to the output line.

        outlines.append(outline);
    return outlines;

############################################################
def outCount(vcf, vcf_passed, vcf_headers, vcf_format, outfile, samples, human_flag):
    debug = False;
    gt_ind = vcf_format.index("GT");
    # Get the index of the genotype field from the FORMAT column.

    out_headers = "CNV key,Pedigree,Individual,Sex,Chromosome,Type,Pos,Length,Qual,";
    for field in vcf_format:
        out_headers += field + ",";
    out_headers += "F1,Transmitted,Denovo,Trio F,Trio M,Trio F1,Paternal age,Maternal age,CNV freq,Ref 95";
    # Output headers.

    if human_flag:
        out_headers += ",Algorithm,Genotyper,Parent AF";
        # Add to the header if we're outputting genes.

    num_calls = float(len(vcf));
    numbars, donepercent, counter = 0, [], 0;

    outfile.write(out_headers + "\n");
    # Write the header line.

    counter = 0;
    for sv_key in vcf:
        result = cnvParse((sv_key, vcf, vcf_passed, vcf_headers, samples, vcf_format, gt_ind, debug, human_flag));
        if result != []:
            for outline in result:
                outfile.write(outline + "\n");
        # counter += 1;
        # if (float(counter) / num_calls) % 10 == 0:
        #     print(str(float(counter) / num_calls) + "%");
    return;
############################################################

