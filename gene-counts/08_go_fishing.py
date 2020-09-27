#!/usr/bin/python3
#######################################################################
# For macaque revisions, 03.2020
# Fisher's test for GO enrichment
#######################################################################

import sys, os, scipy.stats as stats, lib.mqcore as MQ
import argparse

#######################################################################

def optParse():
# Function to parse input options and check for errors.
	parser = argparse.ArgumentParser(description="Fisher's test for GO enrichment");
	parser.add_argument("-q", dest="query_file", help="A subset of GO terms to test for enrichment.");
	parser.add_argument("-b", dest="background_file", help="A set of GO terms used as the background for the enrichment test.");
	parser.add_argument("-a", dest="alpha_input", help="Alpha: the p-value threshold", type=float, default=0.05);
	parser.add_argument("-c", dest="correction_method", help="The multiple test correction method: Bonferroni (b), Dunn-Sidak (ds), or False discovery rate (fdr). Default: None", default="None");
	parser.add_argument("-o", dest="output_file", help="Output file to write enriched GO terms.")
	args = parser.parse_args();

	if None in [args.query_file, args.background_file, args.output_file]:
		sys.exit(MQ.errorOut(1, "All of query (-q), background (-b), and output (-o) files must be specified!"));
	if not os.path.exists(os.path.abspath(args.query_file)):
		sys.exit(MQ.errorOut(2, "Query file not found!"));
	if not os.path.exists(os.path.abspath(args.background_file)):
		sys.exit(MQ.errorOut(3, "Background file not found!"));
	if args.alpha_input <= 0 or args.alpha_input >= 1:
		sys.exit(MQ.errorOut(4, "The p-value threshold (alpha, -a) must be between 0 and 1."))
	if args.correction_method not in ['b', 'ds', 'fdr', "None"]:
		sys.exit(MQ.errorOut(5, "Correction (-c) method must be one of b, ds, and fdr! Alternatively, if -c is unspecified, no correction will be done."));

	return args.query_file, args.background_file, args.alpha_input, args.correction_method, args.output_file;

#########################

def enrichedOut(go_dict, threshold):
	# Function that goes through all test results and reports those below adjusted p-value threshold.
	MQ.PWS("# ----------", outfile);
	outfile.write("# GO Accession\tFamily ID\t# PS w/GO / # PS w/o GO\t# background w/ GO / # background w/o GO\t\tp-value\tOdds ratio\tGO term name\tGO domain\tGO definition\n");
	enriched = 0;
	for go_acc in go_dict:
		pval = go_dict[go_acc][2];
		if pval <= threshold:
			enriched += 1;
			outline = go_acc + "\t";
			for col in go_dict[go_acc]:
				outline += str(col) + "\t";
			outfile.write(outline[:-1] + "\n");

	return enriched

#########################

def correctionStr(c):
# Function to get the nice correction string to report in the log.
	if c == "None":
		cstr = "None"
	if c == "b":
		cstr = "Bonferroni";
	if c == "d":
		cstr = "Dunn-Sidak";
	if c == "fdr":
		cstr = "False Discovery Rate";
	return cstr;

#######################################################################

queryfile, bgfile, alpha, correction, outfilename = optParse();
correction_str = correctionStr(correction);
io_pad = 35;
r_pad = 53;
# Parse input options.

with open(outfilename, "w") as outfile:
	MQ.runTime("# Fisher's test for GO enrichment", outfile);
	MQ.PWS(MQ.spacedOut("# Query file:", io_pad) + queryfile, outfile);
	MQ.PWS(MQ.spacedOut("# Background file:", io_pad) + bgfile, outfile);
	MQ.PWS(MQ.spacedOut("# Alpha (p-value threshold):", io_pad) + str(alpha), outfile);
	MQ.PWS(MQ.spacedOut("# Multiple test correction method:", io_pad) + correction_str, outfile);
	if correction == "None":
		MQ.PWS("# --> WARNING: Not correcting for multiple tests!", outfile);
	MQ.PWS(MQ.spacedOut("# Output file:", io_pad) + outfilename, outfile);
	MQ.PWS("# ----------", outfile);
	# Report I/O information.

	MQ.PWS("# " + MQ.getDateTime() + " Counting total query GO terms...", outfile);
	query_genes = [];
	query_go_count = 0;
	for line in open(queryfile):
		if line[0] == "#":
			continue;
		line = line.strip().split("\t");
		gid = line[0];
		if gid not in query_genes:
			query_genes.append(gid);
		query_go_count += 1;
	# Get count of query GO terms and unique list of features.

	MQ.PWS("# " + MQ.getDateTime() + " Counting total background GO terms...", outfile);
	background_genes = [];
	bg_go_count = 0;
	for line in open(bgfile):
		if line[0] == "#":
			continue;
		line = line.strip().split("\t");
		gid = line[0];
		if gid in query_genes:
			continue;
		if gid not in background_genes:
			background_genes.append(gid);
		bg_go_count += 1;
	# Get count of background GO terms and unique list of features.

	MQ.PWS(MQ.spacedOut("# " + MQ.getDateTime() + " Total query GO terms:", r_pad) + str(query_go_count) + " in " + str(len(query_genes)) + " genes.", outfile);
	MQ.PWS(MQ.spacedOut("# " + MQ.getDateTime() + " Total background GO terms:", r_pad) + str(bg_go_count) + " in " + str(len(background_genes)) + " genes.", outfile);
	MQ.PWS("# " + MQ.getDateTime() + " Running Fisher's Tests...", outfile);

	done, pvals, results_dict, num_tests = [], [], {}, 0;
	numlines = MQ.getFileLen(queryfile);
	i, numbars, donepercent = 0, 0, [];

	print();
	for line in open(queryfile):
		numbars, donepercent, fbar = MQ.loadingBar(i, numlines, donepercent, numbars, disperc=True);
		i += 1;
		if line[0] == "#":
			continue;

		line = line.strip().split("\t");
		if len(line) != 5:
			continue;

		geneid, go_acc = line[0], line[1];
		# Parse the line

		if go_acc in done:
			continue;
		done.append(go_acc);
		# If we've already tested this GO term, skip.

		num_tests += 1;

		w, x, y, z = 0, 0, 0, 0;
		# w = # genes with GO term in query
		# x = # genes with GO term in background
		# y = # genes without GO term in query
		# z = # genes without GO term in background

		for next_line in open(queryfile):
			if next_line[0] == "#":
				continue;
			next_line = next_line.strip().split("\t");
			next_go_acc = next_line[1];
			if next_go_acc == go_acc:
				w += 1;
		y = query_go_count - w;
		# Count the number of genes with the current GO term in the query (w) and
		# the number without it (y).

		for not_line in open(bgfile):
			if not_line[0] == "#":
				continue;
			not_line = not_line.strip().split("\t");
			not_go_acc = not_line[1];
			if not_go_acc == go_acc:
				x += 1;
		z = bg_go_count - x;
		# Count the number of genes with the current GO term in the background (x) and
		# the number without it (z)

		oddsratio, pvalue = stats.fisher_exact([[w,y],[x,z]], alternative='greater');
		try:
			results_dict[go_acc] = [str(w) + "/" + str(y), str(x) + "/" + str(z), pvalue, str(oddsratio), line[2], line[3], line[4]];
		except:
			print(line);
			sys.exit();
		pvals.append(pvalue);
		# Call scipy's stats module to run the Fisher's test and save the output.

	pstring = "100.0% complete.";
	sys.stderr.write('\b' * len(pstring) + pstring);
	sys.stderr.flush();
	print("\n# " + MQ.getDateTime() +  " Done!\n");

	num_tests = float(num_tests);
	MQ.PWS(MQ.spacedOut("# " + MQ.getDateTime() + " Number of tests:", r_pad) + str(num_tests), outfile);
	# Report the number of tests run.

	if correction == "None":
		num_enriched = enrichedOut(results_dict, alpha);

	if correction == "b":
		corrected_alpha = (alpha / num_tests);
		MQ.PWS(MQ.spacedOut("# " + MQ.getDateTime() + " Bonferroni corrected alpha:", r_pad) + str(corrected_alpha), outfile);
		num_enriched = enrichedOut(results_dict, corrected_alpha);

	elif correction == "ds":
		corrected_alpha = (1 - (1 - alpha)**(1/num_tests))
		MQ.PWS(MQ.spacedOut("# " + MQ.getDateTime() + " Dunn-Sidak corrected alpha:", r_pad) + str(corrected_alpha), outfile);
		num_enriched = enrichedOut(results_dict, corrected_alpha);

	elif correction == "fdr":
		pvals = sorted(pvals);
		with open("test2.txt", "w") as pfile:
			pfile.write(",".join([ str(p) for p in pvals ]));
		qvals = [];
		corrected_alpha = 0.0;
		for x in range(len(pvals)):
			pval = pvals[x];
			qval = (float(x)/num_tests) * alpha;
			qvals.append(qval);
		for x in range(len(pvals)):
			if pvals[x] < qvals[x]:
				corrected_alpha = pvals[x];

		MQ.PWS(MQ.spacedOut("# " + MQ.getDateTime() + " FDR corrected alpha:", r_pad) + str(corrected_alpha), outfile);
		num_enriched = enrichedOut(results_dict, corrected_alpha);
	# Report enriched GO terms after using specified multiple testing correction

	outfile.write("# ----------\n")
	MQ.PWS(MQ.spacedOut("# " + MQ.getDateTime() + " Number enriched:", r_pad) + str(num_enriched), outfile);
	# Report total number enriched.