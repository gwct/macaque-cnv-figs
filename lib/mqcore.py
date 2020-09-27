#!/usr/bin/python3
#######################################################################
# For macaque revisions, 03.2020
# Shared functions for the gene counting and annotation steps.
#######################################################################

import sys, datetime, subprocess

#######################################################################

def PWS(o_line, o_stream=False, std_stream=True):
# Function to print a string AND write it to the file.
	if std_stream:
		print(o_line);
	if o_stream:
		o_stream.write(o_line + "\n");

#######################################################################

def getDateTime():
#Function to get the date and time in a certain format.
	return datetime.datetime.now().strftime("%m.%d.%Y | %H:%M:%S");

#######################################################################

def runTime(msg=False, writeout=False):
# Prints relevant system info at the beginning of a script (or whenever it's called...)
	if msg:
		if not msg.startswith("#"):
			msg = "# " + msg;
		PWS(msg, writeout);

	PWS("# PYTHON VERSION: " + ".".join(map(str, sys.version_info[:3])), writeout)
	PWS("# Script call:    " + " ".join(sys.argv), writeout)
	PWS("# Runtime:        " + datetime.datetime.now().strftime("%m/%d/%Y %H:%M:%S"), writeout);
	PWS("# ----------------", writeout);

#######################################################################

def spacedOut(string, totlen, sep=" "):
# Properly adds spaces to the end of a string to make it a given length
	spaces = sep * (totlen - len(string));
	return string + spaces;

#############################################################################

def getFileLen(i_name):
#Calls 'wc -l' to get the number of lines in the file.
	p = subprocess.Popen(['wc', '-l', i_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE);
	result, err = p.communicate();
	if p.returncode != 0:
		raise IOError(err);
	return int(result.strip().split()[0]);

#######################################################################

def loadingBar(counter, length, done, bars, firstbar=False, disperc=False):
#This function serves as a text loading bar for long scripts with counters. The following
#lines must be added within the script to initialize and terminate the script:
#Initilization:
#numlines = core.getFileLen(alnfilename);
#numbars = 0;
#donepercent = [];
#i = 0;
#Termination:
#	pstring = "100.0% complete.";
#	sys.stderr.write('\b' * len(pstring) + pstring);
#	print "\nDone!";
#
#If length is lines in a file use the core.getFileLen function to get the number of lines in the file

	# try:
	# 	if sys.version[0] == '2':
	# 		pchr = u'\u2591'.encode('utf-8');
	# 		lchr = u'\u2588'.encode('utf-8');
	# 	elif sys.version[0] == '3':
	# 		pchr = u'\u2591';
	# 		lchr = u'\u2588';		
	# except:
	# 	pchr, lchr = "*","*";

	try:
		#pchr, lchr=u'\u2591',u'\u2588';
		pchr, lchr=u'\u2588',u'\u2591';
	except:
		pchr, lchr = "=","|";

	percent = float(counter) / float(length) * 100.0;
	percentdone = int(percent);

	p = str(percent)
	pstring = " " + p[:5] + "% complete.";

	if percentdone % 2 == 0 and done != None and percentdone not in done:
		loading = "|";
		j = 0;
		while j < bars:
			loading += pchr;
			j += 1;
		if j <= 49:
			loading += lchr;
		else:
			loading += pchr;
		j += 1;
		if j == 50:
			loading = loading[:-1] + pchr;

		while j < 50:
			loading += "-";
			j += 1;
		loading += "|";

		if disperc:
			loading += "                 ";
		if firstbar:
			sys.stderr.write(loading);
			firstbar = False
		else:
			sys.stderr.write('\b' * len(loading) + loading);

		done.append(percentdone);
		bars = bars + 1;
	if disperc:
		sys.stderr.write('\b' * len(pstring) + pstring);
	sys.stderr.flush();
	
	return bars, done, firstbar;

#######################################################################