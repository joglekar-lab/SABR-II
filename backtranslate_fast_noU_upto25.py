#!/usr/bin/env python
#
# Usage: backtranslate.py IEDB.csv
#
#Modified to accomodate epitopes upto 25 AA
import sys
from subprocess import PIPE, Popen

if len(sys.argv) < 2:
    sys.exit("Usage: backtranslate_fast.py IEDB.csv")

AA_table = {
'A':'GCC', 'C':'TGC', 'D':'GAC', 'E':'GAG', 
'F':'TTC', 'G':'GGC', 'H':'CAC', 'I':'ATC', 
'K':'AAG', 'L':'CTG', 'M':'ATG', 'N':'AAC', 
'P':'CCC', 'Q':'CAG', 'R':'AGG', 'S':'AGC', 
'T':'ACC', 'V':'GTG', 'W':'TGG', 'Y':'TAC',
'X':'NNN'}
    
non_AA = ["B", "J", "O", "U", "Z"]

filename = sys.argv[1]
outfilename = filename + ".nuc.csv"
csv = open(filename, 'r')
outfile = open(outfilename, 'w')

# Set up a counter for progress output
count = 0

# Filter duplicates (after stripping modifications)
AA_unique = []

# Read input csv, directly outputted from IEDB
# Skips first two headers
line = csv.readline()
while line:
    if count%1000 == 0:
        print("Processed " + str(count) + " lines")
    count += 1

    # Read Amino Acid
    line_stripped = line.strip()
    line_array = line_stripped.split(',')
    name = line_array[0]
    AA = line_array[1].upper()

    #Translate AA
    #Strip modifications
    if ' ' in AA:
        AA = AA.split(' ')[0]
    
    #Stringency conditions
    #if ' ' in AA:
    #    line = csv.readline()
    #    continue

    if 'X' in AA:
        line = csv.readline()
        print("\"X\" found: %s" % AA)
        continue

    if '.' in AA:
        line = csv.readline()
        print("\".\" found: %s" % AA)
        continue

    if len(AA) < 7:
        line = csv.readline() 
        print("less than 8 AA: %s" % AA)
        continue

    if len(AA) > 25:
        line = csv.readline()
        print("more than 25 AA: %s" % AA)
        continue

    # Duplicate
    if AA in AA_unique:
        print("Duplicate found: %s" % AA)
        line = csv.readline()
        continue

    AA_unique.append(AA)

    nucleotide = ""
    for current_AA in AA:
	if current_AA in non_AA:
		continue        
	else:
		nucleotide += AA_table[current_AA]
    
    line_out = name + ',' + AA + ',' + nucleotide + '\n'
    
    outfile.write(line_out)

    line = csv.readline()

csv.close()
outfile.close()
