#!/usr/bin/env python
#
# Created by: Alok Joglekar
# Last updated: January 14, 2021
#
# Usage: extract_epitope_fastq.py <SALL_MERGED_R1_fastq> <library_oligos.txt> <>
# Output: a fasta file that extracts the epitope sequence from a SABR library and outputs read counts/epitope
#
#

import sys
import re
import csv

filename = sys.argv[1] #this would be the R1.fastq file 
filename2 = sys.argv[2] #this would be the library oligos file
filename3 = sys.argv[3] #this would be the outputfile
#print (filename, filename2)

###For trial; let's use one output file
#trial_output_filename = "trial_output.csv"


codon_table = {
'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
'TAC':'Y', 'TAT':'Y', 'TAA':'.', 'TAG':'.',
'TGC':'C', 'TGT':'C', 'TGA':'.', 'TGG':'W',
}

epitope_table = {}
with open(filename2,'r') as f:
    for line in f:
       epitope_table[line.strip()] = 0
       #print(epitope_table[line])
       #print(line)

#Reading the fastq file
fastq1 = open(filename, 'r')
line1 = fastq1.readline()
while line1:
    # Get header and sequence
    header1 = line1.strip()
    sequence1 = ""
    line1 = fastq1.readline()
    linecount = 0
    while line1 and (line1[0] != '+'):
        sequence1 = sequence1 + line1.strip()
        line1 = fastq1.readline()
        linecount += 1
        
        # Get quality line
    quality1 = ""
    line1 = fastq1.readline()
    while linecount > 0:
        quality1 = quality1 + line1.strip()
        line1 = fastq1.readline()
        linecount -= 1
        
        
    #print (sequence1) #this works
    #print("\n") #this works 
    #print (barcode1) #this works
    
### This is an example of how to extract subsequence between two strings and output it  -- it works
#s = "abcacbAUG|GAC|UGAfjdalfd"
#start = s.find("AUG|") + len("AUG|")
#end = s.find("|UGA")
#substring = s[start:end]
#print(substring)
# the output is GAC

    # we need to define the start/end substrings from the oligo files
    # entire oligo: CAGGAGGGCTCGGCActgcccagctacgacgaggccgagaggaccaagaccgaggccaccatccccctggtgcccggcagggacgaggacGGAGGTGGTGGCAGCGGCGGCCGCTCATCT
    # oligo_start = 'CAGGAGGGCTCGGCA'
    # oligo_end = 'GGAGGTGGTGGCAG'
    #start = sequence1.find("CAGGAGGGCTCGGCA") + len("CAGGAGGGCTCGGCA")
    #end = sequence1.find("GGAGGTGGTGGCAG")
    #epitope = sequence1[start:end]
    #print(sequence1)
    #print(epitope)
    #print("\n")
    #print(start)
    
    epitope2 = re.search('CAGGAGGGCTCGGCA(.*?)GGAGGTGGTGGCAG', sequence1)
    if epitope2 is None:
        continue
    else:
        epitope = epitope2.group(1)
        #print(epitope2.group(1)) 
        #print(epitope)
        #This is extracting the right portion of the read
        # DNA to AA translation
        amino_acid_sequence = ""
        if "N" in epitope:
            continue
        else:
            for amino_acid in range(len(epitope)//3): #this was changed from / to // -- floor division for updated python
                amino_acid_start = amino_acid*3
                amino_acid_stop = amino_acid*3 + 3
                codon = epitope[amino_acid_start:amino_acid_stop]
                amino_acid_sequence += codon_table[codon]
    
    # Output
            print(amino_acid_sequence) #this works
            if amino_acid_sequence in epitope_table.keys():
                epitope_table[amino_acid_sequence] += 1
                #print (epitope_table[amino_acid_sequence])
            else:
                continue
    #Need to translate the epitope
        

#### Until this, everything is working. The printing is not.. but I
myoutput = open(filename3,"w") ## this is for index 19, so i can test the output
#wr = csv.writer(myoutput,quoting=csv.QUOTE_ALL)
for key, value in epitope_table.items():
    row = str(key) + "," + str(value) + "\n"
    myoutput.write(row)

f.close()
fastq1.close()
myoutput.close()
### THIS SEEMS TO WORK