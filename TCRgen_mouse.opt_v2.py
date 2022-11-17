#!/usr/bin/env python
#
# Created by: Michael Leonard
# Last updated: September 20, 2018
#
# Usage:   TCRgen.py seq_database.fasta species Va Ja CDR3a Vb Jb CDR3b
# Example: TCRgen.py GENE-DB.all.fasta "Homo sapiens" "TRAV29/DV5*01" 
#          "TRAJ52*01" "CAASSSAGGTSYGKLTF" "TRBV6-2*01" "TRBJ2-3*01" 
#          "CASRPRDPVTQYF"
# 
# on 05/13/21: Changed all the variables that say "None" to "0"

import sys
import os
from datetime import datetime
import random
import subprocess

if len(sys.argv) != 9:
    sys.exit("Usage: TCRgen.py seq_database.fasta species Va Ja CDR3a Vb Jb CDR3b")

# Windows line endings in genbank output
windows = '\r'
#windows = ''

# For calculating overlap score
match_multiplier = 2
gap_penalty = 1

fasta_filename = sys.argv[1]
species = sys.argv[2]
Va_gene_name = sys.argv[3]
Ja_gene_name = sys.argv[4]
CDR3a_sequence = sys.argv[5]
Vb_gene_name = sys.argv[6]
Jb_gene_name = sys.argv[7]
CDR3b_sequence = sys.argv[8]
perform_optimization = True

if CDR3a_sequence == "" or CDR3b_sequence == "":
    sys.stderr.write("Error: No CDR3a or CDR3b sequence given\n")
    exit()
CDR3a_sequence = CDR3a_sequence.strip()
CDR3b_sequence = CDR3b_sequence.strip()
CDR3a_sequence = CDR3a_sequence.upper()
CDR3b_sequence = CDR3b_sequence.upper()
    
# Hardcode human constant regions, 
# they contain introns in input file,
# would have to concatenate
# Alpha "TRAC*01"
#For the mouse TCRs, these are taken from IMGT GeneDB
TRAC_AA = "DIQNPEPAVYQLKDPRSQDSTLCLFTDFDSQINVPKTMESGTFITDKTVLDMKAMDSKSNGAIAWSNQTSFTCQDIFKETNATYPSSDVPCDATLTEKSFETDMNLNFQNLSVMGLRILLLKVAGFNLLMTLRL"

TRAC_nucleotide = "GACATCCAGAACCCAGAACCTGCTGTGTACCAGTTAAAAGATCCTCGGTCTCAGGACAGCACCCTCTGCCTGTTCACCGACTTTGACTCCCAAATCAATGTGCCGAAAACCATGGAATCTGGAACGTTCATCACTGACAAAACTGTGCTGGACATGAAAGCTATGGATTCCAAGAGCAATGGGGCCATTGCCTGGAGCAACCAGACAAGCTTCACCTGCCAAGATATCTTCAAAGAGACCAACGCCACCTACCCCAGTTCAGACGTTCCCTGTGATGCCACGTTGACCGAGAAAAGCTTTGAAACAGATATGAACCTAAACTTTCAAAACCTGTCAGTTATGGGACTCCGAATCCTCCTGCTGAAAGTAGCGGGATTTAACCTGCTCATGACGCTGAGGCTG"

# Beta "TRBC1*01"
#In the IMGT Gene DB, mouse TRBC2 starts with X as the first AA. Based on the OTII TCR sequence from Addgene (52112) the first AA used was "E"
TRBC_AA = "EDLRNVTPPKVSLFEPSKAEIANKQKATLVCLARGFFPDHVELSWWVNGKEVHSGVSTDPQAYKESNYSYCLSSRLRVSATFWHNPRNHFRCQVQFHGLSEEDKWPEGSPKPVTQNISAEAWGRADCGITSASYHQGVLSATILYEILLGKATLYAVLVSGLVLMAMVKKKNS"

TRBC_nucleotide = "GAGGATCTGAGAAATGTGACTCCACCCAAGGTCTCCTTGTTTGAGCCATCAAAAGCAGAGATTGCAAACAAACAAAAGGCTACCCTCGTGTGCTTGGCCAGGGGCTTCTTCCCTGACCACGTGGAGCTGAGCTGGTGGGTGAATGGCAAGGAGGTCCACAGTGGGGTCAGCACGGACCCTCAGGCCTACAAGGAGAGCAATTATAGCTACTGCCTGAGCAGCCGCCTGAGGGTCTCTGCTACCTTCTGGCACAATCCTCGAAACCACTTCCGCTGCCAAGTGCAGTTCCATGGGCTTTCAGAGGAGGACAAGTGGCCAGAGGGCTCACCCAAACCTGTCACACAGAACATCAGTGCAGAGGCCTGGGGCCGAGCAGACTGTGGAATCACTTCAGCATCCTATCATCAGGGGGTTCTGTCTGCAACCATCCTCTATGAGATCCTACTGGGGAAGGCCACCCTATATGCTGTGCTGGTCAGTGGCCTGGTGCTGATGGCCATGGTCAAGAAAAAAAATTCC"


# Other sequences for genbank file
In_Fusion_Left = "ATCTCGAATCGAATTC"
LNGFR = "ATGTCCGGAGCGGGCGCAACTGGGCGGGCGATGGATGGACCACGCCTCCTGTTGCTCCTTCTGCTTGGCGTTTCTCTGGGCGGAGCAAAGGAAGCATGCCCTACAGGTCTCTATACACACTCAGGTGAGTGTTGTAAGGCTTGTAACCTGGGAGAAGGCGTCGCCCAGCCCTGCGGAGCTAATCAGACAGTGTGCGAGCCGTGTCTCGATAGCGTGACGTTCTCCGACGTTGTGTCCGCGACTGAACCATGCAAACCTTGTACAGAATGCGTCGGCCTGCAGTCTATGAGCGCTCCCTGCGTGGAGGCTGACGATGCGGTGTGCCGGTGTGCATACGGGTATTATCAAGATGAAACCACAGGGCGATGTGAGGCCTGCCGAGTATGTGAAGCAGGAAGTGGGCTGGTATTCTCCTGTCAGGACAAGCAGAATACCGTTTGCGAGGAGTGTCCCGACGGCACATATAGCGATGAGGCTAACCACGTTGATCCCTGTCTGCCCTGTACAGTCTGCGAAGACACAGAAAGGCAGCTGCGGGAATGTACACGCTGGGCCGACGCGGAGTGTGAGGAGATCCCAGGTAGATGGATCACTCGCTCTACACCACCCGAGGGAAGCGATTCAACAGCACCAAGCACCCAGGAGCCTGAGGCACCGCCGGAGCAGGACCTGATTGCATCAACGGTTGCAGGAGTTGTGACCACTGTAATGGGTAGTTCTCAGCCAGTGGTAACTCGCGGAACCACTGACAACCTCATTCCAGTATACTGTAGCATCCTCGCAGCCGTGGTTGTCGGGCTGGTGGCCTACATAGCCTTCAAACGGTGGAATAGT"
LNGFR_additional = "TCAGGATCCGGT"
P2A = "GCCACCAACTTCAGCCTGCTGAAGCAGGCCGGCGACGTGGAGGAGAACCCCGGCCCCTCT"
HGH_SS_1 = "ATGGCGACGGGTTCAAGAACTTCCCTACTTCTTGCATTTGGCCTGCTTTGTTTGCCGTGGTTACAGGAAGCCTCAGCA"
F2A = "AGGGCAAAACGTTCGGGTTCGGGTGCGCCAGTAAAGCAGACATTAAACTTTGATTTGCTGAAACTTGCAGGTGATGTAGAGTCAAATCCAGGTCCA"
HGH_SS_2 = "ATGGCAACAGGGAGCCGAACCTCTCTGCTCCTTGCTTTCGGGCTCCTTTGCCTACCGTGCCTGCAGGAGGGCTCGGCA"
stop_codon = "TAA"
In_Fusion_Right = "GAATTCGTTAACCTCG"

# Table for translation of germline sequences
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

# Table for reverse translation of CDR3
# Use most abundant codons in human genome
# Use "N" for "X", as occurs in the human constant region
AA_table = {
'A':'GCC', 'C':'TGC', 'D':'GAC', 'E':'GAG', 
'F':'TTC', 'G':'GGC', 'H':'CAC', 'I':'ATC', 
'K':'AAG', 'L':'CTG', 'M':'ATG', 'N':'AAC', 
'P':'CCC', 'Q':'CAG', 'R':'AGG', 'S':'AGC', 
'T':'ACC', 'V':'GTG', 'W':'TGG', 'Y':'TAC',
'X':'AAC'}

germline_database = {}

# Parse IMGT Fasta file
#   Header delimited by pipe "|"
#   The FASTA header contains 15 fields separated by '|':
#
# 0. IMGT/LIGM-DB accession number(s)
# 1. IMGT gene and allele name
# 2. species
# 3. IMGT allele functionality
# 4. exon(s), region name(s), or extracted label(s)
# 5. start and end positions in the IMGT/LIGM-DB accession number(s)
# 6. number of nucleotides in the IMGT/LIGM-DB accession number(s)
# 7. codon start, or 'NR' (not relevant) for non coding labels
# 8. +n: number of nucleotides (nt) added in 5' compared to the 
#     corresponding label extracted from IMGT/LIGM-DB
# 9. +n or -n: number of nucleotides (nt) added or removed in 3' 
#     compared to the corresponding label extracted from IMGT/LIGM-DB
# 10.+n, -n, and/or nS: number of added, deleted, and/or substituted 
#     nucleotides to correct sequencing errors, or 'not corrected' if 
#     non corrected sequencing errors
# 11. number of amino acids (AA): this field indicates that the sequence is in amino acids
# 12. number of characters in the sequence: nt (or AA)+IMGT gaps=total
# 13. partial (if it is)
# 14. reverse complementary (if it is)
# 
fasta = open(fasta_filename, 'r')
line = fasta.readline()
while line:
    # Get header and sequence
    header = line.strip()
    sequence = ""
    line = fasta.readline()
    while line and (line[0] != '>'):
        sequence = sequence + line.strip()
        line = fasta.readline()
    
    # Parse header, extract allele and species
    # for the germline database
    header_array = header.split('|')
    species_trimmed = header_array[2].split('_')[0]
    allele_and_species = header_array[1] + '_' + species_trimmed
    
    # Prepare sequence for database
    # Convert to uppercase 
    # and trim out of frame starting bases
    sequence_start = int(header_array[7]) - 1
    sequence = sequence.upper()
    sequence = sequence[sequence_start:]
    
    # Enter the gene into the database
    germline_entry = {}
    germline_entry["header"] = header
    germline_entry["sequence"] = sequence
    germline_database[allele_and_species] = germline_entry

fasta.close()
    
    
# Get relevant sequences for input TCR
Va_gene = germline_database[Va_gene_name + '_' + species]
Ja_gene = germline_database[Ja_gene_name + '_' + species]
Vb_gene = germline_database[Vb_gene_name + '_' + species]
Jb_gene = germline_database[Jb_gene_name + '_' + species]


# Convert Variable alpha (Va) to AA
Va_nucleotide = Va_gene["sequence"]
Va_AA = ""    
Va_germline_length = len(Va_nucleotide)
for amino_acid in range(int(Va_germline_length/3)):
    amino_acid_start = amino_acid*3
    amino_acid_stop = amino_acid*3 + 3
    codon = Va_nucleotide[amino_acid_start:amino_acid_stop]
    Va_AA += codon_table[codon]
#print(Va_AA)
#print(CDR3a_sequence)

# Convert Joining alpha (Ja) to AA
Ja_nucleotide = Ja_gene["sequence"]
Ja_AA = ""    
Ja_germline_length = len(Ja_nucleotide)
for amino_acid in range(int(Ja_germline_length/3)):
    amino_acid_start = amino_acid*3
    amino_acid_stop = amino_acid*3 + 3
    codon = Ja_nucleotide[amino_acid_start:amino_acid_stop]
    Ja_AA += codon_table[codon]
#print(Ja_AA)


# Convert Variable beta (Vb) to AA
Vb_nucleotide = Vb_gene["sequence"]
Vb_AA = ""    
Vb_germline_length = len(Vb_nucleotide)
for amino_acid in range(int(Vb_germline_length/3)):
    amino_acid_start = amino_acid*3
    amino_acid_stop = amino_acid*3 + 3
    codon = Vb_nucleotide[amino_acid_start:amino_acid_stop]
    Vb_AA += codon_table[codon]
#print(Vb_AA)
#print(CDR3b_sequence)

# Convert Joining beta (Jb) to AA
Jb_nucleotide = Jb_gene["sequence"]
Jb_AA = ""    
Jb_germline_length = len(Jb_nucleotide)
for amino_acid in range(int(Jb_germline_length/3)):
    amino_acid_start = amino_acid*3
    amino_acid_stop = amino_acid*3 + 3
    codon = Jb_nucleotide[amino_acid_start:amino_acid_stop]
    Jb_AA += codon_table[codon]
#print(Jb_AA)


# Find overlap of CDR3 with variable and joining (alpha)
# Start by finding the first amino acid of CDR3 in the Variable sequence
overlap_Va_position = Va_AA.rfind(CDR3a_sequence[:1])
#sys.stderr.write("Initial overlap position: %s\n" % overlap_Va_position)
#sys.stderr.write("Length of Va: %s\n" % len(Va_AA))
#print(Va_AA[overlap_Va_position:])
# Increase overlap until the position changes, no longer near the end of the sequence,
# or the new overlap doesn't exist
size_of_Va_overlap = 1
best_size_of_Va_overlap = 0
best_overlap_Va_position = 0
score_of_Va_overlap = 0
while ( Va_AA.rfind(CDR3a_sequence[:size_of_Va_overlap]) >= 0):
    overlap_Va_position = Va_AA.rfind(CDR3a_sequence[:size_of_Va_overlap])
    #sys.stderr.write("CurrPos: %s, BestPos: %s\n" % (overlap_Va_position, best_overlap_Va_position))
    #sys.stderr.write(Va_AA[overlap_Va_position : overlap_Va_position + size_of_Va_overlap] + "\n")
    number_of_gaps = len(Va_AA) - overlap_Va_position - size_of_Va_overlap
    current_score = size_of_Va_overlap*match_multiplier - number_of_gaps*gap_penalty
    #sys.stderr.write("CurrScore: %s, BestScore: %s\n" % (current_score, score_of_Va_overlap))
    if current_score > score_of_Va_overlap:
        best_overlap_Va_position = overlap_Va_position
        best_size_of_Va_overlap = size_of_Va_overlap
        score_of_Va_overlap = current_score
    size_of_Va_overlap += 1
#sys.stderr.write("BestScore: %s, BestPos: %s\n" % (score_of_Va_overlap, best_overlap_Va_position))
if best_overlap_Va_position is 0:
    sys.stderr.write("Warning: No overlap found between CDR3a and Va\n")
    best_size_of_Va_overlap = 0
    best_overlap_Va_position = len(Va_AA)
size_of_Va_overlap = best_size_of_Va_overlap
overlap_Va_position = best_overlap_Va_position
#print(Va_AA[overlap_Va_position : overlap_Va_position + size_of_Va_overlap] 
#      + ',' + str(size_of_Va_overlap))

      
# Find overlap of CDR3 with variable and joining (alpha)
# Next find the first amino acid of joining in the CDR3 sequence
CDR3a_position = len(CDR3a_sequence) - 1
overlap_Ja_position = Ja_AA.find(CDR3a_sequence[CDR3a_position:])
#print(overlap_Ja_position)
#print(Ja_AA[:overlap_Ja_position + 1])
# Increase overlap until no match is found
size_of_Ja_overlap = 1
best_size_of_Ja_overlap = 0
score_of_Ja_overlap = 0
best_overlap_Ja_position = 0
while Ja_AA.find(CDR3a_sequence[CDR3a_position:]) >= 0:
    overlap_Ja_position = Ja_AA.find(CDR3a_sequence[CDR3a_position:])
    #print("%s,%s" % (best_size_of_Ja_overlap, overlap_Ja_position))
    #print(Ja_AA[overlap_Ja_position : overlap_Ja_position + size_of_Ja_overlap])
    current_score = size_of_Ja_overlap*match_multiplier - overlap_Ja_position*gap_penalty
    #print("%s,%s" % (score_of_Ja_overlap,current_score))
    if current_score > score_of_Ja_overlap:
        best_overlap_Ja_position = overlap_Ja_position
        best_size_of_Ja_overlap = size_of_Ja_overlap
        score_of_Ja_overlap = current_score
    size_of_Ja_overlap += 1
    CDR3a_position = len(CDR3a_sequence) - size_of_Ja_overlap
if score_of_Ja_overlap is 0:
    sys.stderr.write("Warning: No overlap found between CDR3a and Ja\n")
    best_size_of_Ja_overlap = 0
    best_overlap_Ja_position = 0
size_of_Ja_overlap = best_size_of_Ja_overlap
CDR3a_position = len(CDR3a_sequence) - size_of_Ja_overlap
overlap_Ja_position = best_overlap_Ja_position
#print(Ja_AA[overlap_Ja_position : overlap_Ja_position + size_of_Ja_overlap] 
#      + ',' + str(size_of_Ja_overlap))
    
# Find overlap of CDR3 with variable and joining (beta)
# Start by finding the first amino acid of CDR3 in the Variable sequence
overlap_Vb_position = Vb_AA.rfind(CDR3b_sequence[:1])
#print(overlap_Vb_position)
#print(Vb_AA[overlap_Vb_position:])
# Increase overlap until the position changes, no longer near the end of the sequence,
# or the new overlap doesn't exist
size_of_Vb_overlap = 1
best_size_of_Vb_overlap = 0
score_of_Vb_overlap = 0
best_overlap_Vb_position = 0
while ( Vb_AA.rfind(CDR3b_sequence[:size_of_Vb_overlap]) >= 0):
    overlap_Vb_position = Vb_AA.rfind(CDR3b_sequence[:size_of_Vb_overlap])
    #print("CurrPos: %s, BestPos: %s" % (overlap_Vb_position, best_overlap_Vb_position))
    #print(Vb_AA[overlap_Vb_position : overlap_Vb_position + size_of_Vb_overlap])
    number_of_gaps = len(Vb_AA) - overlap_Vb_position - size_of_Vb_overlap
    current_score = size_of_Vb_overlap*match_multiplier - number_of_gaps*gap_penalty
    #print("Len_Vb: %s, OverlapVbPos: %s, OverlapSize: %s" % (len(Vb_AA), overlap_Vb_position, size_of_Vb_overlap))
    #print("NumMatches: %s, NumGaps: %s" % (size_of_Vb_overlap, number_of_gaps))
    #print("CurrScore: %s, BestScore: %s\n" % (current_score, score_of_Vb_overlap))
    if current_score > score_of_Vb_overlap:
        best_overlap_Vb_position = overlap_Vb_position
        best_size_of_Vb_overlap = size_of_Vb_overlap
        score_of_Vb_overlap = current_score
    size_of_Vb_overlap += 1
if best_overlap_Vb_position is 0:
    sys.stderr.write("Warning: No overlap found between CDR3b and Vb\n")
    best_size_of_Vb_overlap = 0
    best_overlap_Vb_position = len(Vb_AA)
size_of_Vb_overlap = best_size_of_Vb_overlap
overlap_Vb_position = best_overlap_Vb_position
#print(Vb_AA[overlap_Vb_position : overlap_Vb_position + size_of_Vb_overlap] 
#      + ',' + str(size_of_Vb_overlap))

      
# Find overlap of CDR3 with variable and joining (beta)
# Next find the first amino acid of joining in the CDR3 sequence
CDR3b_position = len(CDR3b_sequence) - 1
overlap_Jb_position = Jb_AA.find(CDR3b_sequence[CDR3b_position:])
#print(overlap_Jb_position)
#print(Jb_AA[:overlap_Jb_position + 1])
# Increase overlap until no match is found
size_of_Jb_overlap = 1
best_size_of_Jb_overlap = 0
score_of_Jb_overlap = 0
best_overlap_Jb_position = 0
while Jb_AA.find(CDR3b_sequence[CDR3b_position:]) >= 0:
    overlap_Jb_position = Jb_AA.find(CDR3b_sequence[CDR3b_position:])
    #print("%s,%s" % (best_size_of_Jb_overlap, overlap_Jb_position))
    #print(Jb_AA[overlap_Jb_position : overlap_Jb_position + size_of_Jb_overlap])
    current_score = size_of_Jb_overlap*match_multiplier - overlap_Jb_position*gap_penalty
    #print("%s,%s" % (score_of_Jb_overlap,current_score))
    if current_score > score_of_Jb_overlap:
        best_overlap_Jb_position = overlap_Jb_position
        best_size_of_Jb_overlap = size_of_Jb_overlap
        score_of_Jb_overlap = current_score
    size_of_Jb_overlap += 1
    CDR3b_position = len(CDR3b_sequence) - size_of_Jb_overlap
if score_of_Jb_overlap is 0:
    sys.stderr.write("Warning: No overlap found between CDR3b and Jb\n")
    best_size_of_Jb_overlap = 0
    best_overlap_Jb_position = 0
size_of_Jb_overlap = best_size_of_Jb_overlap
CDR3b_position = len(CDR3b_sequence) - size_of_Jb_overlap
overlap_Jb_position = best_overlap_Jb_position
#print(Jb_AA[overlap_Jb_position : overlap_Jb_position + size_of_Jb_overlap] 
#      + ',' + str(size_of_Jb_overlap))


# Back-translate CDR3, only non-overlapping region,
# keep all natural genomic nucleotides possible
if ('X' in CDR3a_sequence) or ('X' in CDR3b_sequence):
    sys.stderr.write("Warning: \"X\" in CDR3 will be converted to \"N\" after alignment\n")
CDR3a_nucleotide = ""
for current_AA in CDR3a_sequence:
    if current_AA not in AA_table:
        sys.stderr.write("Error: Amino Acid \"%s\" in CDR3a is invalid\n" % current_AA)
        exit()
    CDR3a_nucleotide += AA_table[current_AA]
#print(CDR3a_sequence + ": " + CDR3a_nucleotide)
    
CDR3b_nucleotide = ""
for current_AA in CDR3b_sequence:
    if current_AA not in AA_table:
        sys.stderr.write("Error: Amino Acid \"%s\" in CDR3b is invalid\n" % current_AA)
        exit()
    CDR3b_nucleotide += AA_table[current_AA]
#print(CDR3b_sequence + ": " + CDR3b_nucleotide)


# Print Alpha overlaps
if len(Va_AA) - best_overlap_Va_position > 25:
    sys.stderr.write("Warning: Overlap occurs over 25 AA from the end of Va\n")
    Va_end_sequence = Va_AA[best_overlap_Va_position:]
    Va_gap = ""
    CDR3a_gap_stop = 0
    CDR3a_gap = ""
else:
    Va_gap = "                         "
    Va_end_sequence = Va_AA[-25:]
    Va_gap_stop = 25 - len(Va_end_sequence)
    Va_gap = Va_gap[:Va_gap_stop]
    CDR3a_gap = "                                                  "
    CDR3a_gap_stop = 25 - len(Va_AA) + overlap_Va_position
    CDR3a_gap = CDR3a_gap[:CDR3a_gap_stop]
if overlap_Ja_position > 25:
    sys.stderr.write("Warning: Overlap occurs over 25 AA from the beginning of Ja\n")
    Ja_beginning_sequence = Ja_AA
    Ja_gap = ""
else:
    Ja_gap = "                                                                                "
    Ja_beginning_sequence = Ja_AA
    if len(Ja_AA) > 25:
        Ja_beginning_sequence = Ja_AA[:25]
    Ja_gap_stop = CDR3a_gap_stop + CDR3a_position - overlap_Ja_position
    Ja_gap = Ja_gap[:Ja_gap_stop]
sys.stderr.write("Alpha Overlap:\n")
sys.stderr.write("Va:    " + Va_end_sequence + "\n")
sys.stderr.write("CDR3a: " + CDR3a_gap + CDR3a_sequence.replace("X", "N") + "\n")
sys.stderr.write("Ja:    " + Ja_gap + Ja_beginning_sequence + "\n\n")

# Print Beta overlaps
if len(Vb_AA) - best_overlap_Vb_position > 25:
    sys.stderr.write("Warning: Overlap occurs over 25 AA from the end of Vb\n")
    Vb_end_sequence = Vb_AA[best_overlap_Vb_position:]
    Vb_gap = ""
    CDR3b_gap_stop = 0
    CDR3b_gap = ""
else:
    Vb_gap = "                         "
    Vb_end_sequence = Vb_AA[-25:]
    Vb_gap_stop = 25 - len(Vb_end_sequence)
    Vb_gap = Vb_gap[:Vb_gap_stop]
    CDR3b_gap = "                                                  "
    CDR3b_gap_stop = 25 - len(Vb_AA) + overlap_Vb_position
    CDR3b_gap = CDR3b_gap[:CDR3b_gap_stop]
if overlap_Jb_position > 25:
    sys.stderr.write("Warning: Overlap occurs over 25 AA from the beginning of Jb\n")
    Jb_beginning_sequence = Jb_AA
    Jb_gap = ""
else:
    Jb_gap = "                                                                                "
    Jb_beginning_sequence = Jb_AA
    if len(Jb_AA) > 25:
        Jb_beginning_sequence = Jb_AA[:25]
    Jb_gap_stop = CDR3b_gap_stop + CDR3b_position - overlap_Jb_position
    Jb_gap = Jb_gap[:Jb_gap_stop]
sys.stderr.write("Beta Overlap:\n")
sys.stderr.write("Vb:    " + Vb_end_sequence + "\n")
sys.stderr.write("CDR3b: " + CDR3b_gap + CDR3b_sequence.replace("X", "N") + "\n")
sys.stderr.write("Jb:    " + Jb_gap + Jb_beginning_sequence + "\n\n")
 
    
# Trim out-of-frame bases and concatenate sequences 
Va_nucleotide_start = 0
Va_nucleotide_stop = overlap_Va_position*3 + size_of_Va_overlap*3
Va_nucleotide_trimmed = Va_nucleotide[Va_nucleotide_start:Va_nucleotide_stop]

Ja_nucleotide_start = overlap_Ja_position*3
Ja_nucleotide_stop = len(Ja_AA)*3
Ja_nucleotide_trimmed = Ja_nucleotide[Ja_nucleotide_start:Ja_nucleotide_stop]

CDR3a_nucleotide_start = size_of_Va_overlap*3
CDR3a_nucleotide_stop = len(CDR3a_sequence)*3 - size_of_Ja_overlap*3
CDR3a_nucleotide_trimmed = CDR3a_nucleotide[CDR3a_nucleotide_start:CDR3a_nucleotide_stop]

Vb_nucleotide_start = 0
Vb_nucleotide_stop = overlap_Vb_position*3 + size_of_Vb_overlap*3
Vb_nucleotide_trimmed = Vb_nucleotide[Vb_nucleotide_start:Vb_nucleotide_stop]

Jb_nucleotide_start = overlap_Jb_position*3
Jb_nucleotide_stop = len(Jb_AA)*3
Jb_nucleotide_trimmed = Jb_nucleotide[Jb_nucleotide_start:Jb_nucleotide_stop]

CDR3b_nucleotide_start = size_of_Vb_overlap*3
CDR3b_nucleotide_stop = len(CDR3b_sequence)*3 - size_of_Jb_overlap*3
CDR3b_nucleotide_trimmed = CDR3b_nucleotide[CDR3b_nucleotide_start:CDR3b_nucleotide_stop]

#print(Va_nucleotide_trimmed)
#print(CDR3a_nucleotide_trimmed)
#print(Ja_nucleotide_trimmed)
#print(Vb_nucleotide_trimmed)
#print(CDR3b_nucleotide_trimmed)
#print(Jb_nucleotide_trimmed)

alpha_nucleotide = (Va_nucleotide_trimmed + CDR3a_nucleotide_trimmed + 
                    Ja_nucleotide_trimmed + TRAC_nucleotide)
beta_nucleotide = (Vb_nucleotide_trimmed + CDR3b_nucleotide_trimmed + 
                   Jb_nucleotide_trimmed + TRBC_nucleotide)
sys.stderr.write(">TCR_Alpha_nuc\n" + alpha_nucleotide + "\n\n")
sys.stderr.write(">TCR_Beta_nuc\n" + beta_nucleotide + "\n\n")








###############################################################################
# Optional: Codon optimize sequence for mouse or human
# 

# Nucleotide to amino acid table
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

# Most frequent codon for each amino acid
# Use "N" for "X", as occurs in the human constant region on IMGT
AA_table = {
'A':'GCC', 'C':'TGC', 'D':'GAC', 'E':'GAG', 
'F':'TTC', 'G':'GGC', 'H':'CAC', 'I':'ATC', 
'K':'AAG', 'L':'CTG', 'M':'ATG', 'N':'AAC', 
'P':'CCC', 'Q':'CAG', 'R':'AGG', 'S':'AGC', 
'T':'ACC', 'V':'GTG', 'W':'TGG', 'Y':'TAC',
'X':'AAC'}

# Human codon optimization table
# IDT removes any codons less than 10% frequent 
human_opt_table = {
'*':[(1,"TGA")],
'A':[(0.11,"GCG"),(0.34,"GCA"),(0.6,"GCT"),(1,"GCC")],
'C':[(0.45,"TGT"),(1,"TGC")],
'D':[(0.46,"GAT"),(1,"GAC")],
'E':[(0.42,"GAA"),(1,"GAG")],
'F':[(0.45,"TTT"),(1,"TTC")],
'G':[(0.16,"GGT"),(0.41,"GGA"),(0.66,"GGG"),(1,"GGC")],
'H':[(0.41,"CAT"),(1,"CAC")],
'I':[(0.16,"ATA"),(0.52,"ATT"),(1,"ATC")],
'K':[(0.42,"AAA"),(1,"AAG")],
'L':[(0.13,"TTG"),(0.26,"CTT"),(0.46,"CTC"),(1,"CTG")],
'M':[(1,"ATG")],
'N':[(0.46,"AAT"),(1,"AAC")],
'P':[(0.11,"CCG"),(0.38,"CCA"),(0.66,"CCT"),(1,"CCC")],
'Q':[(0.25,"CAA"),(1,"CAG")],
'R':[(0.11,"CGA"),(0.3,"CGC"),(0.5,"AGA"),(0.7,"AGG"),(1,"CGG")],
'S':[(0.15,"TCA"),(0.3,"AGT"),(0.48,"TCT"),(0.7,"TCC"),(1,"AGC")],
'T':[(0.12,"ACG"),(0.36,"ACT"),(0.64,"ACA"),(1,"ACC")],
'V':[(0.11,"GTA"),(0.29,"GTT"),(0.53,"GTC"),(1,"GTG")],
'W':[(1,"TGG")],
'Y':[(0.43,"TAT"),(1,"TAC")],
'X':[(1,"AAC")],
'.':[(1,"TGA")]
}

# Mouse codon optimization table
# IDT removes any codons less than 10% frequent 
mouse_opt_table = {
'*':[(1,"TGA")],
'A':[(0.1,"GCG"),(0.33,"GCA"),(0.62,"GCT"),(1,"GCC")],
'C':[(0.48,"TGT"),(1,"TGC")],
'D':[(0.44,"GAT"),(1,"GAC")],
'E':[(0.4,"GAA"),(1,"GAG")],
'F':[(0.43,"TTT"),(1,"TTC")],
'G':[(0.18,"GGT"),(0.41,"GGG"),(0.67,"GGA"),(1,"GGC")],
'H':[(0.4,"CAT"),(1,"CAC")],
'I':[(0.16,"ATA"),(0.5,"ATT"),(1,"ATC")],
'K':[(0.39,"AAA"),(1,"AAG")],
'L':[(0.13,"TTG"),(0.26,"CTT"),(0.46,"CTC"),(1,"CTG")],
'M':[(1,"ATG")],
'N':[(0.43,"AAT"),(1,"AAC")],
'P':[(0.1,"CCG"),(0.38,"CCA"),(0.68,"CCT"),(1,"CCC")],
'Q':[(0.25,"CAA"),(1,"CAG")],
'R':[(0.12,"CGA"),(0.3,"CGC"),(0.49,"CGG"),(0.7,"AGA"),(1,"AGG")],
'S':[(0.14,"TCA"),(0.29,"AGT"),(0.48,"TCT"),(0.7,"TCC"),(1,"AGC")],
'T':[(0.11,"ACG"),(0.36,"ACT"),(0.65,"ACA"),(1,"ACC")],
'V':[(0.12,"GTA"),(0.29,"GTT"),(0.54,"GTC"),(1,"GTG")],
'W':[(1,"TGG")],
'Y':[(0.43,"TAT"),(1,"TAC")],
'X':[(1,"AAC")],
'.':[(1,"TGA")]
}

# Translate to amino acid
alpha_AA = "" 
alpha_nucleotide_length = len(alpha_nucleotide)   
for amino_acid in range(int(alpha_nucleotide_length/3)):
    amino_acid_start = amino_acid*3
    amino_acid_stop = amino_acid*3 + 3
    codon = alpha_nucleotide[amino_acid_start:amino_acid_stop]
    alpha_AA += codon_table[codon]
#sys.stderr.write("alpha Amino Acid Translation: " + alpha_AA + "\n")

# Allow flag for codon optimization
if perform_optimization:
    # Reconvert to optimized nucleotide
    alpha_optimized = ""
    for current_amino_acid in alpha_AA:
        # Weighted codon selection:
        # Assign each codon a segment of the interval 0 to 1 according to the codon frequency
        # https://www.genscript.com/tools/codon-frequency-table
        # Also ignore any codons that have a frequency of <0.1 (rare codons)
        #
        # Pick a random number, select the codon that corresponds to that segment
        #
        # A: 0[GCG]0.1[-----GCA-----]0.33[-------GCT-------]0.62[---*---GCC-------]1.0
        #                                     random_number = 0.72  ^  codon is "GCC"
        #
        random_number = random.random()
        segment_start = 0
        segment_stop = 0
        optimized_codon = "NNN"
        for tuple in human_opt_table[current_amino_acid]:
            segment_start = segment_stop
            segment_stop = tuple[0]
            current_codon = tuple[1]
            if (random_number > segment_start) and (random_number <= segment_stop):
                optimized_codon = current_codon
        #sys.stderr.write("Amino Acid: " + current_amino_acid + ", Selected codon: " + optimized_codon +  ", Random number: " + str(random_number) + "\n")
        alpha_optimized += optimized_codon
    #sys.stderr.write("alpha Codon Optimized : " + alpha_optimized + "\n")

    alpha_nucleotide = alpha_optimized

sys.stderr.write("Optimized alpha nucleotide sequence: \n" + alpha_nucleotide + "\n")

###############################################################################
# Optimize repeated strings
# 
min_substring_length = 5
long_substring_length = 8

sys.stderr.write("Performing repeat optimization: " + "\n")
sys.stderr.write("Optimize repeats of " + str(min_substring_length) + "bp or longer, " + "\n")
sys.stderr.write("stopping when repeats of length " + str(long_substring_length) + "bp no longer exist" + "\n")

minimal_repeat_count = -1
minimal_repeat_sequence = ""

current_loop = 0
substring_counts = {}
substring_length = min_substring_length
codon_optimize_list = []
long_repeats_found = True

while long_repeats_found and (current_loop < 100):
    long_repeats_found = False
    substring_counts = {}
    substring_length = min_substring_length
    codon_optimize_list = []
    while substring_length < alpha_nucleotide_length:
        for substring_position in range(alpha_nucleotide_length - substring_length + 1):
            substring_position_stop = substring_position + substring_length
            substring = alpha_nucleotide[substring_position:substring_position_stop]
            if substring in substring_counts:
                substring_counts[substring] += 1
                codon_start = int(substring_position/3) * 3
                if codon_start not in codon_optimize_list:
                    codon_optimize_list.append(codon_start)
            else:
                substring_counts[substring] = 1
        substring_length += 1

    #sys.stderr.write("\n" + str(codon_optimize_list) + "\n")
        
    current_repeat_count = 0
    for repeat in substring_counts:
        if (substring_counts[repeat] > 1) and (len(repeat) >= long_substring_length):
            #print("Warning: Repeat found \"" + repeat +  "\" of length " + str(len(repeat)) + " repeated " + str(substring_counts[repeat]) + " times.\n")
            long_repeats_found = True
            current_repeat_count += 1
    
    if (current_repeat_count < minimal_repeat_count) or (minimal_repeat_count < 0):
        minimal_repeat_count = current_repeat_count
        minimal_repeat_sequence = alpha_nucleotide
    
    for current_position in codon_optimize_list:
        current_repeat = alpha_nucleotide[current_position:current_position+3]
        current_repeat_AA = codon_table[current_repeat]
        # Weighted codon selection:
        # Assign each codon a segment of the interval 0 to 1 according to the codon frequency
        # https://www.genscript.com/tools/codon-frequency-table
        # Also ignore any codons that have a frequency of <0.1 (rare codons)
        #
        # Pick a random number, select the codon that corresponds to that segment
        #
        # A: 0[GCG]0.1[-----GCA-----]0.33[-------GCT-------]0.62[---*---GCC-------]1.0
        #                                     random_number = 0.72  ^  codon is "GCC"
        #
        random_number = random.random()
        segment_start = 0
        segment_stop = 0
        optimized_codon = "NNN"
        for tuple in human_opt_table[current_repeat_AA]:
            segment_start = segment_stop
            segment_stop = tuple[0]
            current_codon = tuple[1]
            if (random_number > segment_start) and (random_number <= segment_stop):
                optimized_codon = current_codon
        #sys.stderr.write("Amino Acid: " + current_repeat_AA + ", Selected codon: " + optimized_codon +  ", Random number: " + str(random_number) + "\n")
        alpha_nucleotide = alpha_nucleotide[:current_position] + optimized_codon + alpha_nucleotide[current_position+3:]
    
    #sys.stderr.write("alpha nucleotide sequence optimized to minimize repeats : " + alpha_nucleotide + "\n")
    current_loop += 1
    
    if current_loop == 100:
        sys.stderr.write("Warning: Aborting repeat optimization, " + "\n")
        sys.stderr.write("Reverting to a sequence with " + str(minimal_repeat_count) + " repeats" + "\n")
        alpha_nucleotide = minimal_repeat_sequence
    
sys.stderr.write("alpha nucleotide sequence optimized to minimize repeats: " + alpha_nucleotide + "\n")


###############################################################################
# Optimize strings of a single nucleotide sequence
# 

nucleotide_strings = {"AAAA":0, "TTTT":0, "GGGG":0, "CCCC":0}
position_counts = {}

current_loop = 0
current_string_count = 0
strings_exist = True
while strings_exist and (current_loop < 1000):
    current_string_count = 0
    for current_string in nucleotide_strings:
        current_string_blacklist = nucleotide_strings[current_string]
        # First recognize sites of low complexity  
        if current_string in alpha_nucleotide[(current_string_blacklist + 1):]:
            # Make sure string isn't causing infinite loop ("WG" amino acid sequence, etc.)
            current_string_position = alpha_nucleotide.find(current_string, current_string_blacklist + 1)
            #sys.stderr.write("Warning: Low complexity sequence of \"" + current_string + "\" found in alpha nucleotide sequence at position " + str(current_string_position)  + "\n")
            
            if current_string_position in position_counts:
                position_counts[current_string_position] += 1
            else:
                position_counts[current_string_position] = 1
            
            if position_counts[current_string_position] >= 30:
                sys.stderr.write("Warning: Aborting optimization of low complexity sequence \"" + current_string + "\" found at position " + str(current_string_position) + "\n")
                nucleotide_strings[current_string] = current_string_position
            
            # Then assign replacement codons around the site, assigned by frequency table
            current_string_count += 1
            string_AA_position_start = int(current_string_position / 3)
            string_AA_position_stop = string_AA_position_start + 2
            replaced_codons = ""
            for current_string_AA_position in range(string_AA_position_start, string_AA_position_stop):
                current_string_AA = alpha_AA[current_string_AA_position]
                # Weighted codon selection:
                # Assign each codon a segment of the interval 0 to 1 according to the codon frequency
                # https://www.genscript.com/tools/codon-frequency-table
                # Also ignore any codons that have a frequency of <0.1 (rare codons)
                #
                # Pick a random number, select the codon that corresponds to that segment
                #
                # A: 0[GCG]0.1[-----GCA-----]0.33[-------GCT-------]0.62[---*---GCC-------]1.0
                #                                     random_number = 0.72  ^  codon is "GCC"
                #
                random_number = random.random()
                segment_start = 0
                segment_stop = 0
                optimized_codon = "NNN"
                for tuple in human_opt_table[current_string_AA]:
                    segment_start = segment_stop
                    segment_stop = tuple[0]
                    current_codon = tuple[1]
                    if (random_number > segment_start) and (random_number <= segment_stop):
                        optimized_codon = current_codon
                #sys.stderr.write("Amino Acid: " + current_string_AA + ", Selected codon: " + optimized_codon +  ", Random number: " + str(random_number) + "\n")
                replaced_codons += optimized_codon
            alpha_nucleotide = alpha_nucleotide[:string_AA_position_start*3] + replaced_codons + alpha_nucleotide[string_AA_position_stop*3:]
    #sys.stderr.write("alpha nucleotide sequence edited for low complexity sites : " + alpha_nucleotide + "\n")
    current_loop += 1
    if current_string_count == 0:
        strings_exist = False
        
sys.stderr.write("alpha nucleotide sequence edited for low complexity sites: \n" + alpha_nucleotide + "\n")    


###############################################################################
# Optional: Remove restriction sites
# 

restriction_sites = [
("BsmBI","CGTCTC"),
("BsmBIrc","GAGACG"),
("EcoRI","GAATTC"),
("HpaI","GTTAAC"),
("XhoI","CTCGAG"),
("BamHI","GGATCC"),
("XmaI","CCCGGG"),
("SalI","GTCGAC"),
("KpnI","GGTACC"),
("XbaI","TCTAGA"),
("NaeI","GCCGGC"),
("NdeI","CATATG"),
]

current_loop = 0
current_restriction_count = 0
restriction_sites_exist = True
while restriction_sites_exist and (current_loop < 100):
    current_restriction_count = 0
    for current_enzyme in restriction_sites:
        current_enzyme_name = current_enzyme[0]
        current_enzyme_sequence = current_enzyme[1]
        # First recognize sites 
        if current_enzyme_sequence in alpha_nucleotide:
            current_restriction_count += 1
            sys.stderr.write("Warning: Restriction sequence of " + current_enzyme_name + " found in alpha nucleotide sequence at position " + str(alpha_nucleotide.find(current_enzyme_sequence)) + " and has been edited" + "\n")

            # Then assign replacement codons around the site, assigned by frequency table
            current_enzyme_position = alpha_nucleotide.find(current_enzyme_sequence)
            enzyme_AA_position_start = int(current_enzyme_position / 3)
            enzyme_AA_position_stop = enzyme_AA_position_start + 3
            if current_enzyme_position % 3 == 0:
                enzyme_AA_position_stop = enzyme_AA_position_start + 2
            replaced_codons = ""
            for current_enzyme_AA_position in range(enzyme_AA_position_start, enzyme_AA_position_stop):
                current_enzyme_AA = alpha_AA[current_enzyme_AA_position]
                # Weighted codon selection:
                # Assign each codon a segment of the interval 0 to 1 according to the codon frequency
                # https://www.genscript.com/tools/codon-frequency-table
                # Also ignore any codons that have a frequency of <0.1 (rare codons)
                #
                # Pick a random number, select the codon that corresponds to that segment
                #
                # A: 0[GCG]0.1[-----GCA-----]0.33[-------GCT-------]0.62[---*---GCC-------]1.0
                #                                     random_number = 0.72  ^  codon is "GCC"
                #
                random_number = random.random()
                segment_start = 0
                segment_stop = 0
                optimized_codon = "NNN"
                for tuple in human_opt_table[current_enzyme_AA]:
                    segment_start = segment_stop
                    segment_stop = tuple[0]
                    current_codon = tuple[1]
                    if (random_number > segment_start) and (random_number <= segment_stop):
                        optimized_codon = current_codon
                #sys.stderr.write("Amino Acid: " + current_enzyme_AA + ", Selected codon: " + optimized_codon +  ", Random number: " + str(random_number) + "\n")
                replaced_codons += optimized_codon
            alpha_nucleotide = alpha_nucleotide[:enzyme_AA_position_start*3] + replaced_codons + alpha_nucleotide[enzyme_AA_position_stop*3:]
    #sys.stderr.write("alpha nucleotide sequence edited for restriction sites : " + alpha_nucleotide + "\n")
    current_loop += 1
    if current_restriction_count == 0:
        restriction_sites_exist = False
        
#sys.stderr.write("alpha nucleotide sequence edited for restriction sites: \n" + alpha_nucleotide + "\n")








###############################################################################
# Optional: Codon optimize sequence for mouse or human
# 

# Nucleotide to amino acid table
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

# Most frequent codon for each amino acid
# Use "N" for "X", as occurs in the human constant region on IMGT
AA_table = {
'A':'GCC', 'C':'TGC', 'D':'GAC', 'E':'GAG', 
'F':'TTC', 'G':'GGC', 'H':'CAC', 'I':'ATC', 
'K':'AAG', 'L':'CTG', 'M':'ATG', 'N':'AAC', 
'P':'CCC', 'Q':'CAG', 'R':'AGG', 'S':'AGC', 
'T':'ACC', 'V':'GTG', 'W':'TGG', 'Y':'TAC',
'X':'AAC'}

# Human codon optimization table
# IDT removes any codons less than 10% frequent 
human_opt_table = {
'*':[(1,"TGA")],
'A':[(0.11,"GCG"),(0.34,"GCA"),(0.6,"GCT"),(1,"GCC")],
'C':[(0.45,"TGT"),(1,"TGC")],
'D':[(0.46,"GAT"),(1,"GAC")],
'E':[(0.42,"GAA"),(1,"GAG")],
'F':[(0.45,"TTT"),(1,"TTC")],
'G':[(0.16,"GGT"),(0.41,"GGA"),(0.66,"GGG"),(1,"GGC")],
'H':[(0.41,"CAT"),(1,"CAC")],
'I':[(0.16,"ATA"),(0.52,"ATT"),(1,"ATC")],
'K':[(0.42,"AAA"),(1,"AAG")],
'L':[(0.13,"TTG"),(0.26,"CTT"),(0.46,"CTC"),(1,"CTG")],
'M':[(1,"ATG")],
'N':[(0.46,"AAT"),(1,"AAC")],
'P':[(0.11,"CCG"),(0.38,"CCA"),(0.66,"CCT"),(1,"CCC")],
'Q':[(0.25,"CAA"),(1,"CAG")],
'R':[(0.11,"CGA"),(0.3,"CGC"),(0.5,"AGA"),(0.7,"AGG"),(1,"CGG")],
'S':[(0.15,"TCA"),(0.3,"AGT"),(0.48,"TCT"),(0.7,"TCC"),(1,"AGC")],
'T':[(0.12,"ACG"),(0.36,"ACT"),(0.64,"ACA"),(1,"ACC")],
'V':[(0.11,"GTA"),(0.29,"GTT"),(0.53,"GTC"),(1,"GTG")],
'W':[(1,"TGG")],
'Y':[(0.43,"TAT"),(1,"TAC")],
'X':[(1,"AAC")],
'.':[(1,"TGA")]
}

# Mouse codon optimization table
# IDT removes any codons less than 10% frequent 
mouse_opt_table = {
'*':[(1,"TGA")],
'A':[(0.1,"GCG"),(0.33,"GCA"),(0.62,"GCT"),(1,"GCC")],
'C':[(0.48,"TGT"),(1,"TGC")],
'D':[(0.44,"GAT"),(1,"GAC")],
'E':[(0.4,"GAA"),(1,"GAG")],
'F':[(0.43,"TTT"),(1,"TTC")],
'G':[(0.18,"GGT"),(0.41,"GGG"),(0.67,"GGA"),(1,"GGC")],
'H':[(0.4,"CAT"),(1,"CAC")],
'I':[(0.16,"ATA"),(0.5,"ATT"),(1,"ATC")],
'K':[(0.39,"AAA"),(1,"AAG")],
'L':[(0.13,"TTG"),(0.26,"CTT"),(0.46,"CTC"),(1,"CTG")],
'M':[(1,"ATG")],
'N':[(0.43,"AAT"),(1,"AAC")],
'P':[(0.1,"CCG"),(0.38,"CCA"),(0.68,"CCT"),(1,"CCC")],
'Q':[(0.25,"CAA"),(1,"CAG")],
'R':[(0.12,"CGA"),(0.3,"CGC"),(0.49,"CGG"),(0.7,"AGA"),(1,"AGG")],
'S':[(0.14,"TCA"),(0.29,"AGT"),(0.48,"TCT"),(0.7,"TCC"),(1,"AGC")],
'T':[(0.11,"ACG"),(0.36,"ACT"),(0.65,"ACA"),(1,"ACC")],
'V':[(0.12,"GTA"),(0.29,"GTT"),(0.54,"GTC"),(1,"GTG")],
'W':[(1,"TGG")],
'Y':[(0.43,"TAT"),(1,"TAC")],
'X':[(1,"AAC")],
'.':[(1,"TGA")]
}

# Translate to amino acid
beta_AA = "" 
beta_nucleotide_length = len(beta_nucleotide)   
for amino_acid in range(int(beta_nucleotide_length/3)):
    amino_acid_start = amino_acid*3
    amino_acid_stop = amino_acid*3 + 3
    codon = beta_nucleotide[amino_acid_start:amino_acid_stop]
    beta_AA += codon_table[codon]
#sys.stderr.write("beta Amino Acid Translation: " + beta_AA + "\n")

# Allow flag for codon optimization
if perform_optimization:
    # Reconvert to optimized nucleotide
    beta_optimized = ""
    for current_amino_acid in beta_AA:
        # Weighted codon selection:
        # Assign each codon a segment of the interval 0 to 1 according to the codon frequency
        # https://www.genscript.com/tools/codon-frequency-table
        # Also ignore any codons that have a frequency of <0.1 (rare codons)
        #
        # Pick a random number, select the codon that corresponds to that segment
        #
        # A: 0[GCG]0.1[-----GCA-----]0.33[-------GCT-------]0.62[---*---GCC-------]1.0
        #                                     random_number = 0.72  ^  codon is "GCC"
        #
        random_number = random.random()
        segment_start = 0
        segment_stop = 0
        optimized_codon = "NNN"
        for tuple in human_opt_table[current_amino_acid]:
            segment_start = segment_stop
            segment_stop = tuple[0]
            current_codon = tuple[1]
            if (random_number > segment_start) and (random_number <= segment_stop):
                optimized_codon = current_codon
        #sys.stderr.write("Amino Acid: " + current_amino_acid + ", Selected codon: " + optimized_codon +  ", Random number: " + str(random_number) + "\n")
        beta_optimized += optimized_codon
    #sys.stderr.write("beta Codon Optimized : " + beta_optimized + "\n")

    beta_nucleotide = beta_optimized

sys.stderr.write("Optimized beta nucleotide sequence: \n" + beta_nucleotide + "\n")

###############################################################################
# Optimize repeated strings
# 
min_substring_length = 5
long_substring_length = 8

sys.stderr.write("Performing repeat optimization: " + "\n")
sys.stderr.write("Optimize repeats of " + str(min_substring_length) + "bp or longer, " + "\n")
sys.stderr.write("stopping when repeats of length " + str(long_substring_length) + "bp no longer exist" + "\n")

minimal_repeat_count = -1
minimal_repeat_sequence = ""

current_loop = 0
substring_counts = {}
substring_length = min_substring_length
codon_optimize_list = []
long_repeats_found = True

while long_repeats_found and (current_loop < 100):
    long_repeats_found = False
    substring_counts = {}
    substring_length = min_substring_length
    codon_optimize_list = []
    while substring_length < beta_nucleotide_length:
        for substring_position in range(beta_nucleotide_length - substring_length + 1):
            substring_position_stop = substring_position + substring_length
            substring = beta_nucleotide[substring_position:substring_position_stop]
            if substring in substring_counts:
                substring_counts[substring] += 1
                codon_start = int(substring_position/3) * 3
                if codon_start not in codon_optimize_list:
                    codon_optimize_list.append(codon_start)
            else:
                substring_counts[substring] = 1
        substring_length += 1

    #sys.stderr.write("\n" + str(codon_optimize_list) + "\n")
        
    current_repeat_count = 0
    for repeat in substring_counts:
        if (substring_counts[repeat] > 1) and (len(repeat) >= long_substring_length):
            #print("Warning: Repeat found \"" + repeat +  "\" of length " + str(len(repeat)) + " repeated " + str(substring_counts[repeat]) + " times.\n")
            long_repeats_found = True
            current_repeat_count += 1
    
    if (current_repeat_count < minimal_repeat_count) or (minimal_repeat_count < 0):
        minimal_repeat_count = current_repeat_count
        minimal_repeat_sequence = beta_nucleotide
    
    for current_position in codon_optimize_list:
        current_repeat = beta_nucleotide[current_position:current_position+3]
        current_repeat_AA = codon_table[current_repeat]
        # Weighted codon selection:
        # Assign each codon a segment of the interval 0 to 1 according to the codon frequency
        # https://www.genscript.com/tools/codon-frequency-table
        # Also ignore any codons that have a frequency of <0.1 (rare codons)
        #
        # Pick a random number, select the codon that corresponds to that segment
        #
        # A: 0[GCG]0.1[-----GCA-----]0.33[-------GCT-------]0.62[---*---GCC-------]1.0
        #                                     random_number = 0.72  ^  codon is "GCC"
        #
        random_number = random.random()
        segment_start = 0
        segment_stop = 0
        optimized_codon = "NNN"
        for tuple in human_opt_table[current_repeat_AA]:
            segment_start = segment_stop
            segment_stop = tuple[0]
            current_codon = tuple[1]
            if (random_number > segment_start) and (random_number <= segment_stop):
                optimized_codon = current_codon
        #sys.stderr.write("Amino Acid: " + current_repeat_AA + ", Selected codon: " + optimized_codon +  ", Random number: " + str(random_number) + "\n")
        beta_nucleotide = beta_nucleotide[:current_position] + optimized_codon + beta_nucleotide[current_position+3:]
    
    #sys.stderr.write("beta nucleotide sequence optimized to minimize repeats : " + beta_nucleotide + "\n")
    current_loop += 1
    
    if current_loop == 100:
        sys.stderr.write("Warning: Aborting repeat optimization, " + "\n")
        sys.stderr.write("Reverting to a sequence with " + str(minimal_repeat_count) + " repeats" + "\n")
        beta_nucleotide = minimal_repeat_sequence
    
sys.stderr.write("beta nucleotide sequence optimized to minimize repeats: " + beta_nucleotide + "\n")


###############################################################################
# Optimize strings of a single nucleotide sequence
# 

nucleotide_strings = {"AAAA":0, "TTTT":0, "GGGG":0, "CCCC":0}
position_counts = {}

current_loop = 0
current_string_count = 0
strings_exist = True
while strings_exist and (current_loop < 1000):
    current_string_count = 0
    for current_string in nucleotide_strings:
        current_string_blacklist = nucleotide_strings[current_string]
        # First recognize sites of low complexity  
        if current_string in beta_nucleotide[(current_string_blacklist + 1):]:
            # Make sure string isn't causing infinite loop ("WG" amino acid sequence, etc.)
            current_string_position = beta_nucleotide.find(current_string, current_string_blacklist + 1)
            #sys.stderr.write("Warning: Low complexity sequence of \"" + current_string + "\" found in beta nucleotide sequence at position " + str(current_string_position)  + "\n")
            
            if current_string_position in position_counts:
                position_counts[current_string_position] += 1
            else:
                position_counts[current_string_position] = 1
            
            if position_counts[current_string_position] >= 30:
                sys.stderr.write("Warning: Aborting optimization of low complexity sequence \"" + current_string + "\" found at position " + str(current_string_position) + "\n")
                nucleotide_strings[current_string] = current_string_position
            
            # Then assign replacement codons around the site, assigned by frequency table
            current_string_count += 1
            string_AA_position_start = int(current_string_position / 3)
            string_AA_position_stop = string_AA_position_start + 2
            replaced_codons = ""
            for current_string_AA_position in range(string_AA_position_start, string_AA_position_stop):
                current_string_AA = beta_AA[current_string_AA_position]
                # Weighted codon selection:
                # Assign each codon a segment of the interval 0 to 1 according to the codon frequency
                # https://www.genscript.com/tools/codon-frequency-table
                # Also ignore any codons that have a frequency of <0.1 (rare codons)
                #
                # Pick a random number, select the codon that corresponds to that segment
                #
                # A: 0[GCG]0.1[-----GCA-----]0.33[-------GCT-------]0.62[---*---GCC-------]1.0
                #                                     random_number = 0.72  ^  codon is "GCC"
                #
                random_number = random.random()
                segment_start = 0
                segment_stop = 0
                optimized_codon = "NNN"
                for tuple in human_opt_table[current_string_AA]:
                    segment_start = segment_stop
                    segment_stop = tuple[0]
                    current_codon = tuple[1]
                    if (random_number > segment_start) and (random_number <= segment_stop):
                        optimized_codon = current_codon
                #sys.stderr.write("Amino Acid: " + current_string_AA + ", Selected codon: " + optimized_codon +  ", Random number: " + str(random_number) + "\n")
                replaced_codons += optimized_codon
            beta_nucleotide = beta_nucleotide[:string_AA_position_start*3] + replaced_codons + beta_nucleotide[string_AA_position_stop*3:]
    #sys.stderr.write("beta nucleotide sequence edited for low complexity sites : " + beta_nucleotide + "\n")
    current_loop += 1
    if current_string_count == 0:
        strings_exist = False
        
sys.stderr.write("beta nucleotide sequence edited for low complexity sites: \n" + beta_nucleotide + "\n")    


###############################################################################
# Optional: Remove restriction sites
# 

restriction_sites = [
("BsmBI","CGTCTC"),
("BsmBIrc","GAGACG"),
("EcoRI","GAATTC"),
("HpaI","GTTAAC"),
("XhoI","CTCGAG"),
("BamHI","GGATCC"),
("XmaI","CCCGGG"),
("SalI","GTCGAC"),
("KpnI","GGTACC"),
("XbaI","TCTAGA"),
("NaeI","GCCGGC"),
("NdeI","CATATG"),
]

current_loop = 0
current_restriction_count = 0
restriction_sites_exist = True
while restriction_sites_exist and (current_loop < 100):
    current_restriction_count = 0
    for current_enzyme in restriction_sites:
        current_enzyme_name = current_enzyme[0]
        current_enzyme_sequence = current_enzyme[1]
        # First recognize sites 
        if current_enzyme_sequence in beta_nucleotide:
            current_restriction_count += 1
            sys.stderr.write("Warning: Restriction sequence of " + current_enzyme_name + " found in beta nucleotide sequence at position " + str(beta_nucleotide.find(current_enzyme_sequence)) + " and has been edited" + "\n")

            # Then assign replacement codons around the site, assigned by frequency table
            current_enzyme_position = beta_nucleotide.find(current_enzyme_sequence)
            enzyme_AA_position_start = int(current_enzyme_position / 3)
            enzyme_AA_position_stop = enzyme_AA_position_start + 3
            if current_enzyme_position % 3 == 0:
                enzyme_AA_position_stop = enzyme_AA_position_start + 2
            replaced_codons = ""
            for current_enzyme_AA_position in range(enzyme_AA_position_start, enzyme_AA_position_stop):
                current_enzyme_AA = beta_AA[current_enzyme_AA_position]
                # Weighted codon selection:
                # Assign each codon a segment of the interval 0 to 1 according to the codon frequency
                # https://www.genscript.com/tools/codon-frequency-table
                # Also ignore any codons that have a frequency of <0.1 (rare codons)
                #
                # Pick a random number, select the codon that corresponds to that segment
                #
                # A: 0[GCG]0.1[-----GCA-----]0.33[-------GCT-------]0.62[---*---GCC-------]1.0
                #                                     random_number = 0.72  ^  codon is "GCC"
                #
                random_number = random.random()
                segment_start = 0
                segment_stop = 0
                optimized_codon = "NNN"
                for tuple in human_opt_table[current_enzyme_AA]:
                    segment_start = segment_stop
                    segment_stop = tuple[0]
                    current_codon = tuple[1]
                    if (random_number > segment_start) and (random_number <= segment_stop):
                        optimized_codon = current_codon
                #sys.stderr.write("Amino Acid: " + current_enzyme_AA + ", Selected codon: " + optimized_codon +  ", Random number: " + str(random_number) + "\n")
                replaced_codons += optimized_codon
            beta_nucleotide = beta_nucleotide[:enzyme_AA_position_start*3] + replaced_codons + beta_nucleotide[enzyme_AA_position_stop*3:]
    #sys.stderr.write("beta nucleotide sequence edited for restriction sites : " + beta_nucleotide + "\n")
    current_loop += 1
    if current_restriction_count == 0:
        restriction_sites_exist = False
        
#sys.stderr.write("beta nucleotide sequence edited for restriction sites: \n" + beta_nucleotide + "\n")












final_sequence = ( In_Fusion_Left + LNGFR + LNGFR_additional + P2A + HGH_SS_1 +
                   alpha_nucleotide + F2A + HGH_SS_2 + beta_nucleotide + 
                   stop_codon + In_Fusion_Right )

                   
# Output Genbank Format with annotations 
# Would be nice to automate this, but good enough for now
# Need: species with underscore, today's date in format 27-APR-2018
# Need: Actual regions of V, J, CDR3, overlapping
# ignore vntifkey and make everything CDS except In_Fusion and that region between LNGFR and P2A
# Positions are 1 indexed, inclusive start, inclusive end
#
# Fix "locus tag" to "label"

# Start and stop positions of annotations
In_Fusion_Left_gbk_start = 1
In_Fusion_Left_gbk_stop = In_Fusion_Left_gbk_start + len(In_Fusion_Left) - 1

LNGFR_gbk_start = In_Fusion_Left_gbk_stop + 1
LNGFR_gbk_stop = LNGFR_gbk_start + len(LNGFR) - 1
 
LNGFR_additional_gbk_start = LNGFR_gbk_stop + 1
LNGFR_additional_gbk_stop = LNGFR_additional_gbk_start + len(LNGFR_additional) - 1

P2A_gbk_start = LNGFR_additional_gbk_stop + 1
P2A_gbk_stop = P2A_gbk_start + len(P2A) - 1

HGH_SS_1_gbk_start = P2A_gbk_stop + 1
HGH_SS_1_gbk_stop = HGH_SS_1_gbk_start + len(HGH_SS_1) - 1

Va_gbk_start = HGH_SS_1_gbk_stop + 1
Va_gbk_stop = Va_gbk_start + len(Va_nucleotide_trimmed) - 1

CDR3a_gbk_start = Va_gbk_stop - size_of_Va_overlap*3 + 1
CDR3a_gbk_stop = CDR3a_gbk_start + len(CDR3a_sequence)*3 - 1

Ja_gbk_start = CDR3a_gbk_start + len(CDR3a_sequence)*3 - size_of_Ja_overlap*3
Ja_gbk_stop = Ja_gbk_start + len(Ja_nucleotide_trimmed) - 1

TRAC_gbk_start = Ja_gbk_stop + 1
TRAC_gbk_stop = TRAC_gbk_start + len(TRAC_nucleotide) - 1

F2A_gbk_start = TRAC_gbk_stop + 1
F2A_gbk_stop = F2A_gbk_start + len(F2A) - 1

HGH_SS_2_gbk_start = F2A_gbk_stop + 1
HGH_SS_2_gbk_stop = HGH_SS_2_gbk_start + len(HGH_SS_2) - 1

Vb_gbk_start = HGH_SS_2_gbk_stop + 1
Vb_gbk_stop = Vb_gbk_start + len(Vb_nucleotide_trimmed) - 1

CDR3b_gbk_start = Vb_gbk_stop - size_of_Vb_overlap*3 + 1
CDR3b_gbk_stop = CDR3b_gbk_start + len(CDR3b_sequence)*3 - 1

Jb_gbk_start = CDR3b_gbk_start + len(CDR3b_sequence)*3 - size_of_Jb_overlap*3
Jb_gbk_stop = Jb_gbk_start + len(Jb_nucleotide_trimmed) - 1

TRBC_gbk_start = Jb_gbk_stop + 1
TRBC_gbk_stop = TRBC_gbk_start + len(TRBC_nucleotide) - 1

stop_codon_gbk_start = TRBC_gbk_stop + 1
stop_codon_gbk_stop = stop_codon_gbk_start + len(stop_codon) - 1

In_Fusion_Right_gbk_start = stop_codon_gbk_stop + 1
In_Fusion_Right_gbk_stop = In_Fusion_Right_gbk_start + len(In_Fusion_Right) - 1


Va_short_name = Va_gene_name.split('/')[0]
Va_short_name = Va_short_name.split('-')[0]
Va_short_name = Va_short_name.split('*')[0]
Va_short_name = Va_short_name[2:]
Ja_short_name = Ja_gene_name.split('/')[0]
Ja_short_name = Ja_short_name.split('-')[0]
Ja_short_name = Ja_short_name.split('*')[0]
Ja_short_name = Ja_short_name[2:]
Vb_short_name = Vb_gene_name.split('/')[0]
Vb_short_name = Vb_short_name.split('-')[0]
Vb_short_name = Vb_short_name.split('*')[0]
Vb_short_name = Vb_short_name[2:]
Jb_short_name = Jb_gene_name.split('/')[0]
Jb_short_name = Jb_short_name.split('-')[0]
Jb_short_name = Jb_short_name.split('*')[0]
Jb_short_name = Jb_short_name[2:]


locus = Va_short_name + "_" + Ja_short_name + "__" + Vb_short_name + "_" + Jb_short_name
locus_space = "                                                                "
sequence_length = str(In_Fusion_Right_gbk_stop)
todaydate = datetime.now().strftime("%d-%b-%Y")
todaydate = todaydate.upper()


Va_gene_name = sys.argv[3]
Ja_gene_name = sys.argv[4]
CDR3a_sequence = sys.argv[5]
Vb_gene_name = sys.argv[6]
Jb_gene_name = sys.argv[7]
CDR3b_sequence = sys.argv[8]

print("LOCUS       "+locus+locus_space[:24-len(locus)]+sequence_length+" bp ds-DNA     linear       "+todaydate+windows)
print("SOURCE      "+species+windows)
print("  ORGANISM  "+species+windows)
print("COMMENT     "+species+" "+Va_gene_name+", "+Ja_gene_name+", CDR3a:"+CDR3a_sequence+windows)
print("COMMENT     "+species+" "+Vb_gene_name+", "+Jb_gene_name+", CDR3b:"+CDR3b_sequence+windows)
print("COMMENT     "+windows)
print("COMMENT     ApEinfo:methylated:1"+windows)
print("FEATURES             Location/Qualifiers"+windows)

#In_Fusion_Left
print("     misc_feature    "+str(In_Fusion_Left_gbk_start)+".."+str(In_Fusion_Left_gbk_stop)+windows)
print("                     /label=\"In-Fusion Left EcoRI\""+windows)
print("                     /ApEinfo_fwdcolor=\"#ffff80\""+windows)
print("                     /ApEinfo_revcolor=\"#ffff80\""+windows)
print("                     /ApEinfo_graphicformat=\"arrow_data {{0 1 2 0 0 -1} {} 0}"+windows)
print("                     width 5 offset 0\""+windows)

#LNGFR
print("     CDS             "+str(LNGFR_gbk_start)+".."+str(LNGFR_gbk_stop)+windows)
print("                     /label=\"LNGFR (deltaCD271)\""+windows)
print("                     /ApEinfo_fwdcolor=\"#e9d024\""+windows)
print("                     /ApEinfo_revcolor=\"#e9d024\""+windows)
print("                     /ApEinfo_graphicformat=\"arrow_data {{0 1 2 0 0 -1} {} 0}"+windows)
print("                     width 5 offset 0\""+windows)

#LNGFR_additional
print("     misc_feature    "+str(LNGFR_additional_gbk_start)+".."+str(LNGFR_additional_gbk_stop)+windows)
print("                     /label=\"Additional bases left over from cloning\""+windows)
print("                     /ApEinfo_fwdcolor=\"#ffff80\""+windows)
print("                     /ApEinfo_revcolor=\"#ffff80\""+windows)
print("                     /ApEinfo_graphicformat=\"arrow_data {{0 1 2 0 0 -1} {} 0}"+windows)
print("                     width 5 offset 0\""+windows)

#P2A
print("     CDS             "+str(P2A_gbk_start)+".."+str(P2A_gbk_stop)+windows)
print("                     /label=\"P2A\""+windows)
print("                     /ApEinfo_fwdcolor=\"#7eff74\""+windows)
print("                     /ApEinfo_revcolor=\"#7eff74\""+windows)
print("                     /ApEinfo_graphicformat=\"arrow_data {{0 1 2 0 0 -1} {} 0}"+windows)
print("                     width 5 offset 0\""+windows)

#HGH_SS_1
print("     sig_peptide     "+str(HGH_SS_1_gbk_start)+".."+str(HGH_SS_1_gbk_stop)+windows)
print("                     /label=\"HGH SS 1\""+windows)
print("                     /ApEinfo_fwdcolor=\"pink\""+windows)
print("                     /ApEinfo_revcolor=\"pink\""+windows)
print("                     /ApEinfo_graphicformat=\"arrow_data {{0 1 2 0 0 -1} {} 0}"+windows)
print("                     width 5 offset 0\""+windows)

#Va
print("     CDS             "+str(Va_gbk_start)+".."+str(Va_gbk_stop)+windows)
print("                     /label=\"Va\""+windows)
print("                     /ApEinfo_fwdcolor=\"#ff80ff\""+windows)
print("                     /ApEinfo_revcolor=\"ff80ff\""+windows)
print("                     /ApEinfo_graphicformat=\"arrow_data {{0 1 2 0 0 -1} {} 0}"+windows)
print("                     width 5 offset 0\""+windows)

#CDR3a
print("     CDS             "+str(CDR3a_gbk_start)+".."+str(CDR3a_gbk_stop)+windows)
print("                     /label=\"CDR3a\""+windows)
print("                     /ApEinfo_fwdcolor=\"#ff80ff\""+windows)
print("                     /ApEinfo_revcolor=\"ff80ff\""+windows)
print("                     /ApEinfo_graphicformat=\"arrow_data {{0 1 2 0 0 -1} {} 0}"+windows)
print("                     width 5 offset 0\""+windows)

#Ja_gbk
print("     CDS             "+str(Ja_gbk_start)+".."+str(Ja_gbk_stop)+windows)
print("                     /label=\"Ja\""+windows)
print("                     /ApEinfo_fwdcolor=\"#ff80ff\""+windows)
print("                     /ApEinfo_revcolor=\"ff80ff\""+windows)
print("                     /ApEinfo_graphicformat=\"arrow_data {{0 1 2 0 0 -1} {} 0}"+windows)
print("                     width 5 offset 0\""+windows)

#TRAC
print("     CDS             "+str(TRAC_gbk_start)+".."+str(TRAC_gbk_stop)+windows)
print("                     /label=\"TCR Alpha Constant\""+windows)
print("                     /ApEinfo_fwdcolor=\"#7eff74\""+windows)
print("                     /ApEinfo_revcolor=\"#7eff74\""+windows)
print("                     /ApEinfo_graphicformat=\"arrow_data {{0 1 2 0 0 -1} {} 0}"+windows)
print("                     width 5 offset 0\""+windows)

#F2A
print("     CDS             "+str(F2A_gbk_start)+".."+str(F2A_gbk_stop)+windows)
print("                     /label=\"F2A\""+windows)
print("                     /ApEinfo_fwdcolor=\"#e9d024\""+windows)
print("                     /ApEinfo_revcolor=\"#e9d024\""+windows)
print("                     /ApEinfo_graphicformat=\"arrow_data {{0 1 2 0 0 -1} {} 0}"+windows)
print("                     width 5 offset 0\""+windows)

#HGH_SS_2
print("     sig_peptide     "+str(HGH_SS_2_gbk_start)+".."+str(HGH_SS_2_gbk_stop)+windows)
print("                     /label=\"HGH SS 2\""+windows)
print("                     /ApEinfo_fwdcolor=\"pink\""+windows)
print("                     /ApEinfo_revcolor=\"pink\""+windows)
print("                     /ApEinfo_graphicformat=\"arrow_data {{0 1 2 0 0 -1} {} 0}"+windows)
print("                     width 5 offset 0\""+windows)

#Vb
print("     CDS             "+str(Vb_gbk_start)+".."+str(Vb_gbk_stop)+windows)
print("                     /label=\"Vb\""+windows)
print("                     /ApEinfo_fwdcolor=\"#ff80ff\""+windows)
print("                     /ApEinfo_revcolor=\"ff80ff\""+windows)
print("                     /ApEinfo_graphicformat=\"arrow_data {{0 1 2 0 0 -1} {} 0}"+windows)
print("                     width 5 offset 0\""+windows)

#CDR3b
print("     CDS             "+str(CDR3b_gbk_start)+".."+str(CDR3b_gbk_stop)+windows)
print("                     /label=\"CDR3b\""+windows)
print("                     /ApEinfo_fwdcolor=\"#ff80ff\""+windows)
print("                     /ApEinfo_revcolor=\"ff80ff\""+windows)
print("                     /ApEinfo_graphicformat=\"arrow_data {{0 1 2 0 0 -1} {} 0}"+windows)
print("                     width 5 offset 0\""+windows)

#Jb_gbk
print("     CDS             "+str(Jb_gbk_start)+".."+str(Jb_gbk_stop)+windows)
print("                     /label=\"Jb\""+windows)
print("                     /ApEinfo_fwdcolor=\"#ff80ff\""+windows)
print("                     /ApEinfo_revcolor=\"ff80ff\""+windows)
print("                     /ApEinfo_graphicformat=\"arrow_data {{0 1 2 0 0 -1} {} 0}"+windows)
print("                     width 5 offset 0\""+windows)

#TRBC
print("     CDS             "+str(TRBC_gbk_start)+".."+str(TRBC_gbk_stop)+windows)
print("                     /label=\"TCR Beta Constant\""+windows)
print("                     /ApEinfo_fwdcolor=\"#7eff74\""+windows)
print("                     /ApEinfo_revcolor=\"#7eff74\""+windows)
print("                     /ApEinfo_graphicformat=\"arrow_data {{0 1 2 0 0 -1} {} 0}"+windows)
print("                     width 5 offset 0\""+windows)

#stop_codon
print("     misc_feature    "+str(stop_codon_gbk_start)+".."+str(stop_codon_gbk_stop)+windows)
print("                     /label=\"STOP\""+windows)
print("                     /ApEinfo_fwdcolor=\"#7eff74\""+windows)
print("                     /ApEinfo_revcolor=\"#7eff74\""+windows)
print("                     /ApEinfo_graphicformat=\"arrow_data {{0 1 2 0 0 -1} {} 0}"+windows)
print("                     width 5 offset 0\""+windows)

#In_Fusion_Right
print("     misc_feature    "+str(In_Fusion_Right_gbk_start)+".."+str(In_Fusion_Right_gbk_stop)+windows)
print("                     /label=\"In-Fusion Right EcoRI\""+windows)
print("                     /ApEinfo_fwdcolor=\"#ffff80\""+windows)
print("                     /ApEinfo_revcolor=\"#ffff80\""+windows)
print("                     /ApEinfo_graphicformat=\"arrow_data {{0 1 2 0 0 -1} {} 0}"+windows)
print("                     width 5 offset 0\""+windows)


# Print sequence
sequence_label_spacer = "         "
base_count = 1
current_line =  sequence_label_spacer[len(str(base_count)):] + str(base_count)
print("ORIGIN"+windows)
for base in final_sequence:
    # Each line contains a line number and 60 bases
    if ((base_count - 1) % 60 == 0) and (base_count > 1):
        print(current_line + windows)
        current_line =  sequence_label_spacer[len(str(base_count)):] + str(base_count)
    # Every 10 bases is delimited by a space
    if (base_count - 1) % 10 == 0:
        current_line = current_line + " "
    current_line = current_line + base
    base_count += 1

if (base_count - 1) % 60 != 0:
    print(current_line + windows)

# There was a space before the "//" line endings that ApE didn't like
print("//"+windows)

#print(final_sequence)
