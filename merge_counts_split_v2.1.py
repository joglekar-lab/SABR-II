#!/usr/bin/env python
#
# Created by: Michael Leonard
# Last updated: January 04, 2019
#
# Usage:   merge_counts.py count_file1.txt count_file2.txt ...
# Example: ./merge_counts.py 01.txt 02.txt 03.txt
# 

# This method is slow, a much better method would be to create 
# a 2d list with zeroes each time a new sequence is found, then fill
# the list as counts are read

import sys

filename_list = sys.argv[1:]
master_count_table = {}
master_sequence_list = []

# Get counts
sys.stderr.write("Reading input files: ")
for filename in filename_list:
    sys.stderr.write(filename + " ")
    current_count_file = open(filename, 'r')
    current_count_dictionary = {}
    for line in current_count_file:
        line = line.strip()
        line_array = line.split(",") #modified this because the outputs are csv
        current_sequence = line_array[0]
        current_count = line_array[1]
        current_count_dictionary[current_sequence] = current_count
        if current_sequence not in master_sequence_list:
            master_sequence_list.append(current_sequence)
    master_count_table[filename] = current_count_dictionary
sys.stderr.write("\n")

# Write headers
sys.stderr.write("Outputting compiled table\n")
for filename in filename_list:
    current_index = filename.strip()
    current_index = current_index.split('.')[0]
    #current_index = current_index.split('_')[0]
    sys.stdout.write(current_index + "\t" + current_index + "\t")
sys.stdout.write("\n")

# Write counts
for sequence in master_sequence_list:
    for filename in filename_list:
        sys.stdout.write(sequence + "\t")
        sys.stdout.write(str(master_count_table[filename].get(sequence,0)) + "\t")
    sys.stdout.write("\n")