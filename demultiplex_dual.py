#!/usr/bin/env python
#
# Created by: Michael Leonard
# Last updated: February 14, 2019
#
# Usage:   demultiplex.py input.fastq [input_read2.fastq]
# Example: ./demultiplex.py read1.fastq read2.fastq
# 

import sys

if len(sys.argv) == 2:
    filename = sys.argv[1]
    filename2 = ""
elif len(sys.argv) == 3:
    filename = sys.argv[1]
    filename2 = sys.argv[2]
else:
    sys.exit("Usage: demultiplex.py input.fastq [input_read2.fastq]")

barcode_table = {
'CCGCGGTT+CTAGCGCT':'01', 'TACCGAGG+AGTTCAGG':'33', 'ACACTAAG+ATATGGAT':'65', 
'TTATAACC+TCGATATC':'02', 'CGTTAGAA+GACCTGAA':'34', 'GTGTCGGA+GCGCAAGC':'66', 
'GGACTTGG+CGTCTGCG':'03', 'AGCCTCAT+TCTCTACT':'35', 'TTCCTGTT+AAGATACT':'67', 
'AAGTCCAA+TACTCATA':'04', 'GATTCTGC+CTCTCGTC':'36', 'CCTTCACC+GGAGCGTC':'68', 
'ATCCACTG+ACGCACCT':'05', 'TCGTAGTG+CCAAGTCT':'37', 'GCCACAGG+ATGGCATG':'69', 
'GCTTGTCA+GTATGTTC':'06', 'CTACGACA+TTGGACTC':'38', 'ATTGTGAA+GCAATGCA':'70', 
'CAAGCTAG+CGCTATGT':'07', 'TAAGTGGT+GGCTTAAG':'39', 'ACTCGTGT+GTTCCAAT':'71', 
'TGGATCGA+TATCGCAC':'08', 'CGGACAAC+AATCCGGA':'40', 'GTCTACAC+ACCTTGGC':'72', 
'AGTTCAGG+TCTGTTGG':'09', 'ATATGGAT+TAATACAG':'41', 'CAATTAAC+ATATCTCG':'73', 
'GACCTGAA+CTCACCAA':'10', 'GCGCAAGC+CGGCGTGA':'42', 'TGGCCGGT+GCGCTCTA':'74', 
'TCTCTACT+GAACCGCG':'11', 'AAGATACT+ATGTAAGT':'43', 'AGTACTCC+AACAGGTT':'75', 
'CTCTCGTC+AGGTTATA':'12', 'GGAGCGTC+GCACGGAC':'44', 'GACGTCTT+GGTGAACC':'76', 
'CCAAGTCT+TCATCCTT':'13', 'ATGGCATG+GGTACCTT':'45', 'TGCGAGAC+CAACAATG':'77', 
'TTGGACTC+CTGCTTCC':'14', 'GCAATGCA+AACGTTCC':'46', 'CATAGAGT+TGGTGGCA':'78', 
'GGCTTAAG+GGTCACGA':'15', 'GTTCCAAT+GCAGAATT':'47', 'ACAGGCGC+AGGCAGAG':'79', 
'AATCCGGA+AACTGTAG':'16', 'ACCTTGGC+ATGAGGCC':'48', 'GTGAATAT+GAATGAGA':'80', 
'TAATACAG+GTGAATAT':'17', 'ATATCTCG+ACTAAGAT':'49', 'AACTGTAG+TGCGGCGT':'81', 
'CGGCGTGA+ACAGGCGC':'18', 'GCGCTCTA+GTCGGAGC':'50', 'GGTCACGA+CATAATAC':'82', 
'ATGTAAGT+CATAGAGT':'19', 'AACAGGTT+CTTGGTAT':'51', 'CTGCTTCC+GATCTATC':'83', 
'GCACGGAC+TGCGAGAC':'20', 'GGTGAACC+TCCAACGC':'52', 'TCATCCTT+AGCTCGCT':'84', 
'GGTACCTT+GACGTCTT':'21', 'CAACAATG+CCGTGAAG':'53', 'AGGTTATA+CGGAACTG':'85', 
'AACGTTCC+AGTACTCC':'22', 'TGGTGGCA+TTACAGGA':'54', 'GAACCGCG+TAAGGTCA':'86', 
'GCAGAATT+TGGCCGGT':'23', 'AGGCAGAG+GGCATTCT':'55', 'CTCACCAA+TTGCCTAG':'87', 
'ATGAGGCC+CAATTAAC':'24', 'GAATGAGA+AATGCCTC':'56', 'TCTGTTGG+CCATTCGA':'88', 
'ACTAAGAT+CCGCGGTT':'25', 'TGCGGCGT+TACCGAGG':'57', 'TATCGCAC+ACACTAAG':'89', 
'GTCGGAGC+TTATAACC':'26', 'CATAATAC+CGTTAGAA':'58', 'CGCTATGT+GTGTCGGA':'90', 
'CTTGGTAT+GGACTTGG':'27', 'GATCTATC+AGCCTCAT':'59', 'GTATGTTC+TTCCTGTT':'91', 
'TCCAACGC+AAGTCCAA':'28', 'AGCTCGCT+GATTCTGC':'60', 'ACGCACCT+CCTTCACC':'92', 
'CCGTGAAG+ATCCACTG':'29', 'CGGAACTG+TCGTAGTG':'61', 'TACTCATA+GCCACAGG':'93', 
'TTACAGGA+GCTTGTCA':'30', 'TAAGGTCA+CTACGACA':'62', 'CGTCTGCG+ATTGTGAA':'94', 
'GGCATTCT+CAAGCTAG':'31', 'TTGCCTAG+TAAGTGGT':'63', 'TCGATATC+ACTCGTGT':'95', 
'AATGCCTC+TGGATCGA':'32', 'CCATTCGA+CGGACAAC':'64', 'CTAGCGCT+GTCTACAC':'96'
}

# Perform script for one input file
if not filename2:
    # Open output files
    output_file_table1 = {}
    for current_barcode in barcode_table:
        current_index = barcode_table[current_barcode]
        current_output_filename1 = current_index + "_" + current_barcode + "_R1.fastq"
        output_file_table1[current_index] = open(current_output_filename1 , 'w')
    output_file_table1["Undetermined"] = open("Undetermined_R1.fastq" , 'w')
    
    # Demultiplex
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
        
        # Parse header, extract barcode
        barcode = ""
        header_array1 = header1.split(':')
        barcode1 = header_array1[-1]
        if barcode1 in barcode_table.keys():
            barcode = barcode1
            barcode_index = barcode_table[barcode]
        else:
            barcode_index = "Undetermined"
        
        # Output read to file
        current_output_file1 = output_file_table1[barcode_index]
        current_output_file1.write(header1 + "\n")
        current_output_file1.write(sequence1 + "\n+\n")
        current_output_file1.write(quality1 + "\n")

    # Close files
    for current_file in output_file_table1.values():
        current_file.close()
        

# Perform script for two input files
elif filename2:
    # Open output files
    output_file_table1 = {}
    output_file_table2 = {}
    for current_barcode in barcode_table:
        current_index = barcode_table[current_barcode]
        current_output_filename1 = current_index + "_" + current_barcode + "_R1.fastq"
        current_output_filename2 = current_index + "_" + current_barcode + "_R2.fastq"
        output_file_table1[current_index] = open(current_output_filename1 , 'w')
        output_file_table2[current_index] = open(current_output_filename2 , 'w')
    output_file_table1["Undetermined"] = open("Undetermined_R1.fastq" , 'w')
    output_file_table2["Undetermined"] = open("Undetermined_R2.fastq" , 'w')
    
    # Demultiplex
    fastq1 = open(filename, 'r')
    fastq2 = open(filename2, 'r')
    line1 = fastq1.readline()
    line2 = fastq2.readline()
    while line1:
        # Get header and sequence
        header1 = line1.strip()
        header2 = line2.strip()
        sequence1 = ""
        sequence2 = ""
        line1 = fastq1.readline()
        line2 = fastq2.readline()
        linecount = 0
        while line1 and (line1[0] != '+'):
            sequence1 = sequence1 + line1.strip()
            sequence2 = sequence2 + line2.strip()
            line1 = fastq1.readline()
            line2 = fastq2.readline()
            linecount += 1
        
        # Get quality line
        quality1 = ""
        quality2 = ""
        line1 = fastq1.readline()
        line2 = fastq2.readline()
        while linecount > 0:
            quality1 = quality1 + line1.strip()
            quality2 = quality2 + line2.strip()
            line1 = fastq1.readline()
            line2 = fastq2.readline()
            linecount -= 1
        
        # Parse header, extract barcode
        barcode = ""
        header_array1 = header1.split(':')
        barcode1 = header_array1[-1]
        header_array2 = header2.split(':')
        barcode2 = header_array2[-1]
        if barcode1 in barcode_table.keys():
            barcode = barcode1
            barcode_index = barcode_table[barcode]
        elif barcode2 in barcode_table.keys():
            barcode = barcode2
            barcode_index = barcode_table[barcode]
        else:
            barcode_index = "Undetermined"
        
        # Output read to file
        current_output_file1 = output_file_table1[barcode_index]
        current_output_file1.write(header1 + "\n")
        current_output_file1.write(sequence1 + "\n+\n")
        current_output_file1.write(quality1 + "\n")
        
        current_output_file2 = output_file_table2[barcode_index]
        current_output_file2.write(header2 + "\n")
        current_output_file2.write(sequence2 + "\n+\n")
        current_output_file2.write(quality2 + "\n")

    # Close files
    for current_file in output_file_table1.values():
        current_file.close()
    for current_file in output_file_table2.values():
        current_file.close()
        
