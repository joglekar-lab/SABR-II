# SABR-II

#Description of scripts/files

#Seurat and scRepertoire analyses
#Zdinak-et-al-seurat+screpertoire.R
#This file has all the relevant steps for the analysis of scRNAseq data from the paper

#TCR reconstruction:
#TCRgen_mouse.opt_v2.py
#imgt_tcr_mouse.nuc.fa
#The script TCRgen_mouse.opt_v2.py is used to reconstruct full length TCRs using sequences obtained from the imgt_tcr_mouse.nuc.fa file
#Usage:
#python TCRgen_mouse.opt_v2.py imgt_tcr_mouse.nuc.fa "Mus musculus" "TRAV Allele" "TRAJ Allele" "CDR3Alpha" "TRBV Allele" "TRBJ Allele" "CDR3Beta" 1>TCR.opt.gb 2>TCR.opt.out.txt
#The output is two files: TCR.opt.gb - full length TCR insert; TCR.opt.out.txt - summary and log of errors

#SABR library oligo construction
backtranslate_fast_noU_upto25.py
Usage:
python backtranslate_fast_noU_upto25.py epitope_list.csv
Input file is a list of epitope names and epitopes in two columns.
Output file is a list of epitopes + nucleotide seqeuences in two columns



#SABR screen analysis
#demultiplex_dual.py
#De-multiplexes UDI indexes into individual indexes
#Usage:
#python demultiplex_dual.py R1.fasta R2.fasta
#R1.fasta and R2.fasta are concatenated R1/R2 files from sequencing data

#epitope_extract_fastq_v1.1.py
#Extracts the epitope out from sequencing data and counts reads per epitope in the library
#Usage:
#ls *_R1.fastq | parallel "python epitope_extract_fastq_v1.1.py {} AA.txt {}.v2_SABRcounts.txt"
#AA.txt is a list of all the epitopes in the library
#The output is a file named SAMPLE.v2_SABRcounts.txt

#merge_counts_split_v2.1.py
#Merges all the epitope read count lists into one file
#Usage:
#python merge_counts_split_v2.1.py *.v2_SABRcounts.txt > OUTPUT.txt
#The output is a file named OUTPUT.txt, in which read counts for all epitopes are listed in the same order per sample

