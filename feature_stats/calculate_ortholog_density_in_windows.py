#!/usr/bin/env python3
# Started work on script on 7/11/21
#%%
from audioop import mul
import sys
import os
import argparse
#%%
#%%
def parse_orthologs(orthologs_file):
    with open(orthologs_file, 'r') as ortho_file:
        chr2midpos_single, chr2midpos_multi,chr2midpos_other = {}, {}, {}
        ortho_file.readline()
        for line in ortho_file:
                cols = line.rstrip("\n").split('\t') # get columns
                chr, midpos, status = cols[0], cols[2], cols[3]
                if status == "single_copy":
                    if chr not in chr2midpos_single.keys():
                        chr2midpos_single[chr] = [midpos]
                    else:
                        current_list = chr2midpos_single[chr]
                        current_list.append(midpos)
                        chr2midpos_single[chr] = current_list
                elif status == "multi_copy":
                    if chr not in chr2midpos_multi.keys():
                        chr2midpos_multi[chr] = [midpos]
                    else:
                        current_list = chr2midpos_multi[chr]
                        current_list.append(midpos)
                        chr2midpos_multi[chr] = current_list
                else:
                    if chr not in chr2midpos_other.keys():
                        chr2midpos_other[chr] = [midpos]
                    else:
                        current_list = chr2midpos_other[chr]
                        current_list.append(midpos)
                        chr2midpos_other[chr] = current_list
    return(chr2midpos_single, chr2midpos_multi, chr2midpos_other)

def count_ortholog_types_per_window_and_chr(prefix, chr2midpos_single, chr2midpos_multi, chr2midpos_other, window_size):
    chr_list = set(list(chr2midpos_single.keys()) + list(chr2midpos_multi.keys()) + list(chr2midpos_other.keys())) # make a list of all chr
    with open(str(prefix + '.classified_ortholog_counts_per_chr.tsv'), 'w') as chr_output:
        chr_output.write(("%s\t%s\t%s\t%s\n") % ("scaffold", "num_single_genes", "num_multi_genes", "num_other_genes"))
        with open(str(prefix + '.classified_ortholog_counts_per_' + str(int(window_size/1000)) + 'kb.tsv'), 'w') as window_output:
            window_output.write(("%s\t%s\t%s\t%s\t%s\t%s\n") % ("scaffold", "start", "end", "num_single_genes", "num_multi_genes", "num_other_genes"))
            for chr in chr_list:
                all_pos = [] # list to collect all gene pos on chr
                try:
                    single_genes = chr2midpos_single[chr]
                    num_single_genes = len(single_genes)
                    all_pos = all_pos + single_genes
                except KeyError:
                    num_single_genes, single_genes = 0, []
                try:
                    multi_genes = chr2midpos_multi[chr]
                    num_multi_genes = len(multi_genes)
                    all_pos = all_pos + multi_genes
                except KeyError:
                    num_multi_genes, multi_genes = 0, []
                try:
                    other_genes = chr2midpos_other[chr]
                    num_other_genes = len(other_genes)
                    all_pos = all_pos + other_genes
                except KeyError:
                    num_other_genes, other_genes = 0, []
                chr_output.write(("%s\t%s\t%s\t%s\n") % (chr, num_single_genes, num_multi_genes, num_other_genes))
                max_pos = max(all_pos, key=float) # needed if values aren't stored as numbers in list
                last_window = int(float(max_pos)) + window_size
                window_list = range(0, last_window, window_size) 
                for window in window_list: # First window is 0-window_size
                    single_in_window = [i for i in single_genes if (int(float(i)) >= window) & (int(float(i))  <=(window+window_size))]
                    multi_in_window = [i for i in multi_genes if (int(float(i)) >= window) & (int(float(i))  <=(window+window_size))]
                    other_in_window = [i for i in other_genes if (int(float(i)) >= window) & (int(float(i))  <=(window+window_size))]
                    window_output.write(("%s\t%s\t%s\t%s\t%s\t%s\n") % (chr, window, (window+window_size), len(single_in_window), len(multi_in_window), len(other_in_window)))
    return()
#%%
parser = argparse.ArgumentParser(description='This counts the number of unique/singletons/multicopy genes per window')
parser.add_argument("-o", help="File containing the classified orthologs (single, multi, other)")
parser.add_argument("-w",  help="Window size (in bases)", type=int)
parser.add_argument("-p",  help="Prefix for the output files")

args = parser.parse_args()
orthologs_file = args.o
window_size = args.w
prefix = args.p

print("[+] Parsing input files for " + prefix + "...")
chr2midpos_single, chr2midpos_multi, chr2midpos_other = parse_orthologs(orthologs_file)
print("[+] Counting orthologs per chromosome and per " + str(int(window_size/1000)) + "kb. Writing to output files.")
count_ortholog_types_per_window_and_chr(prefix, chr2midpos_single, chr2midpos_multi, chr2midpos_other, window_size)

#%%
# manual run
#orthologs_file = 'Vanessa_atalanta.classified_orthologs.tsv'
#date = '5678'
#prefix = orthologs_file.split('.')[0]
#window_size = 100000 # Make this adjustable - 100kb suggested by Lewis
