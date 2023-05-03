#!/usr/bin/env python3
#%%
# import sys
import argparse
import re
#%%
# Parse index file and itialise windows
def parse_index(index_file):
    with open(index_file) as table:
        chr2length = {}
        chr_list = []
        for line in table:
            cols = line.rstrip("\n").split()
            chr, length = cols[0], cols[1]
            chr2length[chr] = length
            chr_list.append(chr)
    chr_windows_GC3_dict = {}
    chr_windows_length_dict = {}
    for chr in chr_list:
        chr_windows_GC3_dict[chr] = {}
        chr_windows_length_dict[chr] = {}
        length = chr2length[chr]
        max_window = (int(length)/window_size) + 1
        # was "length" rather than max_window. Think max_window is correct
        list_windows = list(range(0,int(max_window)))
        for window in list_windows:
            chr_windows_GC3_dict[chr][window] = 0
            chr_windows_length_dict[chr][window] = 0
    return chr_list, chr_windows_GC3_dict, chr_windows_length_dict

def calculate_GC3_in_windows(table_file, chr_windows_GC3_dict, chr_windows_length_dict, window_size):
    with open(table_file, 'r') as table:
        # Read in output of gff-stats (fron spliced mode)
        for line in table:
                if not line.startswith('ID'):
                    cols = line.rstrip("\n").split()
                    chr, start, stop, gc3_per = cols[0], cols[2], cols[3], cols[12]
                    int_start = int(int(start)/window_size) # get what window it is in
                    int_stop = int(int(stop)/window_size) # get what window it is in
                    length = abs(int(stop) - int(start)) #  # get absolute length as 'stop' may be higher or lower than 'start' depending on orientation +/-
                    if int_start == int_stop: # if start and stop are in the same window
                        window_name = int_start
                        current_GC3_value = chr_windows_GC3_dict[chr][window_name]
                        added_GC3_value = current_GC3_value + float(gc3_per)*length
                        chr_windows_GC3_dict[chr][window_name] = added_GC3_value
                        current_length_vaue = chr_windows_length_dict[chr][window_name]
                        added_length_value = current_length_vaue + length
                        chr_windows_length_dict[chr][window_name] = added_length_value
    return(chr_windows_GC3_dict, chr_windows_length_dict)

# Write output file
def write_GC3_per_window_file(prefix, chr_list, chr_windows_GC3_dict, chr_windows_length_dict):
    window_size_kb = int(window_size / 1000)
    with open(prefix + "_gc3_per_" + str(window_size_kb) + "kb.tsv", "w") as GC3_windows_file:
            GC3_windows_file.write(("%s\t%s\t%s\t%s\n") % ("fasta_id", "start", "stop", "gc3_per"))
            for chr in chr_list:
                current_dict = chr_windows_GC3_dict[chr]
                for window in current_dict:
                    GC3_length = current_dict[window]
                    total_length = chr_windows_length_dict[chr][window]
                    # Calculate GC3% for window
                    if GC3_length != 0:
                        GC3_per = GC3_length / total_length
                        start = window*window_size
                        stop = (window+1)*window_size
                        GC3_windows_file.write(("%s\t%s\t%s\t%s\n") % (chr, start, stop, GC3_per))
    print("[+]\tSuccessfully written GC3_prop per window to '" + prefix + "_gc3_per_" + str(window_size_kb) + "kb.tsv'")

def write_GC3_per_chr_file(prefix, chr_list, chr_windows_GC3_dict, chr_windows_length_dict):
    with open(prefix + "_gc3_per_chr.tsv", "w") as GC3_per_chr_file:
            GC3_per_chr_file.write(("%s\t%s\t%s\n") % ("fasta_id", "cds_length","gc3_per")) # cds_length is the total length of cds' used to calcualte gc3 per 
            for chr in chr_list:
                current_dict = chr_windows_GC3_dict[chr]
                total_chr_length = 0
                total_gc_length = 0
                for window in current_dict:
                    GC3_length = current_dict[window]
                    total_window_length = chr_windows_length_dict[chr][window]
                    total_chr_length = total_chr_length + total_window_length
                    if total_chr_length != 0:
                        total_gc_length = total_gc_length + GC3_length
                if total_gc_length != 0:
                    GC3_per = total_gc_length / total_chr_length # total GC3_per for chromosome
                    GC3_per_chr_file.write(("%s\t%s\t%s\n") % (chr,total_chr_length, GC3_per))
    print("[+]\tSuccessfully written GC3_prop per chr to '" + prefix + "_gc3_per_chr.tsv'")
#%%
if __name__ == "__main__":
    SCRIPT = "Calculate_GC3_per_window.py"
    # Argument set up
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--index", type=str, help = "Index file for the genome", required=True)
    parser.add_argument("-s", "--stats", type=str, help = "Stats file generated by gff-stats stat", required=True)
    parser.add_argument("-w", "--window_size", type=str, help = "Window size (bases)", required=True)
    parser.add_argument("-o", "--output", type=str, help = "Output filename for the TSV (without extension)", required=True)
    args = parser.parse_args()
    index_file = args.index
    table_file = args.stats
    window_size = args.window_size
    window_size = int(window_size)
    prefix = args.output

    # Run the functions
print("[+] Parsing index file and initialise windows...")
chr_list, chr_windows_GC3_dict, chr_windows_length_dict = parse_index(index_file)
print("[+] Calculate GC3 per window")
chr_windows_GC3_dict, chr_windows_length_dict = calculate_GC3_in_windows(table_file, chr_windows_GC3_dict, chr_windows_length_dict, window_size)
write_GC3_per_window_file(prefix, chr_list, chr_windows_GC3_dict, chr_windows_length_dict)
print("[+] Calculate GC3 per chr")
write_GC3_per_chr_file(prefix, chr_list, chr_windows_GC3_dict, chr_windows_length_dict)
