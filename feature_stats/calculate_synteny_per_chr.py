#!/usr/bin/env python3
#%%
# This script will calculate synteny between two genomes. Started working on 14/10
#%%
import sys
sys.settrace
import argparse
# %%
def parse_busco_table(table_file, table_type, seq2Merian={}):
    with open(table_file, 'r') as table:
        chr2pos, pos2buscoID = {}, {} 
        for line in table:
            if not line.startswith("#"): # ignoring comments at top
                cols = line.rstrip("\n").split() # get columns
                if cols[1] == "Complete":
                    if ':' not in cols[2]: # lines with the ':' bug can also have a bug in coordinates, so best to filter out
                        buscoID, chr, start, end = cols[0], cols[2].split(":")[0], int(cols[3]), int(cols[4])
                        pos = (start + end)/2 # get midpoint coordinate of busco (position)
                        if table_type == "reference": # if parsing the reference busco table..
                            try: 
                                chr2pos[chr].append(pos) # try and add the position to the list of BUSCO positions on that chromosome
                            except KeyError: 
                                chr2pos[chr] = [pos] # or create a new list if chromosome hasn't be seen yet
                            pos2buscoID[str(pos) + "_" + chr] = buscoID # key = position_chr, value = buscoID
                        else: # if parsing the query busco table...
                            try: 
                                Merian = seq2Merian[chr] # Convert sequence to assigned Merian
                            except KeyError:
                                continue # This stops fused/split chr from throwing an error as absent from dict
                            try: 
                                chr2pos[Merian].append(pos) # try and add the position to the list of BUSCO positions on that chromosome
                            except KeyError: 
                                chr2pos[Merian] = [pos] # or create a new list if chromosome hasn't be seen yet
                            pos2buscoID[str(pos) + "_" + Merian] = buscoID # key = position_chr, value = buscoID
    for merian in chr2pos.keys():
        pos_list = chr2pos[merian]
        assert len(pos_list) == len(set(pos_list)), f"expected to have a unique set of positions, this has been violated for Merian {merian}"
    return chr2pos, pos2buscoID

def filter_buscos(pos2busco_dict, chr2pos_dict, buscos_to_remove):
    pos2buscoID_ivd = {v: k for k, v in pos2busco_dict.items()} # create an inverse dict to allow lookups
    for i in buscos_to_remove: # for each busco missing from query/ref
        pos = pos2buscoID_ivd[i] # find the corresponding position
        pos2busco_dict.pop(pos) # remove that pos
        chr_pos, chr = pos.split('_')[0], pos.split('_')[1] # now need to do same to second dict
        current_list_of_pos = chr2pos_dict[chr] # get current list, then remove pos from pos_list
        current_list_of_pos.remove(float(chr_pos)) # have to convert to float, or else value is interpreted as string and isn't found in list
        chr2pos_dict[chr] = current_list_of_pos # updated edited list 
    return(pos2busco_dict, chr2pos_dict)

def filter_chromosomes(chr2pos, seq2Merian_ref, seq2Merian_query): # filter reference set of chr to only keep those that are "ancestral"
    set1 = (list(set(seq2Merian_ref.values()) - set(list(seq2Merian_query.values())))) # merian/merian combos in ref but not in query
    if len(set1) != 0:
        for key, value in seq2Merian_ref.items():
            if value in set1:
                chr2pos.pop(key)
    return(chr2pos)

#def filter_chromosomes(chr2pos_dict, seq2Merian_dict, ref_or_query): # filter reference set of chr to only keep those that are "ancestral"
 #   if ref_or_query == "reference":
  #          rearranged_chr = (list(set(chr2pos_dict.keys()) - set(list(seq2Merian_dict.keys()))))
   # else:
    #    rearranged_chr = (list(set(chr2pos_dict.keys()) - set(list(seq2Merian_dict.values())))) # keys are merians in query chr2pos_dict
    #for chr in rearranged_chr:
     #   chr2pos_dict.pop(chr)
    #return(chr2pos_dict)

def find_pairs(chr2pos, pos2buscoID):
    chr2pos_pair, pos2buscoID_pair = {}, {}
    for seq in chr2pos:
        list_pos = sorted(chr2pos[seq])  # put pos in ascending order within sequence
        total_pairs = len(list_pos) - 1 # Number of pairwise BUSCOs
        for i in list(range(0, total_pairs)):
            first_busco_pos, second_busco_pos = list_pos[i], list_pos[i+1]
            first_pos_seq = str(first_busco_pos) + "_" + str(seq)
            second_pos_seq = str(second_busco_pos) + "_" + str(seq)
            first_busco, second_busco = pos2buscoID[first_pos_seq], pos2buscoID[second_pos_seq]  # look up first_busco_pos & second_busco_pos in pos2buscoID
            if str(first_busco) > str(second_busco): # define 'busco_pair' based on descending lexographic order
                busco_pair = str(first_busco) + ',' + str(second_busco)
            else:
                busco_pair = str(second_busco) + ',' + str(first_busco)
            avg_mid_pos = (first_busco_pos + second_busco_pos)/2
            try: 
                chr2pos_pair[seq].append(avg_mid_pos) # try and add the position to the list of BUSCO positions on that chromosome
            except KeyError: 
                chr2pos_pair[seq] = [avg_mid_pos] # or create a new list if chromosome hasn't be seen yet
            pos2buscoID_pair[str(avg_mid_pos) + "_" + seq] = busco_pair # key = position_chr, value = buscoID
    return chr2pos_pair, pos2buscoID_pair

def parse_assignments_file(assignments_file):
    with open(assignments_file, 'r') as assignments:
        seq2Merian = {}
        for line in assignments:
            if not line.startswith("q"): # ignore header
                cols = line.rstrip("\n").split() # get columns
                #if (cols[1] != "fusion") & (cols[1] != "split"): # Filter out fused chr (for now..) and split chr
                if cols[1] != "split": # Filter out fused chr (for now..) and split chr
                    seq, assigned_M= cols[0], cols[2]
                    seq2Merian[seq] = assigned_M
    return seq2Merian

def calculate_synteny(chr2pos_pair, pos2buscoID_pair, pos2buscoID_query_pair, seq2Merian_ref, window_size, prefix):
    window_size_kb = int(window_size/1000) # convert bases to kb
    with open(prefix + '_per_chr_synteny_results.tsv', "w") as output_chr_file:
        output_chr_file.write(("%s\t%s\t%s\t%s\t%s\t%s\n") % ("ref_seq", "query_seq", "matching_busco_pairs", "mismatching_busco_pairs", "total_buscos", "per_match"))
        with open(prefix + '_' + str(window_size_kb) + 'kb' + '_synteny_results.tsv', "w") as output_file:
            output_file.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") % ("ref_seq", "query_seq", "start","stop","matching_busco_pairs", "mismatching_busco_pairs", "total_buscos", "per_match")) 
            for seq in chr2pos_pair:
                merian = seq2Merian_ref[seq]
                query_chr = '???'
                for key, value in seq2Merian_query.items(): # get cognate chromosome in query species
                    if value == merian:
                        query_chr = key
                total_buscos_chr, num_matches_chr = 0, 0
                pos_list = chr2pos_pair[seq] # Pos_list length is given by len(pos_list))
                last_window = max(pos_list) + window_size # max_pos = max(pos_list)
                window_list = range(window_size, int(last_window), window_size) # was 'range(window_size, int(last_window), window_size)' before
                for window in window_list: # NB: last window will be <100 kb. No explicit adjustment as per_match is reported but worth noting.
                    num_matches, num_mismatches = 0, 0 # number of matching vs mismatching pairs of buscos
                    pos_in_window = [i for i in pos_list if i <= window]   # Find busco_pairs in REF in WINDOW
                    for pos in pos_in_window:
                        pos_list.remove(pos)
                    for position in pos_in_window:
                        busco_pair = pos2buscoID_pair[str(position) + '_' + str(seq)]
                        if busco_pair in pos2buscoID_query_pair.values(): # see if busco_pair is also in query spp
                            num_matches += 1
                            num_matches_chr = num_matches_chr +1
                        else:
                            num_mismatches += 1 # mismatch found
                    total_buscos = num_matches + num_mismatches
                    total_buscos_chr = total_buscos_chr + total_buscos
                    if total_buscos == 0:
                        per_match = 'NA' # by definition as none to match
                    else:
                        per_match = (num_matches/ total_buscos) *100
                    window_start = window - window_size # 'window' is technically the value of the end of the window. Thus to get start of window, need to subtract the window size
                    output_file.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") % (seq, query_chr, window_start, window, num_matches, num_mismatches, total_buscos, per_match))
                    num_mismatches_chr = total_buscos_chr - num_matches_chr
                    if total_buscos_chr == 0:
                        per_match = 'NA'
                    else:
                        per_match_chr = (num_matches_chr / total_buscos_chr)*100
                output_chr_file.write(("%s\t%s\t%s\t%s\t%s\t%s\n") % (seq, query_chr, num_matches_chr, num_mismatches_chr, total_buscos_chr, per_match_chr))
    return()


# %%
if __name__ == "__main__":
    SCRIPT = "calculate_synteny.py"
    # Argument set up
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--reference", type=str, help = "BUSCO table file for reference species", required=True)
    parser.add_argument("-q", "--query", type=str, help = "BUSCO table file for query species", required=True)
    parser.add_argument("-ra", "--ref_assignments", type=str, help = "Assignments of Merian elements to chromosomes for reference species", required=True)
    parser.add_argument("-qa", "--query_assignments", type=str, help = "Assignments of Merian elements to chromosomes for query species", required=True)
    parser.add_argument("-w", "--window_size", type=str, help = "Window size (bases)", required=True)
    parser.add_argument("-o", "--output", type=str, help = "Output filename for the TSV (without extension)", required=True)
    args = parser.parse_args()
    reference_table_file = args.reference
    query_table_file = args.query
    query_assignments_file = args.query_assignments
    reference_assignments_file = args.ref_assignments
    window_size = int(args.window_size) # been using 1000kb previously
    prefix = args.output

print('[+] Processing ', query_table_file)
chr2pos, pos2buscoID = parse_busco_table(reference_table_file, "reference") # read in reference_table
seq2Merian_ref = parse_assignments_file(reference_assignments_file) # Make a dict of assignments to sequence
seq2Merian_query = parse_assignments_file(query_assignments_file) # Make a dict of assignments to sequence
chr2pos_query, pos2buscoID_query = parse_busco_table(query_table_file, "query", seq2Merian_query) # read in reference_table
query_buscos = pos2buscoID_query.values()
ref_buscos = pos2buscoID.values()
missing_query = set(ref_buscos) - set(query_buscos) # buscos in ref missing from query
missing_ref = set(query_buscos) - set(ref_buscos) # buscos in query missing from ref
pos2buscoID,chr2pos = filter_buscos(pos2buscoID, chr2pos, missing_query) # missing_query buscos need to be removed from ref
pos2buscoID_query,chr2pos_query = filter_buscos(pos2buscoID_query, chr2pos_query, missing_ref) # missing_query buscos need to be removed from ref
chr2pos = filter_chromosomes(chr2pos, seq2Merian_ref, seq2Merian_query) # Remove rearranged chr 
if len(seq2Merian_query) == 0:
    print('[+] All chromosomes in query species have undergone fusion/fission events, thus synteny cannot be calculated. Exiting.')
    exit()
chr2pos_pair, pos2buscoID_pair = find_pairs(chr2pos, pos2buscoID)
chr2pos_query_pair, pos2buscoID_query_pair = find_pairs(chr2pos_query, pos2buscoID_query)
calculate_synteny(chr2pos_pair, pos2buscoID_pair, pos2buscoID_query_pair, seq2Merian_ref, window_size, prefix)
exit()
# %%
# Algorithm summary:
#pos2buscoID         # Position_Sequence: busco_ID
#chr2pos             # For each sequence, what positions

# Define a set of pairs of BUSCOs in ref_spp and midpos
#pos2buscoID_pair    # AvgMidPos_sequence: busco_pair_IDs
#chr2pos_pair        # For each sequence, what AvgMidPos'

# Convert sequences to assigned_merians in query_spp
# Define a set of pairs of BUSCOs in query_spp and midpos
#pos2buscoID_query   # AvgMidPos_sequence: busco_pair_IDs
#seq2Merian          # is the dict of sequences to assigned_merian
#chr2pos_query       # For each sequence, what AvgMidPos'

# Then do 100 kb intervals what % of entries from spp1 is in spp2
#%%
# Manual run
reference_table_file ='Melitaea_cinxia.tsv'
query_table_file = 'Bombyx_mori_M23.tsv'
query_assignments_file = 'Bombyx_mori_chromosome_assignments.tsv'
reference_assignments_file = 'Melitaea_cinxia_chromosome_assignments.tsv'
prefix = 'test'
window_size = 1000000
chr2pos, pos2buscoID = parse_busco_table(reference_table_file, "reference") # read in reference_table
seq2Merian_ref = parse_assignments_file(reference_assignments_file) # Make a dict of assignments to sequence
seq2Merian_query = parse_assignments_file(query_assignments_file) # Make a dict of assignments to sequence
chr2pos_query, pos2buscoID_query = parse_busco_table(query_table_file, "query", seq2Merian_query) # read in reference_table
query_buscos = pos2buscoID_query.values()
ref_buscos = pos2buscoID.values()
missing_query = set(ref_buscos) - set(query_buscos) # buscos in ref missing from query
missing_ref = set(query_buscos) - set(ref_buscos) # buscos in query missing from ref
pos2buscoID,chr2pos = filter_buscos(pos2buscoID, chr2pos, missing_query) # missing_query buscos need to be removed from ref
pos2buscoID_query,chr2pos_query = filter_buscos(pos2buscoID_query, chr2pos_query, missing_ref) # missing_query buscos need to be removed from ref
chr2pos = filter_chromosomes(chr2pos, seq2Merian_ref, seq2Merian_query)
if len(seq2Merian_query) == 0:
    print('[+] All chromosomes in query species have undergone fusion/fission events, thus synteny cannot be calculated. Exiting.')
chr2pos_pair, pos2buscoID_pair = find_pairs(chr2pos, pos2buscoID)
chr2pos_query_pair, pos2buscoID_query_pair = find_pairs(chr2pos_query, pos2buscoID_query)
calculate_synteny(chr2pos_pair, pos2buscoID_pair, pos2buscoID_query_pair, seq2Merian_ref, window_size, prefix)
