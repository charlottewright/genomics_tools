#!/usr/bin/env python3
#%%
import argparse

def get_transcript_2_cds_ids(file_path):
    transcript_2_cds_dict = {}
    with open(file_path, 'r') as file:
        for line in file:
            if 'CDS' in line:
                cols = line.strip().split('\t')
                info_list = cols[8].split(';')
                first_value = info_list[0]
                cds_id = first_value.split(':')[1].strip()
                second_value = info_list[1]
                transcript_id = second_value.split(':')[1].strip()
                transcript_2_cds_dict[transcript_id] = cds_id
    return(transcript_2_cds_dict)

def get_gene_2_transcript_ids(file_path):
    gene_2_transcript_id = {}
    with open(file_path, 'r') as file:
        for line in file:
            if 'mRNA' in line:
                cols = line.strip().split('\t')
                feature = cols[2]
                if feature == 'mRNA':
                    info_list = cols[8].split(';')
                    first_value = info_list[0]
                    transcript_id = first_value.split(':')[1].strip()
                    second_value = info_list[1]
                    gene_id = second_value.split(':')[1].strip()
                    gene_2_transcript_id[gene_id] = transcript_id
    return(gene_2_transcript_id)

def read_in_genes(file_path, prefix, transcript_2_cds_dict, gene_2_transcript_id):
    output_filename = prefix + '_gene_locations.bed'
    output_lines = [] 
    with open(file_path, 'r') as file:
        for line in file:
            cols = line.strip().split('\t')
            if len(cols) >= 9 and cols[2] == 'gene':
                chr_val, start_val, stop_val = cols[0], cols[3], cols[4]
                gene_id = cols[8].split(';')[0].split(':')[1].strip()
                # Look up the gene ID in transcript_gene_dict to get the transcript ID
                transcript_id = gene_2_transcript_id[gene_id]
                cds_id = transcript_2_cds_dict[transcript_id] + '.1'
                # Write values of 'chr', 'transcript_id', 'start', and 'stop' to the output file
                output_lines.append(f"{chr_val}\t{cds_id}\t{start_val}\t{stop_val}\n")
    # Write to the output file
    with open(output_filename, 'w') as output_file:
        output_file.writelines(output_lines)
    return()

def process_files(input_file, output_file):
    #read in gff3 and output a bed file of gene locations where gene_ids are the cds_ids. Needed to match pep files.
    print(f"Processing files: {input_file} -> {output_file}_gene_locations.bed")
    transcript_2_cds_dict = get_transcript_2_cds_ids(args.input)
    gene_2_transcript_id = get_gene_2_transcript_ids(args.input)
    read_in_genes(args.input, args.prefix, transcript_2_cds_dict, gene_2_transcript_id)

if __name__ == "__main__":
    # create argument parser
    parser = argparse.ArgumentParser(description="Process gff3 files.")

    # add arguments
    parser.add_argument("--input", type=str, help="Path to input file.")
    parser.add_argument("--prefix", type=str, help="Prefix for output bed file.")

    args = parser.parse_args() # parse arguments

    if args.input and args.prefix: # check both input and output arguments are provided
        process_files(args.input, args.prefix)
    else:
        print("Please provide both input gff3 file path and output prefix using --input and --prefix arguments.")
