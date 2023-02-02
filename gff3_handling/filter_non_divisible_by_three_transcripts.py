#!/usr/bin/env python3

import argparse 
import sys 
import os
# This script parses a gff and checks whether each transcript ID is composed of a set of CDS sequences which total a length which is divisible by three
# For transcripts that are not divisible by three, the sequences are removed from the filtered, output gff3 file


def calculate_spliced_cds_lengths(gff_file):
    transcript2CDSlengths = {}
    with open(gff_file) as file:
        for line in file:
            n=0
            if n <5:
                if not line.startswith('#'):
                    cols = line.split('\t')
                    if cols[2] == "CDS":
                        transcript_ID = cols[8].split(';')[1].replace('Parent=transcript:','')
                        start, stop = int(cols[3]), int(cols[4])
                        cds_length = abs(start - stop) +1
                        if transcript_ID in transcript2CDSlengths:
                            transcript2CDSlengths[transcript_ID].append(cds_length)
                        else:
                            transcript2CDSlengths[transcript_ID] = [cds_length]
    return(transcript2CDSlengths)

def find_indivisible_by_three_sequences(transcript2CDSlengths):
    transcripts_to_remove = []
    for key, value in transcript2CDSlengths.items():
        total_length = sum(value)
        total_length_by_3 = str(total_length/3)
        if (total_length_by_3.split('.')[1]) != str(0):
            print(key,total_length_by_3)
            transcripts_to_remove.append(key)
    return(transcripts_to_remove)

def output_filtered_gff3_file(transcripts_to_remove, gff_file):
        output_file = gff_file.replace('.gff3', '')
        DivByThree_file = output_file + '.DivByThree.gff3'
        NonDivByThree_file = output_file + '.NonDivByThree.gff3'
        with open(NonDivByThree_file, 'w') as NonDiv_file:
                with open(DivByThree_file, 'w') as Div_file:
                        with open(gff_file, 'r') as file:
                                for line in file:
                                        if any(ext in line for ext in transcripts_to_remove):
                                                NonDiv_file.write(line)
                                        else:
                                                Div_file.write(line)

parser = argparse.ArgumentParser(description='This filters a gff3 file to remove transcripts composed of CDS sequences that are indivisible by three.')
parser.add_argument("-g",  help="gff3 file")

args = parser.parse_args()
gff_file = args.g

print("[+] Parsing gff3 file:", gff_file)

transcript2CDSlengths = calculate_spliced_cds_lengths(gff_file)
transcripts_to_remove = find_indivisible_by_three_sequences(transcript2CDSlengths)
output_filtered_gff3_file(transcripts_to_remove, gff_file)
