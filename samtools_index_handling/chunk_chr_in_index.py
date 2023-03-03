#! /usr/bin/env python
import sys
import argparse

def chunk_chr_and_write_to_file(input_file, number_chunks, prefix):
    with open(prefix + "_chunked_chr_" + str(number_chunks) + ".bed", 'w') as output_file:
        with open(input_file, 'r') as file:
            for line in file:
                chr, length = line.split('\t')[0], line.split('\t')[1]
                portion_size = int(length) / number_chunks
                start = 0
                for i in range(1,number_chunks+1,1):
                    stop = start + int(portion_size)
                    output_file.write(("%s\t%s\t%s\n") % (chr, start, stop))
                    start = stop # update value of start
    return
  
if __name__ == "__main__":
    SCRIPT = "chunk_chr_in_index.py"
    # argument set up
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--index_file", type=str, help = "index file for genome (.fai)", required=True)
    parser.add_argument("-p", "--prefix", type=str, help = "prefix for output file", default="fsf")
    parser.add_argument("-n", "--number_chunks", type=int, help = "number of chunks to divide each chr into, default 100", default=100)
    args = parser.parse_args()
    input_file = args.index_file
    prefix = args.prefix
    number_chunks = args.number_chunks
    # run the functions
    print("[+] Parsing index file", input_file,  "and writing chunked bed file")
    chunk_chr_and_write_to_file(input_file, number_chunks, prefix) # chunk each chr up into n_chunks and write to new bed file
