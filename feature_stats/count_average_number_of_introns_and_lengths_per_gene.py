#!/usr/bin/env python3
import sys
import argparse

#%%
def get_gene_set(file):
    chr2genes, mRNA2gene, gene_list = {}, {}, [] 
    with open(file, 'r') as f:
        for line in f:
            if ('ID=gene:' in line) or ('ID=transcript') in line: # i.e. get gene and mRNA lines
                cols = line.split('\t')
                chr = cols[0]
                feature = cols[2]
                if feature == 'gene': # exclude lncRNA etc
                    gene_id = cols[8].split(';')[0].split(':')[1].strip()
                    gene_list.append(gene_id)
                    if chr in chr2genes.keys():
                        current_genes = chr2genes[chr]
                        current_genes.append(gene_id)
                        chr2genes[chr] = list(set(current_genes))
                    else:
                        chr2genes[chr] = [gene_id]
                if feature == 'mRNA':
                    gene_id = cols[8].split(';')[1].split(':')[1].strip()
                    mRNA_id = cols[8].split(';')[0].split(':')[1].strip()
                    mRNA2gene[mRNA_id] = gene_id
    return(chr2genes, mRNA2gene, gene_list)

def read_in_gff(file, mRNA2gene):
    gene2exons  =[] # gene2coords holds tuples of start and end of each exon for each gene
    with open(file, 'r') as f:
        for line in f:
            if 'exon' in line: # i.e. get gene exon lines, but this includes exons of lnrRNAs
                cols = line.split('\t')
                start, stop = cols[3], cols[4]
                strand = cols[6]
                mRNA_id = cols[8].split(';')[1].split(':')[1].strip()
                if mRNA_id in mRNA2gene.keys(): # check is a protein-coding gene
                    gene_id = mRNA2gene[mRNA_id]
                    if strand == "+":
                        exon_tuple = (gene_id, start, stop)
                        gene2exons.append(exon_tuple)
                    else: # else flip start and stop pos if on -ve strand
                        exon_tuple = (gene_id, stop, start)
                        gene2exons.append(exon_tuple)
    return(gene2exons)

def calculate_average_intron_number_and_length(chr2genes,gene2exons):
    chr2intron_distances, chr2intron_numbers = {}, {}
    for chr, chr_genes in chr2genes.items():
        if len(chr_genes) >= 2:
            total_intron_lengths, total_intron_numbers = [], []
            for gene in chr_genes: # doesn't matter about order of genes
                n = 0
                exons = list(filter(lambda x: x[0].startswith(gene), gene2exons))
                sorted_exons = sorted(exons, key=lambda x: x[1]) # sort genes by start pos
                total_intron_numbers.append(len(exons)-1) # number introns is number exons - 1
                if len(exons) >=2: # need at least two exons
                    for i in sorted_exons:
                        if n== 0:
                            intron_start = sorted_exons[n][2]
                            n += 1
                        else:
                            intron_stop = sorted_exons[n][1]
                            intron_length = int(intron_stop) - int(intron_start)
                            intron_start = sorted_exons[n][2]
                            n += 1
                            if intron_length > 0: # just a sanity check
                                total_intron_lengths.append(intron_length)
        avg_distance = sum(total_intron_lengths)/len(total_intron_lengths)
        avg_intron_number = sum(total_intron_numbers)/len(total_intron_numbers)
        chr2intron_distances[chr] = round(avg_distance, 1)
        chr2intron_numbers[chr] = round(avg_intron_number, 1)
    return(chr2intron_distances, chr2intron_numbers)

# %%
if __name__ == "__main__":
    SCRIPT = "count_number_introns_and_lengths_per_gene.py"
    # argument set up
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--gff_file", type=str, help = "gff3 file containing gene information", required=True)
    parser.add_argument("-p", "--prefix", type=str, help = "prefix for output file names", default="output")
    args = parser.parse_args()
    input_file = args.gff_file
    prefix = args.prefix
    chr2genes, mRNA2gene, gene_list = get_gene_set(input_file)
    gene2exons = read_in_gff(input_file, mRNA2gene)
    chr2intron_distances, chr2intron_numbers = calculate_average_intron_number_and_length(chr2genes,gene2exons)
    lengths_file = prefix +'_average_intron_length_per_chr.tsv'
    with open(lengths_file, 'w') as output:
        for chr, length in chr2intron_distances.items():
            output.write("%s\t%s" % (chr, length) + "\n")
    numbers_file = prefix +'_average_number_introns_per_chr.tsv'
    with open(numbers_file, 'w') as output:
        for chr, number in chr2intron_numbers.items():
            output.write("%s\t%s" % (chr, number) + "\n")
    print("[+] Written output files sueccessfully:")
    print( "    [+] ", lengths_file)
    print( "    [+] ", numbers_file)

# Example input
#  prefix = 'Melitaea_cinxia'
#  input_file = 'Melitaea_cinxia.relabelled.genes.filtered_completeCDS.DivByThree.gff3'
