#!/usr/bin/env python3

#%%
import sys
import argparse
import pandas as pd
import seaborn as sns
    
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

def get_gene_lengths(file):
    chr2gene_length  = {} # holds total length per gene
    with open(file, 'r') as f:
        for line in f:
            if ('exon' in line) or ('five_prime_UTR' in line) or ('three_prime_UTR' in line): # i.e. get genic DNA
                cols = line.split('\t')
                chr, start, stop = cols[0], int(cols[3]), int(cols[4])
                if chr in chr2gene_length.keys(): # check is a protein-coding gene
                    current_length = chr2gene_length[chr]
                    updated_length = current_length + (abs(start - stop)) # get absolute length regardless of strand
                    chr2gene_length[chr] = updated_length
                else:
                    chr2gene_length[chr] = abs(start - stop)
    return(chr2gene_length)

def get_intron_span(chr2genes, gene2exons):
    chr2intron_span = {}
    for chr, chr_genes in chr2genes.items():
        total_intron_lengths = 0
        for gene in chr_genes: # doesn't matter about order of genes
            n = 0
            exons = list(filter(lambda x: x[0].startswith(gene), gene2exons))
            sorted_exons = sorted(exons, key=lambda x: x[1]) # sort genes by start pos
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
                            total_intron_lengths = total_intron_lengths + intron_length
        chr2intron_span[chr] = total_intron_lengths
    return(chr2intron_span)

def get_chr_lengths_and_repeat_span(repeats_bed):
    chr2length, chr2repeat_span = {}, {}
    with open(repeats_bed, 'r') as f:
        for line in f:
            cols = line.split('\t')
            chr, repeat_bases, chr_length = cols[0], int(cols[4]), int(cols[5])
            chr2repeat_span[chr] = repeat_bases
            chr2length[chr] = chr_length
    return(chr2length, chr2repeat_span)

def get_repeat_and_gene_count(repeat_gff, chr2genes):
    chr2gene_count, chr2repeat_count = {}, {} # count number of genes and repeats per chr
    for i in chr2genes:
        genes = chr2genes[i]
        chr2gene_count[i] = len(genes)

    with open(repeat_gff, 'r') as f:
        for line in f:
            chr = line.split('\t')[0]
            if chr in chr2repeat_count:
                repeat_count = chr2repeat_count[chr]
                updated_repeat_count = repeat_count + 1
                chr2repeat_count[chr] = updated_repeat_count
            else:
                chr2repeat_count[chr] = 1
    return(chr2gene_count, chr2repeat_count)

def write_to_output(output_file, chr2length, chr2intron_span, chr2gene_span, chr2repeat_span, chr2repeat_count, chr2gene_count):
    with open(output_file, 'w') as output:
        output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % ("chr", "length", "gene_span", "intron_span", "non_genic_length", "repeat_span", "gene_count", "repeat_count") + "\n")
        for chr, length in chr2length.items():
            if chr in chr2intron_span.keys():
                intron_span = chr2intron_span[chr]
                gene_span = chr2gene_span[chr]
                repeat_span = chr2repeat_span[chr]
                non_genic_length = length - intron_span - gene_span # not including repeats due to overlap
                repeat_count = chr2repeat_count[chr]
                gene_count = chr2gene_count[chr]
                output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chr, length, gene_span, intron_span, non_genic_length, repeat_span, gene_count, repeat_count) + "\n")
    return()
#%%
# Main script


if __name__ == "__main__":
    SCRIPT = "get_gene_and_repeats_per_chr_stats.py"
    # argument set up
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--genes_gff", type=str, help = "gff3 file containing gene information", required=True)
    parser.add_argument("-r", "--repeats_gff", type=str, help = "gff3 file containing repeat information", default="output")
    parser.add_argument("-b", "--repeats_bed", type=str, help = "bed file containing average repeat density per chr", default="output")
    parser.add_argument("-p", "--prefix", type=str, help = "prefix for output file names", default="output")
    args = parser.parse_args()
    gene_file = args.genes_gff
    repeat_gff = args.repeats_gff
    repeats_bed = args.repeats_bed # script can be improved by replacing this requirement with a genome fai file to get chr lengths. Get repeat span from repeat gff3.
    prefix = args.prefix
    chr2gene_span = get_gene_lengths(gene_file)
    chr2genes, mRNA2gene, gene_list = get_gene_set(gene_file)
    gene2exons = read_in_gff(gene_file, mRNA2gene)
    chr2intron_span = get_intron_span(chr2genes, gene2exons)
    chr2length, chr2repeat_span = get_chr_lengths_and_repeat_span(repeats_bed)
    chr2gene_count, chr2repeat_count = get_repeat_and_gene_count(repeat_gff, chr2genes)
    output_file =  str(prefix) + '_gene_repeat_span_and_counts_per_chr.tsv'
    write_to_output(output_file, chr2length, chr2intron_span, chr2gene_span, chr2repeat_span, chr2repeat_count, chr2gene_count)
#%%
# Example input:
#   gene_file = 'Melitaea_cinxia.relabelled.genes.filtered_completeCDS.DivByThree.gff3'
#   repeat_gff = 'Melitaea_cinxia.filteredRepeats.gff'
#   repeats_bed = 'Melitaea_cinxia_repeat_counts_per_chr.tsv'
#%%
## Plotting output:
#   df = pd.read_csv("prefix_gene_repeat_span_and_counts_per_chr.tsv", sep='\t')
#   sns.scatterplot(data=df, x="length", y="gene_count")
#   sns.scatterplot(data=df, x="length", y="gene_length")
#   sns.scatterplot(data=df, x="length", y="repeat_count")
#   sns.scatterplot(data=df, x="length", y="intron_length"
#   sns.scatterplot(data=df, x="length", y="non_genic_length")
#   sns.scatterplot(data=df, x="length", y="repeat_length")
