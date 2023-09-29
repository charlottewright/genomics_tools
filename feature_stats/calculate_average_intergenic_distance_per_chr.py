#!/usr/bin/env python3
import argparse
#%%
def read_in_gff(file):
    chr2genes, gene2start, gene2stop  = {}, {}, {}
    with open(file, 'r') as f:
        for line in f:
            if 'ID=gene:' in line: # i.e. get gene lines
                cols = line.split('\t')
                chr = cols[0]
                feature = cols[2]
                start, end = cols[3], cols[4]
                strand = cols[6]
                gene_id = cols[8].split(';')[0].split(':')[1].strip()
                if feature == 'gene': # exclude lncRNA etc
                    if chr in chr2genes.keys():
                        current_genes = chr2genes[chr]
                        current_genes.append(gene_id)
                        chr2genes[chr] = list(set(current_genes))
                    else:
                        chr2genes[chr] = [gene_id]
                    gene2stop[gene_id] = end
                    gene2start[gene_id] = start
    return(chr2genes, gene2start, gene2stop)

def calculate_distance_per_chr(chr2genes, gene2start, gene2stop):
    n = 0
    chr2intergenic_distance = {}
    for chr, chr_genes in chr2genes.items():
        total_intergenic_distances = []
        gene_starts = [(k, v) for k, v in gene2start.items()]
        sorted_genes = sorted(gene_starts, key=lambda x: x[1]) # sort genes by start pos
        genes = [x[0] for x in sorted_genes]
        sorted_chr_genes = [g for g in genes if g in chr_genes]
        if len(sorted_chr_genes) >=2 : # need at least one distance
            for gene in sorted_chr_genes:
                if n == 0:
                    intergenic_start = gene2stop[gene]
                    n =+ 1
                else:
                    intergenic_stop = gene2start[gene]
                    intergenic_distance = int(intergenic_stop) - int(intergenic_start)
                # print(gene, ':' ,intergenic_stop, intergenic_start, ' gives ', intergenic_distance)
                    intergenic_start = gene2stop[gene]
                    n =+ 1
                    if intergenic_distance > 0: # i.e. if genes are no overlapping (which would give a negative value)
                        total_intergenic_distances.append(intergenic_distance)
            avg_ditance = sum(total_intergenic_distances)/len(total_intergenic_distances)
            chr2intergenic_distance[chr] = int(avg_ditance)
    return(chr2intergenic_distance)

# %%
if __name__ == "__main__":
    SCRIPT = "calculate_intergenic_distance.py"
    # argument set up
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--gff_file", type=str, help = "gff3 file containing gene information", required=True)
    parser.add_argument("-p", "--prefix", type=str, help = "prefix for output file names", default="output")
    args = parser.parse_args()
    input_file = args.gff_file
    prefix = args.prefix
    chr2genes, gene2start, gene2stop = read_in_gff(input_file)
    chr2intergenic_distance = calculate_distance_per_chr(chr2genes, gene2start, gene2stop)
    output_file = prefix +'_average_intergenic_distance_per_chr.tsv'
    with open(output_file, 'w') as output:
        for chr, distance in chr2intergenic_distance.items():
            output.write("%s\t%s" % (chr, distance) + "\n")
    print("[+] Written output file", output_file, "successfully!")

#  Example input
#  prefix = 'Melitaea_cinxia'
#  input_file = 'Melitaea_cinxia.relabelled.genes.filtered_completeCDS.DivByThree.gff3'
