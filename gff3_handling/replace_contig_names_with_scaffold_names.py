#!/usr/bin/python3
#%%
import argparse
# %%
def make_dict_of_chr_2_scaffold_name(genome_fasta):
    id_2_name = {}
    with open(genome_fasta, 'r') as genome:
        for line in genome:
            if line.startswith('>'):
                if ('chromosome' in line) or ('mitochondrion' in line) or ('linkage' in line): # only chr and mitochondria sequences get re-labelled in gff, other contigs do not
                    fields = line.strip().split(" ")
                    fields[0] = fields[0].replace('>', '')
                    if ('linkage' in line) or ('Danaus' in line) or ('Dryas' in line) or ('marsaeus' in line):
                        chr_name = fields[7] # this is to account for differnet formatting in Dendrolimus_kikuchii.fasta, Dryas_iulia_moderata and Danaus_plexippus_plexippus.fasta headers
                    elif 'menophilus' in line: # this captures Melinaea_menophilus.gff3
                        chr_name = fields[9]
                    else:
                        chr_name = fields[6]
                    getVals = list([val for val in chr_name
                        if val.isalpha() or val.isnumeric()])
                    reform_chr_name = "".join(getVals)
                    if reform_chr_name in id_2_name:
                        print(id_2_name[reform_chr_name], fields[0])
                        assert id_2_name[reform_chr_name] == fields[0], "A chr already exists in dict of id_2_name, its about to be reassigned to a different chr!"
                   # print("chr was:", fields[6], 'now is', reform_chr_name)
                    id_2_name[reform_chr_name] = fields[0]
    return(id_2_name)
    
def convert_gff_chr_2_scaffold_name(input_data, id_2_name):
    prefix = genome_fasta.split("/")[-1].replace('.fasta', '')
    output_file = prefix + '.relabelled.genes.gff3'
    with open(output_file, 'w') as test_file:
        with open(input_data, 'r') as gff:
            for line in gff:
                if line.startswith('#'):
                    test_file.write(line)
                else:
                    fields = line.strip().split('\t')
                    chr_name = fields[0]
                    try:
                        scaffold_name = str(id_2_name[chr_name])
                    except KeyError:
                        scaffold_name = chr_name # used in cases such as mitochondrion where scaffold id is used in gff3 rather than 'mitochondrion''
                    fields[0] = scaffold_name
                    test_file.write('\t'.join(fields) + '\n')
    return()

def make_dict_of_chr_2_scaffold_name_ignoring_set_of_scaffolds(genome_fasta, scaffolds_to_ignore):
    id_2_name = {}
    with open(genome_fasta, 'r') as genome:
        for line in genome:
            if line.startswith('>'):
                if ('chromosome' in line) or ('mitochondrion' in line) or ('linkage' in line): # only chr and mitochondria sequences get re-labelled in gff, other contigs do not
                    fields = line.strip().split(" ")
                    fields[0] = fields[0].replace('>', '')
                    if ('linkage' in line) or ('Danaus' in line) or ('Dryas' in line):
                        chr_name = fields[7] # this is to account for differnet formatting in Dendrolimus_kikuchii.fasta, Dryas_iulia_moderata and Danaus_plexippus_plexippus.fasta headers
                    else:
                        chr_name = fields[6]
                    getVals = list([val for val in chr_name
                        if val.isalpha() or val.isnumeric()])
                    reform_chr_name = "".join(getVals)
                    if fields[0] not in scaffolds_to_ignore:
                        print(fields[0]) # its a scaffold to ignore - skip rest of function!
                        if reform_chr_name in id_2_name:
                            assert id_2_name[reform_chr_name] == fields[0], "A chr already exists in dict of id_2_name, its about to be reassigned to a different chr!"
                        id_2_name[reform_chr_name] = fields[0]
    return(id_2_name)
#%%
if __name__ == "__main__":
    SCRIPT = "replace_contig_names_with_scaffold_names.py"
    # argument set up
    parser = argparse.ArgumentParser(description='This replaces scaffold names with contig/chr names!')
    parser.add_argument("-i", help="This is the input gff3 file")
    parser.add_argument("-f", help="This is the genome fasta file")
    args = parser.parse_args()
    input_data = args.i # gff3 file to have chr_IDs (e.g. '1') to be replaced with unique IDs (e.g. HG996486.1) (e.g. 'Abrostola_tripartita_subset.genes.gff3)
    genome_fasta = args.f # genome fasta

# run functions
print('[+] Processsing ', input_data)
id_2_name = make_dict_of_chr_2_scaffold_name(genome_fasta)
print(id_2_name)
convert_gff_chr_2_scaffold_name(input_data, id_2_name)

exit()
#%%
# manual run for Dryas_iulia_moderata
#input_data = 'Dryas_iulia_moderata.genes.gff3'
#genome_fasta = 'Dryas_iulia_moderata.fasta'
#Dryas_iulia_moderata_scaffolds_to_ignore = ['JAHESG010000005.1', 'JAHESG010000019.1', 'JAHESG010000013.1', 'JAHESG010000014.1', 'JAHESG010000007.1','JAHESG010000015.1']

#id_2_name = make_dict_of_chr_2_scaffold_name_ignoring_set_of_scaffolds(genome_fasta, Dryas_iulia_moderata_scaffolds_to_ignore)
#convert_gff_chr_2_scaffold_name(input_data, id_2_name)
