import sys
def parse_orthogroups_file(orthogroups_file):
        with open(orthogroups_file, 'r') as orthogroups:
                all_tolIDs, orthodict = [], {}
                for line in orthogroups: # for every orthogroup
                        orthoID, seqlist = line.rstrip("\n").split(": ")[0], line.rstrip("\n").split(": ")[1].split(" ") # get list of sequences in orthogroup
                        tolID_list = [] # make a list to store tolIDs in orthogroup
                        for seq in seqlist: # for every seq
                                tolID = tolID = ".".join(seq.split(".")[0:1]) # get the tolID
                                tolID_list.append(tolID) # append this to list of tolIDs for this orthogroup
                                all_tolIDs.append(tolID) # and to a list of tolIDs seen in ALL orthogroups
                        orthodict[orthoID] = tolID_list # store in dict where orthoID is key and tolID list is value
                tolIDs = set(sorted(all_tolIDs)) # get unique list of every tolID/species in orthogroups file
        return tolIDs, orthodict
#%%
def find_fuzzy_singletons(orthodict, tolIDs, threshold):
        suffix = 'fuzzy_singletons.'
        with open(suffix + str(threshold) + ".tsv", "w") as outfile:
                for ortho, spp_list in orthodict.items():
                        #print('ortho is:', ortho, 'spp list is:', len(spp_list))
                      #  print('number of spp:', len(spp_list), 'set of spp:', len(set(spp_list)))
                        if (len(spp_list) <= (len(tolIDs) + threshold)) & (len(set(spp_list)) >= (len(tolIDs) - threshold)):
                                outfile.write("%s\t%s\t%s\n" % (ortho, set(spp_list), len(spp_list)))
        return()

def find_multicopy_genes(orthodict, tolIDs, threshold): # this finds all genes where they have >=2 copies in a number of species (where number is "all spp - threshold") e.g. a threshold of 10 means duplicated or more in all but 10 species
        suffix = 'multi_copy_genes.' 
        with open(suffix + str(threshold) + ".tsv", "w") as outfile:    
                for ortho, spp_list in orthodict.items():
                        dup_spp_count =0
                        for spp in set(spp_list):
                                count = spp_list.count(spp)
                                if count >= 2:
                                        dup_spp_count = dup_spp_count+1
                        if dup_spp_count >= (len(tolIDs)- threshold):
                                outfile.write("%s\t%s\t%s\n" % (ortho, set(spp_list), dup_spp_count))
        return()
#%%
orthogroups_file, threshold = sys.argv[1], int(sys.argv[2]) # parse args
#orthogroups_file, threshold = 'Orthogroups_test.tsv', 20
#orthogroups_file = 'Orthogroups.test.txt'
#threshold = 10
tolIDs, orthodict = parse_orthogroups_file(orthogroups_file)
find_fuzzy_singletons(orthodict, tolIDs, threshold)
find_multicopy_genes(orthodict, tolIDs, threshold)
