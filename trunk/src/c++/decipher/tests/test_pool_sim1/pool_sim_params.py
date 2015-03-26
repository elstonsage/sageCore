#============================================================================
# File:      pool_sim_params.py
#                                                                          
# Author:    Dan Baechle
#                                                                          
# History:   Created 1/25/6
#                                                                          
# Notes:     Configuration file for pool simulation program which samples a 
#            population of haplotypes w. the given frequencies and produces
#            a list of pool genotypes corresponding to the sample.
#
#============================================================================

pool_count = 2000
pool_size = 4                # Number of haplotypes per pool.
fraction_missing_data = 0

# - haps[x][0] = haplotype frequency.
#   haps[x][1] = list of alleles defining the haplotype.
#   constraints
#     + haplotype frequencies must sum to 1.
#     + haplotypes must be of the same length.
#
haps = [
        [.135, ['A', 'A', 'A']],
        [.135, ['a', 'A', 'A']],
        [.135, ['A', 'a', 'A']],
        [.135, ['a', 'a', 'A']], 
        [.115, ['A', 'A', 'a']],
        [.115, ['a', 'A', 'a']],
        [.115, ['A', 'a', 'a']],
        [.115, ['a', 'a', 'a']] 
       ]
