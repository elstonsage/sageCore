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

pool_count = 100
pool_size = 4                # Number of haplotypes per pool.
fraction_missing_data = .05
fraction_ambiguous_data = 0

# - haps[x][0] = haplotype frequency.
#   haps[x][1] = list of alleles defining the haplotype.
#   constraints
#     + haplotype frequencies must sum to 1.
#     + haplotypes must be of the same length.
#
haps = [
        [.135, ['0', '0', '0']],
        [.135, ['1', '0', '0']],
        [.135, ['0', '1', '0']],
        [.135, ['1', '1', '0']], 
        [.115, ['0', '0', '1']],
        [.115, ['1', '0', '1']],
        [.115, ['0', '1', '1']],
        [.115, ['1', '1', '1']] 
       ]
