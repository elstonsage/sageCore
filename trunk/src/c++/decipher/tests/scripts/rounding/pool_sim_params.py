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

pool_count = 500
pool_size = 4                # Number of haplotypes per pool.
fraction_missing_data = 0

# - Ambiguity refers to the counts that the allele weights correspond to.
#
fraction_ambiguous_data = .35     

# - haps[x][0] = haplotype frequency.
#   haps[x][1] = list of alleles defining the haplotype.
#   constraints
#     + haplotype frequencies must sum to 1.
#     + haplotypes must be of the same length.
#
haps = [
        [.135, ['x', 'x', 'x']],
        [.135, ['y', 'x', 'z']],
        [.135, ['z', 'z', 'x']],
        [.135, ['z', 'y', 'x']], 
        [.115, ['x', 'x', 'y']],
        [.115, ['y', 'z', 'y']],
        [.115, ['x', 'y', 'z']],
        [.115, ['y', 'y', 'y']] 
       ]
