#============================================================================
# File:      sim_params.py
#                                                                          
# Author:    Dan Baechle
#                                                                          
# History:   Created 11/14/3
#                                                                          
# Notes:     Configuration file for simulation program which samples a 
#            population of haplotypes w. the given frequencies and produces
#            a list of genotypes corresponding to the sample.
#
#============================================================================

SAGE =     1
HAPFREQS = 2
ARLEQUIN = 3

format = SAGE


ind_count = 200
fraction_missing_data = 0

# - haps[x][0] = haplotype frequency.
#   haps[x][1] = list of alleles defining the haplotype.
#   constraints
#     + haplotype frequencies must sum to 1.
#     + haplotypes must be of the same length.
#
haps = [[.27, ['A', 'A']],
        [.27, ['a', 'A']],
        [.23, ['A', 'a']],
        [.23, ['a', 'a']], ]
