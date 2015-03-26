To determine best way to resolve situation in which allele weights yeild the
wrong total number of alleles at a locus, 50 replicates of data were
analyzed by four different versions of Decipher.  Each of the Decipher
versions utilized a different method of resolving the descrepancy between
allele count and pool size.

Results were analyzed by calculating i sub f (Excoffier and Slatkin. 1995)
using the .../pool_sim_script/calc_if script.  Simulation parameters are given in
pool_sim_params.py.  The script driver was used to run the
.../pool_sim_script/pool_simulation script, the various Decipher versions
(see reconcile_counts() in pool.cpp) and calc_lf.

results were as follows,

  Method        if coef.
  ------        --------
  RANDOM        .832017
  GREATEST      .799599
  LEAST         .849941
  OVER_UNDER    .862514