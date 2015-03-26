#include <iostream>
#include <math.h>
#include <limits>
#include <functional>
#include "numerics/functions.h"
#include "tdtex/Tests.h"
#include "tdtex/Numerics.h"
#include "numerics/functions.h"
#include "tdtex/Tests.h"
#include "tdtex/Numerics.h"


namespace SAGE  {
namespace TDTEX {

// McNemar Chi squared
chi_squared mcnemar_statistic(const TransmissionTable::Matrix & table, bool continuity_correction)
{
  size_t df = 0;
  double xx = 0.0;

  assert( table.rows() == table.cols() );

  for(size_t i = 0; i < table.rows(); ++i)
  {
    for(size_t j = 0; j < i; ++j)
    {
      if(table(i, j) || table(j, i))
      {
        double diff = (double)table(i, j) - (double)table(j, i);

        if(continuity_correction)
          diff = abs(diff) - 1.0;

        xx += (diff * diff) / (table(i, j) + table(j, i));

        ++df;
      }
    }
  }

  return chi_squared(df, xx);
}

// Marginal homogeneity chi square
chi_squared pearson_marginal_homogeneity_statistic(const TransmissionTable::Matrix& table, bool no_diagonal)
{
  assert( table.rows() == table.cols() );

  std::vector<size_t> margins1(table.rows(), 0),
                      margins2(table.rows(), 0);

  for(size_t i = 0; i < table.rows(); ++i)
  {
    for(size_t j = 0; j < table.rows(); ++j)
    {
      if(!no_diagonal || i!=j)
      {
        margins1[i] += table(i, j);
        margins2[j] += table(i, j);
      }
    }
  }

  return pearson_marginal_homogeneity_statistic(margins1, margins2, no_diagonal);
}

// Marginal homogeneity chi square
chi_squared pearson_marginal_homogeneity_statistic(const std::vector<size_t>& margins1,
                                                   const std::vector<size_t>& margins2, 
                                                   bool no_diagonal)
{
  size_t df = 0;
  double xx = 0.0;

  size_t n = margins1.size();

  for(size_t i=0; i < n; ++i)
  {
    if(margins1[i] || margins2[i])
    {
      int diff = margins1[i] - margins2[i],
          sum  = margins1[i] + margins2[i];

      xx += double(diff*diff)/sum;

      ++df;
    }
  }

  if(df && no_diagonal)
    xx *= (df-1.0)/df;

  if(df)
    --df;

  if(!df)
    return chi_squared(0, std::numeric_limits<double>::quiet_NaN());

  return chi_squared(df, xx);
}

double mcn_prob(const std::vector<int>& vector1, const std::vector<int>& vector2)
{
/* Compute McNamara probability function for a given table */
  size_t n = vector1.size();
  double p = 1.0;
  for(size_t i=0; i<n; i++)
    if(vector1[i] || vector2[i])
      p *=   fact(vector1[i]+vector2[i])
        /  ( fact(vector1[i]) * fact(vector2[i]) )
        * pow(0.5,(double)vector1[i]+vector2[i]);
  return p;
}

double mcn_prob2(const std::vector<int>& count, const std::vector<int>& total)
{
/* Compute McNemar probability function for a given table */
  size_t n = count.size();

  double p = 1.0;
  for(size_t i=0; i<n; i++)
    if(total[i])
      p *=   fact(total[i])
        /  ( fact(total[i]-count[i]) * fact(count[i]) )
        * pow(0.5,(double)total[i]);
  return p;
}

double mcn_log_prob(const std::vector<int>& vector1, const std::vector<int>& vector2)
{
  /* Compute log McNamara probability function for a given table */
  size_t n = vector1.size();
  double p = 0.0;
  for(size_t i=0; i<n; i++)
    if(vector1[i] || vector2[i])
    {
      p += log_fact(vector1[i]+vector2[i])
        -  log_fact(vector1[i])
        -  log_fact(vector2[i])
        +  log(0.5)*(vector1[i]+vector2[i]);
    }
  return p;
}

double mcn_log_prob2(const std::vector<int>& count, const std::vector<int>& total)
{
  /* Compute log McNemar probability function for a given table */
  size_t n = count.size();
  double p = 0.0;
  for(size_t i=0; i<n; i++)
    if(total[i])
      p += log_fact(total[i])
        -  log_fact(total[i]-count[i])
        -  log_fact(count[i])
        +  log(0.5)*total[i];
  return p;
}

int 
McNemarExact::multiplicity(const std::vector<int>& count, const std::vector<int>& total)
{
  int n = count.size(),
      m = 1;

  for(int j = 0; j < n; j++)
    if(count[j] * 2 != total[j])
      m *= 2;

  return m;
}

double 
McNemarExact::mcnemar_exact(const FortranMatrix<size_t>& table)
{
  /* Kick off the permutations for exact test for McNemar prob. tables */

  int n = table.rows();
  std::vector<int> count(n * (n + 1) / 2),
                   total(n * (n + 1) / 2);

  n = 0;

  for(int i = 0; i < (int)table.rows(); ++i)
    for(int j = 0; j < i; ++j)
      if(table(i, j) || table(j, i))
      {
        total[n] = table(i, j) + table(j, i);
        count[n] = std::min(table(i, j), table(j, i));
        n++;
      }

  if(!n)
    return std::numeric_limits<double>::quiet_NaN();

  count.resize(n);
  total.resize(n);

  double obs   = exp(mcn_log_prob2(count, total)) * 1.000001,
         asy_p = mcnemar_statistic(table, true).pvalue(),
         p     = std::numeric_limits<double>::quiet_NaN();

  if(asy_p < 0.95)
  {
    for(int i = 0; i < n; ++i)
      count[i] = total[i]/2;

    p = 1.0 - mcn_perm_down(count, total, obs, 0.0, n, 0, 2);
  }
  else
  {
    for(int i = 0; i < n; ++i)
      count[i] = 0;

    p = mcn_perm_up(count, total, obs, 0.0, n, 0, 2);
  }

  return std::max(p, 0.0);
}

double 
McNemarExact::mcn_perm_down(
  const std::vector<int> & count,
  const std::vector<int> & total, 
  const double             obs,
        double             current, 
        int                df, 
        int                level, 
        int                first)
{
  /* actual work done for permutations done here */

  int    n   = count.size(),
         f   = first;
  double sum = 0.0;

  /* Check if we're done (ie, everything is constrained) */
  if(df < 1)
    return sum;

  int i = n - df;

  /* Create a new working vector */
  std::vector<int> new_v = count;

  /* Adjust vector if it is not at symmetry */
  if(new_v[i] && new_v[i] != total[i] / 2)
  {
    --new_v[i];
    current *= double(new_v[i] + 1.0) / (total[i] - new_v[i]);
  }

  for(first = 1; new_v[i] >= 0; --new_v[i], first = 0, f = 0)
  {
    /* Check if this one counts */
    if( f==2 || (!f && !first) )
    {
      /* For accuracy calculate the prob. for a given table once   */
      /* and then transform to reach other tables.  This preserves */
      /* numerical accuracy                                        */

      if(current == 0.0)
        current = exp(mcn_log_prob2(new_v, total));
      else
        current *= (new_v[i] + 1.0) / (total[i] - new_v[i]);

      if(current <= obs)
      {
        break;
      }

      // Compute multiplicity
      int m = multiplicity(new_v,total);

      sum += m*current;
    }

    sum += mcn_perm_down(new_v, total, obs, current, df - 1, level + 1, first);
  }

  return sum;
}

double 
McNemarExact::mcn_perm_up(
  const std::vector<int> & count,
  const std::vector<int> & total, 
  const double             obs,
        double             current, 
        int                df, 
        int                level, 
        int                first)
{
  /* actual work done for permutations done here */

  int    n   = count.size(),
         f   = first;
  double sum = 0.0;

  /* Check if we're done (ie, everything is constrained) */
  if(df < 1)
    return sum;

  int i = n - df;

  /* Create a new working vector */
  std::vector<int> new_v = count;

  for(first = 1; new_v[i] <= total[i] / 2; ++new_v[i], first = 0, f = 0)
  {
    /* Check if this one counts */
    if(f == 2 || (!f && !first))
    {
      /* For accuracy calculate the prob. for a given table once   */
      /* and then transform to reach other tables.  This preserves */
      /* numerical accuracy                                        */

      if(current == 0.0)
        current = exp(mcn_log_prob2(new_v, total));
      else
        current *= (total[i] - new_v[i] + 1.0) / new_v[i];

      if(current > obs)
      {
        break;
      }

      // Compute multiplicity
      int m = multiplicity(new_v, total);

      sum += m*current;
    }

    sum += mcn_perm_up(new_v, total, obs, current, df - 1, level + 1, first);
  }

  return sum;
}

#define BINOM_MAX 800

//================================================================
//
//  generate_binom_table(...)
//
//================================================================
McNemarMonteCarlo::vector2D 
McNemarMonteCarlo::generate_binom_table(const std::vector<int>& total)
{
  size_t n = total.size();
  vector2D table(n);

  for(int i = 0; i < (int)n; i++)
  {
    size_t size = std::min(BINOM_MAX, total[i]);

    if(!size)
      continue;

    table[i].resize(size);

    double tfact = log_fact(total[i]),
           halft = log(0.5) * total[i];

    double sum = 0;

    for(int j = 0; j < (int)size; j++)
    {
      double prob = exp(tfact - log_fact(j) - log_fact(total[i] - j) + halft);
      sum += prob;
      table[i][j]=sum;
    }
  }
  return table;
}

//==================================================================
//
//  cmc_xx_permute(...)
//
// Permute the total counts total[i] into a matched vector[i] randomly
// chosen from [0,total[i] ]
//
//==================================================================
void 
McNemarMonteCarlo::cmc_xx_permute(const vector<int>& total, vector<int>& vector, const vector2D& table)
{
  int n = vector.size();

  for(int i = 0; i<n; i++)
  {
    int size = std::min(BINOM_MAX, total[i]);

    if(!size || !table[i].size())
    {
      vector[i] = 0;
      continue;
    }

    double r = mt.uniform_real();

    int j;

    for(j = 0; j < size; ++j)
      if(r < table[i][j])
        break;
    vector[i] = j;
  }
}

//============================================================
//
//  mcnemar_monte_carlo(...)
//
//============================================================
simulated_pvalue
McNemarMonteCarlo::mcnemar_monte_carlo(const FortranMatrix<size_t>& table, int nrep, int nbatch)
{
  int n = table.rows();
  std::vector<int> count(n*(n+1)/2);
  std::vector<int> total(n*(n+1)/2);

  n = 0;

  for(int i = 0; i < (int)table.rows(); ++i)
  {
    for(int j = 0; j < i; ++j)
    {
      if(table(i, j) || table(j, i))
      {
        total[n] = table(i, j) + table(j, i);
        count[n] = std::min(table(i, j), table(j, i));
        n++;
      }
    }
  }

  if(!n)
  {
    return simulated_pvalue(std::numeric_limits<double>::quiet_NaN(),
                            std::numeric_limits<double>::quiet_NaN());
  }

  count.resize(n);
  total.resize(n);

  vector2D binom_table = generate_binom_table(total);

  double pobs = exp(mcn_log_prob2(count, total)) * 1.00001;

  int    sig  = 0;
  double pb   = 0,
         pbsq = 0;

  for(int batch = 0; batch < nbatch; ++batch)
  {
    int sigb = 0;

    for(int rep = 0; rep < nrep; ++rep)
    {
      cmc_xx_permute(total, count, binom_table);

      double pcurrent = exp( mcn_log_prob2(count, total));

      if(pcurrent <= pobs)
        ++sigb;
    }
    double pt = (sigb + 1.0) / (nbatch + 1.0);

    sig  += sigb;
    pb   += pt;
    pbsq += pt*pt;
  }

  double p  = (sig + 1.0) / (nrep * nbatch + 1.0),
         se = (double)pbsq / nbatch - (pb / nbatch) * (pb / nbatch);

  return simulated_pvalue(p, se);
}


MargHomoMonteCarlo::binom_table_t 
MargHomoMonteCarlo::generate_binom_table(const FortranMatrix<size_t>& table)
{
  size_t n = table.rows();
  binom_table_t binom_table(n);

  for(int i = 0; i < (int)n; i++)
  {
    for(int j = 0; j < (int)i; ++j)  
    {
      size_t total = table(i,j)+table(j,i),
             size  = std::min<int>(BINOM_MAX, total);

      if(!size)
        continue;

      binom_table(i,j).resize(size);

      double tfact = log_fact(total),
             halft = log(0.5) * total;

      double sum = 0;

      for(int k = 0; k < (int)size; k++)
      {
        double prob = exp(tfact - log_fact(k) - log_fact(total - k) + halft);
        sum += prob;
        binom_table(i,j)[k]=sum;

      }
    }
  }
  return binom_table;
}

void 
MargHomoMonteCarlo::cmc_xx_permute(
  const FortranMatrix<size_t> & table,
        vector<size_t>        & margins1, 
        vector<size_t>        & margins2, 
  const binom_table_t         & binom_table)
{
  int n = table.rows();

  for(int i = 0; i < n; i++)
    margins1[i] = margins2[i] = 0;

  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j < i; j++)
    {
      int size = std::min<int>(BINOM_MAX, table(i,j)+table(j,i));

      if(!size || !binom_table(i,j).size())
        continue;

      double r = mt.uniform_real();

      int k;

      for(k = 0; k < size; ++k)
        if(r < binom_table(i,j)[k])
          break;

      margins1[i] += k;
      margins2[j] += k;
      margins1[j] += size - k;
      margins2[i] += size - k;
    }
  }
}

simulated_pvalue
MargHomoMonteCarlo::mh_monte_carlo(const FortranMatrix<size_t>& table, int nrep, int nbatch)
{
  int n = table.rows();
  std::vector<size_t> margins1(n),
                      margins2(n);

  binom_table_t binom_table = generate_binom_table(table);

  double xxobs = pearson_marginal_homogeneity_statistic(table).xx;

  if(!finite(xxobs))
    return simulated_pvalue(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());

  int    sig  = 0;
  double pb   = 0,
         pbsq = 0;

  for(int batch = 0; batch < nbatch; ++batch)
  {
    int sigb = 0;

    for(int rep = 0; rep < nrep; ++rep)
    {
      cmc_xx_permute(table, margins1, margins2, binom_table);

      double xxcurrent = pearson_marginal_homogeneity_statistic(margins1,margins2).xx;

      if(xxcurrent >= xxobs)
        ++sigb;
    }

    double pt = (sigb + 1.0) / (nbatch + 1.0);

    sig  += sigb;
    pb   += pt;
    pbsq += pt*pt;
  }

  double p  = (sig + 1.0) / (nrep * nbatch + 1.0),
         se = (double)pbsq / nbatch - (pb / nbatch) * (pb / nbatch);

  return simulated_pvalue(p, se);
}

} // End namespace TDTEX
} // End namespace SAGE
