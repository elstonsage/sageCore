#ifndef TDTEX_TESTS_H
#define TDTEX_TESTS_H

#include <vector>
#include "globals/SAGEConstants.h"
#include "numerics/fmatrix.h"
#include "numerics/functions.h"
#include "numerics/mt.h"
#include "tdtex/TransmissionTables.h"

namespace SAGE  {
namespace TDTEX {

extern MersenneTwister mt;

/// Represents a chi-square statistic and its associated p-value.
struct chi_squared
{
  chi_squared() : df((size_t)-1), xx(SAGE::QNAN) { }
  chi_squared(size_t _df, double _xx) : df(_df), xx(_df ? _xx : SAGE::QNAN) { }
  chi_squared(const chi_squared & other) : df(other.df), xx(other.xx) { }

  chi_squared& operator=(const chi_squared & other) 
  {
    if(this != &other)
    {
      df = other.df;
      xx = other.xx;
    }
    
    return *this;
  }
  
  size_t df;
  double xx;
  double value()  const { return xx; }
  double pvalue() const { return xx < 10e-6 && df ? 1.0 : 1.0 - chi_square_cdf(df, xx); }

  void dump() const
  {
    std::cout << "xx = " << xx << ", df = " << df << ", pvalue = " << pvalue() << std::endl;
  }
};

/// Represents a simulated p-value.
struct simulated_pvalue
{
  simulated_pvalue(double p, double se) : p(p), se(se) { }
  double p;
  double se;
  double pvalue() const { return p; }
  double standard_error() const { return se; }
};

/// McNemar Chi squared
chi_squared mcnemar_statistic(const FortranMatrix<size_t>& table, bool continuity_correction = false);

/// Marginal homogeneity chi square
chi_squared pearson_marginal_homogeneity_statistic(const FortranMatrix<size_t>& table, bool no_diagonal = true);

chi_squared pearson_marginal_homogeneity_statistic(const std::vector<size_t>& margins1, const std::vector<size_t>& margins2, bool no_diagonal = true);
                                                   
/// Compute McNamara probability function for a given table
double mcn_prob      (const std::vector<int>& vector1, const std::vector<int>& vector2);
double mcn_prob2     (const std::vector<int>& count,   const std::vector<int>& total);
double mcn_log_prob  (const std::vector<int>& vector1, const std::vector<int>& vector2);
double mcn_log_prob2 (const std::vector<int>& count,   const std::vector<int>& total);

class McNemarExact
{
  public:
  
    static double mcnemar_exact(const FortranMatrix<size_t>& table);

  private:

    static double mcn_perm_down(const std::vector<int>& count, const std::vector<int>& total, const double, double, int, int, int);
    static double mcn_perm_up(const std::vector<int>& count, const std::vector<int>& total, const double, double, int, int, int);
    static int multiplicity(const std::vector<int>& count, const std::vector<int>& total);
};

class McNemarMonteCarlo
{
  public:

    typedef std::vector< std::vector<double> > vector2D;

    static simulated_pvalue mcnemar_monte_carlo(const FortranMatrix<size_t>& table, int nrep, int nbatch);
    
  private:

    static vector2D generate_binom_table(const vector<int>& total);

    static void cmc_xx_permute(const vector<int>& total, vector<int>& vector, const vector2D& table);
};

class MargHomoMonteCarlo
{
  public:

    typedef FortranMatrix< vector<double> > binom_table_t;

    static simulated_pvalue mh_monte_carlo(const FortranMatrix<size_t>& table, int nrep, int nbatch);

  private:
  
    static binom_table_t generate_binom_table(const FortranMatrix<size_t>& table);
  
    static void cmc_xx_permute(
      const FortranMatrix<size_t> & table,
            vector<size_t>        & margins1,
            vector<size_t>        & margins2,
      const binom_table_t         & binom_table);
                    

};

} // End namespace TDTEX
} // End namespace SAGE

#endif
