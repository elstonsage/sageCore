#ifndef SIBPAL_REGRESSION_RESULT_H
#define SIBPAL_REGRESSION_RESULT_H

//=============================================================================
// File:    regress_result.h
//
// Author:  Kevin Jacobs
//
// History: Version 0.0 Initial implementation
//
// Notes:   This file contains definition for following data structures.
//            class parameter_estimate
//            struct result
//            struct F_result_type
//
// Copyright (c) 2001   R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "sibpal/regress_params.h"

namespace SAGE   {
namespace SIBPAL {

class parameter_estimate
{
  public:

    enum parameter_type { tNONE     = 0, tINTERCEPT, tMARKER, tCOVARIATE };
    enum direction_type { dTWOSIDED = 0, dLEFTSIDE, dRIGHTSIDE };

    parameter_estimate();

    void            clear();
    void            clear_simulation();

    void            set_test_direction(direction_type d);
    void            set_value(double v);
    void            set_variance(double ss);
    void            set_expected_value(double e);
    void            set_degrees_of_freedom(size_t df);
    void            set_adjusted_tvalue(double t);
    void            set_adjusted_degrees_of_freedom(size_t df);
    void            set_adjusted_pvalue(double p);

    void            add_simulation_result(double v, double ss);

    direction_type  test_direction()                const;

    double          value()                         const;
    double          variance()                      const;
    double          standard_error()                const;
    double          expected_value()                const;
    size_t          degrees_of_freedom()            const;
    double          tvalue()                        const;
    double          raw_tvalue()                    const;
    double          adjusted_tvalue()               const;
    size_t          adjusted_degrees_of_freedom()   const;

    // P-values

    //      Returns the first available (non-QNAN) in the following list:
    //          adjusted  p-value
    //          adjusted  t-value significance
    //          raw p-value
    //          QNAN (if none of the previous are available)

    double          pvalue()                        const;

    //      Returns the p-value of the raw_tvalue (if df are available);
    //          QNAN otherwise.

    double          raw_pvalue()                    const;

    //      Returns the p-value of the adjusted tvalue (if df are
    //          available); QNAN otherwise.

    double          adjusted_tvalue_significance()  const;

    //      Returns the adjusted pvalue data member

    double          raw_adjusted_pvalue()           const;

    //      Returns the adjusted pvalue data member or if not available, the
    //      significance of the adjusted t-value

    double          adjusted_pvalue()               const;

    //      Returns empirical significance (from simulation) of the
    //          estimate, which is defined to be the ratio of the
    //          significant replicates over the total.  If no replicates
    //          have been completed then it returns QNAN.

    double          empirical_pvalue()              const;

    //      Returns the total number of simulation replicates

    size_t          total_replicates()              const;

    //      Returns the total number of significant simulation replicates

    size_t          significant_replicates()        const;

    //      Returns the total number of same simulation replicates

    size_t          same_replicates()               const;

  private:

    double          compute_pvalue(double t, size_t df, direction_type dir) const;

    double          my_expected_value;
    double          my_value;
    double          my_variance;
    size_t          my_degrees_of_freedom;

    double          my_adjusted_tvalue;
    size_t          my_adjusted_degrees_of_freedom;

    double          my_adjusted_pvalue;

    size_t          my_significant_replicates;
    size_t          my_same_replicates;
    size_t          my_total_replicates;

    direction_type  my_direction;
};

struct reg_result
{
  public:

    enum result_type { none, intercept, param };

    reg_result()
    {
      clear();
    }

    reg_result(const independent_variable& p, size_t i)
    {
      my_type      = param;
      my_pair_type = MIXED;
      my_parameter = p;
      my_index     = i;
    }

    reg_result(const reg_result& r)
    {
      estimate     = r.estimate;
      my_type      = r.my_type;
      my_pair_type = r.my_pair_type;
      my_parameter = r.my_parameter;
      my_index     = r.my_index;
    }

    void clear()
    {
      my_index     = (size_t)-1;
      my_parameter = independent_variable();
      my_type      = none;
      my_pair_type = MIXED;
      estimate     = parameter_estimate();
    }

    void set_intercept()
    {
      clear();
      my_type = intercept; 
    }

    void set_pair_type(pair_type pt)
    {
      my_pair_type = pt; 
    }

    void set_parameter(const independent_variable& p, size_t i)
    {
      clear();
      my_index     = i;
      my_type      = param;
      my_parameter = p;
    }

    result_type                      type() const { return my_type;  }
    pair_type               get_pair_type() const { return my_pair_type;  }
    size_t                          index() const { return my_index; }
    const independent_variable& reg_param() const { return my_parameter; }

    parameter_estimate      estimate;

  protected:

    independent_variable    my_parameter;
    result_type             my_type;
    pair_type               my_pair_type;
    size_t                  my_index;
};

typedef vector<reg_result> result_vector;

struct F_result_type
{
  public:

    F_result_type()
    {
      clear();
    }

    F_result_type(double Fs, double Fp, size_t d1, size_t d2)
    {
      my_F_statistic = Fs;
      my_F_pvalue    = Fp;
      my_df1         = d1;
      my_df2         = d2;
    }

    F_result_type(const F_result_type& r)
    {
      my_F_statistic = r.my_F_statistic;
      my_F_pvalue    = r.my_F_pvalue;
      my_df1         = r.my_df1;
      my_df2         = r.my_df2;
      my_valid       = r.my_valid;
    }

    void clear()
    {
      my_F_statistic = 0.;
      my_F_pvalue    = 1.;
      my_df1         = (size_t)-1;
      my_df2         = (size_t)-1;
      my_valid       = false;

      my_significant_replicates = my_same_replicates = my_total_replicates = 0;
    }

    void set_F_statistic(double Fs)
    {
      my_F_statistic = Fs; 
    }

    void set_F_pvalue(double Fp)
    {
      my_F_pvalue = Fp; 
    }

    void set_df1(size_t d1)
    {
      my_df1 = d1;
    }

    void set_df2(size_t d2)
    {
      my_df2 = d2;
    }

    void validate()
    {
      my_valid = true;
    }

    void invalidate()
    {
      clear();
      my_valid = false;
    }

    void add_simulation_result(double F_pval)
    {
      if( !finite(F_pval) )
        return;

      ++my_total_replicates;

      double curr_p = my_F_pvalue;
      bool sig = false;
      bool same = false;

      if( F_pval < curr_p )
        sig = true;
      else if( F_pval == curr_p )
        same = true;

      if( sig )
        ++my_significant_replicates;
      else if( same )
        ++my_same_replicates;

      return;
    }

    double get_F_statistic() const { return my_F_statistic;  }
    double get_F_pvalue()    const { return my_F_pvalue;  }
    size_t get_df1()         const { return my_df1;  }
    size_t get_df2()         const { return my_df2;  }
    bool   is_valid()        const { return my_valid; }

    double get_F_empirical_pvalue() const
    {
      if( my_significant_replicates + my_same_replicates + my_total_replicates == 0 )
        return std::numeric_limits<double>::quiet_NaN();

      double num   = (double)my_significant_replicates+1.0+(0.5*(double)my_same_replicates);
      double denom = (double)my_total_replicates+1.0;

      return num / denom;
    }

    size_t get_significant_replicates() const { return my_significant_replicates; }
    size_t get_same_replicates()        const { return my_same_replicates; }
    size_t get_total_replicates()       const { return my_total_replicates; }

  protected:

    double                  my_F_statistic;
    double                  my_F_pvalue;
    size_t                  my_df1;
    size_t                  my_df2;
    bool                    my_valid;

    size_t                  my_significant_replicates;
    size_t                  my_same_replicates;
    size_t                  my_total_replicates;
};

struct reg_variances
{
  reg_variances()
   : residual(QNAN), total(QNAN), sum_residual(QNAN), diff_residual(QNAN)
  {}

  double                    residual;       // residule variance from gls
  double                    total;          // total variance from gls
  double                    sum_residual;   // residule variance from gls using sum
  double                    diff_residual;  // residule variance from gls using diff
};

class regression_results
{
  public:

    friend class TraitRegression;

    regression_results();
    ~regression_results();

    void  clear();
    void  clear_full();
    void  clear_half();

    void  set_residual_variance(double rv);
    void  set_total_variance(double tv);
    void  set_sum_residual_variance(double rv);
    void  set_diff_residual_variance(double rv);

    void  set_full_residual_variance(double rv);
    void  set_full_total_variance(double tv);
    void  set_full_sum_residual_variance(double rv);
    void  set_full_diff_residual_variance(double rv);

    void  set_half_residual_variance(double rv);
    void  set_half_total_variance(double tv);
    void  set_half_sum_residual_variance(double rv);
    void  set_half_diff_residual_variance(double rv);

    void  set_full_w(const w_map& ws);
    void  set_half_w(const w_map& ws);

    void  set_residual_info(const SampleInfo& s);
    void  set_residual_square_info(const SampleInfo& s);
    void  set_residual_square_info_reduced(const SampleInfo& s);
    void  set_full_residual_info(const SampleInfo& s);
    void  set_half_residual_info(const SampleInfo& s);

    void  set_results(const result_vector& r);
    void  set_F_result(const F_result_type& f);

    void  add_result(const reg_result& r);
    void  add_full_result(const reg_result& r);
    void  add_half_result(const reg_result& r);

    double            get_residual_variance()            const;
    double            get_total_variance()               const;
    double            get_sum_residual_variance()        const;
    double            get_diff_residual_variance()       const;

    double            get_full_residual_variance()       const;
    double            get_full_total_variance()          const;
    double            get_full_sum_residual_variance()   const;
    double            get_full_diff_residual_variance()  const;

    double            get_half_residual_variance()       const;
    double            get_half_total_variance()          const;
    double            get_half_sum_residual_variance()   const;
    double            get_half_diff_residual_variance()  const;

    const w_map&      get_full_w() const;
    const w_map&      get_half_w() const;

    const SampleInfo& get_residual_info()                const;
    const SampleInfo& get_residual_square_info()         const;
    const SampleInfo& get_residual_square_info_reduced() const;
    const SampleInfo& get_full_residual_info()           const;
    const SampleInfo& get_half_residual_info()           const;

          reg_result& get_result(size_t i);
    const reg_result& get_result(size_t i)      const;

          reg_result& get_full_result(size_t i);
    const reg_result& get_full_result(size_t i) const;

          reg_result& get_half_result(size_t i);
    const reg_result& get_half_result(size_t i) const;

    size_t               get_result_count()     const;
    size_t               get_intercept_count()  const;

          result_vector& get_results();
    const result_vector& get_results()          const;

          F_result_type& get_F_result();
    const F_result_type& get_F_result()         const;

    void  dump(ostream &out) const;

  protected:

    correlation_matrix_map  my_correlation_matrix_map;

    w_map                   my_full_w;
    w_map                   my_half_w;

    SampleInfo              my_residual_info;
    SampleInfo              my_residual_square_info;
    SampleInfo              my_residual_square_info_reduced;
    SampleInfo              my_full_residual_info;
    SampleInfo              my_half_residual_info;

    reg_variances           my_variances;
    reg_variances           my_full_variances;
    reg_variances           my_half_variances;

    result_vector           my_results;
    result_vector           my_full_results;
    result_vector           my_half_results;

    size_t                  my_intercept_count;

    F_result_type           my_F_result;
};

#include "regress_result.ipp"

} //end of namespace SIBPAL
} //end of namespace SAGE

#endif
