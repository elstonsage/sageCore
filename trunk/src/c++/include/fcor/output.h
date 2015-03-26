#ifndef FCOROUTPUT_H
#define FCOROUTPUT_H

//****************************************************************************
//* File:      output.h                                                      *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0                                                   *
//*                                                                          *
//* Notes:     This header file defines the fuctions for fcor results output.*
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "fcor/htest.h"

namespace SAGE {
namespace FCOR {


class FcorView
{
  public:

    FcorView(const PairSetData& pairs);

    void  print_analysis_results(ostream& out, output_type ot,
                                 const pairset_result_vector& sub_results,
                                 const pairset_result_vector& main_results,
                                 const Htest_result_vector&   H_results)    const;

    void  print_var_cov_results (ostream& out,
                                 const var_cov_result_vector& vc_results)   const;

  private:

    void  print_analysis_header(ostream& o, pair_type ptype, output_type otype) const;

    void  print_main_type_title(ostream& o, main_to_sub_type_const_iterator rel_group)   const;
    void  print_pooled_types   (ostream& o, main_to_sub_type_const_iterator rel_group,
                                const pairset_info_vector& sinfos, size_t pcount)        const;
    void  print_table_title    (ostream& o, const pairset_info& pinfo, bool xls = false) const;

    void  print_cor_std_err_result(ostream& o, const pairset_info& pinfo,
                                   const pairset_result& re, output_type otype)             const;

    void  print_cor_std_err_heading(ostream& o, bool self, bool intra)             const;
    void  print_cor_std_err(ostream& o, size_t w, const pairset_result& cor_std, bool self, bool intra)  const;
    void  print_pooled_corr(ostream& o, const pairset_result& cor_std)  const;

    void  print_homogeneity_test_result(ostream& o, size_t subtype_count,
                                        const pairset_info& pinfo, const Htest_result& Hre) const;

    void  print_pair_count (ostream& o, const pairset_result& cor_std, bool self, bool intra)  const;
    void  print_xls_format (ostream& o, const pairset_result& cor_std, bool self, bool intra)  const;

    void  print_ddash_outline        (ostream& o) 			           const;
    void  print_sdash_outline        (ostream& o) 			           const;
    void  print_row_outline          (ostream& o, size_t w) 			   const;
    void  print_column_outline       (ostream& o, size_t w) 		           const;
    void  print_partial_outline_dash (ostream& o, size_t t) 		           const;
    void  print_column_trait_names   (ostream& o, size_t w)			   const;
    void  print_trait_name           (ostream& o, size_t t, size_t w)    	   const;
    void  print_class_type           (ostream& o, bool self, bool intra)           const;
    void  print_xls_row_outline      (ostream& o, size_t w, size_t w1, size_t w2)  const;

    void  print_correlation   (ostream& o,  size_t cnt, size_t w1,
                               string   sp, double cor, size_t w2, size_t f2)         const;
    void  print_standard_error(ostream& o,  double eff, size_t w1, size_t f1,
                               string   sp, double se,  size_t w2, size_t f2, bool r) const;
    void  print_pvalue        (ostream& o, double z, double var, size_t w, size_t f)  const;

    void  print_joint_vc_matrix(ostream& o, const var_cov_matrix& a_vc_result,
                                size_t t1, size_t t2, bool lp)                    const;
    void  print_single_vc_matrix(ostream& o, const var_cov_matrix& a_vc_result,
                                 size_t t, bool lp)                               const;

    void  print_legend(ostream& o, string pos, string cor, string trait1, string trait2) const;
    void  print_joint_legend(ostream& o, size_t c, size_t i, size_t j) const;
    void  print_single_legend(ostream& o, size_t c, size_t i)          const;

    const analysis_option_type& get_analysis_options()           const;

    size_t                      trait_count()                    const;

    string                      get_trait_name(size_t t)         const;
    
    double                      eff_count(double cor, double se) const;

    const PairSetData*  my_pairsetdata;
        
};

#include "fcor/output.ipp"    

} // end of namespace FCOR
} // end of namespace SAGE

#endif
