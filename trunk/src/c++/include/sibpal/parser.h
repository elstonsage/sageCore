#ifndef SIBPAL_PARSER_H
#define SIBPAL_PARSER_H

//=============================================================================
// File:    sibpal_parser.h
//
// Author:  Kevin Jacobs
//
// History: Version 0.0 Initial implementation
//
// Notes:   This file contains definition for following data structures.
//
//            class sibpal_parser : public APP::BasicParser
//
//            class SibMeanTestParser     : public sibpal_parser
//            class TraitRegressionParser : public sibpal_parser
//
// Copyright (c) 2001   R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "app/aparser.h"
#include "sibpal/meantest_params.h"
#include "sibpal/regress_params.h"
#include "sibpal/gls3.h"

namespace SAGE   {
namespace SIBPAL {

class sibpal_parser : public APP::BasicParser
{
  public:

    sibpal_parser(const relative_pairs& pairs, cerrorstream& err = SAGE::sage_cerr);

    virtual void parse_symbols(const SymbolTable* syms) = 0;
    virtual void parse_parameter(const LSFBase* param)  = 0;
    virtual void parse_test_parameter_section(const LSFBase* params) = 0;
    virtual void parse_test_parameter(const LSFBase* param) = 0;

  protected:

    const relative_pairs&   my_pairs;
};

class meantest_parser : public sibpal_parser
{
  public:

    meantest_parser(const relative_pairs& pairs, cerrorstream& err = SAGE::sage_cerr);
    meantest_parser(const relative_pairs& pairs, SymbolTable* syms, cerrorstream& err = SAGE::sage_cerr);
    meantest_parser(const relative_pairs& pairs, LSFBase* params, cerrorstream& err = SAGE::sage_cerr);
    meantest_parser(const relative_pairs& pairs, SymbolTable* syms, const LSFBase* params, cerrorstream& err = SAGE::sage_cerr);

    void clear();
     
    void parse_symbols(const SymbolTable* syms);
    void parse_parameter(const LSFBase* param);
    void parse_test_parameter_section(const LSFBase* params);
    void parse_test_parameter(const LSFBase* param);

    const meantest_parameters& parameters()   const;

    bool wide_output()                    const;
    bool csv_output()                     const;

    void set_wide_output(bool b);
    void set_csv_output(bool b);

  protected:

    void parse_marker(const LSFBase* param);
    void parse_trait(const LSFBase* param);
    void parse_subset(const LSFBase* param);

  private:

    bool                my_wide_output;
    bool                my_csv_output;
    
    meantest_parameters my_parameters;
};

class regression_parser : public sibpal_parser
{
  public:

    regression_parser(const relative_pairs& pairs, cerrorstream& err = SAGE::sage_cerr);
    regression_parser(const relative_pairs& pairs, SymbolTable* syms, cerrorstream& err = SAGE::sage_cerr);
    regression_parser(const relative_pairs& pairs, LSFBase* params, cerrorstream& err = SAGE::sage_cerr);
    regression_parser(const relative_pairs& pairs, SymbolTable* syms, 
                      const LSFBase* params, cerrorstream& err = SAGE::sage_cerr);

    void clear();
     
    void parse_symbols(const SymbolTable* syms);
    void parse_parameter(const LSFBase* param);
    void parse_test_parameter_section(const LSFBase* params);
    void parse_test_parameter(const LSFBase* param);

    void set_regression_type(regression_type rt);
    void check_test_options();

    const vector<LSF_ptr<LSFBase> >&     get_pair_info_file()         const;

    const vector<trait_type>&            get_traits()                 const;
    const vector<trait_type>&            get_subsets()                const;
    const vector<covariate_type>&        get_covariates()             const;
    const vector<marker_type>&           get_markers()                const;
    const vector<independent_variable>&  get_interactions()           const;
    const vector<independent_variable>&  get_batch_interactions()     const;

    const regression_type&               get_regression_type()        const;

    const string&                        get_regression_method_name() const;

    const data_options&                  get_data_options()           const;
    const analysis_options&              get_analysis_options()       const;
    const pvalue_options&                get_pvalue_options()         const;
    const output_options&                get_output_options()         const;

    void dump_parser(ostream &out)                                    const;

  protected:

    void parse_trait(const LSFBase* param);
    void parse_subset(const LSFBase* param);
    void parse_marker(const LSFBase* param);
    void parse_covariate(const LSFBase* param);
    void parse_interaction(const LSFBase* params);
    void parse_marker_type(const LSFBase* param, marker_type& mar, bool inter_marker=false);
    void parse_covariate_type(const LSFBase* param, covariate_type& cov);
    void parse_regression_method(const LSFBase* param);
    void parse_use_pairs(const LSFBase* param);
    void parse_empirical_pvalues(const LSFBase* param);
    void parse_dump_data(const LSFBase* param);
    void parse_design_matrix(const LSFBase* param);
    void parse_correlation_matrix(const LSFBase* param);
    void parse_x_linkage(const LSFBase* param);
    
  private:

    vector<LSF_ptr<LSFBase> >     my_pair_info_file;

    vector<trait_type>            my_traits;
    vector<trait_type>            my_subsets;
    vector<covariate_type>        my_covariates;
    vector<marker_type>           my_markers;
    vector<independent_variable>  my_interactions;
    vector<independent_variable>  my_batch_interactions;

    regression_type               my_reg_type;

    string                        my_reg_method_name;

    data_options                  my_data_opt;
    analysis_options              my_analysis_opt;
    pvalue_options                my_pvalue_opt;
    output_options                my_output_opt;
};

#include "sibpal/parser.ipp"

} // end of namespace SIBPAL
} // end of namespace SAGE

#endif
