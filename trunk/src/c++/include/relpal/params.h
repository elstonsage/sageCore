#ifndef RELPAL_PARAMS_H
#define RELPAL_PARAMS_H

//=============================================================================
// File:    params.h
//
// Author:  Yeunjoo Song
//
// History: Initial implementation                                  yjs Mar. 07
//
// Notes:   This file contains definition for following data structures.
//            struct dependent_variable
//            struct covariate_type
//            struct marker_type
//            struct independent_variable
//
// Copyright (c) 2007   R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "relpal/definitions.h"

namespace SAGE   {
namespace RELPAL {

struct dependent_variable
{
  explicit dependent_variable(size_t t = (size_t)-1, bool bin = false);

  string name(const relative_pairs& pairs) const;

  bool operator==(const dependent_variable& t) const;
  bool operator!=(const dependent_variable& t) const;

  size_t         trait_index;

  bool           binary;
  bool           valid;
};

struct covariate_type
{
  explicit covariate_type(size_t c = (size_t)-1, size_t t = (size_t)-1, bool test = false);

  string name(const relative_pairs& pairs) const;
  string effect_name() const;
  string short_effect_name() const;

  bool operator==(const covariate_type& c) const;
  bool operator!=(const covariate_type& c) const;
  bool operator< (const covariate_type& c) const;

  size_t     covariate_index;
  size_t     adj_trait_index;

  SampleInfo info;

  bool       test_variable;
  bool       valid;
};

struct marker_type
{
  enum eff_type { TOTAL=0, DOMINANCE, ADDITIVE };

  explicit marker_type(size_t m = (size_t)-1, bool test = false);

  string name(const relative_pairs& pairs) const;
  string effect_name() const;
  string short_effect_name() const;

  bool operator==(const marker_type& m) const;
  bool operator!=(const marker_type& m) const;
  bool operator< (const marker_type& m) const;

  size_t     marker_index;

  eff_type   effect;

  bool       test_variable;
  bool       x_linked;
  bool       valid;
};

string covariate_name(const relative_pairs& pairs, const vector<covariate_type>& covs);
string marker_name   (const relative_pairs& pairs, const vector<marker_type>&    mars);
string covariate_effect_name(bool shortname,       const vector<covariate_type>& covs);
string marker_effect_name   (bool shortname,       const vector<marker_type>&    mars);

struct interaction_type
{
  enum inter_type { MM=0, MC, CC };

  interaction_type();

  string name(const relative_pairs& pairs) const;
  string effect_name() const;
  string short_effect_name() const;

  bool operator==(const interaction_type& m) const;
  bool operator!=(const interaction_type& m) const;
  bool operator< (const interaction_type& m) const;

  vector<marker_type>     markers;
  vector<covariate_type>  covariates;

  inter_type              type;

  bool                    batch;
  bool                    test_variable;
  bool                    valid;
};

struct independent_variable
{
  enum indep_var_type { INTERCEPT=0, RANDOM_ERR, COMMON_ENV, POLYGENIC,
                        MARKER, COVARIATE, MM_INTER, MC_INTER, CC_INTER };

  independent_variable();

  string    dump(const relative_pairs& pairs, bool with_t = true) const;

  string    name(const relative_pairs& pairs)             const;
  string    effect_name(bool with_t = true)               const;
  string    short_effect_name(bool with_t = true)         const;

  bool operator==(const independent_variable& c) const;
  bool operator!=(const independent_variable& c) const;
  bool operator< (const independent_variable& c) const;

  void clear();

  vector<marker_type>     markers;
  vector<covariate_type>  covariates;

  indep_var_type          type;

  size_t                  t1;
  size_t                  t2;

  bool                    test_variable;
  bool                    valid;
};

struct filtering_type
{
  explicit filtering_type(size_t t = (size_t)-1, string n = "");

  string name() const;

  bool operator==(const filtering_type& t) const;
  bool operator!=(const filtering_type& t) const;

  size_t subset_index;
  string subset_name;
};

typedef vector<dependent_variable>          trait_vector;
typedef trait_vector::iterator              trait_iterator;
typedef trait_vector::const_iterator        trait_const_iterator;

typedef vector<independent_variable>        parameter_vector;        
typedef parameter_vector::iterator          parameter_iterator;      
typedef parameter_vector::const_iterator    parameter_const_iterator;

typedef vector<filtering_type>              filtering_vector;
typedef filtering_vector::iterator          filtering_iterator;
typedef filtering_vector::const_iterator    filtering_const_iterator;

typedef pair<trait_iterator, bool>          tib_value;
typedef pair<parameter_iterator, bool>      pib_value;
typedef pair<filtering_iterator, bool>      fib_value;

struct data_options
{
  data_options();

  string dump(string pre_space) const;

  filtering_vector  subsets;
  use_member_type   use_members;
  use_pair_type     use_pairs;
};

struct output_options
{
  output_options();

  string dump(string pre_space) const;

  bool              detailed_out;
  bool              export_out;
  bool              data_out;
  bool              debug_out;
  bool              residual_out;

  size_t            param_name_max;
  size_t            param_eff_name_max;
};

struct analysis_options
{
  analysis_options();

  string dump(string pre_space) const;
  string dump_var(string pre_space) const;

  bool              transform_residuals;
  bool              naive_variance;
  bool              sandwich_variance;
  bool              alternative_variance;
  bool              IBD_variance;
  string            state_file_name;
};

struct asymptotic_pvalue_options
{
  asymptotic_pvalue_options();

  string dump(string pre_space) const;

  bool              valid;

  size_t            seed;

  size_t            replicates;
  size_t            min_replicates;
  size_t            max_replicates;

  double            threshold;
  double            width;
  double            confidence;
};

typedef asymptotic_pvalue_options pvalue_options;

#include "relpal/params.ipp"

} // end of namespace RELPAL
} // end of namespace SAGE

#endif
