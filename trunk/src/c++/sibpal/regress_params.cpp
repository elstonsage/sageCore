//=============================================================================
// File:    regress_params.cpp
//
// Author:  Kevin Jacobs 
//
// History: Version 0.0 Initial implementation
//          Took it out from sibregress.cpp                         yjs  Jun.05
//
// Notes:   This file contains implementation for following data structures.
//            struct trait_type
//            struct dependent_variable
//            struct covariate_type
//            struct marker_type
//            struct independent_variable
//            class regression_parameters
//
// Copyright (c) 2001 R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "sibpal/regress_params.h"

using namespace std;

namespace SAGE   {
namespace SIBPAL {

//
//------------------------------------------------------------------------
//

trait_type::trait_type(size_t t, bool bin, double fixed_m, double fixed_p)
{
  trait_index       = t;
  binary            = bin;
  fixed_mean        = fixed_male_mean = fixed_female_mean = fixed_m;
  fixed_correlation = fixed_p;

  trait_used_sibs_info.clear();
  trait_all_sibs_info.clear();
  trait_male_sibs_info.clear();
  trait_female_sibs_info.clear();

  var_b = var_r = correlation
                = fisher_correlation
                = fisher_fsib_correlation
                = fisher_hsib_correlation
                = sum_diff_correlation     = QNAN;

  valid = false;
}

string
trait_type::name(const relative_pairs& pairs) const
{
  string t = pairs.trait_name(trait_index);
  if(!t.size())
    t = "Trait #" + long2str(trait_index);
  return t;
}

//
//------------------------------------------------------------------------
//

dependent_variable::dependent_variable(size_t t, bool bin, double fixed_m, double fixed_p)
                  : trait_type(t, bin, fixed_m, fixed_p)
{
  p_correlation.resize(2, QNAN);
  p_fsib_correlation.resize(2, QNAN);
  p_hsib_correlation.resize(2, QNAN);
  
  p_empirical_correlation.resize(2, QNAN);
  p_fsib_empirical_correlation.resize(2, QNAN);
  p_hsib_empirical_correlation.resize(2, QNAN);
  p_fh_empirical_correlation.resize(2, QNAN);

  p_sum_diff_correlation.resize(3, QNAN);
  p_fsib_sum_diff_correlation.resize(3, QNAN);
  p_hsib_sum_diff_correlation.resize(3, QNAN);
  p_fh_sum_diff_correlation.resize(3, QNAN);

  info.clear();
  mm_pair_info.clear();
  mf_pair_info.clear();
  ff_pair_info.clear();
}

//
//------------------------------------------------------------------------
//

covariate_type::covariate_type()
              : covariate_index((size_t)-1), operation(none), power(1.0), fixed_mean(QNAN)
{}

covariate_type::covariate_type(size_t c, cov_op op, double powr, double f_mean)
              : covariate_index(c), operation(op), power(powr), fixed_mean(f_mean)
{}

string
covariate_type::name(const relative_pairs& pairs) const
{
  if( operation == pair )
  {
    string n = pairs.get_pair_covariate_name(covariate_index);

    if(!n.size())
      n = "Pair Covariate #" + long2str(covariate_index);

    return n;
  }

  string n = pairs.trait_name(covariate_index);
  if(!n.size())
    n = "Covariate #" + long2str(covariate_index);

  if(power != 1)
     n  = "(" + n + ")^" + doub2str(power);

  return n;
}

string
covariate_type::effect_name() const
{
  string n;

  switch( operation )
  {
    case sum:    n = "Cov,sum"; break;
    case diff:   n = "Cov,diff"; break;
    case prod:   n = "Cov,prod"; break;
    case avg:    n = "Cov,avg"; break;
    case single: n = "Cov,single"; break;
    case pair:   n = "Cov,pair"; break;
    case none:  throw std::exception();
    default: break;
  }
  return n;
}

string
covariate_type::short_effect_name() const
{
  string n;

  switch( operation )
  {
    case sum:    n = "C,sum"; break;
    case diff:   n = "C,diff"; break;
    case prod:   n = "C,prod"; break;
    case avg:    n = "C,avg"; break;
    case single: n = "C,single"; break;
    case pair:   n = "C,pair"; break;
    case none: throw std::exception();
    default: break;
  }
  return n;
}

//
//------------------------------------------------------------------------
//

marker_type::marker_type(size_t m, effect_type e, bool x)
           : marker_index(m), effect(e), x_linked(x)
{ }

string
marker_type::name(const relative_pairs& pairs) const
{
  string n = pairs.marker_name(marker_index);
  if(!n.size())
    n = "Marker #" + long2str(marker_index);
  return n;
}

string
marker_type::effect_name() const
{
  string n;

  switch(effect)
  {
    case TOTAL:     n = "(A+D)GenVar"; break;
    case ADDITIVE:  n = "AddGenVar";   break;
    case DOMINANCE: n = "DomGenVar";   break;
    default:        n = "???GenVar";   break;
  }
  return n;
}


string
marker_type::short_effect_name() const
{
  string n;

  switch(effect)
  {
    case TOTAL:     n = "A+D"; break;
    case ADDITIVE:  n = "Add"; break;
    case DOMINANCE: n = "Dom"; break;
    default:        n = "?";   break;
  }
  return n;
}

//
//------------------------------------------------------------------------
//

independent_variable::independent_variable()
                    : valid(false), type(INVALID)
{ }

void
independent_variable::clear()
{
  valid = false;
  markers.resize(0);
  covariates.resize(0);
  info.clear();
  mm_pair_info.clear();
  mf_pair_info.clear();
  ff_pair_info.clear();
}

string
independent_variable::name(const relative_pairs& pairs) const
{
  if( !markers.size() && !covariates.size() )
    return "";

  if(  markers.size() && !covariates.size() )
    return marker_name(pairs);

  if( !markers.size() &&  covariates.size() )
    return covariate_name(pairs);

  return marker_name(pairs) + " x " + covariate_name(pairs);
}

string
independent_variable::effect_name() const
{
  if( !markers.size() && !covariates.size() )
    return "";

  if(  markers.size() && !covariates.size() )
    return marker_effect_name();

  if( !markers.size() &&  covariates.size() )
    return covariate_effect_name();

  return marker_effect_name(true) + " x " + covariate_effect_name(true);
}

string
independent_variable::short_effect_name() const
{
  if( !markers.size() && !covariates.size() )
    return "";

  if(  markers.size() && !covariates.size() )
    return marker_effect_name(true);

  if( !markers.size() &&  covariates.size() )
    return covariate_effect_name(true);

  return marker_effect_name(true) + " x " + covariate_effect_name(true);
}

string
independent_variable::covariate_name(const relative_pairs& pairs) const
{
  string n;

  if( covariates.size() )
    n = covariates[0].name(pairs);

  for(size_t i = 1; i < covariates.size(); ++i)
  {
    n += " x " + covariates[i].name(pairs);
  }
  return n;
}

string
independent_variable::covariate_effect_name(bool shortname) const
{
  if( !covariates.size() )
    return "";

  if( !shortname && covariates.size() == 1)
    return covariates[0].effect_name();

  string n = covariates[0].short_effect_name();
  for(size_t i = 1; i < covariates.size(); ++i)
    n += " x " + covariates[i].short_effect_name();

  return n;
}

string
independent_variable::marker_name(const relative_pairs& pairs) const
{
  if( !markers.size() )
    return "";

  string n = markers[0].name(pairs);
  for(size_t i = 1; i < markers.size(); ++i)
  {
    n += " x " + markers[i].name(pairs);
  }

  return n;
}

string
independent_variable::marker_effect_name(bool shortname) const
{
  if( !markers.size() )
    return "";

  if( !shortname && markers.size() == 1)
    return markers[0].effect_name();

  string n = markers[0].short_effect_name();
  for(size_t i = 1; i < markers.size(); ++i)
    n += " x " + markers[i].short_effect_name();

  return n;
}

//
// ----------------------------------------------------------------------
//

data_options::data_options()
            : use_pairs(FSIB), skip_uninformative_pairs(false),
              use_mm_pair(true), use_mf_pair(true), use_ff_pair(true)
{}

string
data_options::dump(string pre_space) const
{
  string n = pre_space + "use_pair = ";

  if( use_pairs == FSIB )
    n += "full sib pairs only\n";
  else if( use_pairs == HSIB )
    n += "half sib pairs only\n";
  else
    n += "both full and half sib pairs\n";

  n += pre_space + "skip_uninformative_pairs = " + get_bool_to_string(skip_uninformative_pairs) + "\n";
  n += pre_space + "use_mm_pairs = " + get_bool_to_string(use_mm_pair) + "\n";
  n += pre_space + "use_mf_pairs = " + get_bool_to_string(use_mf_pair) + "\n";
  n += pre_space + "use_ff_pairs = " + get_bool_to_string(use_ff_pair) + "\n";

  return n;
}

//
// ----------------------------------------------------------------------
//

analysis_options::analysis_options()
                : standardize_parameter(true),
                  robust_variance(false), leverage_adjustment(false),
                  identity_weight(false), pool_correlation(false),
                  blup_mean(true), sibship_mean(false),
                  sibship_mean_threshold(3), max_iterations(10), w1(0.5)
{}

string
analysis_options::dump(string pre_space) const
{
  string n = pre_space + "standardize_parameters = " + get_bool_to_string(standardize_parameter) + "\n";

  n += pre_space + "robust_variance = " + get_bool_to_string(robust_variance) + "\n";
  n += pre_space + "leverage_adjustment = " + get_bool_to_string(leverage_adjustment) + "\n";
  n += pre_space + "identity_weights = " + get_bool_to_string(identity_weight) + "\n";
  n += pre_space + "pool_correlations = " + get_bool_to_string(pool_correlation) + "\n";
  n += pre_space + "blup_mean = " + get_bool_to_string(blup_mean) + "\n";
  n += pre_space + "sibship_mean = " + get_bool_to_string(sibship_mean) + "\n";
  n += pre_space + "w1 = " + doub2str(w1) + "\n";

  return n;
}

//
// ----------------------------------------------------------------------
//

empirical_pvalue_options::empirical_pvalue_options()
                        : is_on(false), seed(0), width(), confidence(0.95), threshold(0.05),
                          replicates(0), min_replicates(50), max_replicates(10000)
{}

string
empirical_pvalue_options::dump(string pre_space) const
{
  string n = pre_space + "seed = " + long2str(seed) + "\n";

  n += pre_space + "replicates = " + long2str(replicates) + "\n";
  n += pre_space + "min_replicates = " + long2str(min_replicates) + "\n";
  n += pre_space + "max_replicates = " + long2str(max_replicates) + "\n";

  n += pre_space + "threshold = " + doub2str(threshold) + "\n";
  n += pre_space + "width = " + doub2str(width) + "\n";
  n += pre_space + "confidence = " + doub2str(confidence) + "\n";

  return n;
}

//
// ----------------------------------------------------------------------
//

output_options::output_options()
              : wide_out(false), detailed_out(true), export_out(false), print_design_matrix(false),
                print_correl_matrix(false), pvalues_scientific_notation(false),
                design_matrix_rows(10), design_matrix_file(""), correl_matrix_file(""),
                param_name_max_size(11), param_effect_name_max_size(11),
                dump_data(false), dump_debug(false), dump_data_file("")
{}

string
output_options::dump(string pre_space) const
{
  string n = pre_space + "wide_out = " + get_bool_to_string(wide_out) + "\n";

  n += pre_space + "detailed_out = " + get_bool_to_string(detailed_out) + "\n";
  n += pre_space + "export_out = " + get_bool_to_string(export_out) + "\n";
  n += pre_space + "print_design_matrix = " + get_bool_to_string(print_design_matrix) + "\n";
  n += pre_space + "print_correlation_matrix = " + get_bool_to_string(print_correl_matrix) + "\n";
  n += pre_space + "pvalues_scientific_notation = " + get_bool_to_string(pvalues_scientific_notation) + "\n";
  n += pre_space + "dump_data = " + get_bool_to_string(dump_data) + "\n";
  n += pre_space + "dump_debug = " + get_bool_to_string(dump_debug) + "\n";

  return n;
}

//
// ----------------------------------------------------------------------
//

regression_model::regression_model()
{
  clear();
}

void
regression_model::clear()
{
  clear_subsets();
  clear_parameters();

  invalidate();
  my_x_linked = false;
}

void
regression_model::invalidate()
{
  my_valid = false;
  my_analysis_opt.robust_variance = false;
  my_pvalue_opt.is_on = false;
}

void
regression_model::validate()
{
  my_valid = true;
}
/*
double
regression_model::trait_mean(const trait_parameter& t) const
{
  double mean = t.fixed_mean;
  if(!finite(mean))
    mean = t.trait_all_sibs_info.mean();
  return mean;
}

double
regression_model::trait_correlation(const trait_parameter& t) const
{
  double rho = t.fixed_correlation;
  if(!finite(rho))
    rho = t.correlation;
  if(!finite(rho))
    rho = t.fisher_correlation;
  if(!finite(rho))
    rho = my_trait_correlation;
  return rho;
}
*/
pair<regression_model::trait_iterator, bool>
regression_model::add_subset(const trait_type& t)
{
  trait_iterator i = std::find(my_subsets.begin(),
                               my_subsets.end(),   t);

  bool inserted = true;
  if( i == my_subsets.end() )
  {
    invalidate();
    my_subsets.push_back( t );
    i = my_subsets.end();
    --i;
  }
  else
  {
    inserted = false;
  }
  return std::pair<trait_iterator, bool>(i,inserted);
}

void
regression_model::clear_markers()
{
  parameter_vector  new_parameters;

  parameter_iterator i = my_parameters.begin();
  for( ; i != my_parameters.end(); ++i )
  {
    if( i->type != independent_variable::MARKER )
      new_parameters.push_back(*i);      
  }

  my_parameters.clear();
  
  i = new_parameters.begin();
  for( ; i != new_parameters.end(); ++i )
    my_parameters.push_back(*i);

  my_x_linked = false;

  return;
}

void
regression_model::clear_invalid_parameters()
{
  bool x_linked = false;

  parameter_vector  new_parameters;

  parameter_iterator i = my_parameters.begin();
  for( ; i != my_parameters.end(); ++i )
  {
    if( i->valid )
    {
      new_parameters.push_back(*i);

      if( i->type == independent_variable::MARKER && i->markers[0].x_linked )
        x_linked = true;
    }
  }

  my_parameters.clear();
  
  i = new_parameters.begin();
  for( ; i != new_parameters.end(); ++i )
    my_parameters.push_back(*i);

  my_x_linked = x_linked;

  return;
}

regression_model::pib_value
regression_model::add_parameter(const independent_variable& p)
{
  // Check to see if parameter is already exists
  parameter_iterator i = std::find(my_parameters.begin(),
                                   my_parameters.end(),   p);

  bool inserted = true;
  if( i == my_parameters.end() )
  {
    invalidate();
    my_parameters.push_back( p );
    i = my_parameters.end();
    --i;

    if( p.type == independent_variable::MARKER && p.markers[0].x_linked )
      my_x_linked = true;
  }
  else
  {
    inserted = false;
  }
  return std::pair<parameter_iterator, bool>(i,inserted);
}

void
regression_model::dump_model(const relative_pairs& relpairs, ostream &out) const
{
  out << endl
      << "==========================" << endl
      << "  Regression model dump  " << endl
      << "==========================" << endl << endl;
                                  
  if( valid() )
    out << "Model valid!!";
  else
    out << "Model invalid!";
  out << endl << endl;

  if( my_x_linked )
    out << "X_linked" << endl << endl;    

  out << "regression type:   " << get_regression_type_to_string(my_reg_type)
      << endl;
  out << "regression method: " << my_reg_method_name
      << endl << endl;

  out << "Dependent variables: " << my_trait.name(relpairs)
      << endl << endl;

  out << "Independent variables:" << endl;
  for( size_t i = 0; i < my_parameters.size(); ++i )
    out << " " << i+1 << " : " << my_parameters[i].name(relpairs)
        << "(" << my_parameters[i].short_effect_name() << ")"
        << endl;
  out << endl;

  out << "Data options:";
  if( my_subsets.size() )
  {
    out << endl;
    out << " subset:" << endl;
    for( size_t i = 0; i < my_subsets.size(); ++i )
      out << " " << my_subsets[i].name(relpairs) << endl;
  }

  out << endl << my_data_opt.dump(" ");

  out << endl << "Analysis options:" << endl;
  out << my_analysis_opt.dump(" ");

  out << endl << "P-value options:" << endl;
  out << my_pvalue_opt.dump(" ");

  //out << endl << "Output options:" << endl;
  //out << my_output_opt.dump(" ");

  out << endl;

  return;
}

} // end of namespace SIBPAL
} // end of namespace SAGE
