//============================================================================
// File:      parser.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   5/15/01 created        -djb
//                                                                          
// Notes:     Inline implementation of struct, model, and class, parser.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#ifndef SEGREG_PARSER_H
#include "segreg/parser.h"
#endif

namespace SAGE
{

namespace SEGREG
{

inline std::ostream&
operator<<(std::ostream& out, const parser::base_ptrs& ptrs)
{
  out << "title             " << ptrs.title             << "\n"
      << "trait             " << ptrs.trait             << "\n"
      << "composite_trait   " << ptrs.composite_trait   << "\n"
      << "type_mean         " << ptrs.type_mean         << "\n"
      << "type_var          " << ptrs.type_var          << "\n"
      << "type_suscept      " << ptrs.type_suscept      << "\n"
      << "mean_cov          " << ptrs.mean_cov          << "\n"
      << "var_cov           " << ptrs.var_cov           << "\n"
      << "suscept_cov       " << ptrs.suscept_cov       << "\n"
      << "class             " << ptrs.class_            << "\n"
      << "fpmm              " << ptrs.fpmm              << "\n"
      << "resid             " << ptrs.resid             << "\n"
      << "transformation    " << ptrs.transformation    << "\n"
      << "geno_freq         " << ptrs.geno_freq         << "\n"
      << "transmission      " << ptrs.transmission      << "\n"
      << "ascertainment     " << ptrs.ascertainment     << "\n"
      << "prev_constraint   " << ptrs.prev_constraint   << "\n"
      << "prev_estimate     " << ptrs.prev_estimate     << "\n"
      << "output_options    " << ptrs.output_options    << "\n"
      << "output            " << ptrs.output            << endl;
      
  return out;
}

//============================================================================
// IMPLEMENTATION:  parser::base_ptrs
//============================================================================
//
inline
parser::base_ptrs::base_ptrs()
    : title(0),
      trait(0),
      composite_trait(0),
      type_mean(0),
      type_var(0),
      type_suscept(0),
      mean_cov(0),
      var_cov(0),
      suscept_cov(0),
      class_(0),          
      fpmm(0),
      resid(0),
      transformation(0),
      geno_freq(0),
      transmission(0),
      ascertainment(0),
      prev_constraint(0),
      prev_constraints(0),
      prev_estimate(0),
      output_options(0),
      output(0),
      maxfun_options(0)
{}


//============================================================================
// IMPLEMENTATION:  parser
//============================================================================
//
inline const model&
parser::get_model() const
{
  return my_model;
}

inline long
parser::analysis_number() const
{
  return idnum;
}

inline void
parser::set_parameter
   (const string& name, const LSFBase** ptr, const LSFBase* param)
{
  if(!*ptr)
  {
    *ptr = param;
  }
  else
  {
    parameter_duplicated(name);
  }
}

inline bool
parser::model_valid() const
{
  return my_model.m_class != model_INVALID;
}

inline bool 
parser::process_block( const LSFBase* ptr) const
{
  return model_valid() && ptr != NULL;
}


}}
