//============================================================================
// File:      freq_sub_model.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   7/2/01  - created.                        djb
//                                                                          
// Notes:     inlines for genotype_frequency_sub_model class.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#ifndef SEGREG_FREQ_SUB_MODEL_H
#include "segreg/freq_sub_model.h"
#endif

namespace SAGE
{

namespace SEGREG
{

//============================================================================
// IMPLEMENTATION:  genotype_frequency_sub_model
//============================================================================
//
inline  
genotype_frequency_sub_model::genotype_frequency_sub_model
      (cerrorstream& errors)
    : SegregSubmodel(errors)
{
  model_input  correlation(CORR_DEFAULT_VALUE, CORR_DEFAULT_FIXED);

  // **** per gcw 10-8-1 FIX TO 0 BECAUSE ASSORTATIVE MATING IS BROKEN.
  //lint --e{506} --e{774}
  if(! TYPE_CORR_AVAILABLE)
  {
    correlation.value = 0;
    correlation.fixed = true;
  }

  //lint --e{534}
  set(hwe, FREQ_A_DEFAULT_VALUE, QNAN, QNAN, QNAN, correlation, PROBS_FIXED_DEFAULT);

  my_default = true;
}

inline
genotype_frequency_sub_model::genotype_frequency_sub_model
      (const genotype_frequency_sub_model& other)
   : SegregSubmodel(other)
{
  my_option = other.my_option;
  for(int i = 0; i < NUM_OF_PROBS; ++i)
  {
    my_probs[i] = other.my_probs[i];
  }
  
  my_freq_A  = other.my_freq_A;
  my_default = other.my_default;
}

inline genotype_frequency_sub_model&
genotype_frequency_sub_model::operator=
        (const genotype_frequency_sub_model& other)
{
  if(this != &other)
  {
    SegregSubmodel::operator=(other);

    for(int i = 0; i < NUM_OF_PROBS; ++i)
    {
      my_probs[i] = other.my_probs[i];
    }

    my_option  = other.my_option;    
    my_freq_A  = other.my_freq_A;
    my_default = other.my_default;
  }
  
  return *this;
}

inline
genotype_frequency_sub_model::~genotype_frequency_sub_model()
{}

inline genotype_frequency_sub_model::sm_option
genotype_frequency_sub_model::option() const
{
  return my_option;
}

inline string  
genotype_frequency_sub_model::option_description() const
{
  return option_2_description(my_option);
}

inline string  
genotype_frequency_sub_model::name() const
{
  return GENO_FREQ_NAME;
}

inline string
genotype_frequency_sub_model::option_2_description(sm_option opt)
{
  switch(opt)
  {
    case hwe:  return "Hardy-Weinberg equilibrium";
    case nhwe: return "non Hardy-Weinberg equilibrium";
    case NONE: return "";

    default:   SAGE_internal_error();
  }

  return "";
}

inline string
genotype_frequency_sub_model::option_2_parameter(sm_option opt)
{
  switch(opt)
  {
    case hwe:  return "hwe";
    case nhwe: return "nhwe";
    case NONE: return "none";

    default:   SAGE_internal_error();
  }

  return "";
}

inline double
genotype_frequency_sub_model::freq_A() const
{
  return my_freq_A;
}

inline double
genotype_frequency_sub_model::prob(genotype_index i) const
{
  return my_probs[i];
}

inline double          
genotype_frequency_sub_model::prob(genotype_index i, genotype_index s) const
{
  return my_probs[i];
}

inline bool
genotype_frequency_sub_model::default_() const
{
  return my_default;
}

// - Write sub-model values in LSF readable format.
//
inline void
genotype_frequency_sub_model::dump(std::ostream& out) const
{
  // No printing empty options!
  if(my_option == NONE) return;

  //lint --e{666,506,40}

  int  old_precision = out.precision();
  out.precision(DUMP_PRECISION);

  bool  probs_fixed = my_parameters[index_freq_A].initial_type == MAXFUN::Parameter::FIXED;
  
  out << "# " << name() << "\n"
      << "geno_freq\n" 
      << "{\n"
      << "  # " << option_2_description(my_option) << "\n"
      << "  option=" << option_2_parameter(my_option) << "\n"
      << "  probs_fixed=" << std::boolalpha << probs_fixed << "\n"; 
      
      if(my_option == hwe)
      {
        assert(! SAGE::isnan(my_freq_A));
        
        out << "  allele_freq_A, val=" << my_freq_A << "\n}\n";
      } 
      else if(my_option == nhwe)
      {
        assert(! SAGE::isnan(my_probs[index_AA]));
        assert(! SAGE::isnan(my_probs[index_AB]));
        assert(! SAGE::isnan(my_probs[index_BB]));
      
        out << "  prob=AA, val=" << my_probs[index_AA] << "\n"
            << "  prob=AB, val=" << my_probs[index_AB] << "\n"
            << "  prob=BB, val=" << my_probs[index_BB] << "\n}\n";
      }
      
  out.precision(old_precision);
}


inline bool genotype_frequency_sub_model::is_complete() const
{
  bool target_freq_complete = false;

  switch(option())
  {
    case NONE : target_freq_complete = true;
                break;

    case hwe  : //lint -e{734} <- really only a bool function anyway
                target_freq_complete = finite(freq_A());
                break;

    case nhwe : target_freq_complete = finite(prob(index_AA)) &&
                                       finite(prob(index_AB)) &&
                                       finite(prob(index_BB));
                break;

    default   : target_freq_complete = false;

  }                                 

  return target_freq_complete;
}

}}
