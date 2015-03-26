//============================================================================
// File:      fpmm_sub_model.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   8/22/01  - created.                        djb
//                                                                          
// Notes:     inlines for finite_polygenic_mixed_model_sub_model class.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#ifndef SEGREG_FPMM_SUB_MODEL_H
#include "segreg/fpmm_sub_model.h"
#endif

namespace SAGE
{

namespace SEGREG
{

//============================================================================
// IMPLEMENTATION:  finite_polygenic_mixed_model_sub_model
//============================================================================
//
inline  
finite_polygenic_mixed_model_sub_model::finite_polygenic_mixed_model_sub_model
      (cerrorstream& errors)
    : SegregSubmodel(errors)
{
  //lint -e{534}
  set(model_input(FPMM_DEFAULT_VALUE, FPMM_DEFAULT_FIXED), FPMM_DEFAULT_FREQ,
                  FPMM_DEFAULT_LOCI);
}

inline
finite_polygenic_mixed_model_sub_model::finite_polygenic_mixed_model_sub_model
      (const finite_polygenic_mixed_model_sub_model& other)
   : SegregSubmodel(other)
{
  my_variance        = other.my_variance;
  my_variance_fixed  = other.my_variance_fixed;
  my_frequency       = other.my_frequency;
  my_max_pgt         = other.my_max_pgt;
  my_means           = other.my_means;
}

inline finite_polygenic_mixed_model_sub_model&
finite_polygenic_mixed_model_sub_model::operator=
        (const finite_polygenic_mixed_model_sub_model& other)
{
  if(this != &other)
  {
    SegregSubmodel::operator=(other);

    my_variance        = other.my_variance;
    my_variance_fixed  = other.my_variance_fixed;
    my_frequency       = other.my_frequency;
    my_max_pgt         = other.my_max_pgt;
    my_means           = other.my_means;
  }
  
  return *this;
}

inline
finite_polygenic_mixed_model_sub_model::~finite_polygenic_mixed_model_sub_model()
{}


inline string
finite_polygenic_mixed_model_sub_model::option_description() const
{
  return //lint -e(713)
         long2str(loci()) + " loci, " +
         "frequency = " + doub2str(my_frequency);
}

inline string  
finite_polygenic_mixed_model_sub_model::name() const
{
  return FPMM_NAME; 
}

inline double
finite_polygenic_mixed_model_sub_model::variance() const
{
  return my_variance;
}

inline double
finite_polygenic_mixed_model_sub_model::mean(size_t p) const
{
  return my_means[p];
}

inline double
finite_polygenic_mixed_model_sub_model::frequency() const
{
  return my_frequency;
}

inline size_t  
finite_polygenic_mixed_model_sub_model::max_pgt() const
{
  return my_max_pgt;
}

inline size_t  
finite_polygenic_mixed_model_sub_model::loci() const
{
  return (my_max_pgt - 1) / 2;
}

inline bool
finite_polygenic_mixed_model_sub_model::is_complete() const
{
  return finite(my_variance) && my_variance > FPMM_LB + FPMM_EPSILON;
}

inline bool // due to JA
finite_polygenic_mixed_model_sub_model::variance_fixed() 
{
  return my_variance_fixed;
}

// - Write sub-model values in LSF readable format.
//
inline void
finite_polygenic_mixed_model_sub_model::dump(std::ostream& out) const
{
  int  old_precision = out.precision();
  out.precision(DUMP_PRECISION);

  assert(! SAGE::isnan(my_variance));

  out << "# " << name() << "\n"
      << "fpmm\n" 
      << "{\n"
      << "  loci=" << loci() << "\n"
      << "  freq=" << my_frequency << "\n"
      << "  var, val=" << my_variance << ", fixed=" << std::boolalpha << my_variance_fixed << "\n"
      << "}" << std::noboolalpha << std::endl;
      
  out.precision(old_precision);
}

inline double          
finite_polygenic_mixed_model_sub_model::pop_freq(size_t s) const
{
  double p     = frequency();
  size_t l     = loci() * 2;

  return NUMERICS::bin_prob(l,p,s);
}

inline double          
finite_polygenic_mixed_model_sub_model::pop_freq(size_t i2, size_t s2, double corr) const
{
  double p          = frequency();
  size_t l          = 2 * loci();

  double sfrac      = (double) s2 / (double) l;

  double p_bar = p + corr * (sfrac - p);

  return NUMERICS::bin_prob(l, p_bar, i2);
}


}}
