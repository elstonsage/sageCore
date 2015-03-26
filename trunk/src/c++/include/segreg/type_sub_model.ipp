//============================================================================
// File:      type_sub_model.ipp
//                                                                          
// Author:    Geoff Wedig (wedig@darwin.cwru.edu)
//                                                                          
// History:   0.1 gcw Initial Implementation                         May 2001
//                Baechle  reformatted, continued development        Jun 2001
//                djb renamed.  
//                    added genotype_specific_variance_sub_model.    7-30-01
//                djb added genotype_specific_susceptibility model.  1-16-02
//                                                          
// Notes:     inlines for genotype_specific_sub_models.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#ifndef SEGREG_TYPE_SUB_MODEL_H
#include "type_sub_model.h"
#endif

namespace SAGE
{

namespace SEGREG {


//============================================================================
// IMPLEMENTATION:  genotype_specific_sub_model
//============================================================================
//
inline
genotype_specific_sub_model::genotype_specific_sub_model
    (derived_type type_word, cerrorstream& errors)
  : SegregSubmodel(errors),
    my_derived_type(type_word),
    my_option(one),
    my_type_count(1),
    my_default(true),
    my_two_is_dom(true)
{
  my_types[0] = QNAN;
  my_types[1] = QNAN;
  my_types[2] = QNAN;

  my_types_fixed[0] = false;
  my_types_fixed[1] = false;
  my_types_fixed[2] = false;
}

inline
genotype_specific_sub_model::genotype_specific_sub_model
      (const genotype_specific_sub_model& other)
    : SegregSubmodel(other),
      my_derived_type(other.my_derived_type),
      my_option      (other.my_option),
      my_type_count  (other.my_type_count),
      my_default     (other.my_default),
      my_two_is_dom  (other.my_two_is_dom)
{
  for(int i = 0; i < NUM_OF_TYPES; ++i)
  {
    my_types      [i] = other.my_types      [i];
    my_types_fixed[i] = other.my_types_fixed[i];
  }
}

inline genotype_specific_sub_model&
genotype_specific_sub_model::operator=
        (const genotype_specific_sub_model& other)
{
  //lint -esym(1539, genotype_specific_sub_model::my_derived_type)

  if(this != &other)
  {
    MAXFUN::Submodel::operator=(other);
    
    my_option     = other.my_option;
    my_type_count = other.my_type_count;
    my_default    = other.my_default;
    my_two_is_dom = other.my_two_is_dom;
    
    for(int i = 0; i < NUM_OF_TYPES; ++i)
    {
      my_types      [i] = other.my_types      [i];
      my_types_fixed[i] = other.my_types_fixed[i];
    }
  }
  
  return *this;
}

inline size_t 
genotype_specific_sub_model::get_type_count() const
{
  return my_type_count;
}

inline bool
genotype_specific_sub_model::is_two_dom() const
{
  return my_two_is_dom || my_option == two_dom;
}

inline bool
genotype_specific_sub_model::is_two_option() const
{
  return my_option == two     ||
         my_option == two_dom ||
         my_option == two_rec;
}

inline bool
genotype_specific_sub_model::is_three_option() const
{
  return my_option == three     ||
         my_option == three_add ||
         my_option == three_dec ||
         my_option == three_inc;
}

inline string  
genotype_specific_sub_model::option_description() const
{
  return static_option_description(my_option, my_derived_type);
}

inline string  
genotype_specific_sub_model::name() const
{
  switch(my_derived_type)
  {
    case mean     : return TYPE_MEAN_NAME;
    case variance : return TYPE_VAR_NAME;
    case susc     : return TYPE_SUSCEPT_NAME;
    case invalid  :
    default       : break;
  }

  SAGE_internal_error();

  //lint --e{527} <- Statement unreachable, but we don't care.
  return "";
}

inline string
genotype_specific_sub_model::static_option_description
  (sm_option opt, derived_type t)
{
  string type;

  switch(opt)
  {
    case one:       type = "one "   + long_sing_type(t);                     break;
    case two:       type = "two "   + long_plural_type(t);                   break;
    case three:     type = "three " + long_plural_type(t);                   break;
    case two_dom:   type = "two "   + long_plural_type(t) + ", A dominant";  break;
    case two_rec:   type = "two "   + long_plural_type(t) + ", A recessive"; break;
    case three_add: type = "three " + long_plural_type(t) + ", additive";    break;
    case three_dec: type = "three " + long_plural_type(t) + ", decreasing";  break;
    case three_inc: type = "three " + long_plural_type(t) + ", increasing";  break;

    default:                                                                 break;
  }

  if(!type.size())
  {
    SAGE_internal_error();
  }

  return type;
}

inline string
genotype_specific_sub_model::option_2_parameter(sm_option opt)
{
  switch(opt)
  {
    case one           : return "one";
    case two           : return "two";
    case three         : return "three";
    case two_dom       : return "two_dom";
    case two_rec       : return "two_rec";
    case three_add     : return "three_add";
    case three_dec     : return "three_dec";
    case three_inc     : return "three_inc";
    default            : return "";
  }
}

inline genotype_specific_sub_model::sm_option
genotype_specific_sub_model::option() const
{
  return my_option;
}

inline double
genotype_specific_sub_model::parameter(genotype_index gt) const
{
  return my_types[gt];
}

inline string
genotype_specific_sub_model::long_sing_type(derived_type t)
{
  switch(t)
  {
    case mean     : return "mean";
    case variance : return "variance";
    case susc     : return "susceptibility";
    case invalid  :
    default       : return "";
  }
}

inline string
genotype_specific_sub_model::long_sing_type() const
{
  return long_sing_type(my_derived_type);
}

inline string
genotype_specific_sub_model::long_plural_type(derived_type t)
{
  switch(t)
  {
    case mean     : return "means";
    case variance : return "variances";
    case susc     : return "susceptibilities";
    case invalid  :
    default       : return "";
  }
}

inline bool genotype_specific_sub_model::is_complete() const
{
  return finite(parameter(index_AA)) &&
         finite(parameter(index_AB)) &&
         finite(parameter(index_BB));
}

//============================================================================
// IMPLEMENTATION:  genotype_specific_mean_susc_sub_model
//============================================================================
//
inline  
genotype_specific_mean_susc_sub_model::genotype_specific_mean_susc_sub_model
      (derived_type t, cerrorstream& errors)
    : genotype_specific_sub_model(t, errors) 
{
  //lint -e{534}
  set(one,
      model_input(MEAN_DEFAULT_VALUE, MEAN_DEFAULT_FIXED),
      model_input(QNAN, false),
      model_input(QNAN, false));
  
  my_default = true;
}

inline
genotype_specific_mean_susc_sub_model::genotype_specific_mean_susc_sub_model
      (const genotype_specific_mean_susc_sub_model& other)
    : genotype_specific_sub_model(other)
{ }

inline genotype_specific_mean_susc_sub_model&
genotype_specific_mean_susc_sub_model::operator=
        (const genotype_specific_mean_susc_sub_model& other)
{
  if(this != &other)
  {
    genotype_specific_sub_model::operator=(other);
  }
  
  return *this;
}


inline
genotype_specific_mean_susc_sub_model::~genotype_specific_mean_susc_sub_model()
{}

// - Determine 'fixity' of each type.
//
inline void
genotype_specific_mean_susc_sub_model::type_fixity(bool types[]) const
{
  types[index_AA] = my_types_fixed[index_AA];
  types[index_AB] = my_types_fixed[index_AB];
  types[index_BB] = my_types_fixed[index_BB];
}

// - Write sub-model values in LSF readable format.
//
inline void
genotype_specific_mean_susc_sub_model::dump(std::ostream& out) const
{
  int  old_precision = out.precision();
  out.precision(DUMP_PRECISION);

  bool  types[NUM_OF_TYPES];
  type_fixity(types); 

  bool b = !SAGE::isnan(my_types[0]);

  //lint --e{666,506,40}

  assert(b);
  assert(! SAGE::isnan(my_types[1]));
  assert(! SAGE::isnan(my_types[2]));

  string m = short_sing_type();
  if(m == "susc") m = "suscept";

  out << "# " << name() << endl
      << type_param() << endl 
      << "{" << endl
      << "  # " << static_option_description(my_option, my_derived_type) << endl
      << "  option=" << option_2_parameter(my_option) << endl
      << std::boolalpha
      << "  " << m << "=AA, val=" << my_types[0] << ", fixed=" << types[index_AA] << endl 
      << "  " << m << "=AB, val=" << my_types[1] << ", fixed=" << types[index_AB] << endl 
      << "  " << m << "=BB, val=" << my_types[2] << ", fixed=" << types[index_BB] << endl 
      << std::noboolalpha
      << "}" << endl;
      
  out.precision(old_precision);
}

//============================================================================
// IMPLEMENTATION:  genotype_specific_susceptibility_sub_model
//============================================================================
//
inline  
genotype_specific_susceptibility_sub_model::genotype_specific_susceptibility_sub_model
      (cerrorstream& errors)
    : genotype_specific_mean_susc_sub_model(susc, errors) 
{ }

inline
genotype_specific_susceptibility_sub_model::genotype_specific_susceptibility_sub_model
      (const genotype_specific_susceptibility_sub_model& other)
    : genotype_specific_mean_susc_sub_model(other)
{}

inline genotype_specific_susceptibility_sub_model&
genotype_specific_susceptibility_sub_model::operator=
        (const genotype_specific_susceptibility_sub_model& other)
{
  if(this != &other)
  {
    genotype_specific_mean_susc_sub_model::operator=(other);
  }
  
  return *this;
}

inline
genotype_specific_susceptibility_sub_model::~genotype_specific_susceptibility_sub_model()
{}

inline string
genotype_specific_susceptibility_sub_model::option_2_description(sm_option opt)
{
  return genotype_specific_sub_model::static_option_description(opt, susc);
}

//============================================================================
// IMPLEMENTATION:  genotype_specific_mean_sub_model
//============================================================================
//
inline  
genotype_specific_mean_sub_model::genotype_specific_mean_sub_model
      (cerrorstream& errors)
    : genotype_specific_mean_susc_sub_model(mean, errors) 
{ }

inline
genotype_specific_mean_sub_model::genotype_specific_mean_sub_model
      (const genotype_specific_mean_sub_model& other)
    : genotype_specific_mean_susc_sub_model(other)
{}

inline genotype_specific_mean_sub_model&
genotype_specific_mean_sub_model::operator=
        (const genotype_specific_mean_sub_model& other)
{
  if(this != &other)
  {
    genotype_specific_mean_susc_sub_model::operator=(other);
  }
  
  return *this;
}

inline
genotype_specific_mean_sub_model::~genotype_specific_mean_sub_model()
{}

inline string
genotype_specific_mean_sub_model::option_2_description(sm_option opt)
{
  return genotype_specific_sub_model::static_option_description(opt, susc);
}

//============================================================================
// IMPLEMENTATION:  genotype_specific_variance_sub_model
//============================================================================
//
inline  
genotype_specific_variance_sub_model::genotype_specific_variance_sub_model
      (genotype_specific_mean_susc_sub_model* m_ptr, cerrorstream& errors)
    : genotype_specific_sub_model(variance, errors), my_m_ptr(m_ptr) 
{
  //lint -e{534}
  set(one,
      model_input(VAR_DEFAULT_VALUE, VAR_DEFAULT_FIXED),
      model_input(QNAN, false),
      model_input(QNAN, false), false, false);

  my_default = true;
}

inline
genotype_specific_variance_sub_model::genotype_specific_variance_sub_model
      (const genotype_specific_variance_sub_model& other)
    : genotype_specific_sub_model(other)
{
  //lint -e{1554} <- Copy is ok
  my_m_ptr = other.my_m_ptr;
}

inline genotype_specific_variance_sub_model&
genotype_specific_variance_sub_model::operator=
        (const genotype_specific_variance_sub_model& other)
{
  if(this != &other)
  {
    genotype_specific_sub_model::operator=(other);

    //lint -e{1555} <- Copy is ok
    my_m_ptr = other.my_m_ptr;
  }
  
  return *this;
}


inline
genotype_specific_variance_sub_model::~genotype_specific_variance_sub_model()
{
  my_m_ptr = NULL;
}

inline string
genotype_specific_variance_sub_model::option_2_description(sm_option opt)
{
  return genotype_specific_sub_model::static_option_description(opt, variance);
}

// - Write sub-model values in LSF readable format.
//
inline void
genotype_specific_variance_sub_model::dump(std::ostream& out) const
{
  int  old_precision = out.precision();
  out.precision(DUMP_PRECISION);

  bool  var_AA_fixed = my_types_fixed[index_AA];
  bool  var_AB_fixed = my_types_fixed[index_AB];
  bool  var_BB_fixed = my_types_fixed[index_BB];

  //lint --e{666,506,40}

  assert(! SAGE::isnan(my_types[0]));
  assert(! SAGE::isnan(my_types[1]));
  assert(! SAGE::isnan(my_types[2]));

  out << "# " << name() << "\n"
      << "type_var\n" 
      << "{\n"
      << "  # " << option_2_description(my_option) << "\n"
      << "  option=" << option_2_parameter(my_option) << "\n"
      << "  var=AA, val=" << my_types[0] << ", fixed=" << std::boolalpha << var_AA_fixed << "\n" 
      << "  var=AB, val=" << my_types[1] << ", fixed=" << var_AB_fixed << "\n" 
      << "  var=BB, val=" << my_types[2] << ", fixed=" << var_BB_fixed << "\n" 
      << "}" << std::noboolalpha << std::endl;
      
  out.precision(old_precision);
}

}}
