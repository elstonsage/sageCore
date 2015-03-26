//============================================================================
// File:      sub_model_base.ipp
//                                                                          
// Author:    Geoff Wedig
//                                                                          
// History:   gcw      Initial Implementation                    Apr 2001 
//            Baechle  reformatted.  Added various functions.    Jun 2001
//                                                                          
// Notes:     
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================
#ifndef SEGREG_SUB_MODEL_BASE_H
#include "segreg/sub_model_base.h"
#endif

namespace SAGE
{

namespace SEGREG
{

inline
SegregSubmodel::SegregSubmodel (cerrorstream &errors)
  : MAXFUN::Submodel(errors)
{ }

inline
SegregSubmodel::SegregSubmodel(const SegregSubmodel &other)
  : MAXFUN::Submodel(other)
{ }

inline
SegregSubmodel& SegregSubmodel::operator= (const SegregSubmodel &other)
{
  MAXFUN::Submodel::operator=(other);
  
  return *this;
}    

inline
SegregSubmodel::~SegregSubmodel()
{ }

inline ostream& operator<< (ostream& out, const SegregSubmodel& sm)
{
  out << "\n" << sm.name() << " values: \n";
  out << "Option: " << sm.option_description() << std::endl;
  
  // Make sure the my_parameters vector is up to date, and print it out
  // It is often a more concise way of displaying what's in the model
  
  const_cast<SegregSubmodel&>(sm).finalizeConfiguration();

  OUTPUT::Table t;

  vector<MAXFUN::ParameterInput>::const_iterator  iter;
  for(iter = sm.my_parameters.begin(); iter != sm.my_parameters.end(); ++iter)
  {
    out << (OUTPUT::Table()
        << (OUTPUT::TableRow() << "maxfun parameter" << iter->param_name)
        << (OUTPUT::TableRow() << "value"            << iter->initial_estimate)
        << (OUTPUT::TableRow() << "status"           << MAXFUN::ParamTypeEnum2str(iter->initial_type))
        << (OUTPUT::TableRow() << "lower bound"      << iter->lower_bound)
        << (OUTPUT::TableRow() << "upper bound"      << iter->upper_bound));
  }

  return out;
}


inline genotype_index::genotype_index(unsigned int i)
  : my_type((gi_type) i)
{ }

inline genotype_index::genotype_index(gi_type i)
  : my_type(i)
{ }

inline genotype_index::genotype_index(const genotype_index& i)
  : my_type(i.my_type)
{ }

inline genotype_index& genotype_index::operator=(unsigned int i)
{
  my_type = (gi_type) i;

  return *this;
}

inline genotype_index& genotype_index::operator=(gi_type i)
{
  my_type = i;

  return *this;
}

inline genotype_index& genotype_index::operator=(const genotype_index& i)
{
  my_type = i.my_type;

  return *this;
}

inline genotype_index::operator unsigned int() const
{
  return my_type;
}

inline genotype_index& genotype_index::operator++()
{
  my_type = (gi_type) ((unsigned int) my_type + 1);

  return *this;
}

inline genotype_index genotype_index::operator++(int)
{
  genotype_index temp = *this;

  ++*this;

  return temp;
}

inline bool genotype_index::operator==(const genotype_index& g) const
{
  return my_type == g.my_type;
}

inline bool genotype_index::operator!=(const genotype_index& g) const
{
  return my_type != g.my_type;
}

inline bool genotype_index::operator==(gi_type g) const
{
  return my_type == g;
}

inline bool genotype_index::operator!=(gi_type g) const
{
  return my_type != g;
}

inline genotype_info  
index_2_info(genotype_index index)
{
  switch(index)
  {
    case index_AA:
      return AA;
    case index_AB:
      return AB;
    case index_BB:
      return BB;
    default:
      bool  bad_call_to_index_2_info = true;
      assert(! bad_call_to_index_2_info);
      return no_geno;
  }
}

// - Return a number which indicates what combination of model_inputs
//   have values.
//
inline genotype_info
total_info(const model_input& mi_AA, const model_input& mi_AB,
           const model_input& mi_BB)
{  
  return static_cast<genotype_info>( (SAGE::isnan(mi_AA.value) ? no_geno : AA) |
                                     (SAGE::isnan(mi_AB.value) ? no_geno : AB) |
                                     (SAGE::isnan(mi_BB.value) ? no_geno : BB)  );
}

// - Return the index which is not an argument.
//
inline genotype_index
third_index(genotype_index index_one, genotype_index index_two)
{
  genotype_index  return_value = index_AA;
  switch(index_one + index_two)
  {
    case 1:
      return_value = index_BB;
      break;
    case 2:
      return_value = index_AB;
      break;
    case 3:
      return_value = index_AA;
      break;
    default:
      assert(false);
  }
  
  return return_value;
}

// - For debugging.
//
inline std::ostream&  
operator<<(std::ostream& out, const model_input& mi)
{
  out << "value : " << mi.value << "\n"
      << std::boolalpha
      << "fixed:  " << mi.fixed << std::endl;
  
  return out;
}


}}
