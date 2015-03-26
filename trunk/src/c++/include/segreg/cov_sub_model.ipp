//============================================================================
// File:      cov_sub_model.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   7/20/01 - created.                            djb
//                                                                          
// Notes:     inline implementations of mean covariate sub-models for SEGREG.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#ifndef SEGREG_COV_SUB_MODEL_H
#include "segreg/cov_sub_model.h"
#endif

namespace SAGE
{

namespace SEGREG
{


inline std::ostream&  
operator<<(std::ostream& out, const CovariateSubmodel::covariate& cov)
{
  out << "\n" << "Name: " << cov.trait_name << "\n"
      << "Coefficient: " << cov.coefficient << "\n"
      << "Interaction: " << boolalpha << cov.has_interaction << "\n"
      << "Fixed: "       << boolalpha << cov.is_fixed << "\n"
      << "Interaction AA: " << cov.i_taus[index_AA] << "\n"
      << "Interaction AB: " << cov.i_taus[index_AB] << "\n"
      << "Interaction BB: " << cov.i_taus[index_BB] << "\n"
      << "Maxfun Index:   " << cov.maxfun_index << "\n"
      << std::endl;
  
  return out;
}

//============================================================================
// IMPLEMENTATION:  CovariateSubmodel::covariate
//============================================================================
//

inline
CovariateSubmodel::covariate::covariate
    (const std::string& trait, size_t indx, double coeff, bool inter, bool fixed)
    : trait_name(trait),
      trait_index(indx),
      coefficient(coeff),
      has_interaction(inter),
      is_fixed(fixed),
      maxfun_index(size_t(-1))
{
  for(int i = 0; i < NUM_OF_TYPES; i++)
  {
    i_taus[i] = INTERACTION_INIT_VALUE;
  }
}

//============================================================================
// IMPLEMENTATION:  CovariateSubmodel
//============================================================================
//

inline
CovariateSubmodel::CovariateSubmodel
  (const genotype_specific_mean_susc_sub_model* m_ptr,
   cerrorstream& errors)
    : SegregSubmodel(errors), my_m_ptr(m_ptr)
{ }

inline
CovariateSubmodel::CovariateSubmodel(const CovariateSubmodel& other)
  : SegregSubmodel(other)
{
  my_covariates = other.my_covariates;
  my_covariate_means = other.my_covariate_means;
  //lint -e{1554} <- Copy is ok
  my_m_ptr = other.my_m_ptr;
}

inline CovariateSubmodel&  
CovariateSubmodel::operator=(const CovariateSubmodel& other)
{
  if(this != &other)
  {
    SegregSubmodel::operator=(other);
    
    my_covariates = other.my_covariates;
    my_covariate_means = other.my_covariate_means;
    //lint -e{1555} <- Copy is ok
    my_m_ptr = other.my_m_ptr;
  }
  
  return *this;
}

inline
CovariateSubmodel::~CovariateSubmodel()
{
  my_m_ptr = NULL;
}

inline
size_t CovariateSubmodel::get_covariate_count() const
{
  return my_covariates.size();
}

inline
double CovariateSubmodel::get_covariate_mean(size_t i) const
{
  return my_covariate_means[i];
}

inline
void CovariateSubmodel::set_covariate_mean(size_t i, double d)
{
  my_covariate_means[i] = d;
}

inline const std::vector<CovariateSubmodel::covariate>&  
CovariateSubmodel::covariates() const
{
  return my_covariates;
}

inline string CovariateSubmodel::option_description() const
{
  return "";
}

// - Write sub-model values in LSF readable format.
//
inline void
CovariateSubmodel::dump(std::ostream& out) const
{
  int  old_precision = out.precision();
  out.precision(DUMP_PRECISION);

  out << std::boolalpha;
  out << "# " << name() << "\n"
      << get_subblock_name() << "\n" 
      << "{\n";

  vector<covariate>::const_iterator  cov_iter;
  for(cov_iter = my_covariates.begin(); cov_iter != my_covariates.end(); ++cov_iter)
  {
    assert(! SAGE::isnan(cov_iter->coefficient));
    
    out << "  covariate=" << cov_iter->trait_name 
        << ", val="       << cov_iter->coefficient
        << ", fixed="     << cov_iter->is_fixed;
        
    if(name() != COMPOSITE_TRAIT_NAME)
    {
      out << ", interaction=" << cov_iter->has_interaction;
    }
    
    out << "\n";
  }      
  
  out << "}" << std::noboolalpha << std::endl;
  out.precision(old_precision);  
}


inline std::string
CovariateSubmodel::get_cov_type_subblock_name(CovariateSubmodel::CovariateTypeEnum type)
{
  switch(type)
  {
    case ct_MEAN : return "mean_cov";
    case ct_VAR  : return "var_cov";
    case ct_SUSC : return "suscept_cov";
    case ct_COMP : return "composite_trait";
    default:        return "";
  }
}

template <CovariateSubmodel::CovariateTypeEnum T>
TypedCovariateSubmodel<T>::TypedCovariateSubmodel
    (const genotype_specific_mean_susc_sub_model* m_ptr,
     cerrorstream& errors)
  : CovariateSubmodel(m_ptr, errors)
{ }

template <CovariateSubmodel::CovariateTypeEnum T>
TypedCovariateSubmodel<T>::TypedCovariateSubmodel
    (const TypedCovariateSubmodel& other)
  : CovariateSubmodel(other)
{ }

template <CovariateSubmodel::CovariateTypeEnum T>
TypedCovariateSubmodel<T>&  TypedCovariateSubmodel<T>::operator=
    (const TypedCovariateSubmodel& other)
{
  if(this != &other)
  {
    CovariateSubmodel::operator=(other);
  }
  return *this;
}

template <CovariateSubmodel::CovariateTypeEnum T>
TypedCovariateSubmodel<T>::~TypedCovariateSubmodel()
{ }

template <CovariateSubmodel::CovariateTypeEnum T>
inline std::string
TypedCovariateSubmodel<T>::name() const
{
  switch(T)
  {
    case ct_MEAN : return MEAN_COV_NAME;
    case ct_VAR  : return VAR_COV_NAME;
    case ct_SUSC : return SUSCEPT_COV_NAME;
    case ct_COMP : return COMPOSITE_TRAIT_NAME;
    default:        return "";
  }
}
// - Return sub-block name i.e. parameter name in user documentation.
//
template <CovariateSubmodel::CovariateTypeEnum T>
inline std::string
TypedCovariateSubmodel<T>::get_subblock_name() const
{
  return get_cov_type_subblock_name(T);
}

template <CovariateSubmodel::CovariateTypeEnum T>
bool TypedCovariateSubmodel<T>::add_covariate
      (const RPED::MultiPedigree*    mp,
       const string&                 trait_name,
       double                        coefficient,
       bool                          interaction,
       bool                          fixed, 
       bool                          type_missing)
{
  BOOST_STATIC_ASSERT(T != ct_COMP && "This function can only be called on covariate submodels that have interactions");
  
  return base_add_covariate(mp, trait_name, coefficient, interaction, fixed, type_missing);
}

template <CovariateSubmodel::CovariateTypeEnum T>
bool TypedCovariateSubmodel<T>::add_covariate
      (const RPED::MultiPedigree*    mp,
       const string&                 trait_name,
       double                        coefficient,
       bool                          fixed, 
       bool                          type_missing)
{
  BOOST_STATIC_ASSERT(T == ct_COMP && "This function can only be called on covariate submodels that don't have interactions");

  return base_add_covariate(mp, trait_name, coefficient, false, fixed, type_missing);
}

}}
