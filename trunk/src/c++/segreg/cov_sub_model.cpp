//============================================================================
// File:      cov_sub_model.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   7/20/01 - created.                            djb
//                                                                          
// Notes:     implementation of covariate sub_model classes.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "segreg/cov_sub_model.h"
#include "boost/bind.hpp"

using namespace std;
namespace SAGE
{
namespace SEGREG
{

const std::string  MEAN_COV_NAME           = "Mean covariates";
const std::string  VAR_COV_NAME            = "Variance covariates";
const std::string  SUSCEPT_COV_NAME        = "Susceptibility covariates";
const std::string  COMPOSITE_TRAIT_NAME    = "Composite trait covariates";
const double       COV_DEFAULT_VALUE       = 0;
const bool         COV_DEFAULT_FIXED       = false;
const bool         COV_DEFAULT_INTERACTION = false;
const double       INTERACTION_INIT_VALUE  = 0;

//============================================================================
// IMPLEMENTATION:  CovariateSubmodel
//============================================================================
//
void  CovariateSubmodel::print_covariate_table
    (std::ostream& out, size_t char_indent) const
{
  char fill_char = out.fill();

  for(size_t i = 0; i < get_covariate_count(); ++i)
  {
    //lint --e{713}
    if(char_indent)
      out << setw(char_indent) << fill_char;

    out << setw(22) << my_covariates[i].trait_name << fill_char << ":" << fill_char;
    out << "(Mean = "
        << doub2str(my_covariate_means[i], 12, 5, ios::showpoint | ios::fixed)
        << ")";

    if(my_covariates[i].has_interaction)
      out << fill_char << fill_char << "(interaction)";

    out << endl;
  }
}

void CovariateSubmodel::calculate_covariate_means(const FPED::Multipedigree& rmp)
{
  for(size_t i = 0; i < get_covariate_count(); ++i)
    calculate_covariate_mean(rmp, i);
}

void CovariateSubmodel::calculate_covariate_mean(const FPED::Multipedigree& rmp, size_t i)
{
  double cov_mean   = 0.0;
  size_t mean_count = 0;

  const covariate& cov = my_covariates[i];

  FPED::PedigreeConstIterator ped = rmp.pedigree_begin();

  for( ; ped != rmp.pedigree_end(); ++ped)
  {
    //lint --e{115} This is a spurious error from lint

    FPED::MemberConstIterator mem = ped->member_begin();

    for( ; mem != ped->member_end(); ++mem)
    {
      size_t tr_index = rmp.info().trait_find(cov.trait_name);
      double value    = ped->info().trait(mem->index(),tr_index);

      if(finite(value))
      {
        cov_mean += value;

        mean_count++;
      }
    }
  }

  if(mean_count)
  {
    cov_mean /= mean_count;

    my_covariate_means[i] = cov_mean;
  }
}

int  
CovariateSubmodel::update()
{
  assert(my_m_ptr != 0);

  uint  num_of_means = my_m_ptr->get_type_count();

  // Update each covariate
  
  vector<covariate>::iterator  c_iter;
  for(c_iter = my_covariates.begin(); c_iter != my_covariates.end(); c_iter++)
  {
    // Get the covariate's value
    
    size_t  c_index = c_iter->maxfun_index;

    c_iter->coefficient = getParam(c_index);
    
    // If there are interactions, retrieve those too.

    if(c_iter->has_interaction && num_of_means > 1)
    {
      if(num_of_means == 2)
      {
        // In the case of two means, the second tau value is equal to the
        // negative value of the first.

        getParam(c_index + 2) = - getParam(c_index + 1);
        
        // Copy the taus into the tau matrix.  The tau for the AB value depends
        // on whether we're two dominant or two recessive.

        c_iter->i_taus[index_AA] = getParam(c_index+1);
             
        if(my_m_ptr->is_two_dom())
        {
          c_iter->i_taus[index_AB] = getParam(c_index+1);
        }
        else
        {
          c_iter->i_taus[index_AB] = getParam(c_index+2);
        }
      
        c_iter->i_taus[index_BB] = getParam(c_index+2);
      }
      else // num_of_means == 3
      {
        // In the case of three means, the third tau value is equal to the
        // negative sum of the first two.
        
        getParam(c_index + 3) = - (getParam(c_index + 1) + getParam(c_index + 2) );
        
        // Copy the taus into the tau matrix.

        c_iter->i_taus[index_AA] = getParam(c_index+1);
        c_iter->i_taus[index_AB] = getParam(c_index+2);
        c_iter->i_taus[index_BB] = getParam(c_index+3);
      }
    }
  }
  
  return 0;
}

bool
CovariateSubmodel::has_covariate(const std::string& cov) const
{
  return covariate_index(cov) != (size_t) -1;
}

/// \internal 
///
/// Returns true if the covariates share the
/// same trait_name.
///
/// \param lhs The first covariate
/// \param rhs The second covariate

inline
bool  compare_covariate_trait_names
  (const CovariateSubmodel::covariate& lhs,
   const CovariateSubmodel::covariate& rhs)
{
  return toUpper(lhs.trait_name) == toUpper(rhs.trait_name);
}

/// finds the covariate and returns its index.  Returns -1 if not found
///
size_t
CovariateSubmodel::covariate_index(const std::string& cov) const
{
  vector<covariate>::const_iterator i = 
      find_if(my_covariates.begin(), my_covariates.end(),
              boost::bind(compare_covariate_trait_names,
                covariate(string(cov), (size_t) -1, QNAN, false, false), _1));

  if( i == my_covariates.end() ) return (size_t) -1;

  return (size_t) (i - my_covariates.begin());
}

std::vector<CovariateSubmodel::covariate>::iterator 
CovariateSubmodel::find(const covariate& cov)
{
  vector<covariate>::iterator  begin = my_covariates.begin();
  vector<covariate>::iterator  end   = my_covariates.end();
  
  return find_if(begin, end, boost::bind(compare_covariate_trait_names, cov, _1));

}

int
CovariateSubmodel::finalizeConfiguration()
{
  my_parameters.resize(0);
  
  vector<covariate>::iterator  c_iter;
  for(c_iter = my_covariates.begin(); c_iter != my_covariates.end(); ++c_iter)
  {
    c_iter->maxfun_index = my_parameters.size();
    
    my_parameters.push_back(
      MAXFUN::ParameterInput
        ( name(),
          c_iter->trait_name,
          c_iter->is_fixed ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT, 
          c_iter->coefficient,
          NEGATIVE_INF,
          POSITIVE_INF) );

    if(c_iter->has_interaction)
    {
      initialize_interactions(c_iter);
    }
  }

  return 0;
}

void
CovariateSubmodel::initialize_interactions(vector<covariate>::const_iterator c_iter)
{
  assert(my_m_ptr != 0);

  // Initialize trait names, used for two mean models

  string tname1;
  string tname2;

  switch(my_m_ptr->get_type_count())
  {
    case 1:
      break;
  
    case 2:
      //lint -e{788} We don't care about options other than two*
      if(my_m_ptr->is_two_dom())
      {
        tname1 = c_iter->trait_name + " * AA_AB";
        tname2 = c_iter->trait_name + " * BB";
      }
      else
      {
        tname1 = c_iter->trait_name + " * AA";
        tname2 = c_iter->trait_name + " * AB_BB";
      }
        
      my_parameters.push_back(  
        MAXFUN::ParameterInput
          ( name(), tname1,
            MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL, 
            c_iter->i_taus[0],
            NEGATIVE_INF,
            POSITIVE_INF));
      my_parameters.push_back(  
        MAXFUN::ParameterInput
          ( name(), tname2,
            MAXFUN::Parameter::DEPENDENT, 
            c_iter->i_taus[1],
            NEGATIVE_INF,
            POSITIVE_INF));
      break;
      
    case 3:
      my_parameters.push_back(  
        MAXFUN::ParameterInput
          ( name(), c_iter->trait_name + " * AA",
            MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL, 
            c_iter->i_taus[0],
            NEGATIVE_INF,
            POSITIVE_INF));
      my_parameters.push_back(  
        MAXFUN::ParameterInput
          ( name(), c_iter->trait_name + " * AB",
            MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL, 
            c_iter->i_taus[1],
            NEGATIVE_INF,
            POSITIVE_INF));
      my_parameters.push_back(  
        MAXFUN::ParameterInput
          ( name(), c_iter->trait_name + " * BB",
            MAXFUN::Parameter::DEPENDENT, 
            c_iter->i_taus[2],
            NEGATIVE_INF,
            POSITIVE_INF));
      break;
      
    default: 
      assert(false);
  }
}

/// Add an element to the member vector of covariates.  If this covariate 
/// is already part of the sub-model, overwrite it; otherwise, add it.
void  
CovariateSubmodel::insert_covariate(const covariate& cov)
{
  vector<covariate>::iterator  c_iter = find(cov); 
  
  if(c_iter == my_covariates.end())
  {
    my_covariates.push_back(cov);
    
    my_covariate_means.resize(my_covariates.size(), QNAN);
  }
  else     // - Covariate in question already part of sub-model.
  {
    *c_iter = cov;    
  }
}


bool  
CovariateSubmodel::base_add_covariate
    (const RPED::MultiPedigree* mp,
     const string& trait_name,
     double coefficient,
     bool interaction,
     bool fixed, 
     bool type_missing)
{
  bool  return_value = false;

  assert(mp != 0);
  if(mp == 0)
  {
    return return_value;
  }
  
  size_t  trait_index = mp->info().trait_find(trait_name);
  
  if(trait_index != (size_t) -1)
  {
    // - Create 'default' covariate
    //
    covariate  cov = covariate(trait_name, trait_index, COV_DEFAULT_VALUE, COV_DEFAULT_INTERACTION, COV_DEFAULT_FIXED);
  
    // Test the covariate.  The covariate is valid provided it isn't nan while
    // the fixed flag is set.
    if(! (SAGE::isnan(coefficient) && fixed) )
    {
      // - Set interaction, if any
      //
      if(interaction)
      {
        assert(my_m_ptr != 0);
        uint  means = my_m_ptr->get_type_count();
      
        if(means == 2 || means == 3 || type_missing)
        {
          cov.has_interaction = interaction;
        }
        else
        {
          my_errors << priority(error) << "interactions not allowed "
                    << "for covariate '" << trait_name << "' in " 
                    << get_subblock_name() << " sub-block.  " 
                    << "There must be more than one "
                    << my_m_ptr->long_sing_type() << ".  Ignoring ..." << endl;
          cov.has_interaction = false;
        }
      }
      
      // - Set coefficient parameter member values.
      //
      if(! SAGE::isnan(coefficient)) 
      {
        cov.coefficient = coefficient;
      }
    
      // Set fixed status
      cov.is_fixed = fixed;

      insert_covariate(cov);
      
      return_value = true;
    }
    else
    {
      my_errors << priority(critical) << "No coefficient value given for fixed "
                << "covariate '" << trait_name << "' in " 
                << get_subblock_name() << " sub-block.  "  
                << "Skipping analysis ..." << endl;
    }
  }
  else
  {
    my_errors << priority(critical) << "Covariate '" << trait_name << "' specified "
              << "in " << get_subblock_name() << " sub-block "  
              << "is not a valid trait.  Skipping analysis ..." << endl;
  }

  return return_value;
}                                  

std::list<std::string> get_shared_covariates(const CovariateSubmodel& sm1,
                                             const CovariateSubmodel& sm2)
{
  std::list<std::string> shared_covs;
  
  typedef std::vector<CovariateSubmodel::covariate> CovVectorType;
  
  const CovVectorType& covs1 = sm1.covariates();

  for(CovVectorType::const_iterator citer = covs1.begin(); citer != covs1.end();
      ++citer)
  {
    if(sm2.has_covariate(citer->trait_name))
    {
      shared_covs.push_back(citer->trait_name);
    }
  }
  
  return shared_covs;
}

bool are_covariates_exclusive(const CovariateSubmodel& sm1,
                              const CovariateSubmodel& sm2)
{
  typedef std::vector<CovariateSubmodel::covariate> CovVectorType;
  
  const CovVectorType& covs1 = sm1.covariates();

  for(CovVectorType::const_iterator citer = covs1.begin(); citer != covs1.end();
      ++citer)
  {
    if(sm2.has_covariate(citer->trait_name)) return false;
  }
  
  return true;
}

}
}
