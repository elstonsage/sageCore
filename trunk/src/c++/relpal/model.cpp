//=============================================================================
// File:    model.cpp
//
// Author:  Yeunjoo Song
//
// History: Initial implementation                                  yjs Mar. 07
//
// Notes:   This file contains implementation for following data structures.
//            class regression_model
//
// Copyright (c) 2007 R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "relpal/model.h"

using namespace std;

namespace SAGE   {
namespace RELPAL {

regression_model::regression_model()
{
  my_name = "";
  clear();
}

void
regression_model::clear()
{
  clear_traits();
  clear_ind_parameters();
  clear_ped_parameters();

  invalidate();
}

void
regression_model::invalidate()
{
  my_valid = false;
}

void
regression_model::validate()
{
  my_valid = true;
}

tib_value
regression_model::add_trait(const dependent_variable& t)
{
  trait_iterator i = std::find(my_traits.begin(),
                               my_traits.end(),   t);

  bool inserted = true;
  if( i == my_traits.end() )
  {
    invalidate();
    my_traits.push_back( t );
    i = my_traits.end();
    --i;
  }
  else
  {
    inserted = false;
  }
  return std::pair<trait_iterator, bool>(i,inserted);
}

pib_value
regression_model::add_ind_parameter(const independent_variable& p)
{
  // Check to see if parameter is already exists

  parameter_iterator i = std::find(my_ind_parameters.begin(),
                                   my_ind_parameters.end(),   p);

  bool inserted = true;
  if( i == my_ind_parameters.end() )
  {
    invalidate();
    my_ind_parameters.push_back( p );
    i = my_ind_parameters.end();
    --i;
  }
  else
  {
    inserted = false;
  }
  return std::pair<parameter_iterator, bool>(i,inserted);
}

pib_value
regression_model::add_ped_parameter(const independent_variable& p)
{
  // Check to see if parameter is already exists

  parameter_iterator i = std::find(my_ped_parameters.begin(),
                                   my_ped_parameters.end(),   p);

  bool inserted = true;
  if( i == my_ped_parameters.end() )
  {
    invalidate();
    my_ped_parameters.push_back( p );
    i = my_ped_parameters.end();
    --i;
  }
  else
  {
    inserted = false;
  }
  return std::pair<parameter_iterator, bool>(i,inserted);
}

void
regression_model::clear_markers()
{
  parameter_vector  new_parameters;

  parameter_iterator i = my_ped_parameters.begin();
  for( ; i != my_ped_parameters.end(); ++i )
  {
    if( i->type != independent_variable::MARKER )
      new_parameters.push_back(*i);      
  }

  my_ped_parameters.clear();
  
  i = new_parameters.begin();
  for( ; i != new_parameters.end(); ++i )
    my_ped_parameters.push_back(*i);

  return;
}

void
regression_model::dump_model(const relative_pairs& relpairs, ostream &out) const
{
  out << endl
      << "==========================" << endl
      << "  Regression model dump  " << endl
      << "==========================" << endl << endl;

  out << "Model name: " << get_name() << endl;

  if( valid() )
    out << "  Model valid!!";
  else
    out << "  Model invalid!";
  out << endl << endl;

  out << "Dependent variables:" << endl;
  for( size_t i = 0; i < my_traits.size(); ++i )
  {
    out << " " << i+1 << " : " << my_traits[i].name(relpairs) << " ";
  }
  out << endl << endl;

  out << "Independent variables - individual level:" << endl;
  for( size_t i = 0; i < my_ind_parameters.size(); ++i )
  {
    out << " " << i+1 << " : " << my_ind_parameters[i].dump(relpairs)
        << endl;
  }
  out << endl;

  out << "Independent variables - pedigree level:" << endl;
  for( size_t i = 0; i < my_ped_parameters.size(); ++i )
  {
    out << " " << i+1 << " : " << my_ped_parameters[i].dump(relpairs)
        << endl;
  }
  out << endl;

  out << "P-value options:" << endl;
  out << my_pvalue_opt.dump(" ") << endl;

  if( my_analysis_opt.transform_residuals )
    out << "Transform residuals from 1st level to 2nd level" << endl;

  if(    my_analysis_opt.naive_variance  || my_analysis_opt.sandwich_variance
      || my_analysis_opt.alternative_variance || my_analysis_opt.IBD_variance )
  {
    out << "Variance options:" << endl;
    out << my_analysis_opt.dump_var(" ") << endl;
  }

  return;
}

} // end of namespace RELPAL
} // end of namespace SAGE
