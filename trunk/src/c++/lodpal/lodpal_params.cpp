//****************************************************************************
//* File:      lodpal_params.cpp                                             *
//*                                                                          *
//* Author:    Kevin Jacobs & Yeunjoo Song                                   *
//*                                                                          *
//* History:   Version 0.0 Initial implementation                kbj         *
//*                    1.0 one-parameter model added.            yjs Nov. 00 *
//*                    1.1 covariate added.                      yjs Nov. 00 *
//*                    1.2 bad_sib_pair storage & func added.    yjs Mar. 01 *
//*                    1.3 rel_pair_map storage & func added.    yjs Mar. 01 *
//*                    1.3 multipoint, singlepoint seperated.    yjs Apr. 01 *
//*                    1.4 dsp, re-parameterization added.       yjs Apr. 01 *
//*                    1.5 evaluate & update_bound for dsp added.yjs May. 01 *
//*                    1.6 diagnostic option added.              yjs Jun. 01 *
//*                    1.7 parameters seperated from arp.cpp.    yjs Jul. 01 *
//*                    2.0 weight added.                         yjs Jan. 02 *
//*                                                                          *
//* Notes:     This file implements ARPTest class - parameters only.         *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "lodpal/lodpal_params.h"

namespace SAGE   {
namespace LODPAL {

////////////////////////////////////////////////////////////////////////////
//          Implementation of lodpal_params (non-Inline)                  //
////////////////////////////////////////////////////////////////////////////

//
// ----------------------------------------------------------------------
//

parameter_estimate::parameter_estimate()
{
  clear();
}

parameter_estimate::parameter_estimate(double d)
{
  clear();
  set_initial_value(d);  
}

void
parameter_estimate::clear()
{
  const double qNaN = std::numeric_limits<double>::quiet_NaN();

  set_value(qNaN);
  set_initial_value(qNaN);
  set_first_derivative(qNaN);
  release_value();
}

//
// ----------------------------------------------------------------------
//

trait_parameter::trait_parameter(size_t t, trait_type tt, double cp)
{
  trait = t;
  pair_select = tt;
  cutpoint = cp;
  info.clear();
  valid = false;
}

string
trait_parameter::name(const RelativePairs& pairs) const
{
  string t = pairs.trait_name(trait);
  if(!t.size())
    t = "Trait #" + long2str(trait);
  return t;
}

//
// ----------------------------------------------------------------------
//

covariate_type::covariate_type() 
              : covariate((size_t)-1), power(1.0), operation(both),
                adjust(mean), adjust_value(std::numeric_limits<double>::quiet_NaN())
{ }

covariate_type::covariate_type(size_t c, cov_ad ad, cov_op op, double ad_val,
                               double powr, double d1, double d2)
              : covariate(c), power(powr), operation(op), adjust(ad), adjust_value(ad_val),
                delta1(d1), delta2(d2)
{ }

string
covariate_type::name(const RelativePairs& pairs) const
{
  string n;

  if( operation == pair )
    n = pairs.get_pair_covariate_name(covariate);
  else
    n = pairs.trait_name(covariate);

  if(!n.size())
    n = "Covariate #" + long2str(covariate);

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
    case none:    n = "CovC"; break;
    case sum:     n = "Cov+"; break;
    case diff:    n = "Cov-"; break;
    case prod:    n = "Cov*"; break;
    case avg:     n = "Cov+/"; break;
    case single:  n = "CovF"; break;
    case pair:    n = "CovP"; break;
    default:      break;
  }

  switch( adjust )
  {
    case mean:    n += " mean"; break;
    case minimum: n += " minimum = 0.0"; break;
    case dsp:     n += " concordance"; break;
    case prop:    n += " proportion"; break;
  }
  return n;
}

string
covariate_type::short_effect_name() const
{
  string n;

  switch( operation )
  {
    case none:    n = "CC"; break;
    case sum:     n = "C+"; break;
    case diff:    n = "C-"; break;
    case prod:    n = "C*"; break;
    case avg:     n = "C+/"; break;
    case single:  n = "CF"; break;
    case pair:    n = "CP"; break;
    default:      break;
  }

  switch( adjust )
  {
    case mean:    n += " mean"; break;
    case minimum: n += " mini"; break;
    case dsp:     n += " conc"; break;
    case prop:    n += " prop"; break;
  }
  return n;
}

//
// ----------------------------------------------------------------------
//

autosomal_model_type::autosomal_model_type(model_type one,
                                           constraint_type con,
                                           poo_fixed_type po,
                                           bool poo, double al)
                    : model(one), constraint(con), fixed(po), parent_of_origin(poo), alpha(al)
{ }

string
autosomal_model_type::name() const
{
  string m = "one-parameter model";

  if( model == two_parameter )
    m = "two-parameter model";

  if( constraint == constrained )
    m = m + ", constrained";
  else
    m = m + ", unconstrained";

  if( model == one_parameter )
    m = m + ", alpha = " + doub2str(alpha);

  return m;
}

string
autosomal_model_type::poo_name() const
{
  string model_description;

  model_description += "\n               parent_of_origin model";

  if( fixed == maternal )
    model_description += ", lambda1m fixed to 1";
  else if( fixed == paternal )
    model_description += ", lambda1p fixed to 1";

  return model_description;
}

//
// ----------------------------------------------------------------------
//

x_linkage_model_type::x_linkage_model_type(bool l1, bool l2, double al)
                    : lambda1_equal(l1), lambda2_fixed(l2), alpha(al)
{ }

string
x_linkage_model_type::name() const
{
  string l1 = "lambda1 equal";
  string l2 = "lambda2 fixed, alpha = ";

  if( !lambda1_equal )
    l1 = "lambda1 not equal";

  if( !lambda2_fixed )
    l2 = "lambda2 not fixed";
  else
    l2 = l2 + doub2str(alpha);

  return l1 + ", " + l2;
}

//
// ----------------------------------------------------------------------
//

marker_type::marker_type(size_t m, inheritance_type in, double b1, double b2)
           : marker(m), inheritance(in), beta1(b1), beta2(b2),
             beta1mm(0.1), beta1mf(0.1), beta1ff(0.1), beta2ff(0.1)
{ }

string
marker_type::name(const RelativePairs& pairs) const
{
  string n = pairs.marker_name(marker);
  if(!n.size())
    n = "Marker #" + long2str(marker);
  return n;
}

string
marker_type::effect_name() const
{
  if( inheritance == x_linked )
    return "autosomal";

  return "autosomal";
}

string
marker_type::short_effect_name() const
{
  if( inheritance == x_linked )
    return "X";

  return "A";
}

//
// ----------------------------------------------------------------------
//

independent_variable::independent_variable() 
                    : valid(false)
{ }

void
independent_variable::clear()
{
  valid = false;
  markers.resize(0);
  covariates.resize(0);
}

string
independent_variable::covariate_name(const RelativePairs& pairs) const
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
independent_variable::covariate_effect_name() const
{
  string n;

  if( !covariates.size() )
    return "";

  if( covariates.size() == 1 )
    return covariates[0].effect_name();

  n = covariates[0].short_effect_name();
  for(size_t i = 1; i < covariates.size(); ++i)
  {
    n += " x " + covariates[i].short_effect_name();
  }
  return n;
}

string
independent_variable::marker_name(const RelativePairs& pairs) const
{
  if( !markers.size() )
    return "";

  string n = markers[0].name(pairs);
  for(size_t i = 1; i < markers.size(); ++i)
    n += " x " + markers[i].name(pairs);

  return n;
}

string
independent_variable::marker_effect_name() const
{
  string n;

  if( !markers.size() )
    return "";

  if( markers.size() == 1 )
    return markers[0].effect_name();

  n = markers[0].short_effect_name();
  for(size_t i = 1; i < covariates.size(); ++i)
  {
    n += " x " + markers[i].short_effect_name();
  }
  return n;
}

string
independent_variable::name(const RelativePairs& pairs) const
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

  return marker_effect_name() + " x " + covariate_effect_name();
}

//
// ----------------------------------------------------------------------
//

lodpal_weight_type::lodpal_weight_type(size_t c, weight_op op)
                  : weight(c), operation(op)
{ }

string
lodpal_weight_type::name(const RelativePairs& pairs) const
{
  string n;

  if( operation == pair )
    n = pairs.get_pair_weight_name(weight);
  else
    n = pairs.trait_name(weight);

  if(!n.size())
    n = "Weight #" + long2str(weight);

  return n;
}

//
// ----------------------------------------------------------------------
//

lodpal_parameters::lodpal_parameters()
{
  clear();
}

void
lodpal_parameters::clear()
{
  set_skip_uninformative_pairs(false);
  set_diagnostic_marker((size_t)-1);
  set_turn_off_default(false);
  set_print_lambda(false);
  set_sib_pairs_only(false);
  set_autosomal_marker_exist(false);
  set_x_linked_marker_exist(false);
  set_x_linkage_pair_type(true, true, true);

  clear_traits();
  clear_subsets();
  clear_parameters();
  clear_weight();
  clear_autosomal_model();
  clear_x_linkage_model();

  invalidate();
}

void
lodpal_parameters::invalidate()
{
  my_valid = false;
}

void
lodpal_parameters::dump_parameters(const RelativePairs &pairs) const
{
  // Markers name
  cout << "  List of markers : " << endl;
  lodpal_parameters::marker_const_iterator mi;
  for(mi = marker_begin(); mi != marker_end(); ++mi)
  {
    cout << "        " << mi->name(pairs) << " : " << mi->effect_name() << endl;
  }

  // Covariates name
  lodpal_parameters::covariate_const_iterator ci;
  cout << "  List of covariates : " << endl;
  for(ci = covariate_begin(); ci != covariate_end(); ++ci)
  {
    cout << "        " << ci->name(pairs) << " : " << ci->effect_name() << endl;
  }
  cout << endl;
}

lodpal_parameters::tib_value
lodpal_parameters::add_trait(const trait_parameter& t)
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
//    *i = t;   This is a bad idea
  }
  return std::pair<trait_iterator, bool>(i,inserted);
}

lodpal_parameters::tib_value
lodpal_parameters::add_subset(const trait_parameter& t)
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
//    *i = t;   This is a bad idea
  }
  return std::pair<trait_iterator, bool>(i,inserted);
}

lodpal_parameters::mib_value
lodpal_parameters::add_parameter(const marker_type& m)
{
  // Check to see if marker is already exists

  marker_iterator i = std::find(my_parameters.markers.begin(),
                                my_parameters.markers.end(),   m);

  bool inserted = true;
  if( i == my_parameters.markers.end() )
  {
    invalidate();
    my_parameters.markers.push_back( m );
    i = my_parameters.markers.end();
    --i;
  }
  else
  {
    inserted = false;
//    *i = p;  // WHY?  This looks to be very bad
  }
  return std::pair<marker_iterator, bool>(i,inserted);
}

lodpal_parameters::cib_value
lodpal_parameters::add_parameter(const covariate_type& c)
{
  // Check to see if covariate is already exists

  covariate_iterator i = std::find(my_parameters.covariates.begin(),
                                   my_parameters.covariates.end(),   c);

  bool inserted = true;
  if( i == my_parameters.covariates.end() )
  {
    invalidate();
    my_parameters.covariates.push_back( c );
    i = my_parameters.covariates.end();
    --i;
  }
  else
  {
    inserted = false;
//    *i = p;  // WHY?  This looks to be very bad
  }
  return std::pair<covariate_iterator, bool>(i,inserted);
}

bool
lodpal_parameters::add_weight(const lodpal_weight_type& w)
{
  if( my_weight == w )
    return false;

  my_weight = w;

  return true;
}

bool
lodpal_parameters::add_autosomal_model(const autosomal_model_type & m)
{
  if( my_autosomal_model == m )
    return false;

  my_autosomal_model = m;

  return true;
}

bool
lodpal_parameters::add_x_linkage_model(const x_linkage_model_type & m)
{
  if( my_x_linkage_model == m )
    return false;

  my_x_linkage_model = m;

  return true;
}

} // end of namespace LODPAL
} // end of namespace SAGE
