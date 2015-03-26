//****************************************************************************
//* File:      lodpal_parser.cpp                                             *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0 Initial implementation.                yjs Jan 01 *
//*                    1.0 x-linkage added.                       yjs May 02 *
//*                                                                          *
//* Notes:     This header file implements lodpal_parser class.              *
//*                                                                          *
//* Copyright (c) 2001 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "lodpal/lodpal_parser.h"

namespace SAGE   {
namespace LODPAL {

void lodpal_parser::parse_symbols(const SymbolTable* syms)
{
  if(!syms)
    return;

  AttrVal v = syms->query("wide_out");
  if(v.has_value())
  {
    if( toUpper(v.String()) == "TRUE" )
      set_wide_output(true);
    else if( toUpper(v.String()) == "FALSE" )
      set_wide_output(false);
  }

  v = syms->query("csv_out");
  if(v.has_value())
  {
    if( toUpper(v.String()) == "TRUE" )
      set_csv_output(true);
    else if( toUpper(v.String()) == "FALSE" )
      set_csv_output(false);
  }
}

void lodpal_parser::
     parse_parameter(const LSFBase* param)
{
  if( !param || !param->name().size()) return;
  AttrVal a;

  string name = toUpper( param->name() );

  if     (name == "WIDE_OUT")
    parse_boolean(param, my_wide_output);
  else if(name == "CSV_OUT")
    parse_boolean(param, my_csv_output);

  if( has_attr(param,"lambda") )
    my_parameters.set_print_lambda(true);
}

void lodpal_parser::
parse_test_parameter_section(const LSFBase* params)
{
  if(!params)
    return;

  if(!params->List())
    return;

  my_pair_info_file.resize(0);

  LSFList::const_iterator i;
  for(i=params->List()->begin(); i!=params->List()->end(); ++i)
  {
    if( !(*i) || !(*i)->name().size())
      continue;

    if( toUpper((*i)->name()) == "PAIR_INFO_FILE" )
      my_pair_info_file.push_back(*i);
    else
      parse_test_parameter(*i);
  }

  // Set autosomal_model. - backward compatibility
  //
  if( !my_autosomal_model_exist )
  {
    autosomal_model_type::model_type      m       = autosomal_model_type::one_parameter;
    autosomal_model_type::constraint_type cont    = autosomal_model_type::constrained;
    autosomal_model_type::poo_fixed_type  poo_fix = autosomal_model_type::none;

    bool poo_test = false;
    double al     = alpha();

    if( !one_parameter_model() )
      if( al == -std::numeric_limits<double>::infinity() )
        m = autosomal_model_type::two_parameter;
      else if( al == std::numeric_limits<double>::infinity() )
        cont  = autosomal_model_type::unconstrained;

    my_parameters.add_autosomal_model(m, cont, poo_fix, poo_test, al);
  }
}

void lodpal_parser::
     parse_test_parameter(const LSFBase* param)
{
  if( !param || !param->name().size()) return;
  AttrVal a;

  string name = toUpper( param->name() );

  if(name == "ALPHA")
    parse_alpha(param);
  else if(name == "TWO_PARAMETER_MODEL" || name == "TWO_PARAMETER" || name == "TWO")
    parse_two_model(param);
  else if(name == "UNCONSTRAINED" || name == "UNCON")
    parse_unconstrained(param);
  else if(name == "MARKER")
    parse_marker(param);
  else if(name == "TRAIT")
    parse_trait(param);
  else if(name == "SUBSET")
    parse_subset(param);
  else if(name == "COVARIATE")
    parse_covariate(param);
  else if(name == "WEIGHT")
    parse_weight(param);
  else if(name == "DIAGNOSTIC")
    parse_diagnostic(param);
  else if(name == "WIDE_OUT")
    parse_boolean(param, my_wide_output);
  else if(name == "CSV_OUT")
    parse_boolean(param, my_csv_output);
  else if(name == "TURN_OFF_DEFAULT")
    my_parameters.set_turn_off_default(true);
  else if(name == "SIB_PAIRS_ONLY" || name == "SIB")
    my_parameters.set_sib_pairs_only(true);
  else if(name == "X_LINKAGE_MODEL" || name == "X_LINKAGE" )
    parse_x_linkage_model(param);
  else if(name == "AUTOSOMAL_MODEL" || name == "AUTOSOMAL" )
    parse_autosomal_model(param);
  else if(name == "PVAL_SCIENTIFIC_NOTATION" || name == "PVALUES_SCIENTIFIC_NOTATION" )
    parse_boolean(param, my_pval_scientific_notation);
  else
    errors << priority(information)
           << "Invalid parameter for lodpal_analysis : "
           << param->name() << "\n             Skipped..." << endl;
}

void lodpal_parser::
     parse_alpha(const LSFBase* param)
{
  errors << priority(warning)
         << "This option 'alpha' has been revised.  Please refer to the manual." << endl;

  my_one_parameter_model = true;
  my_alpha = 2.634;

  AttrVal a=attr_value(param,0);
  if( a.has_value() )
  {
    if( finite(a.Real()) && a.Real() >= 1 )
      my_alpha = a.Real();
  }
}

void lodpal_parser::
     parse_two_model(const LSFBase* param)
{
  errors << priority(warning)
         << "This option 'two_parameter_model' has been revised.  Please refer to the manual." << endl;

  // Do two_parameter_model only if discordant pairs not pooled with arp.
  if( !discordant() )
  {
    my_one_parameter_model = false;
    my_alpha = -std::numeric_limits<double>::infinity();
  }
}

void lodpal_parser::
     parse_unconstrained(const LSFBase* param)
{
  errors << priority(warning)
         << "This option 'uncontrained' has been revised.  Please refer to the manual." << endl;

  my_one_parameter_model = false;
  my_alpha = std::numeric_limits<double>::infinity();
}

void lodpal_parser::
     parse_marker(const LSFBase* param)
{
  AttrVal a=attr_value(param,0);
  if( a.has_value() )
  {
    string mname = strip_ws(a.String());
    size_t m = pairs.marker_find(mname);
    if(m < pairs.marker_count())
    {
      double init_b1 = 0.1;
      double init_b2 = 0.1;
#if 0
      if( has_attr(param,"initial_beta1"))
      {
        init_b1 = param->attrs()->FloatAttr("initial_beta1");

        if( !finite(init_b1) )
          init_b1 = 0.1;
      }

      if( has_attr(param,"initial_beta2"))
      {
        init_b2 = param->attrs()->FloatAttr("initial_beta2");

        if( !finite(init_b2) )
          init_b2 = 0.1;
      }
#endif

      marker_type::inheritance_type in = marker_type::autosomal;
      if( pairs.is_x_linked(m) )
      {
        in = marker_type::x_linked;
        my_parameters.set_x_linked_marker_exist(true);
      }
      else
        my_parameters.set_autosomal_marker_exist(true);

      my_parameters.add_marker(m, in, init_b1, init_b2);
    }
    else
      errors << priority(error) << "Marker '" << mname << "' not found." << endl;
  }
  else
    errors << priority(error)
           << "Marker parameter is missing name.  Skipping..." << endl;
}

void lodpal_parser::
     parse_trait(const LSFBase* param)
{
  AttrVal a=attr_value(param,0);
  if( a.has_value() )
  {
    string tname = strip_ws(a.String());
    size_t t = pairs.trait_find(tname);
    if(t < pairs.trait_count())
    {
      double cp = 0.;
      if(pairs.fped_info().trait_info(t).type() != RPED::RefTraitInfo::binary_trait)
      {
        if( param->attrs()->has_attr("cutpoint") )
          cp = param->attrs()->FloatAttr("cutpoint");

        errors << priority(warning) 
               << "Trait '" << tname << "' is not a binary trait!  For this "
                  "analysis it will be dichotomized at "
               << cp
               << " (trait values <= "
               << cp
               << " will be treated as unaffected and values > "
               << cp
               << " will be treated as affected)."
               << endl;
      }

      trait_parameter::trait_type tt = trait_parameter::conaff;

      if( param->attrs()->has_attr("condisc") )
      {
        tt = trait_parameter::condisc;

        set_discordant(true);
        set_one_parameter_model(true);
        set_alpha(2.634);

        //cout << "condisc" << endl;
      }
      else if( param->attrs()->has_attr("noconunaff") )
      {
        tt = trait_parameter::noconunaff;

        set_discordant(true);
        set_one_parameter_model(true);
        set_alpha(2.634);

        //cout << "noconunaff" << endl;
      }
      else if( param->attrs()->has_attr("contrast") )
      {
        tt = trait_parameter::contrast;

        //cout << "contrast" << endl;
      }

      my_parameters.add_trait(t, tt, cp);
    }
    else
      errors << priority(error) << "Trait '" << tname << "' not found." << endl;
  }
  else
    errors << priority(error)
           << "Trait parameter is missing name.  Skipping..." << endl;
}

void lodpal_parser::
     parse_subset(const LSFBase* param)
{
  AttrVal a=attr_value(param,0);
  if( a.has_value() )
  {
    string tname = strip_ws(a.String());
    size_t s = pairs.trait_find(tname);
    if(s < pairs.trait_count())
    {
      if(pairs.fped_info().trait_info(s).type() != RPED::RefTraitInfo::binary_trait)
        errors << priority(warning) 
               << "Trait subset '" << tname << "' is not a binary trait!  For this "
                  "analysis it will be dichotomized at 0 (trait values <= "
                  "0 will not be included into subset and values > 0 will "
                  "be included into subset)."
               << endl;
      my_parameters.add_subset(s);
    }
    else
      errors << priority(error) << "Trait subset '" << tname << "' not found." << endl;
  }
  else
    errors << priority(error)
           << "Trait subset parameter is missing name.  Skipping..." << endl;
}

void lodpal_parser::
     parse_covariate(const LSFBase* param)
{
  AttrVal a=attr_value(param,0);
  if( a.has_value() )
  {
    string cname = strip_ws(a.String());
    size_t c = pairs.trait_find(cname);
    if(c < pairs.trait_count())
    {
      double cpower = param->attrs()->FloatAttr("power");

      if( !finite(cpower) || cpower == 0.0 )
        cpower = 1.0;

      double init_d1 = 0.1;
      double init_d2 = 0.1;

#if 0
      if( has_attr(param,"initial_delta1"))
      {
        init_d1 = param->attrs()->FloatAttr("initial_delta1");

        if( !finite(init_d1) )
          init_d1 = 0.1;
      }

      if( has_attr(param,"initial_delta2"))
      {
        init_d2 = param->attrs()->FloatAttr("initial_delta2");

        if( !finite(init_d2) )
          init_d2 = 0.1;
      }
#endif

      covariate_type::cov_op op = covariate_type::sum;
      covariate_type::cov_ad ad = covariate_type::mean;
      double ad_val = std::numeric_limits<double>::quiet_NaN();

      bool set = false;
      if( param->attrs()->has_attr("mean") )
      {
        ad_val = param->attrs()->FloatAttr("mean");        
      }
      else if( param->attrs()->has_attr("minimum") )
      {
        ad     = covariate_type::minimum;
        ad_val = param->attrs()->FloatAttr("minimum");
      }
      else if( param->attrs()->has_attr("prop") )
      {
        if( pairs.fped_info().trait_info(c).type() == RPED::RefTraitInfo::continuous_trait )
        {
          ad     = covariate_type::prop;
          ad_val = param->attrs()->FloatAttr("prop");

          set_one_parameter_model(true);
          set_alpha(2.634);
        }
        else
        {
          ad     = covariate_type::minimum;
          ad_val = param->attrs()->FloatAttr("minimum");
        }
      }
        
      if( param->attrs()->has_attr("sum") )
      {
        set = true;
        op  = covariate_type::sum;
        my_parameters.add_covariate(c, ad, op, ad_val, cpower, init_d1, init_d2);
      }
      if( param->attrs()->has_attr("diff") )
      {
        set = true;
        op  = covariate_type::diff;
        my_parameters.add_covariate(c, ad, op, ad_val, cpower, init_d1, init_d2);
      }
      if( param->attrs()->has_attr("avg") )
      {
        set = true;
        op  = covariate_type::avg;
        my_parameters.add_covariate(c, ad, op, ad_val, cpower, init_d1, init_d2);
      }
      if( param->attrs()->has_attr("prod") )
      {
        set = true;
        op  = covariate_type::prod;
        my_parameters.add_covariate(c, ad, op, ad_val, cpower, init_d1, init_d2);
      }
      if( param->attrs()->has_attr("single") )
      {
        set = true;
        op  = covariate_type::single;
        my_parameters.add_covariate(c, ad, op, ad_val, cpower, init_d1, init_d2);
      }
      if( param->attrs()->has_attr("both") )
      {
        set = true;
        op  = covariate_type::sum;
        my_parameters.add_covariate(c, ad, op, ad_val, cpower, init_d1, init_d2);
        op  = covariate_type::diff;
        my_parameters.add_covariate(c, ad, op, ad_val, cpower, init_d1, init_d2);
      }

      if(!set)
        my_parameters.add_covariate(c, ad, op, ad_val, cpower, init_d1, init_d2);
    }
    else
      errors << priority(error) << "Covariate '" << cname << "' not found." << endl;
  }
  else
    errors << priority(error)
           << "Covariate parameter is missing name.  Skipping..." << endl;
}

void lodpal_parser::
     parse_weight(const LSFBase* param)
{
  AttrVal a=attr_value(param,0);
  if( a.has_value() )
  {
    string wname = strip_ws(a.String());
    size_t w = pairs.trait_find(wname);
    if(w < pairs.trait_count())
    {
      lodpal_weight_type::weight_op op = lodpal_weight_type::single;

      bool set = my_parameters.add_weight(w, op);

      if(!set)
        errors << priority(error) << "Only one weight is allowed.  Skipping..." << endl;
    }
    else
      errors << priority(error) << "Weight '" << wname << "' not found." << endl;
  }
  else
    errors << priority(error)
           << "Weight parameter is missing name.  Skipping..." << endl;
}

void lodpal_parser::
     parse_diagnostic(const LSFBase* param)
{
  AttrVal a=attr_value(param,0);
  if( a.has_value() )
  {
    string mname = strip_ws(a.String());
    size_t m = pairs.marker_find(mname);
    if(m < pairs.marker_count())
    {
      my_parameters.set_diagnostic_marker(m);
      set_diagnostic(true);
    }
    else
      errors << priority(error) << "Marker '" << mname << "' not found." << endl;
  }
  else
    errors << priority(error)
           << "Diagnostic marker parameter is missing name.  Skipping..." << endl;
}

void lodpal_parser::
parse_x_linkage_model(const LSFBase* params)
{
  if(!params)
    return;

  if(!params->List())
    return;

  bool pair_type_chosen = false;
  bool mm_pair = false;
  bool mf_pair = false;
  bool ff_pair = false;

  bool lambda1_equal = true;
  bool lambda2_fixed = true;    
  double alpha       = 2.634;

  LSFList::const_iterator i;
  for(i=params->List()->begin(); i!=params->List()->end(); ++i)
  {
    if( !(*i) || !(*i)->name().size())
      continue;

    const LSFBase* param = *i;

    if( toUpper((*i)->name()) == "PAIR_TYPE" )
    {
      AttrVal a=attr_value(param,0);
      if( a.has_value() )
      {
        pair_type_chosen = true;

        if( toUpper(a.String()) == "M_M" || toUpper(a.String()) == "M-M" || toUpper(a.String()) == "MM")
          mm_pair = true;
        else if( toUpper(a.String()) == "M_F" || toUpper(a.String()) == "M-F" || toUpper(a.String()) == "MF")
          mf_pair = true;
        else if( toUpper(a.String()) == "F_F" || toUpper(a.String()) == "F-F" || toUpper(a.String()) == "FF")
          ff_pair = true;
        else
        {
          mm_pair = true;
          mf_pair = true;
          ff_pair = true;
        }
      }
    }
    else if( toUpper((*i)->name()) == "LAMBDA1_EQUAL" )
    {
      AttrVal a=attr_value(param,0);
      if( a.has_value() )
      {
        if( toUpper(a.String()) == "FALSE" || toUpper(a.String()) == "NO" )
          lambda1_equal = false;
      }
    }
    else if( toUpper((*i)->name()) == "LAMBDA2_FIXED" )
    {
      AttrVal a=attr_value(param,0);
      if( a.has_value() )
      {
        if( toUpper(a.String()) == "FALSE" || toUpper(a.String()) == "NO" )
          lambda2_fixed = false;
      }
      else
      {
        if( param->attrs()->has_attr("alpha") )
        {
          alpha = param->attrs()->FloatAttr("alpha");
          if( !finite(alpha) || alpha < 1 )
            alpha = 2.634;
        }
      }
    }
  }

  if( !pair_type_chosen )
  {
    mm_pair = true;
    mf_pair = true;
    ff_pair = true;
  }

  my_parameters.add_x_linkage_model(lambda1_equal, lambda2_fixed, alpha);
  my_parameters.set_x_linkage_pair_type(mm_pair, mf_pair, ff_pair);
}

void lodpal_parser::
parse_autosomal_model(const LSFBase* params)
{
  if(!params)
    return;

  if(!params->List())
    return;

  my_autosomal_model_exist = true;

  autosomal_model_type::model_type      model   = autosomal_model_type::one_parameter;
  autosomal_model_type::constraint_type cont    = autosomal_model_type::constrained;
  autosomal_model_type::poo_fixed_type  poo_fix = autosomal_model_type::none;

  bool poo_test = false;
  double alpha  = 2.634;

  LSFList::const_iterator i;
  for(i=params->List()->begin(); i!=params->List()->end(); ++i)
  {
    if( !(*i) || !(*i)->name().size())
      continue;

    const LSFBase* param = *i;

    if( toUpper(param->name()) == "MODEL" )
    {
      AttrVal a=attr_value(param,0);
      if( a.has_value() )
      {
        if( toUpper(a.String()) == "TWO" || toUpper(a.String()) == "TWO_PARAMETER"
            || toUpper(a.String()) == "TWO_PARAMETER_model" )
          model = autosomal_model_type::two_parameter;
      }

      if( param->attrs()->has_attr("uncon") || param->attrs()->has_attr("unconstrained") )
      {
        cont = autosomal_model_type::unconstrained;
      }

      if( param->attrs()->has_attr("alpha") )
      {
        alpha = param->attrs()->FloatAttr("alpha");
        if( !finite(alpha) || alpha < 1 )
          alpha = 2.634;
      }
    }
    else if( toUpper(param->name()) == "PARENT_OF_ORIGIN" )
    {
      AttrVal a=attr_value(param,0);
      if( a.has_value() )
      {
        if( toUpper(a.String()) == "TRUE" || toUpper(a.String()) == "YES" )
        {
          poo_test = true;

          if(     param->attrs()->has_attr("all_pairs")
              && !my_parameters.sib_pairs_only() )
            my_parameters.set_sib_pairs_only(false);
          else
            my_parameters.set_sib_pairs_only(true);
        }
      }

      if( param->attrs()->has_attr("fixed") )
      {
        string fix = param->attrs()->StringAttr("fixed");

        if( toUpper(fix) == "MATERNAL" || toUpper(fix) == "FEMALE" )
          poo_fix = autosomal_model_type::maternal;
        else if( toUpper(fix) == "PATERNAL" || toUpper(fix) == "MALE" )
          poo_fix = autosomal_model_type::paternal;
      }
    }
  }

  my_parameters.add_autosomal_model(model, cont, poo_fix, poo_test, alpha);
}

} // end of namespace LODPAL
} // end of namespace SAGE
