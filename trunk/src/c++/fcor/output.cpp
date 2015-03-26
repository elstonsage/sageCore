//****************************************************************************
//* File:      output.cpp                                                    *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0                                                   *
//*                                                                          *
//* Notes:     This source file defines fuctions to view fcor results.       *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "fcor/output.h"

using namespace std;

namespace SAGE {
namespace FCOR {

FcorView::FcorView(const PairSetData& pairs)
{
  my_pairsetdata = &pairs;
}

void
FcorView::print_analysis_results(ostream& o, output_type otype,
                                 const pairset_result_vector& sub_results,
                                 const pairset_result_vector& main_results,
                                 const Htest_result_vector&   H_results) const
{
  if( !o )
    return;

  print_analysis_header(o, get_analysis_options().pairset, otype);

  const pairset_info_vector& sinfos = my_pairsetdata->get_subtype_info();
  const pairset_info_vector& minfos = my_pairsetdata->get_maintype_info();

  const main_to_sub_map& g = my_pairsetdata->get_pairset_group();

  main_to_sub_type_const_iterator rel_group     = g.type_begin();
  main_to_sub_type_const_iterator rel_group_end = g.type_end();

  size_t generation_limit = get_analysis_options().generation_limit;

  for( size_t r = 0; rel_group != rel_group_end; ++rel_group, ++r )
  {
    if( is_invalid_pair_type(minfos[r], generation_limit, trait_count()) )
      continue;

    //print_main_type_title(o, rel_group);

    if( get_analysis_options().pairset != SUBTYPES )
    {
      print_main_type_title(o, rel_group);
      print_pooled_types(o, rel_group, sinfos, minfos[r].total_pair_count);

      switch( otype )
      {
        case SUMMARY    : print_cor_std_err_result(o, minfos[r], main_results[r], otype);
                          break;

        case DETAILED   : print_cor_std_err_result(o, minfos[r], main_results[r], otype);
                          break;

        case PAIR_COUNT : print_pair_count(o, main_results[r], minfos[r].type == SELF, minfos[r].type == INTRA);
                          break;

        case XLS_FORMAT : print_xls_format(o, main_results[r], minfos[r].type == SELF, minfos[r].type == INTRA);
                          break;

        default         : break;
      }
    }

    size_t subtype_count = rel_group->second.size();

    if( get_analysis_options().pairset != MAINTYPES )
    {
      if( get_analysis_options().pairset == BOTH && subtype_count < 2 )
        continue;

      for( size_t s = 0; s < subtype_count; ++s )
      {
        size_t j = rel_group->second[s];

        if( is_invalid_pair_type(sinfos[j], generation_limit, trait_count()) )
          continue;

        switch( otype )
        {
          case SUMMARY    : print_table_title(o, sinfos[j]);
                            print_cor_std_err_result(o, sinfos[j], sub_results[j], otype);
                            break;

          case DETAILED   : print_table_title(o, sinfos[j]);
                            print_cor_std_err_result(o, sinfos[j], sub_results[j], otype);
                            break;

          case PAIR_COUNT : print_table_title(o, sinfos[j]);
                            print_pair_count(o, sub_results[j], sinfos[j].type == SELF, sinfos[j].type == INTRA);
                            break;

          case XLS_FORMAT : print_table_title(o, sinfos[j], true);
                            print_xls_format(o, sub_results[j], sinfos[j].type == SELF, sinfos[j].type == INTRA);
                            break;

          default         : break;
        }
      }
    }

    if(    (otype == SUMMARY || otype == DETAILED)
        && get_analysis_options().homogeneity_test
        && subtype_count > 1
        && minfos[r].type != SELF )
      print_homogeneity_test_result(o, subtype_count, minfos[r], H_results[r]);

    o << endl;

    if( get_analysis_options().pairset != SUBTYPES && minfos[r].type != MATE )
      print_ddash_outline(o);
  }

  return;
}

void
FcorView::print_var_cov_results(ostream& o, const var_cov_result_vector& vc_results) const
{
  if( !o )
    return;

  const vector<var_cov_param>& vc = get_analysis_options().var_covs;

  for( size_t c = 0; c < vc_results.size(); ++c )
  {
    o << endl;
    o << "===========================================================" << endl;
    o << "  " << vc[c].name() << flush
        << endl;
    o << "  **** : Value is not estimable." << endl;
    if( !get_analysis_options().conservative )
      o << "  [  ] : Calculated by setting at least one nonestimable" << endl
          << "         required correlation to a value of 0." << endl;
    o << "===========================================================" << endl;

    o.setf(ios_base::fixed,ios_base::floatfield);
    o.precision(8);

    if( vc[c].traits.size() == 1 || vc[c].m_type == var_cov_param::SINGLE )
    {
      for( size_t i = 0; i < vc[c].traits.size(); ++i )
      {
        size_t t = 0;
        for( size_t anal_t = 0; anal_t < trait_count(); ++anal_t )
        {
          if( get_trait_name(anal_t) == vc[c].traits[i].first )
          {
            t = anal_t;
            break;
          } 
        }
        t = t * trait_count() + t;

        print_single_legend(o, c, i);
        print_single_vc_matrix(o, vc_results[c], t, false);

        o << "  The Smallest Number of Pairs Used in" << endl
          << "    Calculating Required Correlations" << endl;

        print_single_vc_matrix(o, vc_results[c], t, true);
      }
    }
    else
    {
      for( size_t i = 0; i < vc[c].traits.size(); ++i )
      {
        size_t t1 = 0;
        for( size_t anal_t = 0; anal_t < trait_count(); ++anal_t )
        {
          if( get_trait_name(anal_t) == vc[c].traits[i].first )
          {
            t1 = anal_t;
            break;
          } 
        }

        for( size_t j = i+1; j < vc[c].traits.size(); ++j )
        {

          size_t t2 = 0;
          for( size_t anal_t = 0; anal_t < trait_count(); ++anal_t )
          {
            if( get_trait_name(anal_t) == vc[c].traits[j].first )
            {
              t2 = anal_t;
              break;
            } 
          }

          print_joint_legend(o, c, i, j);
          print_joint_vc_matrix(o, vc_results[c], t1, t2, false);

          o << "  The Smallest Number of Pairs Used in" << endl
            << "    Calculating Required Correlations" << endl;

          print_joint_vc_matrix(o, vc_results[c], t1, t2, true);
        }
      }
    }
  }

  return;
}
//
//-----------------------------------------------------------------------------------
//

void
FcorView::print_analysis_header(ostream& o, pair_type ptype, output_type otype) const
{
  if( ! o )
    return;

  o << endl;
  print_ddash_outline(o);

  // print header title
  //
  if( otype == PAIR_COUNT )
    o << "  Tables of the Smallest Number of Pairs Used in" << endl
      << "    Calculating Required Correlations for ";
  else
  {
    o << "  Tables of Correlations ";
    if( get_analysis_options().standard_error )
      o << "+/- Asymptotic Standard Errors ";
  }

  if( ptype == SUBTYPES )
    o << "for Subtypes" << endl;
  else if( ptype == MAINTYPES )
    o << "for Maintypes" << endl;
  else
    o << "for Sub/Maintypes" << endl;

  print_ddash_outline(o);

  // print ped info
  //
  o << endl;
  o << "    Number of pedigrees = " << my_pairsetdata->get_parser()->get_pedigree_count() << endl;
  o << "    Number of traits    = " << trait_count() << endl;

  // print trait info
  //
  o << endl;
  o << "    Trait(s) :";
  for( size_t t = 0; t < trait_count(); ++t )
    o << " " << get_trait_name(t);
  o << endl;

  // print weight info
  //
  weight_type w = get_analysis_options().class_weight;

  if( w != WEIGHT_COUNT )
  {
    o << endl;
    o << "    Weight method used :";

    switch( w )
    {
      case PAIR_WISE : o << " PAIR_WISE " << endl;  break;
      case UNIFORM   : o << " UNIFORM " << endl;    break;
      case MEAN      : o << " AVERAGE of PAIR_WISE & UNIFORM " << endl;  break;
      default        : break;
    }
  }

  // print legend info
  //
  if( otype == PAIR_COUNT )
  {
    o << endl;
    o << "    [      ] : Excluded the number of pairs for nonestimable required" << endl
      << "               correlations." << endl << endl;
  }
  else
  {
    o << endl;
    o << "    Legend :" << endl;
    o << "      ------ : Value is not estimable." << endl;
    o << "      &&&&&& : Pair count is greater than or equal to 100000." << endl;
    
    if( get_analysis_options().standard_error )
    {
      o << "      @@@@@@ : Standard error is greater than or equal to 10.0." << endl;
      o << "      ###### : Equivalent pair count is greater than or equal to 10000." << endl;

      if( !get_analysis_options().conservative )
        o << "     [StdErr]: Calculated by setting at least one nonestimable required" << endl
          << "               correlation to a value of 0." << endl;
    }

    if( !get_analysis_options().gender_name )
    {
      o << "      [M]    : Male" << endl;
      o << "      [F]    : Female" << endl;
    }
    o << endl;
  }

  print_ddash_outline(o);

  return;
}

void
FcorView::print_main_type_title(ostream& o, main_to_sub_type_const_iterator rel_group) const
{
  if( ! o )
    return;

  o << "Main Relationship Type : " << long_relationship_name(rel_group->first) << endl; 
  print_ddash_outline(o);

  return;
}

void
FcorView::print_pooled_types(ostream&                        o,
                             main_to_sub_type_const_iterator rel_group,
                             const pairset_info_vector&      sinfos,
                             size_t                          total_size) const
{
  if( ! o )
    return;


  for( size_t s = 0; s < rel_group->second.size(); ++s )
  {
    size_t j = rel_group->second[s];

    if( sinfos[j].type != MATE )
    {
      if( s == 0 )
      {
        o << "Subtypes Pooled   : "
          << sinfos[j].gname << endl;
      }
      else
        o << "                    " << sinfos[j].gname << endl;
    }
  }

  o << "Total Pairs Found = " << total_size << endl;

  return;
}

void
FcorView::print_table_title(ostream& o, const pairset_info& pinfo, bool xls) const
{
  if( ! o )
    return;

  o << endl << endl;

  if( pinfo.type != SELF && !xls )
  {
    size_t m_pos = pinfo.gname.find_last_of(':', pinfo.gname.size());

    string name1 = pinfo.gname.substr(0, m_pos);
    string name2 = pinfo.gname.substr(m_pos+1, pinfo.gname.size());

    o << "Relationship Type : " << name1 << "(Row):" << name2 << "(Column)" << endl;
  }
  else
    o << "Relationship Type : " << pinfo.gname << endl;

    o << "      Pairs Found = " << pinfo.total_pair_count << endl;

  return;
}

//------------------------------
// Correlation/StdErr output
//------------------------------
void
FcorView::print_cor_std_err_result(ostream&              o,
                                   const pairset_info&   pinfo,
                                   const pairset_result& result,
                                   output_type           otype) const
{
  bool is_self  = pinfo.type == SELF;
  bool is_intra = pinfo.type == INTRA;
  bool is_mate  = pinfo.type == MATE;

  print_cor_std_err_heading(o, is_self, is_intra);

  weight_type wt = get_analysis_options().class_weight;

  if( otype == SUMMARY || wt != WEIGHT_COUNT )
    print_cor_std_err(o, wt, result, is_self, is_intra);
  else
  {
    for( size_t w = 0; w <= WEIGHT_COUNT; ++w )
    {
      switch( (weight_type)w )
      {
        case PAIR_WISE    : o << " PAIR_WISE WEIGHT" << endl
                              << " ----------------" << endl;
                            print_cor_std_err(o, w, result, is_self, is_intra);
                            break;

        case UNIFORM      : o << " UNIFORM WEIGHT" << endl
                              << " --------------" << endl;
                            print_cor_std_err(o, w, result, is_self, is_intra);
                            break;

        case MEAN         : break;

        case WEIGHT_COUNT : o << " OPTIMAL WEIGHT" << endl
                              << " --------------" << endl;
                            print_cor_std_err(o, w, result, is_self, is_intra);
                            break;

        default           : break;
      }
    }
  }

  if( trait_count() > 1 && !(is_self  || is_intra || is_mate) )
    print_pooled_corr(o, result);

  return;
}

void
FcorView::print_cor_std_err_heading(ostream& o, bool selfclass, bool intraclass) const
{
  if( ! o )
    return;

  size_t c_width = 21;
  size_t t_width = 19;
  if( get_analysis_options().standard_error )
  {
    c_width = 21+11;
    t_width = 19+11;
  }

  print_row_outline(o, c_width);

  o << "             ";
  print_column_trait_names(o, t_width);
  print_class_type(o, selfclass, intraclass);
  print_column_outline(o, c_width);  

  o << "             ";
  for( size_t t = 0; t < trait_count(); ++t )
  {
    o << "  Count  Correlation ";
    if( get_analysis_options().standard_error )
      o << " P-value   ";
  }
  o << endl;

  if( get_analysis_options().standard_error )
  {
    o << "             ";
    for( size_t t = 0; t < trait_count(); ++t )
      o << " EqvCnt   +/- StdErr            ";
    o << endl;
  }

  print_row_outline(o, c_width);

  return;
}

void
FcorView::print_cor_std_err(ostream& o, size_t w,
                            const pairset_result& cor_std,
                            bool selfclass, bool intraclass) const
{
  size_t c_width = 21;
  if( get_analysis_options().standard_error )
    c_width += 11;

  for( size_t t1 = 0; t1 < trait_count(); ++t1 )
  {
    print_trait_name(o, t1, 11);
    
    size_t t2 = 0;

    if( intraclass || selfclass )
    {
      print_partial_outline_dash(o, t1);
      t2 = t1;
    }

    for( ; t2 < trait_count(); ++t2 )
    {
      size_t count = cor_std.corr(t1, t2).count;
      if( intraclass )
        count /= 2;
      double correlation = cor_std.corr(t1, t2).correlation[w];

      print_correlation(o, count, 5, "        ", correlation, 7, 4);

      // print p-value.
      //
      if( get_analysis_options().standard_error )
      {
        double standard_error = sqrt(cor_std.std_err(t1, t2).standard_error[w]);
        double eff            = eff_count(correlation, standard_error);

        double z   = 0.5 * (log((1.0 + correlation) / (1.0 - correlation)));
        double var = 1.0 / (eff - 3.0);

        if( (intraclass || selfclass) && (t1 == t2) )
          var = 1.0 / (eff - 1.5);

        if( selfclass && (t1 == t2) )
          o << "  ------    ";
        else
          print_pvalue(o, z, var, 12, 4);
      }
    }
          
    if( get_analysis_options().standard_error )
    {
      o << endl << "             ";

      t2 = 0;

      if( selfclass )
      {
        print_partial_outline_dash(o, t1);

        print_standard_error(o,        std::numeric_limits<double>::quiet_NaN(), 6, 1,
                             "   +/-", std::numeric_limits<double>::quiet_NaN(), 7, 4, false);
        o << "           ";

        t2 = t1 + 1;
      }
      else if( intraclass )
      {
        print_partial_outline_dash(o, t1);
        t2 = t1;
      }

      for( ; t2 < trait_count(); ++t2 )
      {
        bool   replaced       = cor_std.std_err(t1, t2).replaced;
        double standard_error = sqrt(cor_std.std_err(t1, t2).standard_error[w]);
        double eff            = eff_count(cor_std.corr(t1, t2).correlation[w], standard_error);

        print_standard_error(o, eff, 7, 1, "   +/-", standard_error, 6, 4, replaced);

        o << "           ";
      }
    }
    o << endl;
    print_row_outline(o, c_width);
  }

  return;
}

void
FcorView::print_pooled_corr(ostream& o, const pairset_result& cor_std) const
{
  size_t c_width = 21;
  if( get_analysis_options().standard_error )
    c_width += 11;

  o << " Pooled Cross-Correlations" << endl
    << " -------------------------" << endl;

  for( size_t t1 = 0; t1 < trait_count()-1; ++t1 )
  {
    print_trait_name(o, t1, 11);
    
    print_partial_outline_dash(o, t1+1);

    vector< pair<double, double> > std_err_eff;
    for( size_t t2 = t1+1; t2 < trait_count(); ++t2 )
    {
      size_t count       = cor_std.corr(t2, t1).count;
      double correlation = cor_std.pooled_cross_corr(t2, t1).correlation;

      print_correlation(o, count, 5, "        ", correlation, 7, 4);

      // print p-value.
      //
      double standard_error = sqrt(cor_std.pooled_cross_corr(t2, t1).variance);
      double eff            = eff_count(correlation, standard_error);

      std_err_eff.push_back(make_pair(standard_error, eff));

      double z   = 0.5 * (log((1.0 + correlation) / (1.0 - correlation)));
      double var = 1.0 / (eff - 3.0);

      print_pvalue(o, z, var, 12, 4);
    }
          
    // print standard_error
    //
    o << endl << "             ";

    print_partial_outline_dash(o, t1+1);

    for( size_t t2 = t1+1, t = 0; t < std_err_eff.size(); ++t2, ++t )
    {
      bool   replaced       = cor_std.std_err(t2, t1).replaced;
      double standard_error = std_err_eff[t].first;
      double eff            = std_err_eff[t].second;

      print_standard_error(o, eff, 7, 1, "   +/-", standard_error, 6, 4, replaced);

      o << "           ";
    }
    o << endl;
    print_row_outline(o, c_width);
  }

  size_t df   = (trait_count()*(trait_count()-1)) / 2;
  double chsq = cor_std.pooled_cross_corr_chi_square;
  double pval = cor_std.pooled_cross_corr_p_value;

  o << " Test for Homogeneity of Cross-Correlations";

  if( !isnan(chsq) && !isnan(pval) )
  {
    o << endl
      << " ------------------------------------------" << endl
      << "  Chi-Square = " << fp(chsq, 8, 6) << " with " << df << " degree(s) of freedom" << endl
      << "  P-Value    = " << fp(pval, 8, 6) << endl;
  }
  else
    o << " is not available!" << endl;

  print_row_outline(o, c_width);

  return;
}

void
FcorView::print_homogeneity_test_result(ostream& o, size_t subtype_count,
                                        const pairset_info& pinfo,
                                        const Htest_result& H_result) const
{
  if( ! o )
    return;

  size_t df = (subtype_count - 1) * trait_count() * trait_count();

  if( pinfo.type == INTRA || pinfo.type == SELF )
    df = ((subtype_count - 1) * (trait_count()+1) * trait_count()) / 2;

  if( get_analysis_options().individual_homog )
    df = subtype_count - 1;

  o << endl;
  print_sdash_outline(o);
  o << "Test for Homogeneity of Correlations among Subtypes - ";
  if( trait_count() == 1 || get_analysis_options().individual_homog )
    o << "Single Trait";
  else
    o << "All Traits";
  o << endl;

  //if( b->second.chi_square == -1.0 )
  //  o << "Matrix is ill-conditioned!!!" << endl;
  if( !H_result.chi_square.size() )
  {
    o << " - is not available!" << endl;
    print_sdash_outline(o);
    return;
  }

  print_sdash_outline(o);

  for( size_t t = 0; t < H_result.chi_square.size(); ++t )
  {
    if( trait_count() > 1 && get_analysis_options().individual_homog )
    {
      o << get_trait_name(t) << endl;
      o << "  ";
    }

    o << "  Chi-Square = ";
    if( SAGE::isnan(H_result.chi_square[t]) )
      o << "------" << endl;
    else
    {
      if( H_result.replaced )
        o << "[" << fp(H_result.chi_square[t], 8, 6) << "]";
      else
        o << fp(H_result.chi_square[t], 8, 6);
      o << " with " << df << " degree(s) of freedom" << endl;
    }

    if( trait_count() > 1 && get_analysis_options().individual_homog )
      o << "  ";

    o << "  P-Value    = ";
    if( SAGE::isnan(H_result.p_value[t]) )
      o << "------" << endl;
    else
    {
      if( H_result.replaced )
        o << "[" << fp(H_result.p_value[t], 8, 6) << "]" << endl;
      else
        o << fp(H_result.p_value[t], 8, 6) << endl;
    }

    o << endl;
  }

  return;
}

//------------------------------
// Least Pair Count output
//------------------------------

void
FcorView::print_pair_count(ostream& o,
                           const pairset_result& cor_std,
                           bool selfclass, bool intraclass) const
{
  if( ! o )
    return;

  print_row_outline(o, 12);
  print_class_type(o, selfclass, intraclass);
  print_column_trait_names(o, 10);
  print_row_outline(o, 12);
  
  for( size_t t1 = 0; t1 < trait_count(); ++t1 )
  {
    print_trait_name(o, t1, 11);
    
    size_t t2 = 0;

    if( selfclass )
    {
      size_t local_t1 = t1;

      for( size_t i = 0; i < local_t1; ++i )
        o << "            ";
        //o << " ---------- ";

      o << "         *  ";

      t2 = t1 + 1;
    }
    else if( intraclass )
    {
      size_t local_t1 = t1;

      for( size_t i = 0; i < local_t1; ++i )
        o << "            ";
        //o << " ---------- ";

      t2 = t1;
    }

    for( ; t2 < trait_count(); ++t2 )
    {
      size_t pair_count = cor_std.std_err(t1, t2).least_pair_count;

      if( cor_std.std_err(t1, t2).replaced )
        o << " [" << right << setw(8) << pair_count << "] ";
      else
        o << "  " << right << setw(8) << pair_count << "  ";
    }
          
    o << endl;
    print_row_outline(o, 12);
  }
  o << endl;

  return;
}

//------------------------------
// Alternative output
//------------------------------

void
FcorView::print_xls_format(ostream& o,
                           const pairset_result& cor_std,
                           bool selfclass, bool intraclass) const
{
  if( ! o )
    return;

  size_t width = 0;

  for( size_t t = 0; t < trait_count(); ++t )
    width = max(width, get_trait_name(t).size());

  size_t w1 = 6+13;
  size_t w2 = 9+11+10+4;

  // line 1
  print_xls_row_outline(o, width, w1, w2);

  // line 2
  for( size_t t = 0; t < width*2+5; ++t )
    o << " ";
  o << " Count  Correlation";
  if( get_analysis_options().standard_error )
    o << "   EqvCnt  StdError   P-value     ";
  o << endl;

  // line 3 == 1
  print_xls_row_outline(o, width, w1, w2);

  // line 4...
  for( size_t t1 = 0; t1 < trait_count(); ++t1 )
  {
    size_t t2 = 0;

    if( intraclass || selfclass )
      t2 = t1;

    for( ; t2 < trait_count(); ++t2 )
    {
      if( selfclass && t1 == t2 )
        continue;

      print_trait_name(o, t1, width);
      o << "-";
      print_trait_name(o, t2, width);

      size_t count    = cor_std.corr(t1, t2).count;
      if( intraclass )
        count /= 2;

      size_t w = get_analysis_options().class_weight;;
      double correlation = cor_std.corr(t1, t2).correlation[w];

      print_correlation(o, count, 6, "  ", correlation, 11, 7);

      if( get_analysis_options().standard_error )
      {
        bool   replaced       = cor_std.std_err(t1, t2).replaced;
        double standard_error = sqrt(cor_std.std_err(t1, t2).standard_error[w]);
        double eff            = eff_count(correlation, standard_error);

        print_standard_error(o, eff, 9, 1, " ", standard_error, 8, 6, replaced);

        // print p-value.
        //
        double z   = 0.5 * (log((1.0 + correlation) / (1.0 - correlation)));
        double var = 1.0 / (eff - 3.0);
        if( (intraclass || selfclass) && (t1 == t2) )
          var = 1.0 / (eff - 1.5);

        print_pvalue(o, z, var, 13, 6);
      }
      o << endl;
    }
  }

  // last line == 1
  print_xls_row_outline(o, width, w1, w2);

  return;
}

//------------------------------
// Misc.
//------------------------------

void
FcorView::print_ddash_outline(ostream& o) const
{
  if( ! o )
    return;

  o << "============================================================================" << endl;

  return;
}

void
FcorView::print_sdash_outline(ostream& o) const
{
  if( ! o )
    return;

  o << "----------------------------------------------------------------------------" << endl;

  return;
}

void
FcorView::print_row_outline(ostream& o, size_t width) const
{
  if( ! o )
    return;

  o << "-------------";

  print_column_outline(o, width);
}

void
FcorView::print_column_outline(ostream& o, size_t width) const
{
  if( ! o )
    return;

  for( size_t t = 0; t < trait_count(); ++t )
    o << string(width, '-');
  o << endl;

  return;
}

void
FcorView::print_partial_outline_dash(ostream& o, size_t t1) const
{
  if( ! o )
    return;

  if( t1 >= trait_count() )
    t1 -= trait_count();

  for( size_t i = 0; i < t1; ++i )
  {
    o << "                     ";
    if( get_analysis_options().standard_error )
      o << "           ";
  }

  return;
}

void
FcorView::print_column_trait_names(ostream& o, size_t width) const
{
  if( ! o )
    return;

  for( size_t t = 0; t < trait_count(); ++t )
  {
    print_trait_name(o, t, width);
  }
  o << endl;

  return;
}

void
FcorView::print_trait_name(ostream& o, size_t t, size_t width) const
{
  if( ! o )
    return;

  // by default
  //  column width = 19+11(p-value)
  //  row width    = 11
  //
  if( get_trait_name(t).size() >= width )
  {
    string temp = get_trait_name(t).substr(0, width);
    o << " " << left << setw(width) << temp << " ";
  }
  else
  {
    size_t empty = size_t(floor( (width - get_trait_name(t).size()) / 2.));
    string temp(empty, ' ');
    temp.append(get_trait_name(t));
    o << " " << left << setw(width) << temp << " ";
  }

  return;
}

void
FcorView::print_class_type(ostream& o, bool selfclass, bool intraclass) const
{
  if( ! o )
    return;

  if( selfclass )
    o << "   POOLED    ";
  else if( intraclass )
    o << " INTRACLASS  ";
  else
    o << " INTERCLASS  ";

  return;
}

void
FcorView::print_xls_row_outline(ostream& o, size_t width, size_t w1, size_t w2) const
{
  if( ! o )
    return;

  o << string((width*2+5) + w1, '-');
  if( get_analysis_options().standard_error )
    o << string(w2, '-');
  o << endl;

  return;
}

void
FcorView::print_correlation(ostream& o,  size_t count,       size_t w1,
                            string   sp, double correlation, size_t w2, size_t f2) const
{
  if( ! o )
    return;

  if( count >= 100000 )
    o << " " << string(w1-1, '&');
  else
    o << right << setw(w1) << count;

  o << sp << fp(correlation, w2, f2);

  return;
}

void
FcorView::print_standard_error(ostream& o,  double eff,     size_t w1, size_t f1,
                               string   sp, double std_err, size_t w2, size_t f2, bool replaced) const
{
  if( ! o )
    return;

  if( eff < 10000 )
    o << fp(eff, w1, f1);
  else  if( SAGE::isnan(eff) )
    o << " " << string(w1-1, '-');
  else
    o << " " << string(w1-1, '#');

  o << sp;

  if( std_err < 10.0 )
  {
    if( replaced )
      o << "[";
    else
      o << " ";

    o << fp(std_err, w2, f2);

    if( replaced )
      o << "]";
    else
      o << " ";
  }
  else if( SAGE::isnan(std_err) )
    o << " " << string(w2, '-') << " ";
  else
    o << " " << string(w2, '@') << " ";

  return;
}

void
FcorView::print_pvalue(ostream& o, double z, double var, size_t w, size_t f) const
{
  if( ! o )
    return;

  if( !SAGE::isnan(z) && !SAGE::isnan(var) )
  {
    double p_value = (1. - normal_cdf(0.0, sqrt(var), fabs(z))) * 2.;
    o << pval(p_value, w-1, f, 3);
  }
  else
    o << string(w, ' ');

  return;
}

//------------------------------
// Variance-Covariance matrix
//------------------------------

void
FcorView::print_joint_vc_matrix(ostream& o,
                                const var_cov_matrix& a_vc_result,
                                size_t t1, size_t t2, bool lp) const
{
  if( ! o )
    return;

  o << "  " << string(8+13*4, '-') << endl;

  o << "     \\    ";
  for( size_t t = 0; t < 4; ++t )
    o << "   [Col" << t + 1 << "]    ";
  o << "" << endl;

  o << "  " << string(8+13*4, '-') << endl;

  for( size_t anal_t1 = 0, R = 1; anal_t1 < a_vc_result.rows(); ++anal_t1 )
  {
    size_t divR = anal_t1 / trait_count();
    size_t modR = anal_t1 % trait_count();

    if(    (divR == t1 || divR == t2)
        && (modR == t1 || modR == t2) )
    {
      o << "  [Row" << R << "]  ";
      for( size_t anal_t2 = 0; anal_t2 < a_vc_result.cols(); ++anal_t2 )
      {              
        size_t divC = anal_t2 / trait_count();
        size_t modC = anal_t2 % trait_count();

        if(    (divC == t1 || divC == t2)
            && (modC == t1 || modC == t2) )
        {
          double num_val = a_vc_result(anal_t1, anal_t2).var_cov;
          string str_val = fp(num_val, 10, 7, '*');

          if( lp )
          {
            num_val = (double)a_vc_result(anal_t1, anal_t2).least_pair_count;
            str_val = fp(num_val, 11, 0, '*').substr(0, 10);
          }

          if( !SAGE::isnan(num_val) )
          {
            if( a_vc_result(anal_t1, anal_t2).replaced )
              o << "[" << str_val << "] ";
            else
              o << " " << str_val << "  ";
          }
          else
            o << "************ ";
        }
      }
      o << "" << endl;
      ++R;
    }
  }

  o << "  " << string(8+13*4, '-') << endl;
  o << endl;

  return;
}

void
FcorView::print_single_vc_matrix(ostream& o,
                                 const var_cov_matrix& a_vc_result,
                                 size_t t, bool lp) const
{
  if( ! o )
    return;

  double num_val = a_vc_result(t, t).var_cov;
  string str_val = fp(num_val, 10, 7, '*');

  if( lp )
  {
    num_val = (double)a_vc_result(t, t).least_pair_count;
    str_val = fp(num_val, 11, 0, '*').substr(0, 10);
  }

  o << "  " << string(8+13, '-') << endl;
  o << "     \\       [Col1]    " << endl;
  o << "  " << string(8+13, '-') << endl;

  o << "  [Row1]  ";
  if( !SAGE::isnan(num_val) )
  {
    if( a_vc_result(t, t).replaced )
      o << "[" << str_val << "] ";
    else
      o << " " << str_val << "  ";
  }
  else
    o << "************ ";
  o << "" << endl;

  o << "  " << string(8+13, '-') << endl;
  o << endl;

  return;
}

void
FcorView::print_legend(ostream& o, string pos, string cor, string trait1, string trait2) const
{
  if( ! o )
    return;

  string rel1 = cor;
  string rel2 = cor;

  size_t pos1  = cor.find(':');
  if( pos1 < cor.size() )
  {
    rel1 = cor.substr(0, pos1);
    rel2 = cor.substr(pos1+1, cor.size());
  }

  o << "    [" << pos << "] " 
      << rel1
      << "-"
      << trait1
      << " : "
      << rel2
      << "-"
      << trait2
      << endl;

  return;
}

void
FcorView::print_joint_legend(ostream& o, size_t c, size_t i, size_t j) const
{
  if( ! o )
    return;

  const vector<var_cov_param>& vc = get_analysis_options().var_covs;

  o << endl << "  Legend : " << endl;

  print_legend(o, "Row1", vc[c].correlations.first, vc[c].traits[i].first, vc[c].traits[i].first);
  print_legend(o, "Row2", vc[c].correlations.first, vc[c].traits[i].first, vc[c].traits[j].first);
  print_legend(o, "Row3", vc[c].correlations.first, vc[c].traits[j].first, vc[c].traits[i].first);
  print_legend(o, "Row4", vc[c].correlations.first, vc[c].traits[j].first, vc[c].traits[j].first);
  print_legend(o, "Col1", vc[c].correlations.second, vc[c].traits[i].first, vc[c].traits[i].first);
  print_legend(o, "Col2", vc[c].correlations.second, vc[c].traits[i].first, vc[c].traits[j].first);
  print_legend(o, "Col3", vc[c].correlations.second, vc[c].traits[j].first, vc[c].traits[i].first);
  print_legend(o, "Col4", vc[c].correlations.second, vc[c].traits[j].first, vc[c].traits[j].first);

  return;
}

void
FcorView::print_single_legend(ostream& o, size_t c, size_t i) const
{
  if( ! o )
    return;

  const vector<var_cov_param>& vc = get_analysis_options().var_covs;

  o << endl << "  Legend : " << endl;

  print_legend(o, "Row1", vc[c].correlations.first, vc[c].traits[i].first, vc[c].traits[i].first);
  print_legend(o, "Col1", vc[c].correlations.second, vc[c].traits[i].first, vc[c].traits[i].first);

  return;
}

// end of FcorView Implementation

} // end of namespace FCOR
} // end of namespace SAGE
