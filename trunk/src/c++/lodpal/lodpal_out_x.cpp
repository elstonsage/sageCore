//****************************************************************************
//* File:      lodpal_out_x.cpp                                              *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0 Initial implementation.               yjs May. 02 *
//*                                                                          *
//* Notes:     This file implemnets lodpal_test_xfile class.                 *
//*                                                                          *
//* Copyright (c) 2002 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "lodpal/lodpal_out.h"

namespace SAGE   {
namespace LODPAL {

//---------------------------------------------------------------------------------

lodpal_test_xfile::lodpal_test_xfile(ostream& output)
                 : lodpal_test_viewer(output)
{
  my_mm_pair_result.resize(0);
  my_mf_pair_result.resize(0);
  my_ff_pair_result.resize(0);
}

lodpal_test_xfile::lodpal_test_xfile(ostream& output, bool wide_out, bool pval_sci)
                 : lodpal_test_viewer(output, wide_out, pval_sci)
{
  my_mm_pair_result.resize(0);
  my_mf_pair_result.resize(0);
  my_ff_pair_result.resize(0);
}

void
lodpal_test_xfile::print_header(const ARP_base_analysis& test)
{
  ostream& out = output_stream();

  set_format_params(test);

  print_double_line(test);

  print_title_heading(test, true);

  print_double_line(test);

  out << endl;

  print_model_summary(test);

  out << "    Model    : " << test.pairs_info().parameters().x_linkage_model().name();
  out << endl << endl;

  print_double_line(test);

  print_table_heading(test, false);

  print_single_line(test, false);
}

void
lodpal_test_xfile::print_results(const ARP_base_analysis& test)
{
  ostream& out = output_stream();

  lodpal_parameters::marker_const_iterator pi = test.parameters().marker_begin();
  for( ; pi != test.parameters().marker_end(); ++pi)
  {
    double total_lod_score   = 0.;
    size_t total_fsib_count  = 0;
    size_t total_hsib_count  = 0;
    size_t total_other_count = 0;
    size_t total_pair_count  = 0;

    lodpal_parameters::covariate_const_iterator ci = test.parameters().covariate_begin();

    if( test.parameters().use_mm_pair() )
    {
      vector<double> param_est;
      vector<double> fst_deriv;

      param_est.push_back(pi->beta1mm.value());
      fst_deriv.push_back(pi->beta1mm.first_derivative());

      if( !test.parameters().x_linkage_model().lambda2_fixed )
      {
        param_est.push_back(numeric_limits<double>::quiet_NaN());
        fst_deriv.push_back(numeric_limits<double>::quiet_NaN());
      }

      ci = test.parameters().covariate_begin();
      for( ; ci != test.parameters().covariate_end(); ++ci )
      {
        param_est.push_back(ci->delta1mm.value());
        fst_deriv.push_back(ci->delta1mm.first_derivative());
      }

      sub_type_info mm;
      mm.marker_name      = test.relative_pairs().marker_name(pi->marker);
      mm.marker_distance  = test.relative_pairs().marker_distance(pi->marker);
      mm.lod_score        = test.get_lodpal_result().lod_score_mm();
      mm.fsib_pair_count  = test.pairs_info().mm_fsib_pair_count();
      mm.hsib_pair_count  = test.pairs_info().mfm_hsib_pair_count();
      mm.other_pair_count = test.pairs_info().mm_other_pair_count();
      mm.pair_count       = test.pairs_info().mm_pair_count();
      mm.param_estimates  = param_est;
      mm.first_derivative = fst_deriv;
      mm.function_evaluations = test.get_lodpal_result().function_evaluations();
      mm.last_error           = test.get_lodpal_result().last_error();

      my_mm_pair_result.push_back(mm);

      total_lod_score   += test.get_lodpal_result().lod_score_mm();
      total_fsib_count  += test.pairs_info().mm_fsib_pair_count();
      total_hsib_count  += test.pairs_info().mfm_hsib_pair_count();
      total_other_count += test.pairs_info().mm_other_pair_count();
      total_pair_count  += test.pairs_info().mm_pair_count();
    }

    if( test.parameters().use_mf_pair() )
    {
      vector<double> param_est;
      vector<double> fst_deriv;

      param_est.push_back(pi->beta1mf.value());
      fst_deriv.push_back(pi->beta1mf.first_derivative());

      if( !test.parameters().x_linkage_model().lambda2_fixed )
      {
        param_est.push_back(numeric_limits<double>::quiet_NaN());
        fst_deriv.push_back(numeric_limits<double>::quiet_NaN());
      }

      ci = test.parameters().covariate_begin();
      for( ; ci != test.parameters().covariate_end(); ++ci )
      {
        param_est.push_back(ci->delta1mf.value());
        fst_deriv.push_back(ci->delta1mf.first_derivative());
      }

      sub_type_info mf;
      mf.marker_name      = test.relative_pairs().marker_name(pi->marker);
      mf.marker_distance  = test.relative_pairs().marker_distance(pi->marker);
      mf.lod_score        = test.get_lodpal_result().lod_score_mf();
      mf.fsib_pair_count  = test.pairs_info().mf_fsib_pair_count();
      mf.hsib_pair_count  = test.pairs_info().mff_hsib_pair_count();
      mf.other_pair_count = test.pairs_info().mf_other_pair_count();
      mf.pair_count       = test.pairs_info().mf_pair_count();
      mf.param_estimates  = param_est;
      mf.first_derivative = fst_deriv;
      mf.function_evaluations = test.get_lodpal_result().function_evaluations();
      mf.last_error           = test.get_lodpal_result().last_error();

      my_mf_pair_result.push_back(mf);

      total_lod_score   += test.get_lodpal_result().lod_score_mf();
      total_fsib_count  += test.pairs_info().mf_fsib_pair_count();
      total_hsib_count  += test.pairs_info().mff_hsib_pair_count();
      total_other_count += test.pairs_info().mf_other_pair_count();
      total_pair_count  += test.pairs_info().mf_pair_count();
    }

    if( test.parameters().use_ff_pair() )
    {
      vector<double> param_est;
      vector<double> fst_deriv;

      param_est.push_back(pi->beta1ff.value());
      fst_deriv.push_back(pi->beta1ff.first_derivative());

      if( !test.parameters().x_linkage_model().lambda2_fixed )
      {
        param_est.push_back(pi->beta2ff.value());
        fst_deriv.push_back(pi->beta2ff.first_derivative());
      }

      ci = test.parameters().covariate_begin();
      for( ; ci != test.parameters().covariate_end(); ++ci )
      {
        param_est.push_back(ci->delta1ff.value());
        fst_deriv.push_back(ci->delta1ff.first_derivative());
      }

      sub_type_info ff;
      ff.marker_name      = test.relative_pairs().marker_name(pi->marker);
      ff.marker_distance  = test.relative_pairs().marker_distance(pi->marker);
      ff.lod_score        = test.get_lodpal_result().lod_score_ff();
      ff.fsib_pair_count  = test.pairs_info().ff_fsib_pair_count();
      ff.hsib_pair_count  = test.pairs_info().fff_hsib_pair_count();
      ff.other_pair_count = test.pairs_info().ff_other_pair_count();
      ff.pair_count       = test.pairs_info().ff_pair_count();
      ff.param_estimates  = param_est;
      ff.first_derivative = fst_deriv;
      ff.function_evaluations = test.get_lodpal_result().function_evaluations();
      ff.last_error           = test.get_lodpal_result().last_error();

      my_ff_pair_result.push_back(ff);

      total_lod_score   += test.get_lodpal_result().lod_score_ff();
      total_fsib_count  += test.pairs_info().ff_fsib_pair_count();
      total_hsib_count  += test.pairs_info().fff_hsib_pair_count();
      total_other_count += test.pairs_info().ff_other_pair_count();
      total_pair_count  += test.pairs_info().ff_pair_count();
    }

    out << " ";
    out.setf(ios_base::left, ios_base::adjustfield);

    if( test.relative_pairs().valid_distance_exist() )
      out << setw(my_max_marker_name - 6)
          << test.relative_pairs().marker_name(pi->marker)
          << fp(test.relative_pairs().marker_distance(pi->marker), 6, 1) << " ";
    else
      out << setw(my_max_marker_name) << test.relative_pairs().marker_name(pi->marker) << " ";

    out << fp(total_lod_score,9,6) << " ";
    out.setf(ios_base::right, ios_base::adjustfield);
    out << setw(6)   << total_fsib_count << " ";
    if( wide_output() )
      out << setw(6) << total_hsib_count << " "
          << setw(6) << total_other_count << " ";
    
    out << setw(6)   << total_pair_count << "  ";

    if(    test.parameters().x_linkage_model().lambda1_equal
        || (   !test.parameters().x_linkage_model().lambda1_equal
             && test.parameters().used_pair_type() < 2 ) )
    {
      out << fp(pi->beta1mm.value(), 9, 6);

      if( !test.parameters().x_linkage_model().lambda2_fixed && test.parameters().use_ff_pair() )
        out << " " << fp(pi->beta2ff.value(), 9, 6);

      ci = test.parameters().covariate_begin();
      for( ; ci != test.parameters().covariate_end(); ++ci )
        out << " " << fp(ci->delta1mm.value(), 9, 6);

      if( wide_output() )
      {
        out << "  " << fp(pi->beta1mm.first_derivative(), 9, 6);

        if( !test.parameters().x_linkage_model().lambda2_fixed && test.parameters().use_ff_pair() )
          out << " " << fp(pi->beta2ff.first_derivative(), 9, 6);

        ci = test.parameters().covariate_begin();
        for( ; ci != test.parameters().covariate_end(); ++ci )
          out << " " << fp(ci->delta1mm.first_derivative(), 9, 6);
#if 0
        out << " " << setw(4) << test.get_lodpal_result().function_evaluations()
            << " " << setw(3) << test.get_lodpal_result().last_error();
#endif
      }
    }

    out << endl;
  }
}

void
lodpal_test_xfile::print_footer(const ARP_base_analysis& test)
{
  ostream& out = output_stream();

  print_double_line(test);

  out << endl << endl;

  size_t total_pair_type = test.parameters().used_pair_type();

  if( total_pair_type > 1 && !test.parameters().x_linkage_model().lambda1_equal )
  {
    if( test.parameters().use_mm_pair() )
      print_mm_result(test);
    if( test.parameters().use_mf_pair() )
      print_mf_result(test);
    if( test.parameters().use_ff_pair() )
      print_ff_result(test);
  }
}

void
lodpal_test_xfile::print_mm_result(const ARP_base_analysis& test)
{
  ostream& out = output_stream();

  print_double_line(test);

  out << "Male-Male" << endl;

  print_double_line(test);
  print_table_heading(test);
  print_single_line(test);

  for( size_t i = 0; i < my_mm_pair_result.size(); ++i )
  {
    out << " ";
    out.setf(ios_base::left, ios_base::adjustfield);

    if( test.relative_pairs().valid_distance_exist() )
      out << setw(my_max_marker_name - 6)
          << my_mm_pair_result[i].marker_name
          << fp(my_mm_pair_result[i].marker_distance, 6, 1) << " ";
    else
      out << setw(my_max_marker_name) << my_mm_pair_result[i].marker_name << " ";

    out << fp(my_mm_pair_result[i].lod_score,9,6)                   << " ";
    out.setf(ios_base::right, ios_base::adjustfield);
    out << setw(6)   << my_mm_pair_result[i].fsib_pair_count  << " ";
    if( wide_output() )
      out << setw(6) << my_mm_pair_result[i].hsib_pair_count << " "
          << setw(6) << my_mm_pair_result[i].other_pair_count << " ";
    
    out << setw(6)   << my_mm_pair_result[i].pair_count       << "  ";

    for( size_t j = 0; j < my_mm_pair_result[i].param_estimates.size(); ++j )
    {
      out << fp(my_mm_pair_result[i].param_estimates[j], 9, 6) << " ";
    }

    if( wide_output() )
    {
      for( size_t j = 0; j < my_mm_pair_result[i].first_derivative.size(); ++j )
      {
        out << fp(my_mm_pair_result[i].first_derivative[j], 9, 6) << " ";
      }
#if 0
      out << " " << setw(4) << my_mm_pair_result[i].function_evaluations
          << " " << setw(3) << my_mm_pair_result[i].last_error;
#endif
    }
    out << endl;
  }

  print_double_line(test);
}

void
lodpal_test_xfile::print_mf_result(const ARP_base_analysis& test)
{
  ostream& out = output_stream();

  print_double_line(test);

  out << "Male-Female" << endl;

  print_double_line(test);
  print_table_heading(test);
  print_single_line(test);

  for( size_t i = 0; i < my_mf_pair_result.size(); ++i )
  {
    out << " ";
    out.setf(ios_base::left, ios_base::adjustfield);

    if( test.relative_pairs().valid_distance_exist() )
      out << setw(my_max_marker_name - 6)
          << my_mf_pair_result[i].marker_name
          << fp(my_mf_pair_result[i].marker_distance, 6, 1) << " ";
    else
      out << setw(my_max_marker_name) << my_mf_pair_result[i].marker_name << " ";

    out << fp(my_mf_pair_result[i].lod_score,9,6)                   << " ";
    out.setf(ios_base::right, ios_base::adjustfield);
    out << setw(6)   << my_mf_pair_result[i].fsib_pair_count  << " ";
    if( wide_output() )
      out << setw(6) << my_mf_pair_result[i].hsib_pair_count << " "
          << setw(6) << my_mf_pair_result[i].other_pair_count << " ";
    
    out << setw(6)   << my_mf_pair_result[i].pair_count       << "  ";

    for( size_t j = 0; j < my_mf_pair_result[i].param_estimates.size(); ++j )
    {
      out << fp(my_mf_pair_result[i].param_estimates[j], 9, 6) << " ";
    }

    if( wide_output() )
    {
      for( size_t j = 0; j < my_mf_pair_result[i].first_derivative.size(); ++j )
      {
        out << fp(my_mf_pair_result[i].first_derivative[j], 9, 6) << " ";
      }
#if 0
      out << " " << setw(4) << my_mf_pair_result[i].function_evaluations
          << " " << setw(3) << my_mf_pair_result[i].last_error;
#endif
    }
    out << endl;
  }


  print_double_line(test);
}

void
lodpal_test_xfile::print_ff_result(const ARP_base_analysis& test)
{
  ostream& out = output_stream();

  print_double_line(test);

  out << "Female-Female" << endl;

  print_double_line(test);
  print_table_heading(test);
  print_single_line(test);

  for( size_t i = 0; i < my_ff_pair_result.size(); ++i )
  {
    out << " ";
    out.setf(ios_base::left, ios_base::adjustfield);

    if( test.relative_pairs().valid_distance_exist() )
      out << setw(my_max_marker_name - 6)
          << my_ff_pair_result[i].marker_name
          << fp(my_ff_pair_result[i].marker_distance, 6, 1) << " ";
    else
      out << setw(my_max_marker_name) << my_ff_pair_result[i].marker_name << " ";

    out << fp(my_ff_pair_result[i].lod_score,9,6)                   << " ";
    out.setf(ios_base::right, ios_base::adjustfield);
    out << setw(6)   << my_ff_pair_result[i].fsib_pair_count  << " ";
    if( wide_output() )
      out << setw(6) << my_ff_pair_result[i].hsib_pair_count << " "
          << setw(6) << my_ff_pair_result[i].other_pair_count << " ";
    
    out << setw(6)   << my_ff_pair_result[i].pair_count       << "  ";

    for( size_t j = 0; j < my_ff_pair_result[i].param_estimates.size(); ++j )
    {
      out << fp(my_ff_pair_result[i].param_estimates[j], 9, 6) << " ";
    }

    if( wide_output() )
    {
      for( size_t j = 0; j < my_ff_pair_result[i].first_derivative.size(); ++j )
      {
        out << fp(my_ff_pair_result[i].first_derivative[j], 9, 6) << " ";
      }
#if 0
      out << " " << setw(4) << my_ff_pair_result[i].function_evaluations
          << " " << setw(3) << my_ff_pair_result[i].last_error;
#endif
    }
    out << endl;
  }

  print_double_line(test);
}

void
lodpal_test_xfile::print_table_heading(const ARP_base_analysis& test, bool full)
{
  ostream& out = output_stream();

  for( size_t i = 0; i <= my_max_marker_name; ++i )
    out << " ";
  out << "           Full   ";
  if( wide_output() )
    out << "Half          ";

  if( full || ( !full && (    test.parameters().x_linkage_model().lambda1_equal
                           || (   !test.parameters().x_linkage_model().lambda1_equal
                                && test.parameters().used_pair_type() < 2 ) ) ) )
  {
    if( test.pairs_info().parameters().covariate_count() )
      out << "        Parameter Estimates";
    else
    {
      if( !test.pairs_info().parameters().x_linkage_model().lambda2_fixed )
        out << "        Parameter Estimates";
      else
        out << "        Param Est";
    }
    for( size_t ci = 1; ci < test.pairs_info().parameters().covariate_count(); ++ci )
      out << "          ";
    if( wide_output() )
    {
      if( test.pairs_info().parameters().covariate_count() )
        out << "  First Derivatives";
      else
      {
        if( !test.pairs_info().parameters().x_linkage_model().lambda2_fixed )
          out << "  First Derivatives";
        else
          out << "  1st Deriv";
      }
    }
  }
  out << endl;

  for( size_t i = 0; i <= my_max_marker_name; ++i )
    out << " ";
  out << " LOD       Sib    ";
  if( wide_output() )
    out << "Sib    Other  ";
  out << "All     ";

  if( full || ( !full && (    test.parameters().x_linkage_model().lambda1_equal
                           || (   !test.parameters().x_linkage_model().lambda1_equal
                                && test.parameters().used_pair_type() < 2 ) ) ) )
  {
    out << "---------";
    if( !test.pairs_info().parameters().x_linkage_model().lambda2_fixed )
      out << "----------";

    for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
      out << "----------";

    if( wide_output() )
    {
      out << "  ---------";
      if( !test.pairs_info().parameters().x_linkage_model().lambda2_fixed )
        out << "----------";
      for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
        out << "----------";
    }
  }
  out << endl;

  out.setf(ios_base::left, ios_base::adjustfield);
  if( test.relative_pairs().valid_distance_exist() )
    out << setw(my_max_marker_name-4) << "MARKER" << setw(5) << "cM";
  else
    out << setw(my_max_marker_name+1) << "MARKER";

  out << " SCORE     Pairs  ";
  if( wide_output() )
    out << "Pairs  Pairs  ";
  out << "Pairs   ";

  if( full || ( !full && (    test.parameters().x_linkage_model().lambda1_equal
                           || (   !test.parameters().x_linkage_model().lambda1_equal
                                && test.parameters().used_pair_type() < 2 ) ) ) )
  {
    out << "Beta1    ";
    if( !test.pairs_info().parameters().x_linkage_model().lambda2_fixed )
      out << " Beta2    ";

    for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
    {
      out << " ";
      if( test.pairs_info().parameters().covariate_parameters(ci).operation == covariate_type::none )
      {
        string c_name = test.pairs_info().parameters().covariate_parameters(ci).name(test.relative_pairs());
        if( c_name.size() > 2 )
        {
          string temp = c_name.substr(0,2);
          c_name = temp;
        }
        c_name += " discor";
        out.setf(ios_base::left, ios_base::adjustfield);
        out << setw(9) << c_name;
      }
      else
      {
        string c_name = test.pairs_info().parameters().covariate_parameters(ci).name(test.relative_pairs());
        out.setf(ios_base::left, ios_base::adjustfield);
        out << setw(9) << c_name.substr(0,9);
      }
    }
    if( wide_output() )
    {
      out << "  Beta1    ";

      if( !test.pairs_info().parameters().x_linkage_model().lambda2_fixed )
        out << " Beta2    ";

      for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
      {
        out << " ";
        if( test.pairs_info().parameters().covariate_parameters(ci).operation == covariate_type::none )
        {
          string c_name = test.pairs_info().parameters().covariate_parameters(ci).name(test.relative_pairs());
          if( c_name.size() > 2 )
          {
            string temp = c_name.substr(0,2);
            c_name = temp;
          }
          c_name += " discor";
          out << setw(9) << c_name;
        }
        else
        {
          string c_name = test.pairs_info().parameters().covariate_parameters(ci).name(test.relative_pairs());
          out.setf(ios_base::left, ios_base::adjustfield);
          out << setw(9) << c_name.substr(0,9);
        }
      }
#if 0
      out << " NFE  LFL ";
#endif
    }
  }
  out << endl;
}

void
lodpal_test_xfile::print_double_line(const ARP_base_analysis& test)
{
  ostream& out = output_stream();

  for( size_t i = 0; i <= my_max_marker_name; ++i )
    out << "=";
  out << "==================";
  if( wide_output() )
    out << "==============";
  out << "=================";
  if( !test.pairs_info().parameters().x_linkage_model().lambda2_fixed )
    out << "==========";
  for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
    out << "==========";
  if( wide_output() )
  {
    out << "===========";
    if( !test.pairs_info().parameters().x_linkage_model().lambda2_fixed )
      out << "==========";
    for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
      out << "==========";
#if 0
    out << "=========";
#endif
  }
  out << endl;
}

void
lodpal_test_xfile::print_single_line(const ARP_base_analysis& test, bool full)
{
  ostream& out = output_stream();

  for( size_t i = 0; i <= my_max_marker_name; ++i )
    out << "-";
  out << " --------- ------ ";
  if( wide_output() )
    out << "------ ------ ";
  out << "------ ";

  if( full || ( !full && (    test.parameters().x_linkage_model().lambda1_equal
                           || (   !test.parameters().x_linkage_model().lambda1_equal
                                && test.parameters().used_pair_type() < 2 ) ) ) )
  {
    out << " ---------";
    if( !test.pairs_info().parameters().x_linkage_model().lambda2_fixed )
      out << " ---------";

    for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
      out << " ---------";

    if( wide_output() )
    {
      out << "  ---------";
      if( !test.pairs_info().parameters().x_linkage_model().lambda2_fixed )
        out << " ---------";
      for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
        out << " ---------";
#if 0
      out << " ---- ---";
#endif
    }
  }
  out << endl;
}

} // end of namespace LODPAL
} // end of namespace SAGE
