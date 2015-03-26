//==========================================================================
//  File:       reltest_output.cpp
//
//  Author:     Qing Sun & Yeunjoo Song
//
//  History:    Version 1.0
//                      2.0 updated to new libraries.            yjs Jul. 03
//  Notes:
//
//  Copyright (c) 2003 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "reltest/output.h"

namespace SAGE
{

namespace RELTEST
{

reltest_output::reltest_output(const reltest_analysis* p,
                               ostream& out, ostream& fam, ostream& det, cerrorstream& err)
              : my_ra(p),max_ped_name(5),max_ind_name(7)
{
  if( !init() )
    return;
  else
    my_file_ready = true;

  errors = err;

  for( size_t i = 0; i < 5; ++i )
  {
    if(    i == 0
        && !(my_ra->get_parser()->get_analysis_pairtypes()[i])
        && !(my_ra->get_parser()->calculate_cutpoints()) )
      continue;

    if( i && !(my_ra->get_parser()->get_analysis_pairtypes()[i]) )
      continue;

    my_current_pairtype = (putative_type) i;

    if( i == 0 )
    {
      print_header(out, false);

      if( my_ra->get_parser()->generate_detailed_output() )
        print_header(det, true);
    }

    print_summary_file(out);

    if( my_ra->get_parser()->generate_detailed_output() )
      print_detailed_file(det);

    if(    my_current_pairtype == SIB
        && my_ra->get_parser()->generate_nucfam_output() )
      write_sib_info_file(fam);
  }

//  out << "\r\f";
  out << "\n\f";

  for( size_t i = 0; i < 5; ++i )
  {
    if(    i == 0
        && !(my_ra->get_parser()->get_analysis_pairtypes()[i])
        && !(my_ra->get_parser()->calculate_cutpoints()) )
      continue;

    if( i && !(my_ra->get_parser()->get_analysis_pairtypes()[i]) )
      continue;

    my_current_pairtype = (putative_type) i;

    // Print Yj histogram for all pairs.
    //
    size_t yjp_count = print_histogram(out, true);
    assert( yjp_count == 0 );

    // Print Yj* histogram for all pairs.
    //
    if(    my_current_pairtype == SIB
        || my_current_pairtype == PARENT_OFFSPRING )
      yjp_count = print_histogram(out, false, true);

    if(    (    my_current_pairtype == SIB
             || my_current_pairtype == PARENT_OFFSPRING  )
        && yjp_count == my_ra->get_ptt_pairs(my_current_pairtype).size() )
        continue;

    // Print Yj* histogram for reclassified_as_sib pairs.
    //
    yjp_count = print_histogram(out, false);
  }
}

bool
reltest_output::init()
{
  if( my_ra == NULL )
    return false;

  get_format_params();

  return true;
}

void
reltest_output::get_format_params()
{
  for( size_t p = 0; p < 5; ++p )
  {
    const vector<putative_pair>& pairs = my_ra->get_ptt_pairs((putative_type)p);

    for( size_t i = 0; i < pairs.size(); ++i )
    {
      size_t big_ind = max(pairs[i].first()->name().size(),pairs[i].second()->name().size());
      max_ind_name   = max(big_ind, max_ind_name);
      max_ped_name   = max(pairs[i].first()->pedigree()->name().size(), max_ped_name);
    }
  }
} 

bool
reltest_output::print_summary_file(ostream& out)
{
  if( !my_file_ready || !out )
    return false;

  putative_type type = my_current_pairtype;
  string ps =  pair_name(type);

  if( my_ra->get_misc_pairs(type).size() == 0 ) //no misclassified pairs.
    return true;

  out << endl << endl
      << "PUTATIVE " << ps << " PAIRS TO BE RECLASSIFIED :"  << endl;

  print_subtitle(out);
  print_body(out, my_ra->get_misc_pairs(type));
  print_trailer(out);

  return true;
}

void
reltest_output::print_subtitle(ostream& out) 
{
  if( !out )
    return;

  out << endl  << std::left
      << setw(max_ped_name)   << "     "    << setw(max_ind_name-1) << " " 
      << setw(max_ind_name+2) << "    "
      << "Reclassified" << endl
      << left << setw(max_ped_name)   << "  PID"    << setw(max_ind_name-1) << " " 
      << left << setw(max_ind_name+2) << "Pair"
      << "Pair Type         Yj          Yj*          Missing Data" << endl;

  out << "-------------------------------------------------------------------------------"
      << endl;
}

void
reltest_output::print_body(ostream& out, const vector<putative_pair>& pt)
{
  if( !out || pt.size()==0 )
    return;

  putative_pair  p;
  putative_type  ct;
  
  int class_field = 12;
  int numeric_field = 12;

  for( size_t i = 0; i < pt.size(); ++i )
  {
    p = pt[i];

    ct = p.get_pair_type();

    out << right << setw(max_ped_name) << p.first()->pedigree()->name() 
        << right << setw(max_ind_name) << p.first()->name() 
        << '/' 
        << left  << setw(max_ind_name) << p.second()->name();

     if( ct == PARENT_OFFSPRING )
       out << setw(class_field) << "PAR./OFFS." << "    ";
     else if( ct == MARITAL )
       out << setw(class_field) << "MO./FATHER" << "    ";
     else
     {
       string pname = pair_name(ct);
       out << setw(class_field) << pname << "    ";
     }
   
     out << setw(numeric_field) << fp(p.get_Yj(), 7, 4) 
         << setw(numeric_field) << fp(p.get_Yjp(), 7, 4) 
         << "    " 
         << right << setw(2) << floor(p.get_s1_missing_rate() * 100.) <<"% / "
         << right << setw(2) << floor(p.get_s2_missing_rate() * 100.) <<'%'
         << endl;     
  }
}

void
reltest_output::print_trailer(ostream& out)
{
  if( !out )
    return;

  out << "-------------------------------------------------------------------------------"
      << endl << " total putative pairs            : " << my_ra->get_ptt_pairs(my_current_pairtype).size() 
      << endl << " total to be re-classified pairs : " << my_ra->get_misc_pairs(my_current_pairtype).size()
      << endl ;  
}

bool
reltest_output::print_detailed_file(ostream& det)
{
  if( !my_file_ready || !det )
    return false;

  putative_type type = my_current_pairtype;
  string ps =  pair_name(type);

  if( my_ra->get_ptt_pairs(type).size() == 0 ) //no pairs.
    return true;

  det << endl << endl
      << "PUTATIVE " << ps << " PAIRS TO BE RECLASSIFIED :"  << endl << endl ;

  print_subtitle_det(det);
  print_body_det(det, my_ra->get_ptt_pairs(type));
  print_trailer(det);

  return true;
}

void
reltest_output::print_subtitle_det(ostream& out) 
{
  if( !out )
    return;

  out << endl
      << setw(max_ped_name)   << " "    << setw(max_ind_name-1) << " " 
      << setw(max_ind_name+2) << " "
      << "Putative    Reclassified" << endl
      << left << setw(max_ped_name)   << "  PID" << setw(max_ind_name-1) << " " 
      << left << setw(max_ind_name+2) << "Pair"
      << "Pair Type   Pair Type      Yj        Yj*     Missing Data" << endl;

  out << "-------------------------------------------------------------------------------"
      << endl;
}

void
reltest_output::print_body_det(ostream& out, const vector<putative_pair>& pt)
{
  if( !out || pt.size()==0 )
    return;

  putative_type type = my_current_pairtype;
  string ps = pair_name(type);

  putative_pair  p;
  putative_type  ct;
  
  int class_field = 12;
//  int numeric_field = 12;

  for( size_t i = 0; i < pt.size(); ++i )
  {
    p = pt[i];

    ct = p.get_pair_type();

    out << right << setw(max_ped_name) << p.first()->pedigree()->name() 
        << right << setw(max_ind_name) << p.first()->name() 
        << '/' 
        << left  << setw(max_ind_name) << p.second()->name();

     /** output parents of putative sib pairs may be helpful.    
     if( type==SIB )
       out << setw(max_ind_name) 
           << mps->get(p->first(),MP_Base::Mother)->name() 
           << setw(max_ind_name) 
           << mps->get(p->first(),MP_Base::Father)->name();
     **/

     if( ct == PARENT_OFFSPRING )
       out << setw(class_field) << "PAR./OFFS.";
     else if( ct == MARITAL )
       out << setw(class_field) << "MO./FATHER";
     else
       out << setw(class_field) << ps;

     if( ct == PARENT_OFFSPRING )
       out << setw(class_field) << "PAR./OFFS.";
     else if( ct == MARITAL )
       out << setw(class_field) << "MO./FATHER";
     else
     {
       string pname = pair_name(ct);
       out << setw(class_field) << pname;
     }
   
     out << "  "  << fp(p.get_Yj(), 7, 4) 
         << "   " << fp(p.get_Yjp(), 7, 4) 
         << "   " 
         << right << setw(2) << floor(p.get_s1_missing_rate() * 100.) <<"% / "
         << right << setw(2) << floor(p.get_s2_missing_rate() * 100.) <<'%'
         << endl;     
  }
}

void
reltest_output::print_header(ostream& out, bool detailed)
{ 
  if( !out )
    return;

  out << endl
      << "===============================================================================";
  out << endl;

  if( detailed )
    out << "  RELATIONSHIP TEST PROGRAM DETAILED OUTPUT";
  else
    out << "  RELATIONSHIP TEST PROGRAM SUMMARY OUTPUT";

  out << endl
      << "===============================================================================";
  out << endl << endl;

  string analysis_name = my_ra->get_parser()->get_analysis_name();

  out << "    Analysis Name                      : " << analysis_name
      << endl << endl
      << "    Average Marker Information Content : " << my_ra->get_AMIC()
      << endl
      << "    Total Length of Genome             : " << my_ra->get_total_genome_length() << " (cM) "
      << endl << endl;
//option 1
/*
      << "    Cut Points                 :              |  original       adjusted" << endl
      << "                                 -------------|-------------------------" << endl
      << "                                 unrelated    | " << fp(my_ra->get_cutpoints[my_ra->get_Cu],8,5)
                                               << "       " << fp(my_ra->get_adjusted_cutpoints[my_ra->get_Cu],8,5) << endl
      << "                                 half sib     | " << fp(my_ra->get_cutpoints[my_ra->get_Ch],8,5)
                                               << "       " << fp(my_ra->get_adjusted_cutpoints[my_ra->get_Ch],8,5) << endl
      << "                                 MZtwins      | " << fp(my_ra->get_cutpoints[my_ra->get_Cm],8,5)
                                               << "       " << fp(my_ra->get_adjusted_cutpoints[my_ra->get_Cm],8,5) << endl
      << "                                 parent/child | " << fp(my_ra->get_cutpoints[my_ra->get_Cp],8,5)
                                               << "       " << fp(my_ra->get_adjusted_cutpoints[my_ra->get_Cp],8,5) << endl
      << endl;
*/
//option 2

  out << "    Cut-points" << endl
      << "                                      |           |  original       adjusted" << endl
      << "      --------------------------------|-----------|--------------------------" << endl
      << "      Sibling Classification          |unrelated  | " << fp(my_ra->get_cutpoints(Cu),8,5)
                                               << "       " << fp(my_ra->get_adjusted_cutpoints(Cu),8,5) << endl
      << "        Statistics(Yj)                |half sib   | " << fp(my_ra->get_cutpoints(Ch),8,5)
                                               << "       " << fp(my_ra->get_adjusted_cutpoints(Ch),8,5) << endl
      << "                                      |MZtwins    | " << fp(my_ra->get_cutpoints(Cm),8,5)
                                               << "       " << fp(my_ra->get_adjusted_cutpoints(Cm),8,5) << endl
      << "      --------------------------------|-----------|--------------------------" << endl
      << "      Parent/Offspring Classification |parent/    |" << endl
      << "        Statistics(Yj*)               | offspring | " << fp(my_ra->get_cutpoints(Cp),8,5)
                                               << "       " << fp(my_ra->get_adjusted_cutpoints(Cp),8,5) << endl
      << "      -----------------------------------------------------------------------" << endl
      << endl;

  const vector<solution_type>& solution = my_ra->get_solution_Yj();
  size_t picked = my_ra->get_picked_Yj();
  double mu     = my_ra->get_mu_Yj();

  out << "    Sibling Classification Statistics(Yj)" << endl;

  if( finite(mu) )
    out << "       robust (L2) mean     : " << solution[picked].mean << endl
        << "       robust (L2) variance : " << solution[picked].variance << endl;
  else
  {
    for( size_t i = 0; i < solution.size(); ++i )
      out << "       robust (L2) mean     : " << solution[i].mean << endl
          << "       robust (L2) variance : " << solution[i].variance << endl;
  }

  out << endl;

  const vector<solution_type>& solution_p = my_ra->get_solution_Yjp();
  picked     = my_ra->get_picked_Yjp();
  mu         = my_ra->get_mu_Yjp();

  out << "    Parent/Offspring Classification Statistics(Yj*)" << endl;

  if( finite(mu) )
    out << "       robust (L2) mean     : " << solution_p[picked].mean << endl
        << "       robust (L2) variance : " << solution_p[picked].variance << endl;
  else
  {
    for( size_t i = 0; i < solution_p.size(); ++i )
      out << "       robust (L2) mean     : " << solution_p[i].mean << endl
          << "       robust (L2) variance : " << solution_p[i].variance << endl;
  }

  out << endl;

  double mean = my_ra->get_Yj_mean();
  double std2 = 2. * my_ra->get_standard_error();

  out << "    Average Yj of Pairs" << endl
      << "      Reclassified as Full Sibs : " << mean << endl
      << "    Standard Error              : " << my_ra->get_standard_error() << endl
      << "    95% Confidence Interval     : " << mean - std2 << " to " << mean + std2 << endl
      << endl;

  if(    (mean - std2 > 0. && mean + std2 > 0.)
      || (mean - std2 < 0. && mean + std2 < 0.) )
    out << "  ! WARNING : THE MEAN OF THE SIB-PAIR DISTRIBUTION DIFFERS SIGNIFICANTLY FROM" << endl
        << "              ZERO.  YOU MAY HAVE SUBSTANTIAL DATA ERROR OR MISSPECIFICATION OF" << endl
        << "              PARAMETERS SUCH AS ALLELE FREQUENCIES." << endl
        << endl;

  if( my_ra->get_total_genome_length() < 1200 )
    out << "  ! WARNING : HIGH PERCENTAGE OF MISSING DATA MAY RESULT IN MISCLASSIFICATION." << endl
        << endl;

  out << "===============================================================================";
  out << endl;
}

size_t
reltest_output::print_histogram(ostream& out, bool Yj, bool all) 
{
  if( !out )
    return (size_t)-1;

  putative_type  type = my_current_pairtype;

  if( !my_ra->get_ptt_pairs(type).size() )
    return 0;

  putative_pair  p;
  int            Yjp_counter = 0;

  string         pname = pair_name(type);
  
  double Y_max = -numeric_limits<double>::infinity();
  double Y_min =  numeric_limits<double>::infinity();

  for( size_t i = 0; i < my_ra->get_ptt_pairs(type).size(); ++i )
  {
    p = my_ra->get_ptt_pairs(type)[i];

    if( Yj )
    {
      Y_min = min(Y_min, p.get_Yj());
      Y_max = max(Y_max, p.get_Yj());
    }
    else
    {
      if( p.get_pair_type() == SIB )
        Yjp_counter++;            

      if( all )
      {
        Y_min = min(Y_min, p.get_Yjp());
        Y_max = max(Y_max, p.get_Yjp());
      }
      //only pairs classified as sib pairs are meaningful to Yjp*
      else if( p.get_pair_type() == SIB )
      {
        Y_min = min(Y_min, p.get_Yjp());
        Y_max = max(Y_max, p.get_Yjp());
      }
    }
  }

  //no element of Yj* to output when Yj* histogram is required.
  if( !Yj && Yjp_counter==0 && !all )
    return Yjp_counter; 

  int max_bin = (int) (ceil( ((Y_max+0.01) - (Y_min-0.01)) / 0.25 ));  //default bin size=0.25

  Histogram hist_Y(max_bin, Y_min-0.01, (Y_min-0.01)+0.25*max_bin);

  for( size_t i = 0; i < my_ra->get_ptt_pairs(type).size(); ++i )
  {
    putative_pair  p = my_ra->get_ptt_pairs(type)[i];

    if( Yj )
      hist_Y.add(p.get_Yj());
    else
    {
      if( all )
        hist_Y.add(p.get_Yjp());  
      else if( !all && p.get_pair_type() == SIB )
        hist_Y.add(p.get_Yjp());  
    }
  }

  string Y_str;
  if( Yj )
    Y_str = "Yj";
  else
    Y_str = "Yj*";

  out << endl << endl
      << "  ============================================================" << endl
      << "  ==                                                        ==" << endl;

  if( Yj )
    out << "  ==   HISTOGRAM OF SIBLING CLASSIFICATION STATISTIC (Yj)   ==" << endl;
  else
    out << "  ==  HISTOGRAM OF PAR/OFF CLASSIFICATION STATISTIC (Yj*)   ==" << endl;

  out << "  ==   ";
  size_t margin = 50 - (19 + pname.size());
  size_t s = 0;
  for( ; s < margin/2; ++s )
    out << " ";
  out << "FOR PUTATIVE " << pname << " PAIRS";
  for( ; s < margin; ++s )
    out << " ";
  out << "   ==" << endl;
  if( !Yj && !all )
    out << "  ==                RECLASSIFIED AS FULL SIB                ==" << endl;

//  if( !Yj && !all )
//    out << "  ==        FOR PAIRS RECLASSIFIED AS FULL SIB              ==" << endl;
//  else
//    out << "  ==                FOR PUTATIVE PAIRS                      ==" << endl;

  out << "  ==                                                        ==" << endl
//      << "  ==           putative pair type : "
//      <<                       setw(23) << pname    <<              " ==" << endl
      << "  ==                  maximum " << setw(3) << Y_str <<" : "
      <<                          setw(10) << Y_max << "              ==" << endl
      << "  ==                  minimum " << setw(3) << Y_str <<" : "
      <<                          setw(10) << Y_min << "              ==" << endl;
  
  //generate parameter for graphic output
  size_t m_points = hist_Y.maxFrequency().second; 
  double bin_size = hist_Y.binsize();
  double points_per_star = 1;

  if(m_points > 36)
    points_per_star = ceil(m_points/36.0);

  out << "  ==                     bin size : " 
      <<                       setw(10) << bin_size << "              ==" << endl
      << "  ==                                                        ==" << endl;

  out << "  ============================================================" << endl
      << endl << endl
      << "      Interval              Count (one * is equal up to " 
      << points_per_star;

  if( points_per_star > 1 )
    out << " pairs.)" << endl;
  else
    out << " pair.)"  << endl;

  out << "===================================================================="
      << endl;

  //output histogram

  // Calculate the bottom of the first bin.
  double prev_bin = hist_Y.lower_bound();

  for( size_t i = 0; i < hist_Y.max_bins(); ++i )
  {
    out << "  " 
        << fp(prev_bin, 5, 2) << "  to  " ;

    prev_bin += bin_size;

    out << fp(prev_bin, 5, 2) << "\t   " << right << setw(4) << hist_Y[i];

    if( hist_Y[i] == 0 ) 
      out << endl;
    else
      graph_hist(out, (int) (ceil(hist_Y[i]/points_per_star)) );
  }
  
  out << "--------------------------------------------------------------------" << endl
      << "                  Total :" << setw(6) << hist_Y.point_count();

  out << endl << endl;

  return Yjp_counter;
}

void
reltest_output::graph_hist(ostream& out, int count)
{
  if( !out )
    return;

  out << ' ';
  if( count == 0 ) //bins contains less than 'points_per_star'
  {
    out << '*' << endl;
    return;
  }

  for( int i = 0; i < count; ++i )
    out << '*';
  
  out << endl;
}


string
reltest_output::pair_name(putative_type type) const
{
  if( type == SIB )
    return "FULL SIB";
  else if( type == HSIB )
    return "HALF SIB";
  else if( type == MZTWIN )
    return "MZTWINS";
  else if( type == MARITAL )
    return "MOTHER FATHER";
  else if( type == PARENT_OFFSPRING )
    return "PARENT OFFSPRING";
  else
    return "UNRELATED";
}

//generating "nucFam.inf" file which contains all sibling informations of
//to be re-classified putative SIB PAIRS.
void
reltest_output::write_sib_info_file(ostream& fam)
{
  putative_type type = my_current_pairtype;

  if( !fam || my_ra->get_misc_pairs(type).empty() )
    return;

  vector<putative_pair> info_pairs;
 
  size_t current_pid;

  size_t pid = my_ra->get_misc_pairs(type)[0].first()->pedigree()->index();

  for( size_t i = 0; i < my_ra->get_misc_pairs(type).size(); ++i )
  {
    current_pid = my_ra->get_misc_pairs(type)[i].first()->pedigree()->index();

    if( current_pid == pid )
      continue;

    for( size_t j = 0; j < my_ra->get_ptt_pairs(type).size(); ++j )
    {
      if( my_ra->get_ptt_pairs(type)[j].first()->pedigree()->index() == pid )
        info_pairs.push_back(my_ra->get_ptt_pairs(type)[j]);
    }
    pid = current_pid;
  }

  string analysis_name = my_ra->get_parser()->get_analysis_name();

  fam << endl
      << "===============================================================================" << endl
      << "  RELATIONSHIP TEST PROGRAM NUCLEAR FAMILY INFORMATION" << endl
      << "===============================================================================" << endl
      << endl
      << "    Note          : This file contains information about all sib pairs of the" << endl
      << "                    nuclear families in which at least one sib pair should be" << endl
      << "                    reclassified." << endl
      << endl
      << "    Analysis Name : " << analysis_name << endl
      << endl;

  fam << "===============================================================================" << endl;

  print_subtitle(fam);
  print_body(fam, info_pairs);

  fam << "===============================================================================" << endl;
}

} // end of namespace RELTEST

} // end of namespace SAGE
