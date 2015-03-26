//============================================================================
// File:     stats_view.cpp
//                                                                          
// Author:   Dan Baechle
//                                                                          
// History:  Initial version: 10/00 
//                                                                          
// 
// Notes:    Implements the following classes -
//              base_stats_viewer
//              binary_stats_viewer
//              cont_stats_viewer
//              cmpd_stats_viewer
//              general_stats_viewer
//              ped_stats_viewer
//              mp_stats_viewer
//              base_trait_stats_viewer
//              log_histogram
//                                                                          
// Copyright (c) 2000 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "pedinfo/stats_view.h"

using namespace std;

namespace SAGE {
namespace RPED {

//============================================================================
// IMPLEMENTATION:  base_stats_viewer
//============================================================================
//
void
base_stats_viewer::major_header(std::string title) const
{
  outer_line();
  header_body(title);
}

void
base_stats_viewer::header_body(std::string title) const
{
  if(title.length() > LINE_SIZE - 5)
  {
    // Break title into multiple lines if it is too long to fit available space.
    string               remainder = "";
    string               line = title;
    string::size_type    break_pt = 0;
    while(line.length() > LINE_SIZE - 5 && break_pt != string::npos)
    {
      break_pt = line.rfind(' ');  
      if(break_pt != string::npos)
      {
        remainder = line.substr(break_pt) + remainder;
        line = line.substr(0, break_pt);
      }
    }
    
    if(line.length() > LINE_SIZE - 5)            // No breakpoint found.
    {
      // Arbitrarily break the line.
      remainder = line.substr(LINE_SIZE - 6) + remainder;    
      line = line.substr(0, LINE_SIZE - 6);
    }
    else if(line.length() == 0)                   // Breakpoint at position 0.
    {
      // Arbitrarily break the remainder.
      if(remainder.length() > LINE_SIZE - 5)
      {
        line = remainder.substr(0, LINE_SIZE - 6);
        remainder = remainder.substr(LINE_SIZE - 6);
      }
      else             // Should never happen.
      {
        header_line(remainder);
        return;
      }
    }
    else                                         // Line size OK.
    {
      ;
    }
    
    header_line(line);
    header_body(remainder);   
  }
  else
  {
    header_line(title);
  }
}


void
base_stats_viewer::header_line(std::string line) const
{
  // Make sure line has an even number of characters for ease of centering.
  string even_title;
  if(line.length() % 2 == 1)
  {
    even_title = line + " ";
  }
  else
  {
    even_title = line;
  }
  
  // Center title.
  size_t margin = (LINE_SIZE - 2 - even_title.length()) / 2;
  
  my_o << "|" << setw(margin) << "\0" << even_title 
              << setw(margin) << "\0" << "|" << endl;
}

void
base_stats_viewer::display_double(double value, unsigned int width, 
                                  unsigned int prec, unsigned int sig_dig) const
{
  ios_base::fmtflags old_flags = my_o.flags();
  
  if( !finite(value) )
  {
    my_o << setw(width) << right << "--- ";
  }
  else
  {
    if(abs(value) >= 1)
    {
      // Width required to display w/o scientific notation.  '2' allows for sign
      // and decimal point.
      //
      unsigned int  integral_dig  = static_cast<unsigned int>(floor(log10(abs(value)))) + 1;
      unsigned int  req_width     = 2 + integral_dig + prec;
      if(req_width > width)
      {
        my_o.setf(ios_base::scientific, ios_base::floatfield);
        my_o.precision(prec);
      }
      else
      {
        my_o.setf(ios_base::fixed, ios_base::floatfield);
        my_o.precision(prec);
      }
    }
    else
    {
      unsigned int  zeros = 0;
      
      if(value != 0)
      {
        zeros = static_cast<unsigned int>(ceil(abs(log10(abs(value))))) - 1;
      }
      
      // Allow for leading zero as well as sign and decimal pt.
      unsigned int  req_width  = 3 + sig_dig + zeros;
      if(req_width > width)
      {
        my_o.setf(ios_base::scientific, ios_base::floatfield);
        my_o.precision(sig_dig);      
      }
      else
      {
        my_o.setf(ios_base::fixed, ios_base::floatfield);
        my_o.precision(zeros + sig_dig);
      }
    }
    
    my_o << setw(width) << right << value;
  }
  my_o.flags(old_flags);
  my_o.precision(6);                       // Restore default value.
}


//============================================================================
// IMPLEMENTATION:  binary_stats_viewer
//============================================================================
//
// - Display table of individual counts by affection status and gender.
//
void 
binary_stats_viewer::ind_data() const
{
  ind_header();
  single_inner_line();
  
  for(int g = 0; g < 3; ++g)
  {
    ind_gender_data_line(static_cast<bt::gender>(g));
  }
  ind_gender_total_line();
  
  ind_blank_data_line();
  
  for(int f = 0; f < 3; ++f)
  {
    ind_founder_data_line(static_cast<bt::founder_status>(f));
  }
  ind_founder_total_line();
}
    
void     
binary_stats_viewer::ind_gender_data_line(bt::gender g) const
{
  if(0 <= g && g < 3)
  {
    ios_base::fmtflags old_flags = my_o.flags();    
    
    my_o << "|" << setw(FIELD_SIZE) << left << gender_label(g) << "|"
                                    << right;
 
    // Display each count for each affection status and compute total.
    size_t total = 0;
    size_t count;
    for(int st = 0; st < 3; ++st)
    {
      count = my_binary_stats.ind_count_gender(g, static_cast<bt::ind_affection_status>(st));
      total += count;
      my_o << setw(FIELD_SIZE) << count << "|";
    }
    
    my_o << setw(FIELD_SIZE) << total << "|"
         << setw(FIELD_SIZE) << "\0"  << "|"
         << setw(FIELD_SIZE) << "\0"  << "|"
         << endl;
         
    my_o.flags(old_flags);
  }
}

void
binary_stats_viewer::ind_blank_data_line() const
{
  ios_base::fmtflags old_flags = my_o.flags();
  
  my_o << "|";
  for(int i = 0; i < 7; ++i)
  {
    my_o << setw(FIELD_SIZE) << "\0" << "|";
  }
  my_o << endl;
  
  my_o.flags(old_flags);
}

void     
binary_stats_viewer::ind_founder_data_line(bt::founder_status f) const
{
  if(0 <= f && f < 3)
  {
    ios_base::fmtflags old_flags = my_o.flags();    
    
    my_o << "|" << setw(FIELD_SIZE) << left << founder_status_label(f) << "|"
                                    << right;
 
    // Display each count for each affection status and compute total.
    size_t total = 0;
    size_t count;
    for(int st = 0; st < 3; ++st)
    {
      count = my_binary_stats.ind_count_founder(f, static_cast<bt::ind_affection_status>(st));
      total += count;
      my_o << setw(FIELD_SIZE) << count << "|";
    }
    
    my_o << setw(FIELD_SIZE) << total << "|"
         << setw(FIELD_SIZE) << "\0"  << "|"
         << setw(FIELD_SIZE) << "\0"  << "|"
         << endl;
         
    my_o.flags(old_flags);
  }
}

void     
binary_stats_viewer::ind_gender_total_line() const
{
  // Compute the totals.
  size_t totals[4] = {0, 0, 0, 0};
  for(int st = 0; st < 3; ++st)
  {
    for(int g = 0; g < 3; ++g)
    {
      totals[st] += my_binary_stats.ind_count_gender(static_cast<bt::gender>(g), 
                                              static_cast<bt::ind_affection_status>(st)); 
    }
    totals[3] += totals[st];       // Grand total.
  }
  
  // Display total line.
  ios_base::fmtflags old_flags = my_o.flags();
  
  my_o << "|" << setw(FIELD_SIZE) << left  << "  Total" << "|"
                                  << right;
  for(int i = 0; i < 4; ++i)
  {
    my_o << setw(FIELD_SIZE) << totals[i] << "|";
  }
  my_o << setw(FIELD_SIZE) << "\0" << "|";
  my_o << setw(FIELD_SIZE) << "\0" << "|" << endl;
  
  my_o.flags(old_flags);
}

void     
binary_stats_viewer::ind_founder_total_line() const
{
  // Compute the totals.
  size_t totals[4] = {0, 0, 0, 0};
  for(int st = 0; st < 3; ++st)
  {
    for(int f = 0; f < 3; ++f)
    {
      totals[st] += my_binary_stats.ind_count_founder(static_cast<bt::founder_status>(f), 
                                              static_cast<bt::ind_affection_status>(st)); 
    }
    totals[3] += totals[st];       // Grand total.
  }
  
  // Display total line.
  ios_base::fmtflags old_flags = my_o.flags();
  
  my_o << "|" << setw(FIELD_SIZE) << left  << "  Total" << "|"
                                  << right;
  for(int i = 0; i < 4; ++i)
  {
    my_o << setw(FIELD_SIZE) << totals[i] << "|";
  }
  my_o << setw(FIELD_SIZE) << "\0" << "|";
  my_o << setw(FIELD_SIZE) << "\0" << "|" << endl;
  
  my_o.flags(old_flags);
}

// - Display table of pair counts by affection status.
//
void    
binary_stats_viewer::pair_data() const
{
  pair_header();
  single_inner_line();
  for(int t = 0; t < TOTAL_PAIR_TYPES; ++t)
  {
    pair_data_line(static_cast<pg::pair_type>(t + 1));
  }
}

void     
binary_stats_viewer::pair_data_line(pg::pair_type t) const
{
  if(0 < t &&  t <= TOTAL_PAIR_TYPES)
  {
    ios_base::fmtflags old_flags = my_o.flags();
    
    my_o << "|" << setw(FIELD_SIZE) << left << pair_label(t) << "|"
                                    << right;
 
    // Display count for each affection type and compute total.
    size_t total = 0;
    size_t count;
    for(int st = 0; st < 4; ++st)
    {
      count = my_binary_stats.pair_count(t, static_cast<pft::affection_status>(st + 1));
      total += count;
      my_o << setw(FIELD_SIZE) << count << "|";
    }
    
    my_o << setw(FIELD_SIZE) << total << "|";
    double corr = my_binary_stats.correlation(t);
    if(SAGE::isnan(corr) || !finite(corr))
    {
      my_o << setw(FIELD_SIZE) << "--- " << "|" << endl;
    }
    else
    {
      my_o.setf(ios_base::fixed, ios_base::floatfield);
      my_o << setw(FIELD_SIZE) << setprecision(4) << corr << "|" << endl;
    }
    my_o.flags(old_flags);
  }
}

//============================================================================
// IMPLEMENTATION:  cont_stats_viewer
//============================================================================
//
// - Display table of individual trait values by gender.
//
void 
cont_stats_viewer::ind_data() const
{
  ind_header();
  single_inner_line();
  
  for(int g = 0; g < 4; ++g)
  {
    ind_gender_data_line(static_cast<ct::gender>(g));
  }
  
  ind_blank_data_line();
  
  for(int f = 0; f < 4; ++f)
  {
    ind_founder_data_line(static_cast<ct::founder_status>(f));
  }
}

    
void     
cont_stats_viewer::ind_gender_data_line(ct::gender g) const
{
  if(0 <= g && g < 4)
  {
    ios_base::fmtflags old_flags = my_o.flags();
    
    size_t count = my_cont_stats.ind_gender_count(g);
    my_o << "|" << setw(FIELD_SIZE - 1) << left  << gender_label(g)            << "|"
                << setw(FIELD_SIZE - 1) << right << count                      << "|" << " ";
    
    display_double(my_cont_stats.ind_gender_mean(g), FIELD_SIZE - 1);
    my_o << " +/- ";
    display_double(my_cont_stats.ind_gender_std_dev(g), FIELD_SIZE - 1);
    my_o << " (";
    display_double(my_cont_stats.ind_gender_min(g), FIELD_SIZE - 1);
    my_o << ", ";
    display_double(my_cont_stats.ind_gender_max(g), FIELD_SIZE - 1);
    my_o << ") "                          << "|"
         << setw(FIELD_SIZE_CORR) << "\0" << "|"
                                  << endl;
                                  
    my_o.flags(old_flags);
  }
}

void
cont_stats_viewer::ind_blank_data_line() const
{
  ios_base::fmtflags old_flags = my_o.flags();
  
  my_o << "|";
  for(int i = 0; i < 2; ++i)
  {
    my_o << setw(FIELD_SIZE - 1) << "\0" << "|";
  }
  my_o << setw(48) << "\0" << "|";
  my_o << setw(FIELD_SIZE_CORR) << "\0" << "|";
  my_o << endl;
  
  my_o.flags(old_flags);
}

void     
cont_stats_viewer::ind_founder_data_line(ct::founder_status f) const
{
  if(0 <= f && f < 4)
  {
    ios_base::fmtflags old_flags = my_o.flags();
    
    size_t count = my_cont_stats.ind_founder_count(f);
    my_o << "|" << setw(FIELD_SIZE - 1) << left  << founder_status_label(f)    << "|"
                << setw(FIELD_SIZE - 1) << right << count                      << "|" << " ";
    
    display_double(my_cont_stats.ind_founder_mean(f), FIELD_SIZE - 1);
    my_o << " +/- ";
    display_double(my_cont_stats.ind_founder_std_dev(f), FIELD_SIZE - 1);
    my_o << " (";
    display_double(my_cont_stats.ind_founder_min(f), FIELD_SIZE - 1);
    my_o << ", ";
    display_double(my_cont_stats.ind_founder_max(f), FIELD_SIZE - 1);
    my_o << ") "                          << "|"
         << setw(FIELD_SIZE_CORR) << "\0" << "|"
                                  << endl;
                                  
    my_o.flags(old_flags);
  }
}

// - Display table of pair trait values.
//
void    
cont_stats_viewer::pair_data() const      
{
  pair_header();
  single_inner_line();
  for(int i = 0; i < 5; ++i)
  {
    pair_data_line(static_cast<pg::pair_type>(i + 1));
  }
}

// - Display table of pair trait values w. additional pair statistics.
//
void    
cont_stats_viewer::pair_data_alt() const      
{
  pair_header_alt();
  single_inner_line();
  for(int t = 0; t < TOTAL_PAIR_TYPES; ++t)
  {
    pair_data_line_alt(static_cast<pg::pair_type>(t + 1));
  }
}

void     
cont_stats_viewer::pair_data_line(pg::pair_type t) const
{
  ios_base::fmtflags old_flags = my_o.flags();
   
  pair_data_half_line(t);
  my_o << "|";
  pair_data_half_line(static_cast<pg::pair_type>(t + 5));
  my_o << "|" << endl;
                                
  my_o.flags(old_flags);
}

void
cont_stats_viewer::pair_data_half_line(pg::pair_type t) const
{
  ios_base::fmtflags old_flags = my_o.flags();
  
  if(0 < t && t <= TOTAL_PAIR_TYPES)
  {
    size_t count = my_cont_stats.pair_count(t);
    
    my_o << "|" << setw(FIELD_SIZE + 1) << left  << pair_label(t) << "|"
                << setw(FIELD_SIZE - 1) << right << count << "|";
    
    double corr = my_cont_stats.correlation(t);
    if(SAGE::isnan(corr) || !finite(corr))
    {
      my_o << setw(FIELD_SIZE_CORR + FIELD_SIZE - 2) << "--- ";
    }
    else
    {
      my_o.setf(ios_base::fixed, ios_base::floatfield);
      my_o.precision(4);
      my_o << setw(FIELD_SIZE_CORR + FIELD_SIZE - 2) << corr;
    }
  }
  else if(t == TOTAL_PAIR_TYPES + 1)
  {
    my_o <<  "|" << setw(FIELD_SIZE + 1) << "\0" << "|"
                 << setw(FIELD_SIZE - 1) << "\0" << "|"
                 << setw(FIELD_SIZE_CORR + FIELD_SIZE - 2) << "\0";   
  }
  
  my_o.flags(old_flags);
}

void     
cont_stats_viewer::pair_data_line_alt(pg::pair_type t) const
{
  if(0 < t &&  t <= TOTAL_PAIR_TYPES)
  {
    ios_base::fmtflags old_flags = my_o.flags();
    
    size_t count = my_cont_stats.pair_count(t);
    my_o << "|" << setw(FIELD_SIZE - 1) << left  << pair_label(t)       << "|"
                << setw(FIELD_SIZE - 1) << right << count               << "|" << " ";
                //<< setw(48) << "\0" << "|";
    
    // - Per Dr. Elston on 10/9/00, this data not meaningful.
    //
    display_double(my_cont_stats.pair_mean(t), FIELD_SIZE - 1);
    my_o << " +/- ";
    display_double(my_cont_stats.pair_std_dev(t), FIELD_SIZE - 1);
    my_o << " (";
    display_double(my_cont_stats.pair_min(t), FIELD_SIZE - 1);
    my_o << ", ";
    display_double(my_cont_stats.pair_max(t), FIELD_SIZE - 1);
    my_o << ") "                          << "|";
    
    
                                  
    double corr = my_cont_stats.correlation(t);
    if(SAGE::isnan(corr) || !finite(corr))
    {
      my_o << setw(FIELD_SIZE_CORR) << "--- "                                << "|"
                                    << endl;
    }
    else
    {
      my_o.setf(ios_base::fixed, ios_base::floatfield);
      my_o.precision(4);
      my_o << setw(FIELD_SIZE_CORR) << corr                                  << "|"
                                    << endl;
    }
                                  
    my_o.flags(old_flags);
  }
}

//============================================================================
// IMPLEMENTATION:  cmpd_stats_viewer
//============================================================================
//
void          
cmpd_stats_viewer::view() const
{
  my_o << "\n\n";
  std::string title = build_title();
  major_header(title);
  my_base_trait_stats_viewer.view();
  pair_ind_data();
  outer_line();
}

std::string          
cmpd_stats_viewer::build_title() const
{
  std::string temp;
  std::string name;
  
  name = my_cmpd_stats.pedigree_name();
  temp += "Trait Statistics: ";
  if(name == "")
  {
    temp += "All Pedigrees";
  }
  else
  {
    temp += "Pedigree - ";
    temp += name;
  }
  temp += ", Traits - ";
  temp += my_cmpd_stats.traits_name();
  return temp;
}

void          
cmpd_stats_viewer::pair_ind_data() const
{
  pair_ind_header();
  single_inner_line();
  
  size_t total = 0;
  for(int i = 0; i < TOTAL_PAIR_TYPES; ++i)
  {
    pair_half_line(static_cast<pg::pair_type>(i + 1));
    switch(i)
    {
      case 0: case 1: case 2:               // Individual counts male, female, unknown.
      case 4:                               // Blank.
      case 5: case 6: case 7:               // Individual counts founder, non-founder, unconnected.
        total += ind_half_line(i);
        break;
      case 3: case 8:                       // Total count by gender, by founder status.
        total_ind_half_line(total);
        total = 0;
        break;
      default:
      ;
    }
    my_o << endl;
  }
}

void          
cmpd_stats_viewer::pair_half_line(pg::pair_type t) const
{
  if(0 < t &&  t <= TOTAL_PAIR_TYPES)
  {
    ios_base::fmtflags old_flags = my_o.flags();
    
    my_o << "|" << setw(FIELD_SIZE_PED * 2) << left << pair_label(t) << "|"
                                            << right;
 
    size_t count = my_cmpd_stats.pair_count(t);
    my_o << setw(FIELD_SIZE_PED) << count;
    
    my_o.flags(old_flags);
  }
}

size_t        
cmpd_stats_viewer::ind_half_line(int i) const
{
  ios_base::fmtflags old_flags = my_o.flags();
  
  my_o << "||" << setw(FIELD_SIZE_PED * 2);
  size_t count = 0;
  
  switch(i)
  {
    case 0:
      my_o << left << "  Male" << "|"
           << right;
      count = my_cmpd_stats.ind_count_gender(cp::MALE);
      my_o << setw(FIELD_SIZE_PED) << count << "|";
      break;
    case 1:
      my_o << left << "  Female" << "|"
           << right;
      count = my_cmpd_stats.ind_count_gender(cp::FEMALE);
      my_o << setw(FIELD_SIZE_PED) << count << "|";
      break;
    case 2:
      my_o << left << "  Unknown" << "|"
           << right;
      count = my_cmpd_stats.ind_count_gender(cp::UNKNOWN);
      my_o << setw(FIELD_SIZE_PED) << count << "|";
      break;
    case 4: 
      my_o << left << "\0" << "|"
           << right;
      my_o << setw(FIELD_SIZE_PED) << "\0" << "|";
      break;
    case 5:
      my_o << left << "  Founder" << "|"
           << right;
      count = my_cmpd_stats.ind_count_founder(cp::FOUNDER);
      my_o << setw(FIELD_SIZE_PED) << count << "|";
      break;
    case 6:
      my_o << left << "  Non-founder" << "|"
           << right;
      count = my_cmpd_stats.ind_count_founder(cp::NON_FOUNDER);
      my_o << setw(FIELD_SIZE_PED) << count << "|";
      break;
    case 7:
      my_o << left << "  Singleton" << "|"
           << right;
      count = my_cmpd_stats.ind_count_founder(cp::UNCONNECTED);
      my_o << setw(FIELD_SIZE_PED) << count << "|";
      break;
    default:
    ;
  }
  my_o.flags(old_flags);
  return count;
}

void          
cmpd_stats_viewer::total_ind_half_line(size_t total) const
{
  ios_base::fmtflags old_flags = my_o.flags();
  
  my_o << "||" << setw(FIELD_SIZE_PED * 2) << left   << "        Total" << "|"
                                              << right
              << setw(FIELD_SIZE_PED) << total << "|";
  
  my_o.flags(old_flags);
}


//============================================================================
// IMPLEMENTATION:  general_stats_viewer
//============================================================================
//
void
general_stats_viewer::view(std::string title) const
{
  // Non-trait data.
  if(title == "")         // All pedigrees.
  {
    ;
  }
  else                    // Individual pedigrees.
  {
    my_o << "\n\n";
    major_header(title);
    
    // - Write maximum inheritance vector bits for ea. pedigree (scr 505).
    //   Added 6-25-3, djb.
    //
    double_inner_line();
    iv_bit_line();    
  }
  
  sib_header();
  single_inner_line();
  //nuc_line();
  sib_line();
  sub_loop_header();
  sub_loop_line();
  pair_ind_data();
}

void
general_stats_viewer::iv_bit_line() const
{
  my_o << "|Inheritance Vector Bits   " << setw(LINE_SIZE - 28) 
       << left << my_gen_stats.likelihood_bits() << "|" << right << endl;
}

void
general_stats_viewer::sib_line() const
{
  ios_base::fmtflags old_flags = my_o.flags();
  
  my_o << "|" << setw(FIELD_SIZE + FIELD_SIZE_CORR) << left  << "Sibships"  << "|";
  size_t count = my_gen_stats.sibship().count();
  my_o        << setw(FIELD_SIZE - 1)                       << right << count        << "|" << " ";
  
  display_double(count == 0 ? numeric_limits<double>::quiet_NaN() : my_gen_stats.sibship().mean(), 
                 FIELD_SIZE - 1);
  my_o << " +/- ";
  display_double(count == 0 ? numeric_limits<double>::quiet_NaN() : my_gen_stats.sibship().standard_deviation(), 
                 FIELD_SIZE - 1);
  my_o << " (";
  display_double(count == 0 ? numeric_limits<double>::quiet_NaN() : my_gen_stats.sibship().min(), 
                 FIELD_SIZE - 1, 0, 0);
  my_o << ", ";
  display_double(count == 0 ? numeric_limits<double>::quiet_NaN() : my_gen_stats.sibship().max(), 
                 FIELD_SIZE - 1, 0, 0);
  my_o << ") "                          << "|"
                                        << endl;
                                
  my_o.flags(old_flags);
}

void
general_stats_viewer::sub_loop_line() const
{                                           
  ios_base::fmtflags old_flags = my_o.flags();
  
  my_o << "|" << std::setw(FIELD_SIZE_PED + 2) << left   << "Pedigrees"                       << "|"
              << std::setw(FIELD_SIZE)         << right  << my_gen_stats.total_subpedigrees() << "||"
              << std::setw(FIELD_SIZE_PED + 1) << left   << " Rings"                          << "|"
              << std::setw(FIELD_SIZE)         << right  << my_gen_stats.marriage_loops()     << "||"
              << std::setw(FIELD_SIZE_PED + 1) << left   << " Loops"                          << "|"
              << std::setw(FIELD_SIZE - 1)     << right  << my_gen_stats.non_marriage_loops() << "|"
              << std::endl;
              
  my_o.flags(old_flags);
}

// - Display pair and individual counts.
//
void
general_stats_viewer::pair_ind_data() const
{
  pair_ind_header();
  single_inner_line();
  
  size_t total = 0;
  for(int i = 0; i < TOTAL_PAIR_TYPES; ++i)
  {
    pair_half_line(static_cast<pg::pair_type>(i + 1));
    switch(i)
    {
      case 0: case 1: case 2:               // Individual counts male, female, unknown.
      case 4:                               // Blank.
      case 5: case 6: case 7:               // Individual counts founder, non-founder, unconnected.
        total += ind_half_line(i);
        break;
      case 3: case 8:                       // Total count by gender, by founder status.
        total_ind_half_line(total);
        total = 0;
        break;
      default:
      ;
    }
    my_o << endl;
  }
}

void          
general_stats_viewer::pair_half_line(pg::pair_type t) const
{
  if(0 < t &&  t <= TOTAL_PAIR_TYPES)
  {
    ios_base::fmtflags old_flags = my_o.flags();
    
    my_o << "|" << setw(FIELD_SIZE_PED * 2) << left << pair_label(t) << "|"
                                            << right;
 
    size_t count = my_gen_stats.pairs(t);
    my_o << setw(FIELD_SIZE_PED) << count;
    
    my_o.flags(old_flags);
  }
}

size_t          
general_stats_viewer::ind_half_line(int i) const
{
  ios_base::fmtflags old_flags = my_o.flags();
  
  my_o << "||" << setw(FIELD_SIZE_PED * 2);
  size_t count = 0;
  
  switch(i)
  {
    case 0:
      my_o << left << "  Male" << "|"
           << right;
      count = my_gen_stats.male_count();
      my_o << setw(FIELD_SIZE_PED) << count << "|";
      break;
    case 1:
      my_o << left << "  Female" << "|"
           << right;
      count = my_gen_stats.female_count();
      my_o << setw(FIELD_SIZE_PED) << count << "|";
      break;
    case 2:
      my_o << left << "  Unknown" << "|"
           << right;
      count = my_gen_stats.unknown_sex_count();
      my_o << setw(FIELD_SIZE_PED) << count << "|";
      break;
    case 4: 
      my_o << left << "\0" << "|"
           << right;
      my_o << setw(FIELD_SIZE_PED) << "\0" << "|";
      break;
    case 5:
      my_o << left << "  Founder" << "|"
           << right;
      count = my_gen_stats.founder_count();
      my_o << setw(FIELD_SIZE_PED) << count << "|";
      break;
    case 6:
      my_o << left << "  Non-founder" << "|"
           << right;
      count = my_gen_stats.nonfounder_count();
      my_o << setw(FIELD_SIZE_PED) << count << "|";
      break;
    case 7:
      my_o << left << "  Singleton" << "|"
           << right;
      count = my_gen_stats.unconnecteds();
      my_o << setw(FIELD_SIZE_PED) << count << "|";
      break;
    default:
    ;
  }
  
  my_o.flags(old_flags);
  return count;
}

void          
general_stats_viewer::total_ind_half_line(size_t total) const
{
  ios_base::fmtflags old_flags = my_o.flags();
  
  my_o << "||" << setw(FIELD_SIZE_PED * 2) << left   << "        Total" << "|"
                                              << right
              << setw(FIELD_SIZE_PED) << total << "|";
  
  my_o.flags(old_flags);
}

void
general_stats_viewer::multiple_mates() const
{
  ios_base::fmtflags old_flags = my_o.flags();

  header_line("Individuals With Multiple Mates");
  
  if(my_gen_stats.multiple_mates().empty())
  {
    my_o << left << setw(LINE_SIZE - 1) << "|none" << "|" << endl;
  }
  else
  {
    my_o << left << setw(LINE_SIZE - 1) << "|(Pedigree, Individual)   Mates" << "|" << endl; 
    single_inner_line();
         
    mate_map::const_iterator  mm_iter;
    for(mm_iter = my_gen_stats.multiple_mates().begin(); 
        mm_iter != my_gen_stats.multiple_mates().end(); 
        ++mm_iter)
    {
      ostringstream  mate_line;
      mate_line << left << "|" << setw(24) << mm_iter->first.name() << " ";
      
      size_t  mates_in_line = 0;
      size_t  count = mm_iter->second.size();
      for(size_t i = 0; i < count; ++i)
      {
        string  temp = (mm_iter->second)[i] + ((i < count - 1) ? ", " : "");
        if(mates_in_line < 1 || wrap_not_needed(mate_line, temp))
        {
          mate_line << temp;
          ++mates_in_line;
        }
        else    // wrap line
        {
          my_o << setw(LINE_SIZE - 1) << mate_line.str() << "|" << endl;
          mate_line.str(""); 
          mate_line << left << setfill(' ') << setw(26) << "|";
          mate_line << temp;
          mates_in_line = 1; 
        }
      }
      
      my_o << setw(LINE_SIZE - 1) << mate_line.str() << "|" << endl;
    }
  }
  
  my_o.flags(old_flags);  
}

// - Can temp fit on mate_line or should it be wrapped?
//
bool
general_stats_viewer::wrap_not_needed(const ostringstream& mate_line, const string& temp)
{
  return  mate_line.str().size() + temp.size() <= LINE_SIZE - 1;
}

void
general_stats_viewer::cons_pairs() const
{
  ios_base::fmtflags old_flags = my_o.flags();

  pp_map  pairs = my_gen_stats.cons_pairs();

  header_line("Consanguineous Mating Pairs");
  
  no_pairs np = for_each(pairs.begin(), pairs.end(), no_pairs());
  if(np.my_no_pairs)
  {
    my_o << left << setw(LINE_SIZE - 1) << "|none" << "|" << endl;
  }
  else
  {
    my_o << left << setw(LINE_SIZE - 1) << "|Pedigree         Pair" << "|" << endl;
    single_inner_line();
    
    pp_map::const_iterator  pp_iter;
    for(pp_iter = pairs.begin(); pp_iter != pairs.end(); ++pp_iter)
    {
      for(size_t i = 0; i < pp_iter->second.size(); ++i)
      {
        ostringstream  cons_mate_line;
        cons_mate_line << left << "|" << setw(17) << pp_iter->first
                                          << pp_iter->second[i].p1 << ", "
                                          << pp_iter->second[i].p2;
        my_o << setw(LINE_SIZE - 1) << cons_mate_line.str() << "|" << endl;                                   
      }
    }
  }
  
  my_o.flags(old_flags);    
}


//============================================================================
// IMPLEMENTATION:  mp_stats_viewer
//============================================================================
//
void
mp_stats_viewer::view(bool each, bool suppress_general) const
{
  if(! suppress_general)
  {
    
    // Display data specific to MP_stats object.
    major_header(build_title());
    ped_size_header();
    single_inner_line();
    ped_size_line();
    histograms_header();
    histograms();
    
    // - Display data which MP_stats object has in common w. a Ped_stats object.
    //
    summary_stats_viewer   s_view(my_o, my_mp_stats);
    s_view.view();
    double_inner_line();
    s_view.multiple_mates();
    double_inner_line();
    s_view.cons_pairs();
    outer_line();
  }

  // Display MP_stats trait data.
  const vector<Binary_trait_stats>& bt = my_mp_stats.Bt_stats();
  for(size_t stat = 0; stat < bt.size(); ++stat)
  {
    binary_stats_viewer b_view(my_o, bt[stat]);
    b_view.view();
  }
  
  const vector<Cont_trait_stats>& ct = my_mp_stats.Ct_stats();
  for(size_t stat = 0; stat < ct.size(); ++stat)
  {
    cont_stats_viewer ct_view(my_o, ct[stat]);
    ct_view.view();
  }
  
  const vector<Cmpd_trait_stats>& cd = my_mp_stats.Cd_stats();
  for(size_t stat = 0; stat < cd.size(); ++stat)
  {
    cmpd_stats_viewer cd_view(my_o, cd[stat]);
    cd_view.view();
  }

  // Display data for each pedigree.
  if(each)
  {
    for(size_t p = 0; p < my_mp_stats.size(); ++p)
    {
      if(! suppress_general)
      {
        ped_stats_viewer p_view(my_o, *(my_mp_stats[p]), my_mp_stats[p]->pedigree_name());
        p_view.view();
      }
      
      // Display individual trait data.
      const vector<Binary_trait_stats>& bt = my_mp_stats[p]->Bt_stats();
      for(size_t stat = 0; stat < bt.size(); ++stat)
      {
        binary_stats_viewer b_view(my_o, bt[stat]);
        b_view.view();
      }
      
      const vector<Cont_trait_stats>& ct = my_mp_stats[p]->Ct_stats();
      for(size_t stat = 0; stat < ct.size(); ++stat)
      {
        cont_stats_viewer ct_view(my_o, ct[stat]);
        ct_view.view();
      }
      
      const vector<Cmpd_trait_stats>& cd = my_mp_stats[p]->Cd_stats();
      for(size_t stat = 0; stat < cd.size(); ++stat)
      {
        cmpd_stats_viewer cd_view(my_o, cd[stat]);
        cd_view.view();
      }
    }
  }
}

void
mp_stats_viewer::ped_size_line() const
{
  ios_base::fmtflags old_flags = my_o.flags();
  
  my_o << "|" << setw(FIELD_SIZE + FIELD_SIZE_CORR) << left  << "Pedigrees "  << "|";
  size_t count = my_mp_stats.pedigree_size_info().count();
  my_o        << setw(FIELD_SIZE - 1)                       << right << count        << "|" << " ";
  
  display_double(my_mp_stats.pedigree_size_info().mean(), FIELD_SIZE - 1);
  my_o << " +/- ";
  display_double(my_mp_stats.pedigree_size_info().standard_deviation(), FIELD_SIZE - 1);
  my_o << " (";
  display_double(my_mp_stats.pedigree_size_info().min(), FIELD_SIZE - 1, 0, 0);
  my_o << ", ";
  display_double(my_mp_stats.pedigree_size_info().max(), FIELD_SIZE - 1, 0, 0);
  my_o << ") "                          << "|"
                                        << endl;
                                
  my_o.flags(old_flags);
}

// - Display histograms of # of pedigrees by generation size, number of families
//   # of inheritance vector bits.
//
void
mp_stats_viewer::histograms() const
{
  single_inner_line();
  
  const Histogram&     generations = my_mp_stats.generation_freq();
  const log_histogram  families(my_mp_stats.family_count_freq(), 2, 16);
  const log_histogram  vector_bits(my_mp_stats.likelihood_bits_freq(), 2, 16);
  Histogram  test_h;
  
  // - Determine number of lines to be displayed.
  //
  size_t  g_max_size  = generations.size();
  size_t  f_max_size  = families.size();
  size_t  v_max_size  = vector_bits.size();
  size_t  g_offset    = 0;
  size_t  f_offset    = families.offset();
  size_t  v_offset    = vector_bits.offset();
  
  // Determine value of g_offset.
  for(size_t i = 0; i < g_max_size; ++i)
  {
    if(generations[i] > 0)
    {
      g_offset = i;
      break;
    }
  }
  
  size_t  temp        = max(g_max_size - g_offset, f_max_size - f_offset);
  size_t  line_count  = max(temp, v_max_size - v_offset);
  
  // Display.
  ios_base::fmtflags old_flags = my_o.flags();
  my_o << right;
  
  for(size_t line = 0; line < line_count; ++line)
  {
    size_t  g_index = line + g_offset;
  
    // Generational statistics.
    if(g_index < g_max_size)
    {
      if(g_index == 0)
      {
        // - Pedigrees w. non-marriage loops have 0 generations because number
        //   of generations is potentially undefined.  This code added to 
        //   indicate that the number of generations in this case is undetermined.
        //   6-25-3, -djb
        //
        my_o << "|" << setw(FIELD_SIZE)     << "undet."             << "|"; 
      }
      else
      {
        my_o << "|" << setw(FIELD_SIZE)     << g_index              << "|";
      }
      
      my_o << setw(FIELD_SIZE) << generations[g_index]  << "||";  
    }
    else
    {
      my_o << "|" << setw(FIELD_SIZE)     << "\0"                   << "|"
           << setw(FIELD_SIZE) << "\0"                          << "||";  
    }
    
    // Nuclear family statistics.
    size_t  f_index = line + f_offset;
    if(f_index < f_max_size)
    {
      ostringstream bin_label;
      bin_label << right;
      if(f_index == f_max_size - 1 && families.limited())        // Last line of limited log_histogram.
      {
        bin_label << ">= " << setw(5) << families.boundaries(f_index).first;
      }
      else
      {
        bin_label << setw(5) << families.boundaries(f_index).first << " - " 
                  << setw(5) << families.boundaries(f_index).second;
      }
      my_o << setw(FIELD_SIZE_BIN)  << bin_label.str()   << "|"
           << setw(FIELD_SIZE - 1)      << families[f_index]    << "||";
    }
    else
    {
      my_o << setw(FIELD_SIZE_BIN)  << "\0"              << "|"
           << setw(FIELD_SIZE - 1)      << "\0"              << "||";  
    }
    
    // Inheritance vector bit statistics.
    size_t  v_index = line + v_offset;
    if(v_index < v_max_size)
    {
      ostringstream bin_label;
      bin_label << right;
      if(v_index == v_max_size - 1 && vector_bits.limited())     // Last line of limited log_histogram.
      {
        bin_label << ">= " << setw(5) << vector_bits.boundaries(v_index).first;
      }
      else
      {
        bin_label << setw(5) << vector_bits.boundaries(v_index).first << " - " 
                  << setw(5) << vector_bits.boundaries(v_index).second;
      }
      my_o << setw(FIELD_SIZE_BIN)  << bin_label.str()   << "|"
           << setw(FIELD_SIZE)      << vector_bits[v_index] << "|";  
    }
    else
    {
      my_o << setw(FIELD_SIZE_BIN)  << "\0"              << "|"
           << setw(FIELD_SIZE)      << "\0"              << "|";  
    }
    my_o << endl;
  }
  
  my_o.flags(old_flags);
}

//============================================================================
// IMPLEMENTATION:  base_trait_stats_viewer
//============================================================================
//
void
base_trait_stats_viewer::sib_line() const
{
  ios_base::fmtflags old_flags = my_o.flags();
  
  my_o << "|" << setw(FIELD_SIZE + FIELD_SIZE_CORR) << left  << "Sibships"  << "|";
  size_t count = my_base_trait_stats.sibship_size_info().count();
  my_o        << setw(FIELD_SIZE - 1)                       << right << count        << "|" << " ";
  
  display_double(count == 0 ? numeric_limits<double>::quiet_NaN() : my_base_trait_stats.sibship_size_info().mean(), 
                 FIELD_SIZE - 1);
  my_o << " +/- ";
  display_double(count == 0 ? numeric_limits<double>::quiet_NaN() : my_base_trait_stats.sibship_size_info().standard_deviation(), 
                 FIELD_SIZE - 1);
  my_o << " (";
  display_double(count == 0 ? numeric_limits<double>::quiet_NaN() : my_base_trait_stats.sibship_size_info().min(), 
                 FIELD_SIZE - 1, 0, 0);
  my_o << ", ";
  display_double(count == 0 ? numeric_limits<double>::quiet_NaN() : my_base_trait_stats.sibship_size_info().max(), 
                 FIELD_SIZE - 1, 0, 0);
  my_o << ") "                          << "|"
                                        << endl;
                                
  my_o.flags(old_flags);
}

void
base_trait_stats_viewer::ped_size_line() const
{
  ios_base::fmtflags old_flags = my_o.flags();
  
  my_o << "|" << setw(FIELD_SIZE + FIELD_SIZE_CORR) << left  << "Pedigrees "  << "|";
  size_t count = my_base_trait_stats.pedigree_size_info().count();
  my_o        << setw(FIELD_SIZE - 1)                       << right << count        << "|" << " ";
  
  display_double(my_base_trait_stats.pedigree_size_info().mean(), FIELD_SIZE - 1);
  my_o << " +/- ";
  display_double(my_base_trait_stats.pedigree_size_info().standard_deviation(), FIELD_SIZE - 1);
  my_o << " (";
  display_double(my_base_trait_stats.pedigree_size_info().min(), FIELD_SIZE - 1, 0, 0);
  my_o << ", ";
  display_double(my_base_trait_stats.pedigree_size_info().max(), FIELD_SIZE - 1, 0, 0);
  my_o << ") "                          << "|"
                                        << endl;
                                
  my_o.flags(old_flags);
}

void
base_trait_stats_viewer::nuclear_families() const
{
  // - header.
  //
  double_inner_line();
  my_o << "|" << "  Nuc  Family Statistics "   << "  "
              << "                     "       << "  "
              << "                          "  << "|"
              << endl;
  my_o << "|" << " # of Nuc Fams |# of Peds"   << "  "
              << "                     "       << "  " 
              << "                          "  << "|"
              << endl;
              
  single_inner_line();
  
  // - data.
  //
  const log_histogram  families(my_base_trait_stats.family_count_freq(), 2, 16);
  
  ios_base::fmtflags old_flags = my_o.flags();
  my_o << right;
  
  size_t  f_size   = families.size();
  size_t  f_offset = families.offset();
  
  for(size_t line = 0; line < f_size - f_offset; ++line)
  {
    size_t  f_index = line + f_offset;
    
    my_o << "|";
  
    ostringstream bin_label;
    bin_label << right;
    
    if(f_index == f_size - 1 && families.limited())        // Last line of limited log_histogram.
    {
      bin_label << ">= " << setw(5) << families.boundaries(f_index).first;
    }
    else
    {
      bin_label << setw(5) << families.boundaries(f_index).first << " - " 
                << setw(5) << families.boundaries(f_index).second;
    }
    
    my_o << setw(FIELD_SIZE_BIN)  << bin_label.str()   << "|"
         << setw(FIELD_SIZE - 1)  << families[f_index]    << "  " 
         << setw(78 - (4 + FIELD_SIZE_BIN + (FIELD_SIZE - 1))) << "|";
    
    my_o << endl;
  }
  
  my_o.flags(old_flags);
}

//============================================================================
// IMPLEMENTATION:  log_histogram
//============================================================================
//
// - Create a log histogram by populating bins w. ranges 0-10, 11-100, 101-1000, etc., for
//   base 10 for example, w.values obtained by consolidating bin counts from a "conventional" histogram.
//
//   prelim_size() returns the number of bins necessary to match the range of the
//   conventional histogram.  If this exceeds the size limit, the last bin is made
//   as big as necessary to account for all of the values in the conventional histo-
//   gram.
//
log_histogram::log_histogram(const Histogram& histogram, unsigned int base, size_t size_limit)
      : my_base(base),
        my_size(min(prelim_size(histogram), size_limit)), 
        my_limited(prelim_size(histogram) > size_limit)
{
  // Fill bins.
  for(size_t lbin = 0; lbin < my_size; ++lbin)
  {
    my_counts.push_back(0);
    for(size_t hbin = boundaries(lbin).first; 
        hbin <= boundaries(lbin).second && hbin < histogram.size(); ++hbin)
    { 
      my_counts[lbin] += histogram[hbin];
    }
  }
}

// - Calculate number of bins necessary for log_histogram to accomodate the range of 
//   the conventional histogram by successively comparing powers of the base 
//   to the largest count.
//   Note that for the histogram bin[0] -> [0,1), bin[1] -> [1,2), etc. while for
//   the log_histogram bin[0] -> [0,2], bin[1] -> (2,4], etc. 
//
size_t  
log_histogram::prelim_size(const Histogram& histogram) const
{
  size_t  num = histogram.size() - 1;
  
  size_t  bin_count = 1;
  while(static_cast<size_t>(std::pow(static_cast<double>(my_base), static_cast<double>(bin_count))) < num)
  {
    ++bin_count;
  }
         
  return  bin_count;
}

// - Returns number of bins before first non-empty bin.
//
size_t
log_histogram::offset() const
{
  size_t off = 0;
  for(size_t i = 0; i < my_counts.size(); ++i)
  {
    if(my_counts[i] > 0)
    {
      off = i;
      break;
    }
  } 
  return off;
} 

} // End namespace RPED
} // End namespace SAGE


