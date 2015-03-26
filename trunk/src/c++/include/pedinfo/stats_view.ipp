//============================================================================
// File:     stats_view.ipp
//                                                                          
// Author:   Dan Baechle
//                                                                          
// History:  Initial version: 10/00
//                                                                          
// Notes:    With stats_view.h supercedes stat_view.h
// 
//           Inline implementation of the following classes -
//              base_stats_viewer
//              base_trait_stats_viewer
//              binary_stats_viewer
//              cont_stats_viewer
//              cmpd_stats_viewer
//              general_stats_viewer
//              ped_stats_viewer
//              summary_stats_viewer
//              mp_stats_viewer
//              log_histogram
//                                                                          
// Copyright (c) 2000 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

namespace SAGE {
namespace RPED {


//============================================================================
// IMPLEMENTATION:  base_stats_viewer
//============================================================================
//
inline
base_stats_viewer::base_stats_viewer(std::ostream& o)
      : my_o(o)
{}
  
inline void
base_stats_viewer::outer_line() const
{
  my_o << std::setw(LINE_SIZE) << std::setfill('=') << "" << std::endl;
  my_o << std::setfill(' ');
}

inline void
base_stats_viewer::single_inner_line() const
{
  my_o << "|" << std::setw(LINE_SIZE - 2) << std::setfill('-') << "" << "|" << std::endl;
  my_o << std::setfill(' ');
}

inline void
base_stats_viewer::double_inner_line() const
{
  my_o << "|" << std::setw(LINE_SIZE - 2) << std::setfill('=') << "" << "|" << std::endl;
  my_o << std::setfill(' ');
}

inline std::string
base_stats_viewer::pair_label(pg::pair_type type) const
{
  switch(type)
  {
    case pg::PARENTAL:
      return "Parent/Off";
    case pg::SIBSIB:
      return "Sib/Sib";
    case pg::SISSIS:
      return "Sis/Sis";
    case pg::BROBRO:
      return "Bro/Bro";
    case pg::BROSIS:
      return "Bro/Sis";
    case pg::GRANDP:
      return "Grandp.";
    case pg::AVUNC:
      return "Avunc.";
    case pg::HALFSIB:
      return "Half Sib";
    case pg::COUSIN:
      return "Cousin";
    default:
      return "";
  }
}

inline void
base_stats_viewer::sib_header() const
{
  double_inner_line();
  ios_base::fmtflags old_flags = my_o.flags();
  
  my_o << "|" << std::setw(FIELD_SIZE + FIELD_SIZE_CORR) << left  << ""     << "|"
              << std::setw(FIELD_SIZE - 1)                       << right << "Count"  << "|"
              << " Mean Size +/- Std. Dev. (     Min.,      Max.) "                 << "|"
              << std::endl;
  
  my_o.flags(old_flags);
}

//============================================================================
// IMPLEMENTATION:  base_trait_stats_viewer
//============================================================================
//
inline
base_trait_stats_viewer::base_trait_stats_viewer(std::ostream& o, 
                                         const Base_trait_stats& base_trait_data)
      : base_stats_viewer(o), my_base_trait_stats(base_trait_data)
{}
    
inline void          
base_trait_stats_viewer::view() const
{
  double_inner_line();
  header();
  data_line();
  
  sib_header();
  single_inner_line();
  sib_line();
  
  if(my_base_trait_stats.mped_member())
  {
    single_inner_line();
    ped_size_line();
    
    nuclear_families();
  }
}
    
inline void          
base_trait_stats_viewer::header() const
{
  ios_base::fmtflags old_flags = my_o.flags();

  my_o << "|" << std::setw(FIELD_SIZE_LABEL) << left  << ""                 << "|"
              << std::setw(FIELD_SIZE_PED) << right << "0 Parents w. Data"    << "|"
              << std::setw(FIELD_SIZE_PED)          << "1 Parent w. Data"     << "|"
              << std::setw(FIELD_SIZE_PED)          << "2 Parents w. Data"    << "|"
              << std::endl;
  
  my_o.flags(old_flags);
}

inline void          
base_trait_stats_viewer::data_line() const
{
  ios_base::fmtflags old_flags = my_o.flags();
  single_inner_line();
  
  my_o << "|" << std::setw(FIELD_SIZE_LABEL) << left  << "Sibships"               << "|"
              << std::setw(FIELD_SIZE_PED) << right << my_base_trait_stats.sibship(0) << "|"
              << std::setw(FIELD_SIZE_PED) << right << my_base_trait_stats.sibship(1) << "|"
              << std::setw(FIELD_SIZE_PED) << right << my_base_trait_stats.sibship(2) << "|"
              << std::endl;
  
  my_o.flags(old_flags);
}
  

//============================================================================
// IMPLEMENTATION:  binary_stats_viewer
//============================================================================
//
inline
binary_stats_viewer::binary_stats_viewer(std::ostream& o, const Binary_trait_stats& binary_stats)
      : base_stats_viewer(o), my_binary_stats(binary_stats), my_base_trait_stats_viewer(o, binary_stats)
{}

inline void
binary_stats_viewer::view() const
{
  my_o << "\n\n";
  std::string title = build_title();
  major_header(title);
  my_base_trait_stats_viewer.view();
  pair_data();
  ind_data();
  outer_line();
}

inline std::string   
binary_stats_viewer::build_title() const
{
  std::string temp;
  std::string name;
  
  name = my_binary_stats.pedigree_name();
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
  temp += ", Trait - ";
  temp += my_binary_stats.trait_name();
  return temp;
}

inline void     
binary_stats_viewer::ind_header() const
{
  ios_base::fmtflags old_flags = my_o.flags();
  double_inner_line();

  my_o << "|" << std::setw(FIELD_SIZE) << left  << "Indiv.'s"    << "|"
              << std::setw(FIELD_SIZE) << right << "Affected"    << "|"
              << std::setw(FIELD_SIZE)          << "Unaff."      << "|"
              << std::setw(FIELD_SIZE)          << "Missing"     << "|"
              << std::setw(FIELD_SIZE)          << "Total"       << "|"
              << std::setw(FIELD_SIZE)          << ""          << "|"
              << std::setw(FIELD_SIZE)          << ""          << "|"
              << std::endl;
  
  my_o.flags(old_flags);
}

inline std::string
binary_stats_viewer::gender_label(bt::gender g) const
{
  switch(g)
  {
    case bt::MALE:
      return "Male";
    case bt::FEMALE:
      return "Female";
    case bt::UNKNOWN:
      return "Unknown";
    default:
      return "";
  }
}

inline std::string
binary_stats_viewer::founder_status_label(bt::founder_status f) const
{
  switch(f)
  {
    case bt::FOUNDER:
      return "Founder";
    case bt::NON_FOUNDER:
      return "Nonfound.";
    case bt::UNCONNECTED:
      return "Singleton";
    default:
      return "";
  }
}

inline void     
binary_stats_viewer::pair_header() const
{
  double_inner_line();
  ios_base::fmtflags old_flags = my_o.flags();
  
  // Line one.
  my_o << "|" << std::setw(FIELD_SIZE) << left  << "Pairs"         << "|"
              << std::setw(FIELD_SIZE) << right << "Concord."      << "|"
              << std::setw(FIELD_SIZE)          << "Discord."      << "|"
              << std::setw(FIELD_SIZE)          << "Concord."      << "|"
              << std::setw(FIELD_SIZE)          << "Uninform."     << "|"
              << std::setw(FIELD_SIZE)          << "Total"         << "|"
              << std::setw(FIELD_SIZE)          << "Corr."         << "|"
              << std::endl;
              
  // Line two.
  my_o << "|" << std::setw(FIELD_SIZE) << left  << ""            << "|"
              << std::setw(FIELD_SIZE) << right << "Unaff."        << "|"
              << std::setw(FIELD_SIZE)          << ""            << "|"
              << std::setw(FIELD_SIZE)          << "Aff."          << "|"
              << std::setw(FIELD_SIZE)          << ""            << "|"
              << std::setw(FIELD_SIZE)          << ""            << "|"
              << std::setw(FIELD_SIZE)          << ""            << "|"
              << std::endl;
  
  my_o.flags(old_flags);
}

//============================================================================
// IMPLEMENTATION:  cont_stats_viewer
//============================================================================
//
inline
cont_stats_viewer::cont_stats_viewer(std::ostream& o, const Cont_trait_stats& cont_stats)
      : base_stats_viewer(o), my_cont_stats(cont_stats), my_base_trait_stats_viewer(o, cont_stats)
{}

inline void
cont_stats_viewer::view() const
{
  my_o << "\n\n";
  std::string title = build_title();
  major_header(title);
  my_base_trait_stats_viewer.view();
  
  // - Choose normal or detailed version.
  //
  //#define DETAILED
  #ifdef DETAILED
    pair_data_alt();       
  #else
    pair_data(); 
  #endif
                           
  ind_data();
  outer_line();
}

inline std::string   
cont_stats_viewer::build_title() const
{
  std::string temp;
  std::string name;
  
  name = my_cont_stats.pedigree_name();
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
  temp += ", Trait - ";
  temp += my_cont_stats.trait_name();
  return temp;
}

inline void     
cont_stats_viewer::ind_header() const
{
  double_inner_line();
  ios_base::fmtflags old_flags = my_o.flags();
  
  my_o << "|" << std::setw(FIELD_SIZE - 1) << left  << "Indiv.'s"      << "|"
              << std::setw(FIELD_SIZE - 1) << right << "Count"         << "|"
              << "      Mean +/- Std. Dev. (     Min.,      Max.) "  << "|"
              << std::setw(FIELD_SIZE_CORR)     << ""            << "|"
              << std::endl;
              
  my_o.flags(old_flags);
}

inline std::string
cont_stats_viewer::gender_label(ct::gender g) const
{
  switch(g)
  {
    case ct::MALE:
      return "Male";
    case ct::FEMALE:
      return "Female";
    case ct::UNKNOWN:
      return "Unknown";
    case ct::ALL:
      return "    All";
    default:
      return "";
  }
}

inline std::string
cont_stats_viewer::founder_status_label(ct::founder_status f) const
{
  switch(f)
  {
    case ct::FOUNDER:
      return "Founder";
    case ct::NON_FOUNDER:
      return "Nonfound.";
    case ct::UNCONNECTED:
      return "Singleton";
    case ct::ALL:
      return "    All";
    default:
      return "";
  }
}

inline void     
cont_stats_viewer::pair_header() const
{
  double_inner_line();
  ios_base::fmtflags old_flags = my_o.flags();
  
  
  
  my_o << "|" << std::setw(FIELD_SIZE + 1) << left  << "\0"      << "|"
              << std::setw(FIELD_SIZE - 1) << right << "\0"         << "|"
              << std::setw(FIELD_SIZE_CORR + FIELD_SIZE - 2) << "Pair       " << "||"
              << std::setw(FIELD_SIZE + 1) << left  << "\0"      << "|"
              << std::setw(FIELD_SIZE - 1) << right << "\0"         << "|"
              << std::setw(FIELD_SIZE_CORR + FIELD_SIZE - 2) << "Pair       " << "|"
              << endl;
  
  my_o << "|" << std::setw(FIELD_SIZE + 1) << left  << "Pairs"      << "|"
              << std::setw(FIELD_SIZE - 1) << right << "Count"         << "|"
              << std::setw(FIELD_SIZE_CORR + FIELD_SIZE - 2) << "Correlation" << "||"
              << std::setw(FIELD_SIZE + 1) << left  << "Pairs"      << "|"
              << std::setw(FIELD_SIZE - 1) << right << "Count"         << "|"
              << std::setw(FIELD_SIZE_CORR + FIELD_SIZE - 2) << "Correlation" << "|"
              << endl;
              
  my_o.flags(old_flags);
}

inline void     
cont_stats_viewer::pair_header_alt() const
{
  double_inner_line();
  ios_base::fmtflags old_flags = my_o.flags();
  
  my_o << "|" << std::setw(FIELD_SIZE - 1) << left  << "Pairs"      << "|"
              << std::setw(FIELD_SIZE - 1) << right << "Count"         << "|"
              << "      Mean +/- Std. Dev. (     Min.,      Max.) "  << "|"
              << std::setw(FIELD_SIZE_CORR)     << "Corr."            << "|"
              << std::endl;
              
  my_o.flags(old_flags);
}

//============================================================================
// IMPLEMENTATION:  cmpd_stats_viewer
//============================================================================
//
inline
cmpd_stats_viewer::cmpd_stats_viewer(std::ostream& o, const Cmpd_trait_stats& cmpd_stats)
      : base_stats_viewer(o), my_cmpd_stats(cmpd_stats), 
        my_base_trait_stats_viewer(o, cmpd_stats)  
{}

inline void
cmpd_stats_viewer::pair_ind_header() const
{
  double_inner_line();
  ios_base::fmtflags old_flags = my_o.flags();
  
  my_o << "|" << std::setw(FIELD_SIZE_PED * 2)      << left   << "Pairs"          << "|"
              << std::setw(FIELD_SIZE_PED)          << right  << "Count"          << "||"
              << std::setw(FIELD_SIZE_PED * 2    )  << left   << "  Individuals"  << "|"
              << std::setw(FIELD_SIZE_PED)          << right  << "Count"          << "|"
                                                              << std::endl;
                                                            
  my_o.flags(old_flags);
}

//============================================================================
// IMPLEMENTATION:  general_stats_viewer
//============================================================================
//
inline
general_stats_viewer::general_stats_viewer(std::ostream& o, const General_Stats& gen_stats)
      : base_stats_viewer(o), my_gen_stats(gen_stats)
{}



inline void
general_stats_viewer::sub_loop_header() const
{
  double_inner_line();
  ios_base::fmtflags old_flags = my_o.flags();

  my_o << "|" << std::setw(FIELD_SIZE_PED + 2) << left   << "Constituent"  << "|"
              << std::setw(FIELD_SIZE)         << right  << ""             << "||"
              << std::setw(FIELD_SIZE_PED + 1) << left   << " Marriage"    << "|"
              << std::setw(FIELD_SIZE)         << right  << ""             << "||"
              << std::setw(FIELD_SIZE_PED + 1) << left   << ""             << "|"
              << std::setw(FIELD_SIZE - 1)     << right  << ""             << "|"
              << std::endl;
              
  my_o.flags(old_flags);
}

inline void
general_stats_viewer::pair_ind_header() const
{
  double_inner_line();
  ios_base::fmtflags old_flags = my_o.flags();
  
  my_o << "|" << std::setw(FIELD_SIZE_PED * 2)      << left   << "Pairs"          << "|"
              << std::setw(FIELD_SIZE_PED)          << right  << "Count"          << "||"
              << std::setw(FIELD_SIZE_PED * 2)  << left   << "  Individuals"  << "|"
              << std::setw(FIELD_SIZE_PED)          << right  << "Count"          << "|"
                                                              << std::endl;
                                                            
  my_o.flags(old_flags);
}

//============================================================================
// IMPLEMENTATION:  ped_stats_viewer
//============================================================================
//
inline
ped_stats_viewer::ped_stats_viewer(std::ostream& o, const Ped_stats& ped_stats, std::string ped_name)
      : general_stats_viewer(o, ped_stats), my_ped_name(ped_name)
{}

inline void
ped_stats_viewer::view() const
{
  general_stats_viewer::view(build_title());
  outer_line();
}

inline std::string   
ped_stats_viewer::build_title() const
{
  std::string temp;
  temp += "General Statistics: ";
  temp += "Pedigree - ";
  temp += my_ped_name;
  return temp;
}

//============================================================================
// IMPLEMENTATION:  summary_stats_viewer
//============================================================================
//
inline
summary_stats_viewer::summary_stats_viewer(std::ostream& o, const MP_stats& multi_stats)
      : general_stats_viewer(o, multi_stats)
{}

inline void
summary_stats_viewer::view() const
{
  general_stats_viewer::view("");
}



//============================================================================
// IMPLEMENTATION:  mp_stats_viewer
//============================================================================
//
inline
mp_stats_viewer::mp_stats_viewer(std::ostream& o, const MP_stats& multi_stats)
      : base_stats_viewer(o), my_mp_stats(multi_stats)
{}

inline void
mp_stats_viewer::ped_size_header() const
{
  double_inner_line();
  ios_base::fmtflags old_flags = my_o.flags();
  
  my_o << "|" << std::setw(FIELD_SIZE + FIELD_SIZE_CORR) << left  << ""     << "|"
              << std::setw(FIELD_SIZE - 1)                       << right << "Count"  << "|" 
              << " Mean Size +/- Std. Dev. (     Min.,      Max.) "                 << "|"
              << std::endl;
  
  my_o.flags(old_flags);
}

inline void
mp_stats_viewer::histograms_header() const
{
  double_inner_line();
  my_o << "|" << "Generation Statistics"       << "||"
              << "  Nuc  Family Statistics "   << "||"
              << "  Inh  Vector Bit Stats   "  << "|"
              << endl;
  my_o << "|" << "# of Gens | # of Peds"       << "||" 
              << " # of Nuc Fams |# of Peds"   << "||"
              << "   # of Bits   |# of Peds "  << "|"
              << endl;
}

inline std::string   
mp_stats_viewer::build_title() const
{
  std::string temp;
  temp += "General Statistics: ";
  temp += "All Pedigrees";
  return temp;
}

//============================================================================
// IMPLEMENTATION:  log_histogram
//============================================================================
//
inline size_t                     
log_histogram::size() const
{
  return my_counts.size();
}

inline size_t                     
log_histogram::operator[](size_t i) const
{
  if(i < my_counts.size())
  {
    return my_counts[i];
  }
  else
  {
    return 0;
  }
}

// - Calculate boundaries for the bin w. the indicated index.
//
inline std::pair<size_t, size_t>  
log_histogram::boundaries(size_t i) const
{
  std::pair<size_t, size_t> temp;

  if(i < my_size)                           // Valid index.
  {
    if(i == 0)
    {
      temp.first = 0;
    }
    else
    {
      temp.first = static_cast<size_t>(std::pow(static_cast<double>(my_base), static_cast<double>(i))) + 1;
    }
                                           
    if(my_limited && i == (my_size - 1))    // Last bin of "limited" histogram.
    {
      temp.second = (size_t)(-1);
    }
    else
    {
      temp.second = static_cast<size_t>(std::pow(static_cast<double>(my_base), static_cast<double>(i + 1)));
    }
  }
  else
  {
    temp.first   = 0;
    temp.second  = 0;
  }
  
  return temp;
}

inline bool
log_histogram::limited() const
{
  return my_limited;
}

} // End namespace RPED
} // End namespace SAGE
