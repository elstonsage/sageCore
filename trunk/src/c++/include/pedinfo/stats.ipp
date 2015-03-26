//============================================================================
// File:      stats.ipp
//                                                                          
// Author:    
//                                                                          
// History:   8/00 - created from original stats.h file and extended to.
//                   to utilize pair generator class and to include traits.
//                   - Dan Baechle                                                    
//                                                                          
// Notes:     Inline implementation for the following classes:
//              General_Stats
//              Ped_stats
//              MP_stats
//              Base_trait_stats
//              Binary_trait_stats
//              Cont_trait_stats
//              Cmpd_trait_stats
//                                                                          
// Copyright (c) 2000 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

namespace SAGE {
namespace RPED {


inline
full_name::full_name(const string& ped, const string& ind)
      : ped_name(ped), ind_name(ind)
{}

inline string
full_name::name() const
{
  return  "(" + ped_name + ", " + ind_name + ")";
}

inline bool
operator<(const full_name& lhs, const full_name& rhs)
{
  if(lhs.ped_name != rhs.ped_name)
  {
    return  lhs.ped_name < rhs.ped_name; 
  }
  else
  {
    return  lhs.ind_name < rhs.ind_name;
  }
}

inline
parental_pair::parental_pair(const string& p_one, const string& p_two)
      : p1(p_one), p2(p_two)
{}

inline
no_pairs::no_pairs()
      : my_no_pairs(true)
{}

inline void
no_pairs::operator()(const pp_map::value_type& elem)
{
  my_no_pairs = my_no_pairs && elem.second.empty();
}

//============================================================================
// IMPLEMENTATION:  General_Stats
//============================================================================
//
inline double 
General_Stats::mean_sibship_size() const  
{ 
  return my_sibship.mean();    
}

inline double 
General_Stats::var_sibship_size() const  
{ 
  return my_sibship.variance();
}
  
inline size_t 
General_Stats::family_count() const 
{ 
  return my_sibship.count();   
}

inline size_t 
General_Stats::member_count() const  
{ 
  return my_total_inds;        
}
  
inline size_t 
General_Stats::female_count() const 
{ 
  return my_total_females;     
}

inline size_t 
General_Stats::male_count() const  
{ 
  return my_total_males;       
}

inline size_t 
General_Stats::unknown_sex_count() const  
{ 
  return my_total_unknown_sex; 
}

inline size_t 
General_Stats::founder_count() const  
{ 
  return my_total_founders;    
}

inline size_t 
General_Stats::nonfounder_count() const  
{ 
  return my_total_nonfounders; 
}

inline size_t 
General_Stats::generations() const  
{ 
  return my_num_gen;           
}

inline size_t 
General_Stats::sib_pairs() const  
{ 
  return my_pairs[pg::SIBSIB - 1];          
}

inline size_t 
General_Stats::hsib_pairs() const 
{ 
return my_pairs[pg::HALFSIB - 1];          
}

inline size_t 
General_Stats::avuncular_pairs() const  
{ 
  return my_pairs[pg::AVUNC - 1];          
}

inline size_t 
General_Stats::grandp_pairs() const 
{ 
  return my_pairs[pg::GRANDP - 1];          
}

inline size_t 
General_Stats::cousin_pairs() const  
{ 
  return my_pairs[pg::COUSIN - 1];          
}
  
inline size_t 
General_Stats::brother_pairs() const 
{ 
  return my_pairs[pg::BROBRO - 1];          
}

inline size_t 
General_Stats::sister_pairs() const  
{ 
  return my_pairs[pg::SISSIS - 1];          
}
  
inline size_t 
General_Stats::brother_sister_pairs() const 
{ 
  return my_pairs[pg::BROSIS - 1];         
}

inline size_t
General_Stats::parental_pairs() const
{
  return my_pairs[pg::PARENTAL - 1];
}

inline size_t 
General_Stats::marriage_loops() const  
{ 
  return my_marriage_loops;    
}
 
inline size_t 
General_Stats::non_marriage_loops() const  
{ 
  return my_non_marriage_loops;
}

inline size_t 
General_Stats::unconnecteds() const  
{ 
  return my_unconnecteds;      
}

inline size_t 
General_Stats::total_subpedigrees() const  
{ 
  return my_total_subpedigrees;
}

inline size_t 
General_Stats::sib_count() const 
{ 
  return (size_t)my_sibship.sum(); 
}

inline size_t
General_Stats::stats_type() const  
{ 
  return my_stats_type;    
}

inline size_t 
General_Stats::pairs(pg::pair_type t) const  
{
  if(t > TOTAL_PAIR_TYPES)   
  {
   return (size_t)(-1);
  }
  else
  {                           
   return my_pairs[t - 1];
  }
}

inline bool   
General_Stats::loops() const
{
  return ( marriage_loops() || non_marriage_loops() ) ? true : false;
}

inline bool   
General_Stats::valid() const  
{ 
  return my_valid;             
}

// - likelihood bits is maximum over _each_ subpedigree.
//
inline size_t            
General_Stats::likelihood_bits() const  
{ 
  return my_max_bits;          
}

inline const Histogram&  
General_Stats::sib_size_freq() const  
{ 
  return my_sibship_size_freq; 
}

inline const mate_map&
General_Stats::multiple_mates() const
{
  return my_multiple_mates;
}

inline const pp_map&
General_Stats::cons_pairs() const
{
  return my_cons_pairs;
}

inline const SampleInfo& 
General_Stats::sibship() const  
{ 
  return my_sibship;           
}

inline const std::vector<Binary_trait_stats>&  
General_Stats::Bt_stats() const
{
  return my_binary_trait_stats;
}

inline const std::vector<Cont_trait_stats>&    
General_Stats::Ct_stats() const
{
  return my_cont_trait_stats;
}

inline const std::vector<Cmpd_trait_stats>&    
General_Stats::Cd_stats() const
{
  return my_cmpd_trait_stats;
}

inline void   
General_Stats::invalidate()         
{ 
  my_valid = false;            
}

inline
General_Stats::~General_Stats() 
{}

inline
General_Stats::General_Stats(cerrormultistream &e) 
      : errors(e), my_stats_type(PEDIGREE_STATS) 
{ 
  init(); 
}

//============================================================================
// IMPLEMENTATION:  Ped_stats
//============================================================================
//
inline RefMultiPedigree::pedigree_const_pointer   
Ped_stats::ped() const 
{ 
  return my_pedigree; 
}

inline std::string
Ped_stats::pedigree_name() const
{
  if(my_pedigree != 0)
  {
    return my_pedigree->name();
  }
  else
  {
    return "";
  }
}

//============================================================================
// IMPLEMENTATION:  MP_stats
//============================================================================
//
inline const RefMultiPedigree* 
MP_stats::mped() const
{ 
  return my_multipedigree; 
}

inline const SampleInfo&  
MP_stats::pedigree_size_info() const 
{ 
  return my_pedigree_size_info; 
}

inline const SampleInfo&  
MP_stats::family_count_info() const 
{ 
  return my_family_count_info;  
}

inline size_t  
MP_stats::pedigree_count() const  
{ 
  return my_pedigree_size_info.count();   
}

inline size_t  
MP_stats::size() const  
{ 
  return my_ped_stats.size();        
}

inline const Histogram& 
MP_stats::likelihood_bits_freq() const 
{ 
  return my_likelihood_bits_freq;
}

inline const Histogram& 
MP_stats::family_count_freq() const 
{ 
  return my_family_count_freq;   
}

inline const Histogram& 
MP_stats::generation_freq() const 
{ 
  return my_generations_freq;    
}

inline const Ped_stats*  
MP_stats::operator[](size_t i) const
{
  if(i < pedigree_count())  
  {
    return &my_ped_stats[i]; 
  }
  else 
  {
    return NULL;
  }
}

//============================================================================
// IMPLEMENTATION:  Base_trait_stats
//============================================================================
//
inline bool
Base_trait_stats::mped_member() const
{
  return  mped_member_flag;
}

inline size_t 
Base_trait_stats::sibship(int parents) const
{
  if(0 <= parents && parents < 3)
  {
    return my_sibships[parents];
  }
  else
  {
    return (size_t)(-1);
  }
}

inline const SampleInfo&
Base_trait_stats::pedigree_size_info() const
{
  return  my_pedigree_size_info;
}

inline const SampleInfo&
Base_trait_stats::sibship_size_info() const
{
  return  my_sibship_size_info;
}

inline const Histogram&   
Base_trait_stats::family_count_freq() const
{
  return  my_family_count_freq;
}

inline    
Base_trait_stats::Base_trait_stats(const RefPedigree* p, cerrormultistream& e,
                                   bool mped_member)
      :my_pedigree(p), my_errors(e), mped_member_flag(mped_member)
{
  init();
}

//============================================================================
// IMPLEMENTATION:  Binary_trait_stats
//============================================================================
//
inline
Binary_trait_stats::Binary_trait_stats(size_t trait, cerrormultistream& e, std::string name,
                                       bool mped_member)
      : Base_trait_stats(0, e, mped_member), my_valid(false), my_trait(trait), my_trait_name(name)
{
  init();
}

inline    
Binary_trait_stats::Binary_trait_stats(const RefPedigree* p, size_t trait, cerrormultistream& e,
                                       bool mped_member)
      :Base_trait_stats(p, e, mped_member), my_valid(false), my_trait(trait), my_trait_name("")
{
  init();
  compute();
}

inline RefMultiPedigree::pedigree_const_pointer   
Binary_trait_stats::pedigree() const
{
  return my_pedigree;
}

inline size_t                                     
Binary_trait_stats::trait() const
{
  return my_trait;
}

inline std::string
Binary_trait_stats::pedigree_name() const
{
  if(my_pedigree != 0)
  {
    return my_pedigree->name();
  }
  else
  {
    return "";
  }
}

inline const std::string& 
Binary_trait_stats::trait_name() const
{
  return my_trait_name;
}

inline bool                                       
Binary_trait_stats::valid() const
{
  return my_valid;
}

inline size_t                                     
Binary_trait_stats::ind_count_gender(gender g, ind_affection_status st) const
{
  if(g < 0 || g > 2 || st < 0 || st > 2)
  {
    return (size_t)(-1);
  }
  else
  {
    return my_ind_counts_gender[g][st];
  }
}

inline size_t                                     
Binary_trait_stats::ind_count_founder(founder_status f, ind_affection_status st) const
{
  if(f < 0 || f > 2 || st < 0 || st > 2)
  {
    return (size_t)(-1);
  }
  else
  {
    return my_ind_counts_founder[f][st];
  }
}

inline size_t                                     
Binary_trait_stats::pair_count(pg::pair_type type, pft::affection_status st) const
{
  if(type < 1 || type > TOTAL_PAIR_TYPES || (st - 1) < 0 || (st - 1) > 3)
  {
    return (size_t)(-1);
  }
  else
  {
    return my_pair_counts[type - 1][st - 1];
  }
}

inline double    
Binary_trait_stats::correlation(pg::pair_type type) const
{
  if(type < 1 || type > TOTAL_PAIR_TYPES)
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
  else
  {
    return my_corr_infos[type - 1].correlation();
  }
}

//============================================================================
// IMPLEMENTATION:  Cont_trait_stats
//============================================================================
//
inline
Cont_trait_stats::Cont_trait_stats(size_t trait, cerrormultistream& e, std::string name,
                                   bool mped_member)
      :Base_trait_stats(0, e, mped_member), my_valid(false), my_trait(trait), my_trait_name(name)
{
  init();
}

inline    
Cont_trait_stats::Cont_trait_stats(const RefPedigree* p, size_t trait, cerrormultistream& e,
                                   bool mped_member)
      :Base_trait_stats(p, e, mped_member), my_valid(false), my_trait(trait), my_trait_name("")
{
  init();
  compute();
}

inline RefMultiPedigree::pedigree_const_pointer   
Cont_trait_stats::pedigree() const
{
  return my_pedigree;
}

inline size_t                                     
Cont_trait_stats::trait() const
{
  return my_trait;
}

inline std::string
Cont_trait_stats::pedigree_name() const
{
  if(my_pedigree != 0)
  {
    return my_pedigree->name();
  }
  else
  {
    return "";
  }
}

inline const std::string& 
Cont_trait_stats::trait_name() const
{
  return my_trait_name;
}

inline bool                                       
Cont_trait_stats::valid() const
{
  return my_valid;
}

inline size_t                                    
Cont_trait_stats::ind_gender_count(int g) const
{
  if(g >= 0 && g < 4)
  {
    return my_inds_gender[g].count();
  }
  else
  {
    return (size_t)(-1);
  }
}

inline double                                     
Cont_trait_stats::ind_gender_mean(int g) const
{
  if(g >= 0 && g < 4)
  {
    return (my_inds_gender[g].count() == 0 ? std::numeric_limits<double>::quiet_NaN() :
                                             my_inds_gender[g].mean());
  }
  else
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
}

inline double                                     
Cont_trait_stats::ind_gender_std_dev(int g) const
{
  if(g >= 0 && g < 4)
  {
    return (my_inds_gender[g].count() == 0 ? std::numeric_limits<double>::quiet_NaN() :
                                             my_inds_gender[g].standard_deviation());
  }
  else
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
}

inline double                                     
Cont_trait_stats::ind_gender_min(int g) const
{
  if(g >= 0 && g < 4)
  {
    return (my_inds_gender[g].count() == 0 ? std::numeric_limits<double>::quiet_NaN() :
                                             my_inds_gender[g].min());
  }
  else
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
}

inline double 
Cont_trait_stats::ind_gender_max(int g) const
{
  if(g >= 0 && g < 4)
  {
    return (my_inds_gender[g].count() == 0 ? std::numeric_limits<double>::quiet_NaN() :
                                             my_inds_gender[g].max());
  }
  else
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
}

inline size_t                                    
Cont_trait_stats::ind_founder_count(int f) const
{
  if(f >= 0 && f < 4)
  {
    return my_inds_founder[f].count();
  }
  else
  {
    return (size_t)(-1);
  }
}

inline double                                     
Cont_trait_stats::ind_founder_mean(int f) const
{
  if(f >= 0 && f < 4)
  {
    return (my_inds_founder[f].count() == 0 ? std::numeric_limits<double>::quiet_NaN() :
                                             my_inds_founder[f].mean());
  }
  else
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
}

inline double                                     
Cont_trait_stats::ind_founder_std_dev(int f) const
{
  if(f >= 0 && f < 4)
  {
    return (my_inds_founder[f].count() == 0 ? std::numeric_limits<double>::quiet_NaN() :
                                             my_inds_founder[f].standard_deviation());
  }
  else
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
}

inline double                                     
Cont_trait_stats::ind_founder_min(int f) const
{
  if(f >= 0 && f < 4)
  {
    return (my_inds_founder[f].count() == 0 ? std::numeric_limits<double>::quiet_NaN() :
                                             my_inds_founder[f].min());
  }
  else
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
}

inline double 
Cont_trait_stats::ind_founder_max(int f) const
{
  if(f >= 0 && f < 4)
  {
    return (my_inds_founder[f].count() == 0 ? std::numeric_limits<double>::quiet_NaN() :
                                             my_inds_founder[f].max());
  }
  else
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
}

inline size_t                                     
Cont_trait_stats::pair_count(pg::pair_type type) const
{
  if(type > 0 && type <= TOTAL_PAIR_TYPES)
  {
    return my_pairs[type - 1].count() / 2;       // Convert from indv. to pairs.
  }
  else
  {
    return (size_t)(-1);
  }
}

inline double                                     
Cont_trait_stats::pair_mean(pg::pair_type type) const
{
  if(type > 0 && type <= TOTAL_PAIR_TYPES)
  {
    return (pair_count(type) == 0 ? std::numeric_limits<double>::quiet_NaN() :
                                    my_pairs[type - 1].mean());
  }
  else
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
}

inline double                                     
Cont_trait_stats::pair_std_dev(pg::pair_type type) const
{
  if(type > 0 && type <= TOTAL_PAIR_TYPES)
  {
     return (pair_count(type) == 0 ? std::numeric_limits<double>::quiet_NaN() :
                                     my_pairs[type - 1].standard_deviation());
  }
  else
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
}

inline double                                     
Cont_trait_stats::pair_min(pg::pair_type type) const
{
  if(type > 0 && type <= TOTAL_PAIR_TYPES)
  {
    return (pair_count(type) == 0 ? std::numeric_limits<double>::quiet_NaN() :
                                    my_pairs[type - 1].min());
  }
  else
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
}

inline double                                     
Cont_trait_stats::pair_max(pg::pair_type type) const
{
  if(type > 0 && type <= TOTAL_PAIR_TYPES)
  {
    return (pair_count(type) == 0 ? std::numeric_limits<double>::quiet_NaN() :
                                    my_pairs[type - 1].max());
  }
  else
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
}

inline double                                     
Cont_trait_stats::correlation(pg::pair_type type) const
{
  if(type > 0 && type <= TOTAL_PAIR_TYPES)
  {
    return my_corr_infos[type - 1].correlation();
  }
  else
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
}

//============================================================================
// IMPLEMENTATION:  Cmpd_trait_stats
//============================================================================
//
inline
Cmpd_trait_stats::Cmpd_trait_stats(trait_list traits, cerrormultistream& e, std::string name,
                                   bool mped_member)
      :Base_trait_stats(0, e, mped_member), my_valid(false), my_traits(traits), my_traits_name(name)
{
  init();
}

inline    
Cmpd_trait_stats::Cmpd_trait_stats(const RefPedigree* p, trait_list traits, cerrormultistream& e,
                                   bool mped_member)
      :Base_trait_stats(p, e, mped_member), my_valid(false), my_traits(traits), my_traits_name("")
{
  init();
  compute();
}

inline RefMultiPedigree::pedigree_const_pointer   
Cmpd_trait_stats::pedigree() const
{
  return my_pedigree;
}

inline const Cmpd_trait_stats::trait_list&
Cmpd_trait_stats::traits() const
{
  return my_traits;
}

inline std::string
Cmpd_trait_stats::pedigree_name() const
{
  if(my_pedigree != 0)
  {
    return my_pedigree->name();
  }
  else
  {
    return "";
  }
}

inline const std::string& 
Cmpd_trait_stats::traits_name() const
{
  return my_traits_name;
}

inline bool                                       
Cmpd_trait_stats::valid() const
{
  return my_valid;
}

inline size_t                                     
Cmpd_trait_stats::ind_count_gender(gender g) const
{
  if(g < 0 || g > 2)
  {
    return (size_t)(-1);
  }
  else
  {
    return my_ind_counts_gender[g];
  }
}

inline size_t                                     
Cmpd_trait_stats::ind_count_founder(founder_status f) const
{
  if(f < 0 || f > 2)
  {
    return (size_t)(-1);
  }
  else
  {
    return my_ind_counts_founder[f];
  }
}

inline size_t                                     
Cmpd_trait_stats::pair_count(pg::pair_type type) const
{
  if(type < 1 || type > TOTAL_PAIR_TYPES)
  {
    return (size_t)(-1);
  }
  else
  {
    return my_pair_counts[type - 1];
  }
}


} // End namespace RPED
} // End namespace SAGE
