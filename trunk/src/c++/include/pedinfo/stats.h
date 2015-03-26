#ifndef  MP_STATS_H
#define  MP_STATS_H
//============================================================================
// File:     stats.h
//                                                                          
// Author:   
//                                                                          
// History:  8/00 - Modified and extended to utilize pair_generator class
//                  and to include trait information.  - Dan Baechle                                                    
//                                                                          
// Notes:    Declares the following classes -
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


#include <list>
#include <map>
#include <vector>
#include <set>
#include <algorithm>
#include <string>
#include <limits>
#include <ostream>
#include "error/errorstream.h"
#include "mped/sp.h"
#include "mped/mp.h"
#include "numerics/sinfo.h"
#include "numerics/histogram.h"
#include "numerics/corinfo.h"
#include "rped/rped.h"
#include "rped/loop.h"
#include "pairs/relpair.h"

namespace SAGE {
namespace RPED {

/// PEDIGREE_STATS define
#define   PEDIGREE_STATS          0

/// MULTIPEDIGREE_STATS define
#define   MULTIPEDIGREE_STATS     1

/// TOTAL_PAIR_TYPES define
#define   TOTAL_PAIR_TYPES        9

/// Typedef-ed pair_generator
typedef pair_generator             pg;

/// Typedef-ed pair_filter_trait
typedef pair_filter_trait          pft;

class Binary_trait_stats;
class Cont_trait_stats;
class Cmpd_trait_stats;

/** \internal
  * \brief For finding multiple mates
  *
  */
struct full_name
{
  full_name(const string& ped, const string& ind);
  
  string  name() const;

  string  ped_name;
  string  ind_name;
};

bool  operator<(const full_name& lhs, const full_name& rhs);

typedef std::map<full_name, std::vector<string> >  mate_map;


// - For finding consanguineous mating pairs.
//
struct parental_pair
{
  parental_pair(const string& p_one, const string& p_two);
  
  string  p1;
  string  p2;
};

typedef std::vector<parental_pair>  pp_vector;
typedef std::map<string, pp_vector>  pp_map;    // string is pedigree name.

struct no_pairs
{
  no_pairs();
  void  operator()(const pp_map::value_type& elem);

  bool  my_no_pairs;
};

/** \brief Declares a data structure which contains the statistics common to Ped_stats and MP_stats
  *
  * \par Introduction
  *
  * This data structure stores statistics common to Ped_stats and MP_stats. This includes various counts
  * (sibship size, family size, etc.), as well as SampleInfo's describing constituent traits.
  *
  * \par Trait information
  *
  * Although this class provides storage containers for trait statistics, it does not in fact provide
  * functionality for calculating those statistics. Such functionality is accomplished through the two
  * classes dervied from General_States: Ped_stats and MP_stats.
  */
class General_Stats
{
public:
  
  /// @name Operators / destructor
  //@{
  
    ///
    /// Operator+=.
    /// \param other Object to be added
    General_Stats & operator+=(const General_Stats& other);

    ///
    /// Destructor.
    ~General_Stats();

  //@}

  /// @name Statistics
  //@{
  
    ///
    /// Returns the mean size of sibships
    double mean_sibship_size() const;  
    
    ///
    /// Returns the variance of sibship sizes
    double var_sibship_size() const;  
    
    ///
    /// Returns the number of families.
    size_t family_count() const;  
    
    ///
    /// Returns the number of members.
    size_t member_count() const; 
    
    ///
    /// Returns the number of female members.
    size_t female_count() const;  
    
    ///
    /// Returns the number of male members (haha).
    size_t male_count() const;  
    
    ///
    /// Returns the number of individuals of unknown sex.
    size_t unknown_sex_count() const;  
    
    ///
    /// Returns the number of founders.
    size_t founder_count() const; 
    
    ///
    /// Returns the number of nonfounders.
    size_t nonfounder_count() const;  
    
    ///
    /// Returns the number of generations.
    size_t generations() const;  
    
    ///
    /// Returns the number of sib pairs.
    size_t sib_pairs() const;  
    
    ///
    /// Returns the number of half-sib pairs.
    size_t hsib_pairs() const;  
    
    ///
    /// Returns the number of avuncular pairs.
    size_t avuncular_pairs() const;  
    
    ///
    /// Returns the number of grandparental pairs.
    size_t grandp_pairs() const;  
    
    ///
    /// Returns the number of cousin pairs.
    size_t cousin_pairs() const; 
    
    ///
    /// Returns the number of brother-brother pairs.
    size_t brother_pairs() const;  
    
    ///
    /// Returns the number of sister-sister pairs.
    size_t sister_pairs() const;  
    
    ///
    /// Returns the number of brother-sister pairs.
    size_t brother_sister_pairs() const; 
    
    ///
    /// Returns the number of parental pairs.
    size_t parental_pairs() const;
    
    ///
    /// Returns the number of marriage loops.
    size_t marriage_loops() const;  
    
    ///
    /// Returns the number of non-marriage loops.
    size_t non_marriage_loops() const;  

    ///
    /// Returns the number of unconnecteds (singletons).
    size_t unconnecteds() const;  

    ///
    /// Returns the number of subpedigrees.
    size_t total_subpedigrees() const;  

    ///
    /// Returns the number of sibs.
    size_t sib_count() const;  

    ///
    /// Returns the number of ???
    size_t stats_type() const;  

    ///
    /// Returns the number of pairs of the indicated type.
    /// \param t The type of pair
    size_t pairs(pg::pair_type t) const;  

    ///
    /// Returns a boolean value indicating whether or not there are loops present.
    bool loops() const;
    
    ///
    /// Returns a boolean value indicating whether or not this object is valid. ???
    bool valid() const;  

    ///
    /// Returns the number of likelihood bits. ???
    size_t likelihood_bits() const;  
    
    ///
    /// Returns a histogram of sibship size frequencies.
    const Histogram & sib_size_freq() const;  
    
    ///
    /// Returns a mate_map of all the multiple mates.
    const mate_map & multiple_mates() const;

    ///
    /// Returns a pp_map of all consanguinous pairs.
    const pp_map & cons_pairs() const;

    ///
    /// Returns a SampleInfo with information about sibships. ???
    const SampleInfo & sibship() const;
    
    ///
    /// Returns a vector of Binary_trait_stats describing all the binary traits in the dataset.
    const std::vector<Binary_trait_stats> & Bt_stats() const;
    
    ///
    /// Returns a vector of Cont_trait_stats describing all the continuous traits in the dataset.
    const std::vector<Cont_trait_stats> & Ct_stats() const;
    
    ///
    /// Returns a vector of Cmpd_trait_stats describing all the compound traits in the dataset.
    const std::vector<Cmpd_trait_stats> & Cd_stats() const;

  //@}
      
    void invalidate();
    void init();

    void print_traits();                                       // For testing.

  protected: 
    General_Stats(cerrormultistream &e);

    // Data members.
    cerrormultistream errors;    
    
    SampleInfo      my_sibship;             // containing 4 moments of sibship
    size_t          my_num_gen;             // number of generations in the pedigree
    size_t          my_total_inds;          // total individuals in the pedigree
    size_t          my_total_founders;
    size_t          my_total_nonfounders;
    size_t          my_total_females;
    size_t          my_total_males;
    size_t          my_total_unknown_sex;
    size_t          my_max_bits;
    size_t          my_pairs[TOTAL_PAIR_TYPES];          
    size_t          my_marriage_loops;
    size_t          my_non_marriage_loops;
    size_t          my_unconnecteds;
    size_t          my_total_subpedigrees;
    Histogram       my_sibship_size_freq;    // sib_freq[5]= # of families with 6 children.
    mate_map        my_multiple_mates;       // individuals w. multiple mates.
    pp_map          my_cons_pairs;           // consanguineous mating pairs by pedigree.
    
    std::vector<Binary_trait_stats>    my_binary_trait_stats;
    std::vector<Cont_trait_stats>      my_cont_trait_stats;
    std::vector<Cmpd_trait_stats>      my_cmpd_trait_stats;
    
    bool              my_valid;
    size_t            my_stats_type;
};
 
/** \brief Calculate general statistics for an individual pedigree
  *
  * \par Introduction
  *
  * Given a single pedigree, this class calculates all sorts of interesting statistics.
  * These includes general count information (see General_Stats), as well as trait information.
  *
  * \par Trait information
  *
  * Using either the constructor or the various compute() functions, you can ask Ped_stats to calculate
  * statistics on a set of traits. Once calculated, those statistics are available through the parent class
  * accessors (see General_Stats::Bt_stats(), General_Stats::Ct_stats(), and General_Stats::Cd_stats() ).
  */
class Ped_stats : public General_Stats
{
public:

  /// @name Constructor / destructor / operators
  //@{

    ///
    /// Constructor
    /// \param p Pedigree whose stats this object will describe.
    /// \param e Errorstream to which error messages will be directed.
    Ped_stats(const RefPedigree *p, cerrormultistream &e);

    ///
    /// Constructor
    /// \param p Pedigree whose stats this object will describe.
    /// \param trait Id number of the trait you want analyzed
    /// \param e Errorstream to which error messages will be directed.
    Ped_stats(const RefPedigree *p, size_t trait, cerrormultistream &e);

    ///
    /// Constructor
    /// \param p Pedigree whose stats this object will describe.
    /// \param traits Vector of id numbers of the traits you want analyzed
    /// \param e Errorstream to which error messages will be directed.
    Ped_stats(const RefPedigree *p, std::vector<size_t> traits, cerrormultistream &e);

    ///
    /// Constructor
    /// \param e Errorstream to which error messages will be directed.
    Ped_stats(cerrormultistream &e);

    ///
    /// Destructor
    ~Ped_stats();

  //@}
  
  /// @name Basic info
  //@{
  
    ///
    /// Returns the name of the pedigree associated with this object.
    std::string pedigree_name() const;

    ///
    /// Returns a pointer to the pedigree associated with this object.
    RefMultiPedigree::pedigree_const_pointer ped() const;
    
  //@}
  
  /// @name Pedigree analysis
  //@{
  
    ///
    /// Computes basic statistics for the indicated pedigree.
    /// \param p The Pedigree in question
    /// \retval true Computation was successful
    /// \retval false Computation was \b not successful
    bool compute(const RefPedigree *p);

    ///
    /// Computes statistics for the indicated pedigree's trait.
    /// \param p The Pedigree in question
    /// \param trait The id number of the trait in question
    /// \retval true Computation was successful
    /// \retval false Computation was \b not successful
    bool compute(const RefPedigree *p, size_t trait);

    ///
    /// Computes statistics for the indicated pedigree's traits.
    /// \param p The Pedigree in question
    /// \param traits A vector of trait id numbers
    /// \retval true Computation was successful
    /// \retval false Computation was \b not successful
    bool compute(const RefPedigree *p, std::vector<size_t> traits);
    
  //@}
    
  private:

    void get_generation( std::list<member_const_pointer>& );
    void non_trait_computations(const RefPedigree* pedig);
    bool consanguineous(RefPedigree::family_const_iterator& family);
    static void fill_set(std::set<string>& s, RefPedigree::member_const_pointer ind);
    void compute_traits();
    void compute_traits(size_t trait);
    void compute_traits(std::vector<size_t> traits);
    void compute_pairs(const RefPedigree* pedig);
    void compute_loops(const RefPedigree* pedig);
    void compute_bits(const RefPedigree* pedig);
    
    // Data members.
    RefMultiPedigree::pedigree_const_pointer    my_pedigree;
};

/** \brief Calculate general statistics for an multipedigree
  *
  */
class MP_stats : public General_Stats
{
public:

  /// @name Constructor / destructor
  //@{
  
    ///
    /// Constructor.
    /// \param e Errorstream to which error messages will be directed
    MP_stats(cerrormultistream &e);

    ///
    /// Constructor.
    /// \param mp The RefMultiPedigree whose stats this object will calculate
    /// \param e Errorstream to which error messages will be directed
    MP_stats(const RefMultiPedigree *mp, cerrormultistream &e);

    ///
    /// Constructor.
    /// \param mp The RefMultiPedigree whose stats this object will calculate
    /// \param trait The id number of the trait you want analyzed
    /// \param e Errorstream to which error messages will be directed
    MP_stats(const RefMultiPedigree *mp, size_t trait, cerrormultistream &e);

    ///
    /// Constructor.
    /// \param mp The RefMultiPedigree whose stats this object will calculate
    /// \param traits A vector of trait id numbers that you want analyzed
    /// \param e Errorstream to which error messages will be directed
    MP_stats(const RefMultiPedigree *mp, std::vector<size_t> traits, cerrormultistream &e);

    ///
    /// Destructor.
    ~MP_stats();

  //@}
  
  /// @name Basic info
  //@{
  
    ///
    /// Returns a pointer to the RefMultiPedigree used at construction.
    const RefMultiPedigree * mped() const;

    ///
    /// Returns the number of constituent pedigrees in the RefMultiPedigree used at construction.
    size_t pedigree_count() const;

    ///
    /// Returns the number of Ped_stats calculated.
    size_t size() const;
    
    ///
    /// Returns a pointer to the Ped_stats instanced associated with the indiciated pedigree
    /// \param i The id number of the pedigree in question
    const Ped_stats * operator[](size_t i) const;
  
  //@}

  /// @name Multipedigree analysis
  //@{

    ///
    /// Computes non-trait information for the given RefMultiPedigree
    /// \param mp The RefMultiPedigee you want analyzed
    /// \retval true Computation was successful
    /// \retval false Computation was \b not successful
    bool compute(const RefMultiPedigree *mp);

    ///
    /// Computes statistics on a single trait
    /// \param mp The RefMultiPedigee you want analyzed
    /// \param trait The id number of the trait you want analyzed
    /// \retval true Computation was successful
    /// \retval false Computation was \b not successful
    bool compute(const RefMultiPedigree *mp, size_t trait);

    ///
    /// Computes statistics on a number of traits
    /// \param mp The RefMultiPedigee you want analyzed
    /// \param traits A vector of trait id numbers you want analyzed
    /// \retval true Computation was successful
    /// \retval false Computation was \b not successful
    bool compute(const RefMultiPedigree *mp, std::vector<size_t> traits);

  //@}
  
  /// @name Descriptive statistics
  //@{

    ///
    /// Returns a SampleInfo describing pedigree size.
    const SampleInfo & pedigree_size_info() const;

    ///
    /// Returns a SampleInfo describing family counts.
    const SampleInfo & family_count_info() const;

    ///
    /// Returns a SampleInfo describing likelihood bits frequencies.
    const Histogram & likelihood_bits_freq() const;

    ///
    /// Returns a SampleInfo describing family count frequencies.
    const Histogram & family_count_freq() const;

    ///
    /// Returns a SampleInfo describing generation frequencies.
    const Histogram & generation_freq() const;

  //@}    
    
  private:
    void          init();
    bool          trait_valid(size_t trait);
    bool          cmpd_trait_valid(std::vector<size_t> traits);
    void          compute_pedigrees();
    void          compute_pedigrees(size_t trait);
    void          compute_pedigrees(std::vector<size_t> traits);
    std::string   compute_traits_name(std::vector<size_t> traits)   const;
    
    // Data members.
    const RefMultiPedigree*   my_multipedigree;
    SampleInfo                my_pedigree_size_info;    // 4 moments of pedigree size
    SampleInfo                my_family_count_info;     // 4 moments of # of nuclear families 
                                                        // in pedigrees
    Histogram                 my_likelihood_bits_freq;
    Histogram                 my_family_count_freq;
    Histogram                 my_generations_freq;
    vector<Ped_stats>         my_ped_stats;
};

/** \brief Tabulates sibship parental informativity data for all trait classes
  *
  * Contains enumerations and utility functions common to all trait_stats classes.   
  */
class Base_trait_stats
{
  public:
    
    /// Gender enum
    enum  gender { MALE = 0, /*!< \hideinitializer Code for male */ 
                   FEMALE,   /*!< Code for female                */
                   UNKNOWN   /*!< Code for unknown sex           */ };
    

    /// Founder status enum    
    enum  founder_status   { FOUNDER = 0, /*!< \hideinitializer Code for founder */
                             NON_FOUNDER, /*!< Code for nonfounder               */
                             UNCONNECTED  /*!< Code for singleton                */ };
  
    ///
    /// ???
    bool  mped_member() const;
  
    ///
    /// ???
    size_t sibship(int parents) const;

    ///
    /// Returns a SampleInfo object describing pedigree size ???
    const SampleInfo&  pedigree_size_info() const;
    const SampleInfo&  sibship_size_info() const;
    const Histogram&   family_count_freq() const;
    
  private:
    void   init();

  protected:
    // Constructor/destructor.
    Base_trait_stats(const RefPedigree* p, cerrormultistream& e, bool mped_member);
    
    void                 compute_sibships(const ind_filter& i_filter);
    Base_trait_stats&    operator+=(const Base_trait_stats& other); 
    void                 print()   const;                      // For testing only.
    
    gender               determine_gender(MPED::SexCode ind_sex)                       const;
    founder_status       determine_founder_status(member_const_pointer P1,
                                                  member_const_pointer P2,
                                                  size_t offspring_count)    const;
    void                 compute_pedigree_size(const ind_filter& i_filter);
    void                 compute_sibship_size(const ind_filter& i_filter);
    void                 compute_family_count(const ind_filter& i_filter);
    
    // Data members.
    RefMultiPedigree::pedigree_const_pointer  my_pedigree;
    cerrormultistream   my_errors;
    bool  mped_member_flag;            // Is this object a member of an MP_stats object?
    
    size_t         my_sibships[3];          // Counts of sibships w. 0, 1, and 2 informative
                                            // parents.
    SampleInfo     my_pedigree_size_info;   // containing 4 moments of pedigree size
    SampleInfo     my_sibship_size_info;    // containing 4 moments of sibship size  
    Histogram      my_family_count_freq;
};

/** \brief Generates statistics over a pedigree or multipedigree for a binary trait
  *
  */
class Binary_trait_stats : public Base_trait_stats
{
  public:
    /// Individual affection status enum
    enum ind_affection_status { AFF = 0, /*!< \hideinitializer Code for affected */
                                UNAFF,   /*!< Code for unaffected                */
                                UNINF    /*!< Code for uninformative             */ };
    
    // Constructor/destructor.
    Binary_trait_stats(size_t trait, cerrormultistream& e, std::string name = "", bool mped_member = false);
    Binary_trait_stats(const RefPedigree* p, size_t trait, cerrormultistream& e, bool mped_member = false);

    // Gets.
    RefMultiPedigree::pedigree_const_pointer   pedigree()                                           const;
    size_t                                     trait()                                              const;
    std::string                                pedigree_name()                                      const;
    const std::string&                         trait_name()                                         const;
    bool                                       valid()                                              const;
    size_t                                     ind_count_gender(gender g, ind_affection_status st)  const;
    size_t                                     ind_count_founder(founder_status f, 
                                                                  ind_affection_status st)          const;
    size_t                                     pair_count(pg::pair_type type, 
                                                          pft::affection_status st)                 const;
    double                                     correlation(pg::pair_type type)                      const;
    
    Binary_trait_stats&                        operator+=(const Binary_trait_stats& other); 
    void                                       print()     const;           // For test purposes only.   
  
  private:
    void                                       compute();
    void                                       init();         
    void                                       compute_sibships();
    void                                       compute_pedigree_size();
    void                                       compute_sibship_size();
    void                                       compute_family_count();
    void                                       compute_inds_gender();
    void                                       compute_inds_founder();
    void                                       compute_pairs();
    void compute_correlations(RefPedigree::member_pointer member1,
                              RefPedigree::member_pointer member2, pg::pair_type p_type);
    bool member1_is_uncle(size_t index1, RefPedigree::member_pointer member2);
  
    // Data members.
    bool              my_valid;                               // Stats have been computed.    
    size_t            my_trait;
    std::string       my_trait_name;
    size_t            my_ind_counts_gender[3][3];             // Rows:    males, females, unknown
                                                              // Columns: aff., unaff., missing
    size_t            my_ind_counts_founder[3][3];            // Rows:    founders, non-founders, unconnected
                                                              // Columns: aff., unaff., missing
    size_t            my_pair_counts[TOTAL_PAIR_TYPES][4];    // Columns: concord. aff., discord.,
                                                              //          concord. unaff., uninform
    SimpleCorrelationInfo   my_corr_infos[TOTAL_PAIR_TYPES];
};


/** \brief Generates statistics over a pedigree or multipedigree for a continuous trait
  *
  */
class Cont_trait_stats : public Base_trait_stats
{
  public:
    enum {ALL = 3};
    
    // Constructor/destructor.
    Cont_trait_stats(size_t trait, cerrormultistream& e, std::string name = "", bool mped_member = false);
    Cont_trait_stats(const RefPedigree* p, size_t trait, cerrormultistream& e, bool mped_member = false);

    // Gets.
    RefMultiPedigree::pedigree_const_pointer   pedigree()                               const;
    size_t                                     trait()                                  const;
    std::string                                pedigree_name()                          const;
    const std::string&                         trait_name()                             const;
    bool                                       valid()                                  const;
    size_t                                     ind_gender_count(int g)                  const;
    double                                     ind_gender_mean(int g)                   const;
    double                                     ind_gender_std_dev(int g)                const;
    double                                     ind_gender_min(int g)                    const;
    double                                     ind_gender_max(int g)                    const;
    size_t                                     ind_founder_count(int f)                 const;
    double                                     ind_founder_mean(int f)                  const;
    double                                     ind_founder_std_dev(int f)               const;
    double                                     ind_founder_min(int f)                   const;
    double                                     ind_founder_max(int f)                   const;
    size_t                                     pair_count(pg::pair_type type)           const;
    double                                     pair_mean(pg::pair_type type)            const;
    double                                     pair_std_dev(pg::pair_type type)         const;
    double                                     pair_min(pg::pair_type type)             const;
    double                                     pair_max(pg::pair_type type)             const;
    double                                     correlation(pg::pair_type type)          const;
    
    Cont_trait_stats&                          operator+=(const Cont_trait_stats& other); 
    void                                       print()    const;           // For test purposes only.   
  
  private:
    void                                       compute();
    void                                       init();         
    void                                       compute_sibships();
    void                                       compute_pedigree_size();
    void                                       compute_sibship_size();
    void                                       compute_family_count();
    void                                       compute_inds_gender();
    void                                       compute_inds_founder();
    void                                       compute_pairs();
    void compute_correlations(RefPedigree::member_pointer member1,
                              RefPedigree::member_pointer member2, 
                              size_t index1, size_t index2,
                              double value1, double value2, pg::pair_type p_type);
    bool member1_is_uncle(size_t index1, RefPedigree::member_pointer member2);
    
    // Data members.
    bool              my_valid;                               // Stats have been computed.    
    size_t            my_trait;
    std::string       my_trait_name;
    SampleInfo        my_inds_gender[4];                      // Rows:  males, females, unknown, all.
    SampleInfo        my_inds_founder[4];                     // Rows:  founders, non-founders
                                                              //        unconnected.
    SampleInfo        my_pairs[TOTAL_PAIR_TYPES];            
    SimpleCorrelationInfo   my_corr_infos[TOTAL_PAIR_TYPES];        
};

//----------------------------------------------------------------------------
//  Class:    Cmpd_trait_stats
//                                                                          
//  Purpose:  Generates statistics over a pedigree or multipedigree
//            for a compound trait.  To be included an individual must have
//            a value for each trait specified.
//                                                                          
//----------------------------------------------------------------------------
//
class Cmpd_trait_stats : public Base_trait_stats
{
  public:
    typedef std::vector<size_t>  trait_list;
    
    // Constructor/destructor.
    Cmpd_trait_stats(trait_list traits, cerrormultistream& e, std::string name = "", bool mped_member = false);
    Cmpd_trait_stats(const RefPedigree* p, trait_list traits, cerrormultistream& e, bool mped_member = false);

    // Gets.
    RefMultiPedigree::pedigree_const_pointer   pedigree()                                   const;
    const trait_list&                          traits()                                     const;
    std::string                                pedigree_name()                              const;
    const std::string&                         traits_name()                                const;
    bool                                       valid()                                      const;
    size_t                                     ind_count_gender(gender g)                   const;
    size_t                                     ind_count_founder(founder_status f)          const;
    size_t                                     pair_count(pg::pair_type type)               const;
    
    Cmpd_trait_stats&                          operator+=(const Cmpd_trait_stats& other); 
    
    void                                       print()    const;           // For test purposes only.   
  
  private:
    void                                       compute();
    void                                       compute_traits_name();
    void                                       init();         
    void                                       make_ind_filter();
    void                                       compute_sibships();
    void                                       compute_pedigree_size();
    void                                       compute_sibship_size();
    void                                       compute_family_count();
    void                                       compute_inds_gender();
    void                                       compute_inds_founder();
    void                                       compute_pairs();
  
    // Data members.
    bool              my_valid;                               // Stats have been computed.    
    trait_list        my_traits;
    std::string       my_traits_name;
    ind_filter        my_ind_filter;                          // filter for indivs. informative for all
                                                              // traits in the trait list
                                                              
    size_t            my_ind_counts_gender[3];                // males, females, unknown
    size_t            my_ind_counts_founder[3];               // founder, non-founder, unconnected
    size_t            my_pair_counts[TOTAL_PAIR_TYPES];       
};

} // End namespace RPED
} // End namespace SAGE

#include "pedinfo/stats.ipp"

#endif

