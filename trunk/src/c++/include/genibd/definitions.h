#ifndef GENIBD_DEFS_H
#define GENIBD_DEFS_H

//==========================================================================
//  File:      definitions.h
//
//  Author:    Yeunjoo Song
//
//  History:   Initial implementation.                              Nov. 03
//
//  Notes:     This header file contains various type definition statements
//             used through out the genibd files.
//
//  Copyright (c) 2003 R.C. Elston
//    All Rights Reserved
//==========================================================================

#include "numerics/print_util.h"
#include "util/dots.h"
#include "rped/genome_description.h"
#include "fped/fped.h"
#include "pairs/relpair.h"
#include "ibd/basic_storage_ibd.h"
#include "ibd/exact_ibd_analysis.h"
#include "ibd/ibdfile.h"
#include "mcmc/mcmc_parser.h"
#include "mcmc/mcmc_simulator.h"
#include "globals/SAGEConstants.h"

namespace SAGE   {
namespace GENIBD {

typedef RPED::RefMultiPedigree                               RefMultiPedigree;
typedef RPED::RefPedigree                                    RefPedigree;
typedef RPED::RefSubpedigree                                 RefSubpedigree;
typedef RPED::RefFamily                                      RefFamily;
typedef RPED::RefMPedInfo                                    RefMPedInfo;

typedef FPED::FilteredMultipedigree                          filtered_multipedigree;

typedef FPED::Member                                         member_type;             
typedef FPED::Family                                         family_type;
typedef FPED::Subpedigree                                    subped_type;             

typedef FPED::FilteredPedigreeInfo                           filtered_pedigree_info;
typedef FPED::FilteredMemberInfo                             filtered_member_info;

typedef FPED::PedigreeIterator                               fpedigree_iterator;
typedef FPED::PedigreeConstIterator                          fpedigree_const_iterator;
typedef FPED::FamilyConstIterator                            ffamily_const_iterator;
typedef FPED::MemberConstIterator                            fmember_const_iterator;
typedef FPED::SiblingConstIterator                           fsibling_const_iterator;
typedef FPED::OffspringConstIterator                         foffspring_const_iterator;
typedef FPED::MateConstIterator                              fmate_const_iterator;

typedef FPED::MemberConstPointer                             fmember_const_pointer;    

typedef RPED::genome_description                             genome_description;
typedef RPED::genome_description::region_type                region_type;

typedef Likelihood_Vector                                    lvector;
typedef lvector::size_type                                   size_type;

typedef MLOCUS::allele                                       allele;
typedef MLOCUS::child_genotype_set                           child_genotype_set;
typedef MLOCUS::phased_genotype                              phased_genotype;
typedef MLOCUS::inheritance_model                            inheritance_model;
typedef inheritance_model::unphased_penetrance_iterator      unphased_pen_iter;
typedef inheritance_model::phased_penetrance_iterator        phased_pen_iter;

typedef MCMC::mcmc_parameters                                mcmc_parameters;
typedef MCMC::mcmc_parser                                    mcmc_parser;
typedef MCMC::McmcMeiosisMap                                 mcmc_meiosis_map;
typedef MCMC::mcmc_data_accessor                             mcmc_data_accessor;
typedef MCMC::mcmc_simulator                                 mcmc_simulator;

enum mode_type          { SINGLEPOINT, MULTIPOINT };
enum choice_type        { NO = 0, YES = 1, ALWAYS = 2 };
enum pair_category_type { SIB = 1, ALL_SIB = 2, RELATIVE = 3, ALL = 4 };

struct genibd_region_type
{
  genibd_region_type() : name(), output() { }
  genibd_region_type(const string& n, const string& o) : name(n), output(o) { }

  string name;
  string output;
};

typedef std::list<genibd_region_type>             genibd_region_list;
typedef genibd_region_list::const_iterator        genibd_region_iterator;

struct analysis_data
{
  analysis_data()
    : allow_single(false), allow_exact_single(false), allow_exact_multi(false),
      allow_sim(false), intervals(false), max_ped_size(0), max_loci(0)
  {}
  
  bool allow_single;
  bool allow_exact_single;
  bool allow_exact_multi;
  bool allow_sim;
  
  bool intervals;
  
  size_t max_ped_size;
  size_t max_loci;
};

struct relative_pair_less : std::binary_function<const pair_generator::relative_pair&,
                                                 const pair_generator::relative_pair&,
                                                 bool>
{
  relative_pair_less() {}
  
  bool operator()(const pair_generator::relative_pair& r1,
                  const pair_generator::relative_pair& r2) const
  {
    if( r1.type() != r2.type() )
      return r1.type() < r2.type();
    else
    {
      if( r1.member_one()->pedigree()->name() != r2.member_one()->pedigree()->name() )
        return r1.member_one()->pedigree()->name() < r2.member_one()->pedigree()->name();
      else
      {
        string r1_member_one = r1.member_one()->name();
        string r1_member_two = r1.member_two()->name();

        if( r1_member_one > r1_member_two )
        {
          r1_member_one = r1.member_two()->name();
          r1_member_two = r1.member_one()->name();
        }
        
        string r2_member_one = r2.member_one()->name();
        string r2_member_two = r2.member_two()->name();

        if( r2_member_one > r2_member_two )
        {
          r2_member_one = r2.member_two()->name();
          r2_member_two = r2.member_one()->name();
        }
        
        if( r1_member_one != r2_member_one )
          return r1_member_one < r2_member_one;
        else
          if( r1_member_two != r2_member_two )
            return r1_member_two < r2_member_two;
      }
    }

    return false;
  }
};

typedef set<pair_generator::relative_pair, relative_pair_less> relpair_set_type;

struct filtered_relative_pair
{
  filtered_relative_pair()
    : member_one(NULL), member_two(NULL),
      connector_one(NULL), connector_two(NULL),
      type(pair_generator::NULL_TYPE)
  {}

  filtered_relative_pair(fmember_const_pointer m1, fmember_const_pointer m2,
                         fmember_const_pointer c1, fmember_const_pointer c2,
                         pair_generator::pair_type t)
    : member_one(m1), member_two(m2), connector_one(c1), connector_two(c2), type(t)
  {}

  fmember_const_pointer      member_one;
  fmember_const_pointer      member_two;
  fmember_const_pointer      connector_one;
  fmember_const_pointer      connector_two;
  pair_generator::pair_type  type;
};

typedef vector<filtered_relative_pair>   relative_pairs;

struct ind_genotype
{
  ind_genotype(phased_genotype m)   : mg(m) {}

  ind_genotype(const MLOCUS::penetrance_model::phased_penetrance_iterator& upi)
    : mg(upi.phased_geno()) {}

  void    print()     const { cout << mg.name(); } 
  double  frequency() const
  {
    if( mg.allele1().name() == "~Y" )
      return mg.allele2().frequency();
    else if( mg.allele2().name() == "~Y" )
      return mg.allele1().frequency();

    return mg.frequency();
  }
  
  MLOCUS::phased_genotype    mg;  // Phased marker genotype.
};

struct conditional_genotype
{
  child_genotype_set child_geno;
  bool               b[4];
};

inline bool both_members_exist(const pair_generator::relative_pair& rel_pair,
                               const subped_type&                   subped)
{
  bool mem1_exist = false;
  
  for( fmember_const_iterator mem1 = subped.member_begin(); mem1 != subped.member_end(); ++mem1 )
    if( mem1->name() == rel_pair.member_one()->name() )
    {
      mem1_exist = true;
      break;
    }

  bool mem2_exist = false;
  
  for( fmember_const_iterator mem2 = subped.member_begin(); mem2 != subped.member_end(); ++mem2 )
    if( mem2->name() == rel_pair.member_two()->name() )
    {
      mem2_exist = true;
      break;
    }

  return mem1_exist && mem2_exist;
}

inline bool is_Y_genotype(const phased_genotype& pg)
{
  string al1 = pg.allele1().name();
  string al2 = pg.allele2().name();

  if( al1 != "~Y" && al2 != "~Y" )
    return false;

  return true;
}

} // end of namespace GENIBD
} // end of namespace SAGE

#endif
