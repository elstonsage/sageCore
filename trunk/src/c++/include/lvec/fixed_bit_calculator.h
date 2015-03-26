#ifndef FIXED_BIT_CALCULATOR_H
#define FIXED_BIT_CALCULATOR_H

//==========================================================================
//  File:    fixed_bit_calculator.h
//
//  Author:  Geoff Wedig
//
//  History: 0.1 Initial Implementation
//           1.0 Updated to new libraries                        yjs Sep. 04
//
//  Notes:
//
//  Copyright (c) 1998 R. C. Elston
//  All Rights Reserved
//==========================================================================

#include "output/Output.h"
#include "fped/fped.h"

namespace SAGE
{

class fixed_bit_calculator;

class fixed_bit_container
{
  public:

    typedef fixed_bit_container                self;

    enum bit_type { unknown, fixed, no_rel, dont_care };

    fixed_bit_container(const FPED::Subpedigree& subped);
    fixed_bit_container(const self& fbc);

    fixed_bit_container& operator=(const self& fbc);
    
    ~fixed_bit_container();

    // Accessors
    //

    const FPED::Subpedigree& get_pedigree() const;

    bool mother_fixed(const FPED::Member& i) const;
    bool father_fixed(const FPED::Member& i) const;

    bool mother_dont_care(const FPED::Member& i) const;
    bool father_dont_care(const FPED::Member& i) const;

    // Synchronization of bits for mother/father.  Note that it does not say
    // whether they are the same (0,0 <-> 1,1) or different (0,1 <-> 1,0) as
    // that is irrelevant.
    
    bool mother_synchronized(const FPED::Member& i, const FPED::Member& j) const;
    bool father_synchronized(const FPED::Member& i, const FPED::Member& j) const;

    void dump(ostream& o) const;

  protected:

    friend class fixed_bit_calculator;

    // Building functions - for fixed_bit_calculator
    //
    // ptype - 0 mother, 1 father
    void set_parent_bit_type(bool ptype, const FPED::Member&, bit_type t);

    void synchronize_parent(bool ptype, const FPED::Member&, const FPED::Member& j);

    struct individual_data
    {
      individual_data();
      
      bit_type mother_bit_type;
      bit_type father_bit_type;
      
    };
    struct family_data
    {
      family_data();
      
      std::vector<bool> my_mother_sync_bits;
      std::vector<bool> my_father_sync_bits;
    };

    std::vector<individual_data> my_individual_data;
    std::vector<family_data>     my_family_data;

    std::vector<size_t>          my_ind_sib_indices;
    
    const FPED::Subpedigree* my_ped;
};

/// The fixed bit calculator, given a pedigree and an associated marker,
/// creates a data structure (the fixed bit container) which can be queried
/// as to the status of various bits.  Bits can be in any of four states:
/// fixed, no relationship, or don't care, as well as an
/// unknown status (assumed the same as no relationship in algorithms)
///
/// In addition to the individual bits, sibling bits are classified as
/// synchronous or asynchronous.  Synchronous bits may be either the same,
/// indicating that the sibs always share an allele from the given parent,
/// or different, indicating that they cannot share an allele from that
/// parent.  Only the synchronicity status is stored, not the ibd status of
/// synchronous bits.
///
/// The algorithms used to calculate these statistics assume the pedigree is
/// without loops.  There are extremely unlikely fixed bits that can occur
/// due to inbreeding which these algorithms will miss, but as they are
/// extremely unlikely to occur, they can be safely ignored at a small cost
/// to performance.
///
/// NOTE:  This system *relies* on all parents being sexed.
class fixed_bit_calculator
{
  public:

    fixed_bit_calculator(const FPED::Subpedigree&         mmap, 
                         const MLOCUS::inheritance_model& imodel);

    ~fixed_bit_calculator();
    
    const fixed_bit_container& get_fixed_bit_container() const;

    void dump(ostream& out) const;

  protected:

    /// \brief Stores information about potential fixed and synchronous bits
    ///
    struct potential_data
    {
      potential_data() 
        : is_father(false),
          id1(0),
          id2(0),
          bit_pattern(-1)
      { }
    
      potential_data(bool is_f, FPED::MemberConstPointer i, FPED::MemberConstPointer j = 0)
        : is_father(is_f),
          id1(i),
          id2(j),
          bit_pattern(-1)
      { }

      bool is_father;  ///< Stores whether the bit relates to the father or mother
   
      FPED::MemberConstPointer id1;
      FPED::MemberConstPointer id2; // used only for synchronizations
      
      int bit_pattern; ///< Storage of the patterns we've seen already
    };

    typedef list<potential_data>    potential_list;
    typedef potential_list::iterator potential_iterator;

    // Building functions
    //

    void run_nuclear_family(const FPED::Family&              fam,
                            const MLOCUS::inheritance_model& imodel);

    void test_parent(const FPED::Member&              par,
                     const MLOCUS::inheritance_model& imodel);

    void set_all_no_rel            (bool                             is_father);
    void do_basic_fixed_check      (bool                             is_father,
                                    const MLOCUS::inheritance_model& imodel,
                                    const vector<int>&               allele_sources);
    void set_all_homozygous        (bool                             is_father);
    void create_all_potential_sync (bool                             is_father);

    void test_family  (const MLOCUS::inheritance_model& imodel);
    void test_children(const MLOCUS::inheritance_model& imodel,
                       const MLOCUS::phased_genotype&   mg,
                       const MLOCUS::phased_genotype&   fg);
    
    // Members
    //

    fixed_bit_container my_data;
    
    FPED::FamilyConstPointer my_current_family;

    potential_list      my_potential_fixed_bit;
    potential_list      my_potential_synchronized_bits;
};

}

#include "lvec/fixed_bit_calculator.ipp"

#endif
