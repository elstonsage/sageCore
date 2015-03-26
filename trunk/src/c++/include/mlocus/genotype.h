#ifndef GENOTYPE_H
#define GENOTYPE_H

//============================================================================
//  File:       genotype.h
//
//  Author:     Bob Steagall
//
//  History:    Version 1.00.00
//              X, Y-linkage added    - yjs  Mar. 2002
//              Significantly revised - gcw  Nov. 2002
//
//  Notes:    
//
//  Copyright (c) 2002 R.C. Elston
//  All Rights Reserved
//============================================================================
//

/// @file 
/// @brief genotype_model related objects
///
/// genotype.h contains all the genotype related classes, including alleles.
///
/// A genotype is a pair of alleles at a specific marker.  The alleles may
/// be either phased (parent of origin known) or unphased (parent of origin
/// unknown).  Determination of phase information is done by the use of
/// 'separators', characters which define the relationship between the
/// alleles.
///
/// The primary class is the PRIVATE::genotype_model_info.  All the other classes
/// reference it, either directly or indirectly.  Alleles, phased and
/// unphased genotypes, and their iterators are all only valid when speaking
/// about a particular model.  However, though all the classes refer to a
/// particular info class, the info class itself is a mostly invisible
/// letter class, and not to be used by applications.  Instead, the
/// genotype_model provides the primary interface to be used.  From the
/// genotype_model, alleles, genotypes, and iterators can be gotten.

#include <algorithm>
#include <iomanip>
#include <list>
#include <vector>
#include <map>

#include <boost/smart_ptr.hpp>
#include <boost/iterator_adaptors.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/sequenced_index.hpp>
#include <boost/multi_index/key_extractors.hpp>
#include "globals/data_types.h"
#include "mped/mpfwd.h"
#include "util/disambiguator.h"

using namespace std;

namespace SAGE   {
namespace MLOCUS {

#ifndef SAGE_NPOS_DEFINED
    #
    const uint  NPOS = ((uint) -1) / 2;
    #define SAGE_NPOS_DEFINED
    #
#endif

#define SEX_ALLELE_ID ((uint) -1)

enum Ordering {Unphased = 0, PhasedForward, PhasedBackward};

enum GenotypeModelType
{
  AUTOSOMAL,
  X_LINKED,
  Y_LINKED
};

class allele;
class phased_genotype;
class unphased_genotype;

class allele_iterator;
class phased_genotype_iterator;
class unphased_genotype_iterator;

namespace PRIVATE
{
    class genotype_model_info;

    struct allele_info
    {
        enum SexType
        {
          st_NULL_X = 1, ///< Represents a no allele case on the X chromosome (ie, Y-Linked marker)
          st_NULL_Y = 2, ///< Represents a no allele case on the Y chromosome (ie, X-Linked marker)
          st_NORMAL = 3, ///< Represents a non-null allele
        };
      
        allele_info();
        allele_info(const string& name, double freq, uint id);
        allele_info(SexType sex);
        allele_info(const allele_info&);
        allele_info& operator= (const allele_info&);

        bool    operator==(const allele_info& rhs) const;
        bool    operator!=(const allele_info& rhs) const;

        string      name;
        double      frequency;
        uint        id;
        
        SexType     my_sex_type;
    };

    typedef std::pair<const allele_info*, const allele_info*> AlleleInfoPair;
    
    inline bool operator<(const AlleleInfoPair& a1, const AlleleInfoPair& a2)
    {
      return a1.first->id < a2.first->id ||
             a1.first->id == a2.first->id && a1.second->id < a2.second->id;
    }
    
    struct genotype_info
    {
      public:
      
        genotype_info();
        genotype_info(const allele_info*         a1,
                      const allele_info*         a2,
                      uint                       phased_id,
                      uint                       flipped_id,
                      uint                       unphased_id,
                      const genotype_model_info* info);
                             

        genotype_info(const genotype_info&);
        genotype_info& operator=(const genotype_info&);
        
        AlleleInfoPair        my_alleles;
        uint                  my_phased_id;
        uint                  my_flipped_id;
        uint                  my_unphased_id;
        
        const genotype_model_info* my_ginfo;
    };
    
    extern allele_info   invalid_allele_info;
    extern genotype_info invalid_phased_genotype_info;
    extern genotype_info invalid_unphased_genotype_info;

    // Definitions of our multi_index container, for readability.  This isn't
    // bad, since it's in the PRIVATE namespace, so doesn't clutter up MLOCUS
    using namespace boost::multi_index;

    // Structs used for tagging of types
    
    struct sequence{};
    struct id{};
    struct alleles{};
    
    // Orderings
    
    typedef sequenced<tag<sequence> >                        SequencedOrdering;
    typedef ordered_unique<tag<id>,
              member<genotype_info,uint,
                     &genotype_info::my_phased_id> >         PhasedIdOrdering;
    typedef ordered_unique<tag<id>,
              member<genotype_info,uint,
                     &genotype_info::my_unphased_id> >       UnphasedIdOrdering;
    typedef ordered_unique<tag<alleles>,
              member<genotype_info,AlleleInfoPair,
                     &genotype_info::my_alleles> >           AlleleOrdering;
    
    // Phased storage Container
                     
    typedef multi_index_container<
              genotype_info,
              indexed_by<
                SequencedOrdering,
                PhasedIdOrdering,
                AlleleOrdering> >                            PhasedGenotypeStorage;
                     
    typedef PhasedGenotypeStorage::index<sequence>::type     PhasedGenotypeBySequence;
    typedef PhasedGenotypeStorage::index<id>::type           PhasedGenotypeById;
    typedef PhasedGenotypeStorage::index<alleles>::type      PhasedGenotypeByAlleles;

    // Unphased storage Container

    typedef multi_index_container<
              genotype_info,
              indexed_by<
                SequencedOrdering,
                UnphasedIdOrdering,
                AlleleOrdering> >                            UnphasedGenotypeStorage;
                
    typedef UnphasedGenotypeStorage::index<sequence>::type   UnphasedGenotypeBySequence;
    typedef UnphasedGenotypeStorage::index<id>::type         UnphasedGenotypeById;
    typedef UnphasedGenotypeStorage::index<alleles>::type    UnphasedGenotypeByAlleles;

    //----------------------------------------------------------------------------
    //  Class:      genotype_model_info
    //----------------------------------------------------------------------------
    /// The letter class for storing a set of genotypes.  In actual fact, it
    /// does not store genotypes, but creates them at need from the underlying
    /// alleles, which are stored as allele_info.
    class genotype_model_info
    {
      public:
      
        genotype_model_info(GenotypeModelType type = AUTOSOMAL);
        genotype_model_info(const genotype_model_info&);
        genotype_model_info(const string& name, GenotypeModelType type = AUTOSOMAL);
        
        genotype_model_info& operator=(const genotype_model_info&);
    
        ~genotype_model_info();
        
        void configure_sex_type(GenotypeModelType);
    
        void rebuild_genotypes();
        
        void add_autosomal_genotypes();
        void add_x_linked_genotypes();
        void add_y_linked_genotypes();
        
        void add_autosomal_genotypes_for_allele(allele_iterator a);

        void add_x_linked_male_genotype_for_allele(allele_iterator a);

        void add_y_linked_female_genotype();
        void add_y_linked_male_genotypes();
        
        void add_phased_genotype         (const genotype_info& info);

        void add_unsexed_phased_genotype (const genotype_info& info);
        void add_male_phased_genotype    (const genotype_info& info);
        void add_female_phased_genotype  (const genotype_info& info);

        void add_unphased_genotype         (const genotype_info& info);

        void add_unsexed_unphased_genotype (const genotype_info& info);
        void add_male_unphased_genotype    (const genotype_info& info);
        void add_female_unphased_genotype  (const genotype_info& info);

        void initialize_sex_iterators();
        
        phased_genotype      get_phased_genotype(allele a1, allele a2) const;
    
        phased_genotype      get_phased_genotype   (uint id) const;
        unphased_genotype    get_unphased_genotype (uint id) const;

        typedef std::vector<PRIVATE::allele_info> allele_vector;
        typedef std::map<string, uint>            allele_map;
        typedef std::list<string>                 allele_name_buffer;

        string                      name;
        char                        separators[4];
        string                      missing_allele_name;
        allele_vector               alleles;
        allele_map                  allele_names;
        allele_name_buffer          remap_buffer;
    
        bool                        dynamic_alleles;
    
        GenotypeModelType           my_type;
        
        /// Stores information relating to the sex-specific allele
        /// in cases where the marker is not autosomal
        allele_info                 my_sex_allele_info;
     
        PRIVATE::PhasedGenotypeStorage   my_phased_genotypes;
        PRIVATE::UnphasedGenotypeStorage my_unphased_genotypes;
        
        std::vector<const PRIVATE::genotype_info*> my_phased_genotypes_quick_lookup;
        std::vector<const PRIVATE::genotype_info*> my_unphased_genotypes_quick_lookup;
        
        PhasedGenotypeBySequence::iterator my_male_phased_genotype_begin;
        PhasedGenotypeBySequence::iterator my_male_phased_genotype_end;
        PhasedGenotypeBySequence::iterator my_female_phased_genotype_begin;
        PhasedGenotypeBySequence::iterator my_female_phased_genotype_end;

        UnphasedGenotypeBySequence::iterator my_male_unphased_genotype_begin;
        UnphasedGenotypeBySequence::iterator my_male_unphased_genotype_end;
        UnphasedGenotypeBySequence::iterator my_female_unphased_genotype_begin;
        UnphasedGenotypeBySequence::iterator my_female_unphased_genotype_end;
    };

};


//----------------------------------------------------------------------------
//  Class:      allele
//----------------------------------------------------------------------------
/// The allele class is an accessor into the genotype model, similar to a reference.
/// It provides information about specific alleles, their names, ids and frequencies.
///
/// Ids are assumed to count from 0 - (n-1)

class allele
{
  public:
    friend class genotype_model;
    friend class allele_iterator;
    friend class phased_genotype;
    friend class unphased_genotype;
    friend class PRIVATE::genotype_model_info;

    allele();
    allele(const allele& a);
    ~allele();

    allele& operator=(const allele& a);

    bool    operator==(const allele& rhs) const;
    bool    operator!=(const allele& rhs) const;
    bool    operator<  (const allele& rhs) const;
    
    double          frequency() const;
    const string&   name() const;
    uint            id() const;
    
    bool is_valid() const;
    
    bool is_sex_allele() const;
    
    bool is_null_y_allele() const;
    bool is_null_x_allele() const;

  protected:

    explicit allele(const PRIVATE::allele_info* info);

  private:
    const PRIVATE::allele_info*  my_info;
};

// Private constructors desired in this case.
//lint -esym(1704, phased_genotype::*, unphased_genotype::*)

//----------------------------------------------------------------------------
//  Class:      phased_genotype
//----------------------------------------------------------------------------
/// The phased_genotype is a reference class to the PRIVATE::genotype_model_info.  It
/// is assumed that the genotype is phase-known, such that the first allele
/// comes from a specific, known parent (generally the mother) and the
/// second from the other.
///
/// phased_genotype ids are integers and can be both positive and negative. 
/// Two phenotypes are opposite (ex "A<B" and "B<A") when the absolute value
/// of their ids are equal.  Homozygous genotypes are always positive, and
/// the negative id is not used.

class phased_genotype
{
  public:
    friend class allele;
    friend class child_genotype_set;
    friend class child_x_genotype_set;
    friend class genotype_model;
    friend class PRIVATE::genotype_model_info;
    friend class unphased_genotype;
    friend class phased_genotype_iterator;

    phased_genotype ();
    phased_genotype (const phased_genotype& g);
    ~phased_genotype ();

    phased_genotype&    operator=(const phased_genotype& g);

    bool    operator== (const phased_genotype& rhs) const;
    bool    operator!= (const phased_genotype& rhs) const;
    bool    operator<  (const phased_genotype& rhs) const;

    bool    equivalent(const phased_genotype& g) const;
    bool    equivalent(const unphased_genotype& g) const;
    bool    homozygous() const;

    allele          allele1() const;
    allele          allele2() const;
    const string    name() const;
    string          name(char c) const;
    uint            get_id() const;
    double          frequency() const;
    
    bool            is_flippable() const;

    phased_genotype   get_flipped_phased_genotype()      const;
    unphased_genotype get_equivalent_unphased_genotype() const;
    
    bool is_valid() const;
    
    bool is_sex_specific()      const;
    bool is_male_compatible()   const;
    bool is_female_compatible() const;
    
  private:

    const PRIVATE::genotype_info* my_info;
    
    phased_genotype(const PRIVATE::genotype_info* info);
};


//----------------------------------------------------------------------------
//  Class:      unphased_genotype
//----------------------------------------------------------------------------
/// The unphased_genotype is a reference class to the PRIVATE::genotype_model_info. 
/// It is assumed that the genotype is unphased, such that the two alleles
/// do not come from any specific parent.
///
/// unphased_genotype ids are uints and can only be positive.  The two
/// phased_genotype ids equivalent to the unphased_genotype have the
/// absolute values of their ids equal to the unphased_genotype id.
class unphased_genotype
{
  public:
    friend class child_genotype_set;
    friend class child_x_genotype_set;
    friend class genotype_model;
    friend class PRIVATE::genotype_model_info;
    friend class phased_genotype;
    friend class unphased_genotype_iterator;

    unphased_genotype();
    unphased_genotype(const unphased_genotype& g);
    ~unphased_genotype();

    unphased_genotype&  operator=(const unphased_genotype& g);

    bool    operator==(const unphased_genotype& rhs) const;
    bool    operator!=(const unphased_genotype& rhs) const;
    bool    operator< (const unphased_genotype& rhs) const;

    bool    equivalent(const phased_genotype& g) const;
    bool    equivalent(const unphased_genotype& g) const;
    bool    homozygous() const;
    
    size_t  get_equivalent_phased_genotype_count() const;

    allele          allele1() const;
    allele          allele2() const;
    string          name() const;
    string          name(char c) const;
    uint            get_id() const;
    double          frequency() const;

    phased_genotype   get_equivalent_phased_genotype1()      const;
    phased_genotype   get_equivalent_phased_genotype2()      const;
    
    bool is_valid() const;
    
    bool is_sex_specific()      const;
    bool is_male_compatible()   const;
    bool is_female_compatible() const;
    
  private:
    
    const PRIVATE::genotype_info* my_info;
    
    unphased_genotype(const PRIVATE::genotype_info* info);
};

//----------------------------------------------------------------------------
//  Class:      allele_iterator
//----------------------------------------------------------------------------
/// Bi-directional iterator through the alleles contained within a
/// PRIVATE::genotype_model_info.

class allele_iterator : public
    boost::iterator_adaptor
          <allele_iterator,
           std::vector<PRIVATE::allele_info>::const_iterator,
           allele,
           boost::use_default,
           allele>
{
public:
    
  friend class genotype_model;
  friend class PRIVATE::genotype_model_info;

  allele_iterator();
  allele_iterator(const allele_iterator&);

protected:
  friend class boost::iterator_core_access;

  explicit allele_iterator(const std::vector<PRIVATE::allele_info>::const_iterator&);

  allele dereference() const;
};

//----------------------------------------------------------------------------
//  Class:      phased_genotype_iterator
//----------------------------------------------------------------------------
/// Bi-directional iterator through the phased genotypes contained within a
/// PRIVATE::genotype_model_info.

class phased_genotype_iterator : public
    boost::iterator_adaptor
          <phased_genotype_iterator,
           PRIVATE::PhasedGenotypeBySequence::const_iterator,
           phased_genotype,
           boost::use_default,
           phased_genotype>
{
public:
    
  friend class genotype_model;

  phased_genotype_iterator();
  phased_genotype_iterator(const phased_genotype_iterator&);

protected:
  friend class boost::iterator_core_access;

  explicit phased_genotype_iterator(const PRIVATE::PhasedGenotypeBySequence::const_iterator&);

  phased_genotype dereference() const;
};

//----------------------------------------------------------------------------
//  Class:      unphased_genotype_iterator
//----------------------------------------------------------------------------
/// Bi-directional iterator through the unphased genotypes contained within
/// a PRIVATE::genotype_model_info.

class unphased_genotype_iterator : public
    boost::iterator_adaptor
          <unphased_genotype_iterator,
           PRIVATE::UnphasedGenotypeBySequence::const_iterator,
           unphased_genotype,
           boost::use_default,
           unphased_genotype>
{
public:
    
  friend class genotype_model;

  unphased_genotype_iterator();
  unphased_genotype_iterator(const unphased_genotype_iterator&);

protected:
  friend class boost::iterator_core_access;

  explicit unphased_genotype_iterator(const PRIVATE::UnphasedGenotypeBySequence::const_iterator&);

  unphased_genotype dereference() const;
};

//----------------------------------------------------------------------------
//  Class:      genotype_model
//----------------------------------------------------------------------------
/// The genotype model is the primary class in this library.  Primary
/// responsibilities are to provide information about the genotypes
/// contained within.  Counts of alleles, phased and unphased genotypes,
/// iteration through them, their frequencies, etc.  Also included are
/// structural information such as the separators used, what the missing
/// allele is (if any), if dynamic alleles are used, etc.

class genotype_model
{
  public:
    friend class penetrance_model;

    genotype_model();
    genotype_model(const string& name);
    genotype_model(const genotype_model& m);
    ~genotype_model();

    genotype_model&   operator=(const genotype_model& m);

    //- Counts of various contained elements.
    //
    uint    allele_count() const;
    uint    phased_genotype_count() const;
    uint    unphased_genotype_count() const;

    //- Iteration over alleles.
    //
    allele_iterator     allele_begin() const;
    allele_iterator     allele_end() const;

    //- Iteration over genotypes.
    //
    phased_genotype_iterator    phased_genotype_begin() const;
    phased_genotype_iterator    phased_genotype_end() const;
    unphased_genotype_iterator  unphased_genotype_begin() const;
    unphased_genotype_iterator  unphased_genotype_end() const;

    phased_genotype_iterator    phased_genotype_begin(const MPED::SexCode&) const;
    phased_genotype_iterator    phased_genotype_end(const MPED::SexCode&) const;
    unphased_genotype_iterator  unphased_genotype_begin(const MPED::SexCode&) const;
    unphased_genotype_iterator  unphased_genotype_end(const MPED::SexCode&) const;

    //- Element indexing by ID.
    //
    allele               get_allele            (uint id) const;
    phased_genotype      get_phased_genotype   (uint id) const;
    unphased_genotype    get_unphased_genotype (uint id) const;

    phased_genotype      get_phased_genotype   (allele a1, allele a2) const;
    unphased_genotype    get_unphased_genotype (allele a1, allele a2) const;

    //- Direct element lookup by name.
    //
    allele              get_allele(const string& name) const;
    phased_genotype     get_phased_genotype(const string& name) const;
    unphased_genotype   get_unphased_genotype(const string& name) const;

    allele get_sex_specific_allele() const;

    //- Other properties and services.
    //
    genotype_model  clone() const;
    const string&   name() const;
    const string&   missing_allele_name() const;
    bool            dynamic_alleles() const;
    
    char        backward_separator() const;
    char        forward_separator() const;
    char        unphased_separator() const;
    string      separators() const;

    string      parse_genotype_name(const string& gname, char sep='\0') const;
    string      parse_genotype_name(const string& gname, string& a1, string& a2, Ordering& order, char sep='\0') const;

    // Added for X-linkage - yjs Mar. 2002
    //
    bool        is_x_linked()  const;
    bool        is_y_linked()  const;
    bool        is_autosomal() const;
    
    GenotypeModelType get_model_type() const;
    
    void        set_model_type(GenotypeModelType);
    
    // Added for allele_frequency adjustment - yjs Dec. 2002
    void        modify_allele_frequency(allele id, double freq=0.0);
    void        modify_allele_frequency(allele_iterator ai, double freq=0.0);
    void        modify_allele_frequency(const string& name, double freq=0.0);

    //- Modifiers.
    //
    void        add_allele(const string& name, double freq=0.0);
    void        normalize();
    void        clear();

    void        mark_for_remap(const string& name);
    void        remap();
    void        remap(genotype_model& m);

    void        set_backward_separator(char c);
    void        set_forward_separator(char c);
    void        set_unphased_separator(char c);
    void        set_separators(const string& seps);

    void        set_name(const string& name);
    void        set_missing_allele_name(const string& name);
    void        set_dynamic_alleles(bool d);
    
  private:
    typedef PRIVATE::allele_info             allele_info;
    typedef std::list<string>                allele_name_buffer;
    typedef std::vector<allele_info>         allele_vector;
    typedef std::map<string, uint>           allele_map;

    boost::shared_ptr<PRIVATE::genotype_model_info>  my_info;

    //- Misc. methods.
    //
    string      parse_allele_name(const string& aname) const;

    void        uniquify();

    static void remap(const genotype_model& src, genotype_model& dst);

};

//----------------------------------------------------------------------------
//  Class:      child_genotype_set
//----------------------------------------------------------------------------
/// This class is both a data structure and a functional object.  Given two
/// parental genotypes, creates a list of all possible child genotypes. 
/// Child genotypes are phased with respect to their parents, such that the first
/// allele comes from the first parent, and the second allele from the second.
/// Typically, the first parent is the mother.
///
/// In addition, the set can reduce itself to the minimal set of genotypes in
/// the presence of homozygous parents. The 'reduce' flag controls whether 
/// the child gneotypes are reduced to
/// the minimal set of genotypes given the parents (in the presence of
/// homozygous parents, not all genotypes need be unique) Note that the
/// children's genotypes are phased, so genotypes such as A/B and B/A are
/// not considered equivalent and are *not* reduced when reduce is true.

class child_genotype_set
{
  public:
  
    class iterator : public std::forward_iterator_tag
    {
      //lint -e613

      public:
        friend class child_genotype_set;

        typedef const phased_genotype*         pointer;
        typedef const phased_genotype&         reference;

        iterator()                                                  : pdata(0) {}

        bool        operator==(const iterator& i) const            { return pdata == i.pdata; }
        bool        operator!=(const iterator& i) const            { return pdata != i.pdata; }
        reference   operator*() const                              { return *pdata; }
        pointer     operator->() const                             { return pdata;  }

        iterator&   operator++()                                   { ++pdata; return *this; }
        iterator    operator++(int)                                { iterator  tmp(*this); ++pdata; return tmp; }

      private:
        const phased_genotype* pdata;

        //lint --e{1704}
        iterator(const phased_genotype* p)                                 : pdata(p) {}

        //lint +e613
    };


    child_genotype_set();
    child_genotype_set(const phased_genotype& a,
                       const phased_genotype& b,
                       bool reduce = false);
    child_genotype_set(const unphased_genotype& a,
                       const unphased_genotype& b,
                       bool reduce = false);

    uint            size() const;

    iterator        begin() const;
    iterator        end() const ;

    const phased_genotype& operator[](uint i) const;

    // True if one of the possible child genotypes is equivalent to the one given.
    bool contains(const phased_genotype&)   const;
    bool contains(const unphased_genotype&) const;

    void reduce(); // Reduce to minimal genoset if not previously done.

  private:
    phased_genotype my_genotypes[4];
    uint            my_size;
};

void genotype_model_test(ostream& o, const genotype_model& m);
void genotype_model_print(ostream& o, const genotype_model& m);

} // End namespace MLOCUS
} // End namespace SAGE

#include "mlocus/genotype.ipp"

#endif
