#ifndef REL_PAIR_H
#define REL_PAIR_H

//============================================================================
// File:      relpair.h                                                        
//                                                                          
// Author:    Dan Baechle & Kevin Jacobs                                    
//                                                                          
// History:   7/00   created.  - djb
//            4/4/01 modified so that relative_pair class includes connecting 
//                   members.  -djb
//            8/9/04 corrected error where by order of cousin connectors was
//                   - djb
//            1/18/7 added EVERY type which consists of all unique pairs in a subpedigree
//                   regardless of relationship.  -djb
//                                                                          
// Notes:     Declares the following classes -                               
//              pair_generator 
//                relative_pair
//                base_iterator
//                iterator
//                const_iterator 
//              ind_filter_trait
//              pair_filter_trait                                 
//              ind_filter
//              pair_filter    
//              filtering_pair_generator 
//              filtering_pair_generator_rep
//                base_filtering_iterator
//                filtering_iterator
//                const_filtering_iterator
//
//             In order to provide an iterator as well as a const_iterator,
//             the pair_generator class takes a RefPedigree* as a constructor
//             argument instead of a const RefPedigree*.  This means that client
//             code needing a const_iterator must supply a RefPedigree* to create
//             a pair_generator, even if it means an ugly const_cast.
//
//             Maybe this code should be rewritten so that there is a pair_generator
//             which produces an iterator and a const_pair_generator which produces
//             a const_iterator.  This would mean duplicating a lot of code, however.
//
//                                                                             - djb
//             
// Copyright (c) 2000 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include <vector>
#include <limits>
#include <memory>
#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>
#include "rped/rped.h"
#include "pairs/reldefs.h"
#include "pairs/reltype.h"
#include "pairs/relmatrix.h"

namespace SAGE
{

class filtering_pair_generator_rep;

//---------------------------------------------------------------------------
//  Class:     pair_generator                                                         
//                                                                          
//  Purpose:   Create pairs and iterators of specified relationship 
//             type(s) for a given pedigree.          
//                                                                          
//---------------------------------------------------------------------------
//
// ***** If pedigree contains loops, pairs may double counted or overlooked.  This
//       needs to be analysed for each pair type/loop type combination.
//        
//
class pair_generator
{
  public:
  
    // - EVERY pair_type yeilds all unique, unordered pairs without regard to relationship.
    //   If EVERY is set, other pair types are desregarded.
    //
    enum pair_type  { PARENTAL = 1, SIBSIB, SISSIS, BROBRO, BROSIS, GRANDP, AVUNC,
                      HALFSIB, COUSIN, NULL_TYPE, EVERY };
    
    enum mask { PARENTAL_MASK     =   1, SIBSIB_MASK    =    2, SISSIS_MASK =   4,
                  BROBRO_MASK     =   8, BROSIS_MASK    =   16, GRANDP_MASK =  32,
                   AVUNC_MASK     =  64, HALFSIB_MASK   =  128, COUSIN_MASK = 256, 
                   NULL_TYPE_MASK = 512, EVERY_MASK     = 1024                    };
                   
    // - all types in the pair_type enumeration except EVERY and NULL_TYPE.
    //
    static const unsigned int ALL_TYPES = PARENTAL_MASK | SIBSIB_MASK  | SISSIS_MASK |
                                          BROBRO_MASK   | BROSIS_MASK  | GRANDP_MASK |
                                          AVUNC_MASK    | HALFSIB_MASK | COUSIN_MASK  ;
                   
    static std::string   pair_type_to_string(pair_type type);
    static pair_type     mask_to_pair_type(mask m);
    static mask          pair_type_to_mask(pair_type p);
    
    class relative_pair;
    class base_iterator;
    class iterator;
    class const_iterator;
    friend class base_iterator;
    friend class iterator;
    friend class const_iterator;

    // Constructor/destructor.
    pair_generator(unsigned int types = ALL_TYPES);
    pair_generator(RPED::RefPedigree* p, unsigned int types = ALL_TYPES);
    
    // iteration
    iterator         begin();
    iterator         end();
    const_iterator   begin()  const;
    const_iterator   end()    const;
                                               
    // gets
    unsigned int     types()                            const;
    bool             check_type(mask m)                 const;
    pair_type        first_type()                       const;
    pair_type        next_type(pair_type current_type)  const;
    
    // sets
    void             set_pedigree(RPED::RefPedigree* pedigree);
    void             set_types(unsigned int types);
    void             set_type(pair_type p);
    void             add_type(mask m);
    bool             operator ==(const pair_generator& other)  const;
    bool             operator !=(const pair_generator& other)  const;
    
  private:              
  
    // Data members.                          
    RPED::RefPedigree*     my_p;
    unsigned int     my_types;      // flag whose bits represent values of the
                                    // enumeration relative_pair::pair_type 
    pair_type        my_first_type;

  public:                                 
    //---------------------------------------------------------------------------
    //  Class:     pair_generator::relative_pair                                                        
    //                                                                          
    //  Purpose:   Represent a pair of individuals in a pedigree for the        
    //             purpose of calculating general pedigree statistics.
    //
    //  Note:  member_one, member_two, connector_one, and connector two 
    //         are pair-type specific as follows:
    //
    //         pair type    member_one    member_two   connector_one      connector_two
    //         =========    ==========    ==========   =============      =============
    //         parental     parent        child        NULL               NULL
    //         sibsib       sib           sib          parent 1           parent 2
    //         sissis       sib           sib          parent 1           parent 2
    //         brobro       sib           sib          parent 1           parent 2
    //         brosis       sib           sib          parent 1           parent 2
    //         grandp       gparent       gchild       parent             NULL
    //         avunc        uncle*        nephew*      parent of nephew*  NULL
    //         halfsib      hsib          hsib         common parent      NULL
    //         cousin       cous 1        cous 2       parent of cous 1   parent of cous 2
    //
    //    * uncle means uncle or aunt and nephew means nephew or neice.
    //                                                                          
    //---------------------------------------------------------------------------
    //
    class relative_pair
    {
      friend class base_iterator;
      friend class iterator;
      friend class const_iterator;
      friend std::ostream& operator <<(std::ostream& out, relative_pair pr);

      public:
        typedef RPED::RefPedigree::member_pointer member;
                          
        // gets
        pair_type              type()                   const; 
        
        
        member                            member_one()     const;
        member                            member_two()     const;
        member                            connector_one()  const;
        member                            connector_two()  const;
        
        bool   operator ==(const relative_pair& other)  const;
        bool   operator !=(const relative_pair& other)  const;

      private:
        // Constructor/destructor
        relative_pair();
        relative_pair(member member_one, member member_two, 
                      member connector_one, member connector_two, pair_type type);
        
        // Data members.     
        member       my_member_one;
        member       my_member_two;
        member       my_connector_one;
        member       my_connector_two;
        pair_type    my_type;

    }; // end class relative_pair

  
    //---------------------------------------------------------------------------
    //  Class:     pair_generator::base_iterator                                                         
    //                                                                          
    //  Purpose:   Base for iterator and const_iterator classes.
    //                                                                          
    //---------------------------------------------------------------------------
    //
    class base_iterator
    {
        friend class pair_generator;
        friend class filtering_pair_generator_rep;
    
      public:                
        // Constructor/destructor.
         base_iterator(const base_iterator& other);
         base_iterator& operator  =(const base_iterator& other);
        virtual  ~base_iterator();
      
        // Iteration
         bool            operator ==(const base_iterator& other)  const;
         bool            operator !=(const base_iterator& other)  const;
        
      protected:
        enum parent_number { P1, P2 };
        enum grandp_number { GP1, GP2 };
      
        // Constructor/destructor
         base_iterator();
         base_iterator(pair_generator* generator,
                      RPED::RefPedigree::family_iterator family);
                      
         base_iterator(const pair_generator* generator,
                      RPED::RefPedigree::family_iterator family);
                      
        // Functions used in the implementation of operator ++().         
         void           init_family();   // Set set state variables to search a new family.
         void                  seek();   // Find next pair.  

         void            init_every();   // When pair_type == EVERY.
         void        set_every_pair();
         void                 every();
         void            seek_every();        
        
         void     set_parental_pair();
         void              parental();
        
         void             incr_sibs();
         void      set_sibling_pair();      
         void                sibsib();
         void                sissis();
         void                brobro();
         void                brosis();
        
         void          incr_grandps();
         void       set_grandp_pair();
         void                grandp();
        
         void           incr_avuncs();
         void        set_avunc_pair();
         void                 avunc();
        
         void      set_halfsib_pair();
         void incr_parents_mates_hs();
         void      incr_children_hs();
         void      incr_halfsibs_hs();
         void         init_halfsibs();
         bool         other_parent(RPED::RefPedigree::member_pointer mate);
         void               halfsib();
        
         void       set_cousin_pair();
         void         incr_avuncs_c();
         void   incr_avuncs_mates_c();
         void       incr_children_c();
         void        incr_cousins_c();
         void          init_cousins();
         void                cousin();
        
         void             null_type();              
      
      
        // Data members.
        //
        std::auto_ptr<pair_generator>        my_generator;
        relative_pair                        my_pair;             // temporary pair
        
        // iterators
        RPED::RefPedigree::family_iterator         my_family;
        RPED::RefPedigree::offspring_iterator      my_children;
        RPED::RefPedigree::offspring_iterator      my_fullsibs;  
        RPED::RefPedigree::offspring_iterator      my_cousins;
        RPED::RefPedigree::offspring_iterator      my_halfsibs;
        RPED::RefPedigree::mate_iterator           my_parents_mates;
        RPED::RefPedigree::mate_iterator           my_avuncs_mates;
        RPED::RefPedigree::sibling_iterator        my_avuncs;
        
        // For pair_type, EVERY.
        RPED::RefMultiPedigree::subpedigree_iterator    my_subped;
        RPED::RefMultiPedigree::member_iterator         my_member_one;
        RPED::RefMultiPedigree::member_iterator         my_member_two;
        
        // member pointers
        RPED::RefPedigree::member_pointer          my_parent;
        RPED::RefPedigree::member_pointer          my_grandp;
        
        // other state variables
        parent_number                        my_parent_number;
        grandp_number                        my_grandp_number;
        pair_generator::pair_type            my_current_type;
        bool                                 pair_valid;
        bool                                 at_end;
        bool                                 halfsibs_init;      // Has search for halfsibs been initiated?
        bool                                 cousins_init;       // Has search for cousins been initiated?
        
    }; // end class base_iterator                                               
    
  
    //---------------------------------------------------------------------------
    //  Class:     pair_generator::iterator                                                         
    //                                                                          
    //  Purpose:   Access pairs of pedigree members.
    //                                                                          
    //---------------------------------------------------------------------------
    //
    class iterator : public base_iterator
    {
        friend class pair_generator;
        friend class filtering_pair_generator_rep;
    
      public:
        // Constructor/destructor.
         iterator(const iterator& other);
         iterator& operator  =(const iterator& other);
         ~iterator();
      
        // Iteration
               iterator*        operator ++();
         relative_pair&   operator  *();
         relative_pair*   operator ->();
        
      private:
        // Constructor/destructor
         iterator();
        
        iterator(const pair_generator* generator,
                 RPED::RefPedigree::family_iterator family);
        
    }; // end class iterator
                                                   
    
    //---------------------------------------------------------------------------
    //  Class:     pair_generator::const_iterator                                                         
    //                                                                          
    //  Purpose:   Access pairs of pedigree members.
    //                                                                          
    //---------------------------------------------------------------------------
    //
    class const_iterator : public base_iterator
    {
        friend class pair_generator;
        friend class filtering_pair_generator_rep;
    
      public:
        // Constructor/destructor.
         const_iterator(const const_iterator& other);
         const_iterator& operator  =(const const_iterator& other);
         ~const_iterator();
      
        // Iteration
               const_iterator*        operator ++();
         const relative_pair&   operator *();
         const relative_pair*   operator ->();
        
      private:
        // Constructor/destructor
         const_iterator();

        const_iterator(const pair_generator* generator,
                       RPED::RefPedigree::family_iterator family);
        
    }; // end class const_iterator
    
}; // end class pair_generator

//---------------------------------------------------------------------------
//  Class:     ind_filter_trait                                                         
//                                                                          
//  Purpose:   Represent a trait for purposes of filtering a individual.        
//                                                                          
//---------------------------------------------------------------------------
//
class ind_filter_trait
{
  public:
    // Constructor/destructor
    ind_filter_trait(size_t t  =  (size_t)(-1),
                 double min    = -std::numeric_limits<double>::infinity(),
                 double max    =  std::numeric_limits<double>::infinity());

    // gets
    size_t            get_trait()           const;
    double            get_min()             const;
    double            get_max()             const;
    double            get_affected_min()    const;
    double            get_affected_max()    const;
    double            get_unaffected_min()  const;
    double            get_unaffected_max()  const;

    // sets
    void              clear();
    void              set_trait(size_t trait);
    void              set_min(double min);
    void              set_max(double max);
    void              set_affected_min(double affected_min);
    void              set_affected_max(double affected_max);
    void              set_unaffected_min(double unaffected_min);
    void              set_unaffected_max(double unaffected_max);
    void              set_affected_range(double min, double max);   
    void              set_unaffected_range(double min, double max); 
    void              set_threshold(double threshold);
    
    bool              operator ==(const ind_filter_trait& other)  const;   
    bool              operator !=(const ind_filter_trait& other)  const;

    // - Filtration functions.  Function ranges are as follows:
    //     in_range        [min, max]
    //     unaffected      [unaffected_min, unaffected_max]
    //     affected        (affected_min, affected_max)
    //
    //  ***** A more complete implementation would include emin's and emax's for these ranges.
    //  Min's and Max's would be understood to be inclusive.  Functions would then be defined to either
    //  select the largest range based on the two sets of min's and max's or meet all criteria.
    //
    bool              informative (const RPED::RefPedigree::member_pointer member)  const;
    bool              in_range    (const RPED::RefPedigree::member_pointer member)  const;
    bool              affected    (const RPED::RefPedigree::member_pointer member)  const;
    bool              unaffected  (const RPED::RefPedigree::member_pointer member)  const;

  protected:
    // Data members.  
    size_t            my_trait;
    double            my_min;               // Minimum meaningful trait.
    double            my_max;               // Maximum meaningful trait.
    double            my_affected_min;
    double            my_affected_max;
    double            my_unaffected_min;
    double            my_unaffected_max;

}; // end class ind_filter_trait


//---------------------------------------------------------------------------
//  Class:     pair_filter_trait                                                         
//                                                                          
//  Purpose:   Represent a trait for purposes of filtering a pair.        
//                                                                          
//---------------------------------------------------------------------------
//
class pair_filter_trait : private ind_filter_trait
{
  public:
    using ind_filter_trait::get_trait;
    using ind_filter_trait::get_min;
    using ind_filter_trait::get_max;
    using ind_filter_trait::get_affected_min;
    using ind_filter_trait::get_affected_max;
    using ind_filter_trait::get_unaffected_min;
    using ind_filter_trait::get_unaffected_max;
    using ind_filter_trait::set_trait;
    using ind_filter_trait::set_min;
    using ind_filter_trait::set_max;
    using ind_filter_trait::set_affected_min;
    using ind_filter_trait::set_affected_max;
    using ind_filter_trait::set_unaffected_min;
    using ind_filter_trait::set_unaffected_max;
    using ind_filter_trait::set_affected_range;   
    using ind_filter_trait::set_unaffected_range; 
    using ind_filter_trait::set_threshold;
  
    enum affection_status { CONCORD_UNAFF = 1, DISCORD, CONCORD_AFF, UNINFORM, NULL_STATUS };
    enum mask             { CONCORD_UNAFF_MASK = 1, DISCORD_MASK  = 2, 
                              CONCORD_AFF_MASK = 4, UNINFORM_MASK = 8, NULL_STATUS_MASK = 16 };
                              
    static const unsigned int ALL_INFORMATIVE_STATUSES =   CONCORD_UNAFF_MASK | DISCORD_MASK
                                                         | CONCORD_AFF_MASK; 
    static affection_status  mask_to_status(mask m);
    static mask              status_to_mask(affection_status s);
  
    // Constructor/destructor
    pair_filter_trait(size_t t        =  (size_t)(-1),
                      unsigned int a  =  ALL_INFORMATIVE_STATUSES,
                      double min      = -std::numeric_limits<double>::infinity(),
                      double max      =  std::numeric_limits<double>::infinity());

    // gets
    unsigned int      get_status()          const;

    // sets
    void              clear();
    void              set_statuses(unsigned int status);
    void              set_status(affection_status s);
    
    bool              operator ==(const pair_filter_trait& other)  const;   
    bool              operator !=(const pair_filter_trait& other)  const;
    
    bool              valid(const pair_generator::relative_pair& pair)  const;

  private:
    // Data members.  
    unsigned int      my_status;                // Affection status.
    bool              valid_status;             // Valid affection status.

}; // end class pair_filter_trait


//---------------------------------------------------------------------------
//  Class:     ind_filter                                                         
//                                                                          
//  Purpose:   Filter individuals by specified characteristics.
//                                                                          
//---------------------------------------------------------------------------
//
class ind_filter
{
  public:
    typedef std::list<ind_filter_trait> ind_filter_traits;
  
    // Constructor/destructor
    ind_filter();
    ind_filter(ind_filter_trait trait);
  
    // Trait specification
    void   add_trait(ind_filter_trait trait);
    void   set_trait(ind_filter_trait trait);   // Adds trait only if not a duplicate.
    void   clear_traits();
    
    // Filtration
    bool   informative (const RPED::RefPedigree::member_pointer member)  const;
    bool   in_range    (const RPED::RefPedigree::member_pointer member)  const;
    bool   affected    (const RPED::RefPedigree::member_pointer member)  const;
    bool   unaffected  (const RPED::RefPedigree::member_pointer member)  const;

  private:
  
    // Data members.
    ind_filter_traits my_traits;
    
}; // end ind_filter


//---------------------------------------------------------------------------
//  Class:     pair_filter                                                         
//                                                                          
//  Purpose:   Filter pairs by specified characteristics.
//                                                                          
//---------------------------------------------------------------------------
//
class pair_filter
{
  public:
    typedef std::list<pair_filter_trait> pair_filter_traits;
  
    // Constructor/destructor
    pair_filter();
    pair_filter(pair_filter_trait trait);
  
    // Trait specification
    void   add_trait(pair_filter_trait trait);
    void   set_trait(pair_filter_trait trait);
    void   clear_traits();
    
    // Filtration
    bool   valid(const pair_generator::relative_pair& pair) const;    

  private:
  
    // Data members.
    pair_filter_traits my_traits;
    
}; // end pair_filter


class filtering_pair_generator;

//---------------------------------------------------------------------------
//  Class:     filtering_pair_generator_rep
//                                                                          
//  Purpose:   The "real" filtering pair generator.  Contains state info. 
//             COPY-ON-WRITE.
//                                                                         
//---------------------------------------------------------------------------
//
class filtering_pair_generator_rep
{
  typedef boost::shared_ptr<filtering_pair_generator_rep> rep_ptr_type;
  
  public:

    friend class filtering_pair_generator;

    //----------------------------------------------------------------------------
    //  Class:    filtering_pair_generator_rep::cache
    //                                                                          
    //  Purpose:  Store information about validity of pedigree members w. regard
    //            to a specific filter.
    //
    //  Note:     ####### Not implemented. #######
    //                                                                          
    //----------------------------------------------------------------------------
    //
    class cache
    {
      public:
        enum cache_status { UNKNOWN, VALID, INVALID };
        
        cache_status     check_status(RPED::RefPedigree::member_pointer member);
        void             cache_add(RPED::RefPedigree::member_pointer member, cache_status status);        
        
      private:
        // some type of container...
            
    }; // end class cache 

  
  // - Access to following member functions is soley thru filtering_pair_generator class.
  //  
  private:   
    // Constructor/destructor
    filtering_pair_generator_rep();
    filtering_pair_generator_rep(const pair_generator& generator, pair_filter filter);
    filtering_pair_generator_rep(const filtering_pair_generator_rep& other);
    
    // Intentionally not implemented.
    filtering_pair_generator_rep& operator =(const filtering_pair_generator_rep& other);
    void                                  set_generator(const pair_generator& generator);
    void                                  set_filter(const pair_filter& filter);
    
    // Data members.
    pair_generator                        my_generator;
    pair_generator::iterator              my_end_iterator;
    pair_filter                           my_filter;
    cache                                 my_cache;     // ####### Not implemented. #######
    
  public:
    // Accessors.
    pair_generator*                       get_generator();
    pair_generator::iterator*             get_end();
    pair_filter*                          get_filter(); 
    cache*                                get_cache();  // ####### Not implemented. #######
    
  public:    
    //---------------------------------------------------------------------------
    //  Class:     filtering_pair_generator_rep::base_iterator  
    //                                                                          
    //  Purpose:   Base for iterator and const_iterator classes.
    //                                                                          
    //---------------------------------------------------------------------------
    //
    class base_iterator
    {
        friend class filtering_pair_generator;
    
      public:
      
        //Constructor/destructor.
         base_iterator(const base_iterator& other);
         base_iterator& operator =(const base_iterator& other);
        virtual  ~base_iterator();
        
        // Iteration.
         bool            operator ==(const base_iterator& other)   const;
         bool            operator !=(const base_iterator& other)   const;
        
      protected:
        // Constructor/destructor
         base_iterator(rep_ptr_type rep, pair_generator::iterator iter);      
        
        // Data members.
        rep_ptr_type                      my_rep;
        pair_generator::iterator          my_iterator;
        bool                              initializing;
        
    }; // end class base_iterator                                               
    
    
    //---------------------------------------------------------------------------
    //  Class:     filtering_pair_generator_rep::iterator   
    //                                                                          
    //  Purpose:   Access pairs of pedigree members.
    //                                                                          
    //---------------------------------------------------------------------------
    //
    class iterator : public base_iterator
    {
        friend class filtering_pair_generator;
    
      public:
        // Constructor/destructor.
         iterator(const iterator& other);
         iterator& operator =(const iterator& other);
                                                            
         ~iterator();
      
        // Iteration
               iterator*                            operator ++();
         pair_generator::relative_pair&       operator  *();
         pair_generator::relative_pair*       operator ->();
        
        // For testing.
         void                                 print_references();
        
      private:
        // Constructor/destructor
         iterator(rep_ptr_type rep, pair_generator::iterator iter);
        
    }; // end class iterator
                                                   
    
    //---------------------------------------------------------------------------
    //  Class:     filtering_pair_generator_rep::const_iterator   
    //                                                                          
    //  Purpose:   Access pairs of pedigree members.
    //                                                                          
    //---------------------------------------------------------------------------
    //
    class const_iterator : public base_iterator
    {
        friend class filtering_pair_generator;
    
      public:
        // Constructor/destructor.
         const_iterator(const const_iterator& other);
         const_iterator& operator =(const const_iterator& other);
            
         ~const_iterator();
      
        // Iteration
         const_iterator*                            operator ++();
         const pair_generator::relative_pair&       operator  *();
         const pair_generator::relative_pair*       operator ->();
        
        // Testing.
        void                                 print_references();
        
      private:
        // Constructor/destructor
         const_iterator(rep_ptr_type rep, pair_generator::iterator iter);
    }; // end class const_iterator
  
}; // end class filtering_pair_generator_rep


//---------------------------------------------------------------------------
//  Class:     filtering_pair_generator
//                                                                          
//  Purpose:   Create iterators of a specified relationship 
//             type and set of trait values for a given pedigree.          
//             Contains interface to be used by client code.             
//                                        
//---------------------------------------------------------------------------
//
class filtering_pair_generator
{
    typedef filtering_pair_generator_rep::rep_ptr_type     rep_ptr_type;
  
  public:
    typedef filtering_pair_generator_rep::iterator         iterator;
    typedef filtering_pair_generator_rep::const_iterator   const_iterator;
  

    // Constructor/destructor
    filtering_pair_generator();
    filtering_pair_generator(const pair_generator& generator, pair_filter filter);
    filtering_pair_generator(const filtering_pair_generator& other);
    filtering_pair_generator& operator =(const filtering_pair_generator& other);
    ~filtering_pair_generator();
    
    iterator           begin();
    iterator           end();
    const_iterator     begin()             const;
    const_iterator     end()               const;
    
    
    // Gets                                   
    pair_generator  get_pair_generator()   const;        
    pair_filter             get_filter()   const;
    
    // Sets
    void             set_pair_generator(pair_generator& generator);
    void             set_filter(const pair_filter& filter);               
    
    // For testing.
    void            print_references();                                  
    
  private:
    filtering_pair_generator_rep::rep_ptr_type      make_unique();  
    
    filtering_pair_generator_rep::rep_ptr_type      my_rep;

}; // end class filtering_pair_generator


// Helper functions not belonging to one class.
bool operator ==(const pair_generator::iterator& first, const pair_generator::const_iterator& second);
bool operator !=(const pair_generator::iterator& first, const pair_generator::const_iterator& second);
bool operator ==(const pair_generator::const_iterator& first, const pair_generator::iterator& second);
bool operator !=(const pair_generator::const_iterator& first, const pair_generator::iterator& second);

bool operator ==(const filtering_pair_generator::iterator& first, const filtering_pair_generator::const_iterator& second);
bool operator !=(const filtering_pair_generator::iterator& first, const filtering_pair_generator::const_iterator& second);
bool operator ==(const filtering_pair_generator::const_iterator& first, const filtering_pair_generator::iterator& second);
bool operator !=(const filtering_pair_generator::const_iterator& first, const filtering_pair_generator::iterator& second);

#include "pairs/relpair.ipp"

}  // end namespace SAGE

#endif
