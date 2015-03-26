#ifndef IBD_H
#define IBD_H

//==========================================================================
//  IBD sharing storage interface  
//  -- General data structures for allele sharing IBD
//                                    
//
//  Author: Kevin Jacobs (jacobs@darwin.cwru.edu)
//
//  History:   0.1  kbj Initial implementation              Aug 07 98
//             0.2  kbj Added per pedigree interface        Aug 21 98
//             Interface updated for FPED only              yjs Apr 2004
//
//  Copyright (c) 1998  R.C. Elston
// -------------------------------------------------------------------------
// Introduction:
//   This class interface forms the basis for component access to storage
//   classes and algorithm classes which deal with IBD (inheritance by
//   descent) of pairs of individuals at one or more marker loci.
//   A typical storage class will implement the IBD class by supplying
//   several virtual function definitions that expose a subset of pedigree,
//   marker, and IBD storage information.
//
// Pedigree interface:
//   1) The only part of the pedigree data structure exposed
//      are pointers the pedigree individual_type, which are primarily
//      obtained through IBD::get_pair(ped,i1,i2,...).  If no pedigree is
//      present then NULL pointers will be returned, which makes it
//      extremely important that the pointers are checked pedigree
//      information is required.
//
// IBD storage:
//   1) The conventon for IBD is to store both f0 and f2 (where fn is the
//      probability of two individuals sharing 'n' allees IBD at a given
//      marker).  Consequently, f1 = 1 - f0 - f2.
//
//   2) Any IBD value that is not supplied or is in some way invalid should
//      be stored as std::numeric_limits<double>::quiet_NaN() (NaN = Not a
//      Number).  As a result, any IBD sharing value that is returned must
//      be tested with the SAGE::isnan() predicate to determine if it is valid
//      before it may be used.
//==========================================================================

#include "ibd/definitions.h"

namespace SAGE {

class IBD
{
  public:

    enum error_t {no_error=0, bad_pedigree, bad_ind1, bad_ind2, bad_ind_both};

    // Virtual destructor
    virtual ~IBD() { };

    // Perform any necessary actions to finalize a set of pairs.
    // Many torage structures have no build requirement and this function will
    // simply be a NOP.
    // Postcondition: built() == true iff build was successful
    virtual void build() = 0;

    // Return if the pairs have been finalized by a build().  If no build step
    // is required alwats return true.
    virtual bool built() const = 0;

    // Return true if any pedigree data is available (some storage objects and
    // analyses do not need a pedigree; just pairs.)
    virtual bool has_pedigree() = 0;

    // Return true if the given pedigree has data
    bool has_pedigree(const std::string &ped);

    // Inform the storage object that a possibly unknown marker is going
    // to be read in and return an index for it.
    virtual size_t add_marker(const std::string &name, double dis, gmodel_type mt);

    // Return the number of markers that are to be read in or are currently
    // stored.
    virtual size_t marker_count() const;

    // Translate a marker name to a marker index.  A value >= marker_count()
    // is returned if the marker_name is not found.
    virtual size_t marker_index(const std::string &name) const;

    // Translate a marker index (0..marker_count()-1) to its name.
    virtual string marker_name(size_t m) const;

    // Translate a marker index (0..marker_count()-1) to its ibd_marker_info.
    virtual const ibd_marker_info& get_marker_info(size_t m) const;

    // Get the number of pairs currently stored
    // Precondition:  built() == true
    virtual size_t pair_count() const = 0;

    // Retrieve a pair of individual_type pointers for a particular pair
    // specified by either index value or string names.
    // The pair does may not yet be added when get_pair is called.

    // Precondition:  built() == true   [for only get_pair(size_t)]
    virtual const id_pair get_pair(size_t i) const             = 0;
    virtual id_pair       get_pair(size_t i)                   = 0;

    // NB:  This is the main interface to the underlying pedigree
    //      data structure
    virtual const id_pair get_pair(const std::string &ped,
                                   const std::string &i1,
                                   const std::string &i2,
                                   error_t &e) const       = 0;

    virtual id_pair get_pair(const std::string &ped,
                             const std::string &i1,
                             const std::string &i2,
                             error_t &e)             = 0;

    // This variant is an interface between the pair data and pedigree data
    virtual bool get_pair(size_t i, std::string &ped,
                                    std::string &i1,
                                    std::string &i2) const = 0;

    // Returns false if the pair will not be used and all pair specific
    // information can be safely skipped; true otherwise.  The pair
    // _does not_ have to be added for use_pair to return true.
    virtual bool use_pair(size_t i) const = 0;

    virtual bool use_pair(const mem_pointer i1, const mem_pointer i2) const = 0;

    // Helper use_pair

    bool use_pair(const std::string &ped, const std::string &i1,
                  const std::string &i2, error_t &e) const;

    // Returns false if the pair is not stored (or does not contain valid
    // data that can be read or written to); true otherwise.
    // Precondition: built() == true
    virtual bool valid_pair(size_t i) const = 0;

    bool valid_pair(const mem_pointer i1, const mem_pointer i2) const;

    // Returns false if the pair cannot be invalidated; true if it can be or
    // it does not exist.
    // Precondition: built() == true

    virtual bool invalidate_pair(size_t i) const = 0;

    bool invalidate_pair(const mem_pointer i1, const mem_pointer i2) const;

    // Add a pair to a data object.  A pair _must_ be added before it
    // can referred to.  Adding a pair more than once has undefined
    // behavior.
    //
    // Returns an index value for the current pair.  This index is only
    // valid until the next call to add_pair or a build is performed.
    // If a pair cannot be added (ie. some IBD objects will not allow
    // anonymous pairs, those that contain individual_type pointers = NULL)
    // then (size_t)-1 should be returned to signal the error.
    virtual size_t add_pair(mem_pointer i1, mem_pointer i2,
                            pair_type pt = pair_generator::NULL_TYPE) = 0;

    // Add a pair by name of pedigree and individuals
    size_t add_pair(const std::string &ped, const std::string &i1,
                    const std::string &i2, error_t &e);

    // Find and index to a pair.  Returns the index (0..pair_count()-1) if
    // the pair exists; >= pair_count() otherwise.
    // Precondition:  built() == true
    virtual size_t pair_index(const mem_pointer i1, const mem_pointer i2) const  = 0;

    size_t pair_index(const std::string &ped, const std::string &i1,
                      const std::string &i2, error_t &e) const;

    // Store/retrieve pair specific IBD information
    //
    // All IBD accessor functions return a boolean value 'true' if the
    // operation suceeded; 'false' otherwise.  'false' may not be returned if
    // the IBD information is set to quiet_NaN, since that is a valid missing
    // data condition.
    //
    // IBD accessors must be able to accept the following:
    //   1) a pair of individuals,
    //      a marker name or marker id,
    //      the P(pair sharing 0 alleles) and P(pair sharing 2 alleles)
    //      and p(pair sharing mother_1 - father_1 alleles) - addition of parent_of_origin
    //        at the given marker
    //
    //   2) a pair of individuals,
    //      a vector of P(pair sharing 0 alleles) and
    //      a vector of P(pair sharing 2 alleles)
    //        for all markers specified in this object
    //         (see marker_count(), marker_name(), marker_index())

    // IBD storage functions

    // set_ibd returns 'false' when either f0 or f2 values cannot be set or one
    // of the following is true:
    //    1) The size of vector<double>f0, vector<double>f2, and marker_count()
    //       are different.
    //    2) pair index i is out of range for the implementation
    //
    // If any particular (f0,f2) pair is not finite or violate the bounds for
    // IBD (f0,f2 > 0 & f0+f2 <= 1) both values are set to quiet_NaN and no
    // error is returned.  A warning issued via an error stream may be added
    // eventually.
    
    virtual bool set_ibd(size_t i, size_t m, double f0, double f2)    = 0;
    virtual bool set_ibd(size_t i, const std::vector<double> &f0,
                                   const std::vector<double> &f2)     = 0;

    virtual bool set_ibd(size_t i, size_t m, double f0, double f1, double f2)    = 0;
    virtual bool set_ibd(size_t i, const std::vector<double> &f0,
                                   const std::vector<double> &f1,
                                   const std::vector<double> &f2)     = 0;

    // IBD storage helper functions

    bool set_ibd(const mem_pointer i1, const mem_pointer i2,
                 const std::string &marker, double f0, double f2);

    bool set_ibd(size_t i, const std::string &marker, double f0, double f2);

    bool set_ibd(const mem_pointer i1, const mem_pointer i2,
                 size_t m, double f0, double f2);

    bool set_ibd(const mem_pointer i1, const mem_pointer i2,
                 const std::vector<double> &f0, const std::vector<double> &f2);

    // Additional IBD storage helper functions for parent-of-origin

    bool set_ibd(const mem_pointer i1, const mem_pointer i2,
                 const std::string &marker, double f0, double f1, double f2);

    bool set_ibd(size_t i, const std::string &marker, double f0, double f1, double f2);

    bool set_ibd(const mem_pointer i1, const mem_pointer i2,
                 size_t m, double f0, double f1, double f2);

    bool set_ibd(const mem_pointer i1, const mem_pointer i2,
                 const std::vector<double> &f0, const std::vector<double> &f1,
                 const std::vector<double> &f2);

    // IBD retrieval functions

    // get_ibd returns 'false' when either f0 or f2 values cannot be
    // retrieved or one of the following is true:
    //    1) The size of vector<double>f0, vector<double>f2, and marker_count()
    //       are different.
    //    2) the pair index i is out of range for the implementation

    virtual bool get_ibd(size_t i, size_t m, double &f0, double &f2) const = 0;
    virtual bool get_ibd(size_t i, std::vector<double> &f0,
                                   std::vector<double> &f2) const          = 0;

    virtual bool get_ibd(size_t i, size_t m, double &f0, double &f1, double &f2) const = 0;
    virtual bool get_ibd(size_t i, std::vector<double> &f0,
                                   std::vector<double> &f1,
                                   std::vector<double> &f2) const          = 0;

    // Retrieval helper functions:

    bool get_ibd(const mem_pointer i1, const mem_pointer i2,
                 const std::string &marker, double &f0, double &f2) const;

    bool get_ibd(const mem_pointer i1, const mem_pointer i2,
                 size_t m, double &f0, double &f2) const;

    bool get_ibd(size_t i, const std::string &marker, double &f0, double &f2) const;

    bool get_ibd(const mem_pointer i1, const mem_pointer i2,
                 std::vector<double> &f0, std::vector<double> &f2) const;

    // Additional Retrieval helper functions for parent-of-origin

    bool get_ibd(const mem_pointer i1, const mem_pointer i2,
                 const std::string &marker, double &f0, double &f1, double &f2) const;

    bool get_ibd(const mem_pointer i1, const mem_pointer i2,
                 size_t m, double &f0, double &f1, double &f2) const;

    bool get_ibd(size_t i, const std::string &marker, double &f0, double &f1, double &f2) const;

    bool get_ibd(const mem_pointer i1, const mem_pointer i2,
                 std::vector<double> &f0, std::vector<double> &f1, std::vector<double> &f2) const;

    // The following provide conversions and the ability to get
    // fn where n are the alleles shared ibd
    bool get_ibd(size_t i, size_t m, size_t n, double &fn) const;

    bool get_ibd(const mem_pointer i1, const mem_pointer i2,
                 const std::string &marker, size_t n, double &fn) const;

    bool get_ibd(const mem_pointer i1, const mem_pointer i2,
                 size_t m, size_t n, double &fn) const;

    bool get_ibd(size_t i, const std::string &marker, size_t n, double &fn) const;
    bool get_ibd(size_t i, size_t n, std::vector<double> &fn) const;

    bool get_ibd(const mem_pointer i1, const mem_pointer i2,
                 size_t n, std::vector<double> &fn) const;

    // IBD STATE storage functions
    // IBD STATE retrieval functions

    virtual const sped_pointer get_subped(size_t i) const      = 0;
    virtual sped_pointer       get_subped(size_t i)            = 0;

    virtual bool set_ibd_state(size_t m, const a_marker_ibd_state& i_state) = 0;
    virtual bool get_ibd_state(size_t m,       a_marker_ibd_state& i_state) const = 0;

    virtual bool set_ibd_state(const sped_pointer sp, const ibd_state_info& i_info) = 0;
    virtual bool get_ibd_state(const sped_pointer sp, ibd_state_info& i_info) const = 0;

    // non-virtual functions

    void                   set_ibd_option(const ibd_option_type& opt);
    const ibd_option_type& get_ibd_option() const;


  protected:

    ibd_option_type             my_ibd_option;

    vector<ibd_marker_info>     my_markers;
};

#include "ibd/ibd.ipp"

} // end of namespace SAGE

#endif

