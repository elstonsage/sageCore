#ifndef COMMON_ALLELE_SET_H
#define COMMON_ALLELE_SET_H

#include <vector>
#include "mlocus/genotype.h"
#include "error/internal_error.h"

namespace SAGE
{
namespace MCMC
{

/// This class determines the alleles common (shared) between a set of
/// allele pairs (genotypes).  It maintains a count of allele occurances so
/// that allele pairs can be added or removed.
///
/// There are four states the set may be in at any given time.  These are:
///
/// <table>
///   <tr>
///     <td> \b State  </td>
///     <td> \b Reason </td>
///   </tr>
///   <tr>
///     <td>All alleles valid</td><td>No allele pairs have been added. </td>
///   </tr>
///   <tr>
///     <td>Two valid alleles</td><td>At least one genotype and all
///                                   allele pairs share the same two alleles.
///                                   </td>
///   </tr>
///   <tr>
///     <td>One valid allele </td><td>At least one genotype and set of
///                                   alleles shared by allele pairs contains
///                                   onlye one allele.  </td>
///   </tr>
///   <tr>
///     <td>No valid alleles </td><td>More than one genotype and the set of
///                                   alleles shared by all allele pairs is
///                                   empty.  </td>
///   </tr>
/// </table>
///
/// Alleles are added in pairs, and can be removed as well.  It does not
/// maintain a list of the allele pairs but only an occurance count, so care
/// needs to be taken that the same pairs that are added are removed.  If this
/// rule is not followed, behavior is undefined, and may cause SAGE internal
/// errors.
///
/// While adding and removing, the following statistics are maintained:
///
/// - the set of alleles shared,
///
/// - a status variable indicating the state (see above),
///
/// - a count of the number of allele pairs that have been added.
///
/// - for each allele, a count of the number of allele pairs that have had
///   that allele.  This is a count of the pairs.  Pairs which contained two
///   copies of the same allele do not count twice.
///
/// \anchor CAS_Optimization
/// \par Optimization Issues
///
/// The removal of allele pairs sometimes requires an exhaustive search
/// of all alleles to find the alleles that are valid after the removal.
/// Since the algorithms which use this object are extremely time intensive,
/// rather than recompute this at every removal,
/// the actual recomputation of state is delayed until the status or the
/// alleles are requested.  This means that the recomputation only happens
/// once for potentially many additions and removals.
///
/// Each addition and removal returns a status variable indicating a potential
/// change to internal status, again for potential external optimizations.
/// This is used to reflect a 'dirty' state of the CommonAlleleSet.
/// A is_dirty() function is also provided.  When the alleles and/or state is
/// requested, the actual lookups are processed, and the dirty flag is cleared.
///
/// Typically, the is_dirty() test will not be used, using the potential change
/// boolean return from add_allele_pair() and remove_allele_pair() to determine
/// the need to revisit the state of the CommonAlleleSet when writing algorithms.

class CommonAlleleSet
{
  public:

    /// The possible states in which the CommonAlleleSet can be
    ///
    /// \internal
    /// 
    /// The ennumerated values are converted to and from size_t, so
    /// the 0,1,2 values should not change, If they are changed, underlying
    /// code will break.
    enum StateEnum 
    {
      INVALID   = 0,  ///< No alleles common to all pairs
      ONE_VALID = 1,  ///< one allele common to all pairs
      TWO_VALID = 2,  ///< two alleles common to all pairs
      EMPTY     = 3   ///< all alleles valid (no pairs included in set)
    };

    /// The AlleleID is simply a numerical value from 0 to n-1 for
    /// n alleles.
    typedef MLOCUS::allele AlleleID;
    
    /// \name Constructors and Copy Operators
    //@{
  
    /// Default constructor.
    ///
    CommonAlleleSet();

    /// The basic constructor.  This takes the genotype model so
    /// that it can maintain the data structure.
    explicit CommonAlleleSet(const MLOCUS::genotype_model& gmodel);
    
    /// Copy Constructor
    ///
    /// \param cas The CommonAlleleSet to copy
    CommonAlleleSet(const CommonAlleleSet& cas);

    /// Copy Operator
    ///
    /// \param cas The CommonAlleleSet to copy
    CommonAlleleSet& operator=(const CommonAlleleSet& cas);
 
    //@}

    /// \name Status Indicators
    //@{

    /// Get the maximum allele count.
    ///
    size_t get_max_allele_count() const;
    
    /// Get the current status of the object
    ///
    StateEnum get_status() const;
    
    /// Get the n-th allele of the set.  Alleles are always returned sorted, so that
    /// the smaller allele is get_allele(0) and the larger get_allele(1).
    ///
    /// \param a The allele you want.  This must be 0 or 1 or undefined behavior results.
    ///     
    /// \return The allele of the allele requested.  If the allele is invalid (asking
    ///         for the 1 allele when the state is ONE_VALID or INVALID or
    ///         the 0 allele when the state is INVALID) returns NPOS
    AlleleID get_allele_id(size_t a) const;

    /// Returns the number of allele pairs that have been added to the set.
    ///
    size_t get_allele_pair_count() const;

    /// Returns the number of allele pairs added to the set that have
    /// contained at least one copy of the allele a.
    ///
    /// \param a the allele in question.
    size_t get_allele_pair_count_for_allele(AlleleID a) const;

    /// Indicates that there has been potential changes since the last
    /// view of the status.  See \ref CAS_Optimization section for details.

    bool is_dirty() const;
    //@}
    
    /// \name Mutators
    //@{
 
   /// Adds another allele pair
    ///
    /// \param a1 The first  allele of the pair
    /// \param a2 The second allele of the pair
    ///
    /// \return \c true if this results in a \b potential change to the current state,
    ///         \c false otherwise.  See \ref CAS_Optimization "Optimization Issues"
    ///         section for details on potential changes.
    bool add_allele_pair      (AlleleID a1, AlleleID a2);

    /// Removes an allele pair.
    ///
    ///
    /// \param a1 The first  allele of the pair
    /// \param a2 The second allele of the pair
    ///
    /// \return \c true if this results in a \b potential change to the current state,
    ///         \c false otherwise.  See \ref CAS_Optimization "Optimization Issues"
    ///         section for details on potential changes.
    bool remove_allele_pair      (AlleleID a1, AlleleID a2);

    //@}
    
    /// \name Boolean Tests
    //@{

    //lint --e{1739} Lint produces a warning for the boolean tests below,
    //               but the error is spurious since constructor is explicit
    
    /// Returns \c true if and only if the alleles currently present in the
    /// set are the same.  It does not care about the number/types of
    /// genotypes that have been added are, only about the allele state
    bool operator==(const CommonAlleleSet&) const;

    /// Returns \c false if and only if the alleles currently present in the
    /// set are different.
    bool operator!=(const CommonAlleleSet&) const;

    //@} 
    
  private:

    /// Compare input parameters a1 and a2 against existing valid alleles.
    /// Keep only those alleles which still belong to the set and update state
    /// Does not require the alleles to be sorted or non-unique.
    ///
    /// \param a1 The first allele to test
    /// \param a2 The second allele to test
    void determine_and_set_valid_alleles(AlleleID a1, AlleleID a2);

    /// When allele pairs are removed, it is often necessary to search for
    /// alleles that were invalid (not in all pairs) that have been made
    /// valid (in all remaining pairs).  find_and_set_valid_alleles() finds and
    /// sets these alleles
    ///
    /// Note:  Technically, not a const function, this works on the mutable
    ///        members of the class.
    void find_and_set_valid_alleles() const;

    /// Set current state to EMPTY.
    ///
    void set_to_empty();
    
    /// Set the state to TWO_VALID using a1 and a2.  Alleles must
    /// be unique (a1 != a2), but need not be sorted.
    ///
    /// \param a1 The first allele to set
    /// \param a2 The second allele to set
    void set_to_two_valid(AlleleID a1, AlleleID a2);
    
    /// Set the state to ONE_VALID using a
    ///
    /// \param a The allele to set
    void set_to_one_valid(AlleleID a);
    
    /// Set the state to INVALID
    ///
    void set_to_invalid();

    MLOCUS::genotype_model my_model;
    
    /// The current state
    ///
    /// This data member is mutable to allow for lazy evaluation methods
    mutable StateEnum my_state;
    
    /// The alleles that are currently in the subset
    ///
    /// This data member is mutable to allow for lazy evaluation methods
    mutable AlleleID my_valid_alleles[2];

    /// Contains a value for each allele, the number of allele pairs that
    /// have had that allele as a member
    std::vector<size_t> my_allele_genotype_counts;

    /// The number of allele pairs that have been added to this set
    ///
    size_t my_allele_pair_count;

    mutable bool my_is_dirty;
};

} // End MCMC
} // End SAGE

#include "mcmc/common_allele_set.ipp"

#endif
