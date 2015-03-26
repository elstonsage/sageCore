#ifndef TDTEX_CONFIG_H
#define TDTEX_CONFIG_H

#include "LSF/LSF.h"
#include "app/aparser.h"
#include "error/errorstream.h"
#include "LSF/LSFsymbol.h"

namespace SAGE  {
namespace TDTEX {

///
/// Describes a tdtex analysis to be performed.
class Configuration
{
public:

  /// @name Enums
  //@{
  
    ///
    /// Describes which type of transmissions are to be scored (see 'sample' option in manual).
    enum SamplingMethod { ALLELES = 0, GENOTYPES };
  
  //@}

  /// @name Constructors
  //@{
  
    Configuration();
    
    Configuration(const Configuration & other);

    Configuration & operator= (const Configuration & other);
  
  //@}

  /// @name Accessors
  //@{

    ///
    /// Returns the marker requested for analysis.
    size_t get_marker() const { return my_marker; }

    ///
    /// Returns the affection status trait.
    size_t get_trait() const { return my_trait; }

    ///
    /// Returns the parental trait (used as an indicator for selecting subsets of pairs for analysis).
    size_t get_parent_trait() const { return my_parent_trait; }

    ///
    /// Returns the type of transmission is to be scored.
    SamplingMethod get_method() const { return my_method; }

    ///
    /// Returns the max number of children, or -1 if unlimited.
    size_t get_max_children() const { return my_max_children; }

    ///
    /// Returns the max number of sibpairs (per family), or -1 if unlimited.
    size_t get_max_sib_pairs() const { return my_max_sib_pairs; }

    ///
    /// Returns whether or not a set of three tests is to be performed (see user manual).
    bool get_sex_differential() const { return my_sex_differential; }
    
    ///
    /// Returns whether or not to skip the mc test.
    bool get_skip_mc_test() const { return my_skip_mc_test; }
  
    ///
    /// Returns whether or not to skip the mcmh test.
    bool get_skip_mcmh_test() const { return my_skip_mcmh_test;  }
    
    ///
    /// Returns whether or not to skip the permutation test.
    bool get_skip_permutation_test() const { return my_skip_permutation_test; }

    ///
    /// Returns the output filename for analysis.
    const std::string & get_ofilename() const { return my_ofilename; }

  //@}

  /// @name Mutators
  //@{

    ///
    /// Sets the marker to be analyzed.
    void set_marker(size_t m) { my_marker = m; }

    ///
    /// Sets the affection status trait.
    void set_trait(size_t t) { my_trait = t; }

    ///
    /// Sets the parental trait.
    void set_parent_trait(size_t t) { my_parent_trait = t; }

    ///
    /// Sets the sampling method.
    void set_method(SamplingMethod s) { my_method = s; }

    ///
    /// Sets the max number of children (per family) to use; -1 indicates unlimited.
    void set_max_children(size_t m) { my_max_children = m; }

    ///
    /// Sets the max number of sib pairs (per family) to use; -1 indicates unlimited.
    void set_max_sib_pairs(size_t m) { my_max_sib_pairs = m; }

    ///
    /// Sets whether or not to skip the mc test.
    void set_skip_mc_test(bool skip) { my_skip_mc_test = skip; }
  
    ///
    /// Sets whether or not to skip the mcmh test.
    void set_skip_mcmh_test(bool skip) { my_skip_mcmh_test = skip; }

    ///
    /// Sets whether or not to skip the permutation test.
    void set_skip_permutation_test(bool skip) { my_skip_permutation_test = skip; }

    ///
    /// Sets whether or not to perform a set of three tests (see user manual).
    void set_sex_differential(bool skip) { my_sex_differential = skip; }

    ///
    /// Sets the output filename.
    void set_ofilename(const std::string & s) { my_ofilename = s; }

  //@}

  // Debugging
  void dump() const;

  private:

  // Data members

    size_t         my_marker;
    size_t         my_trait;
    size_t         my_parent_trait;
    SamplingMethod my_method;
    size_t         my_max_children;
    size_t         my_max_sib_pairs;
    bool           my_skip_permutation_test;
    bool           my_skip_mc_test;
    bool           my_skip_mcmh_test;
    bool           my_sex_differential;
    std::string    my_ofilename;
};

} // End namespace TDTEX
} // End namespace SAGE


#endif
