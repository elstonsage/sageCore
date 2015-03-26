#ifndef PEDIGREE_ANALYSIS_SAMPLE_H
#define PEDIGREE_ANALYSIS_SAMPLE_H

#include "rped/rped.h"
#include "fped/fped.h"
#include "output/Output.h"
#include "mlod/analysis_parameters.h"

namespace SAGE {
namespace MLOD {

/// \brief Creates and stores a set of pedigrees according to their usability within an analysis
/// 
/// Based upon an analysis, provides a FPED::Multipedigree of all the data
/// that can be analyzed using that analysis.  Provides a classification 
/// of each member of the source RPED::MultiPedigree according to their status
/// (included or excluded) and the rationale (informative, informative due 
/// to structure, uninformative, informative singleton) for their 
/// inclusion/exclusion.  Provides several basic stats, such as number of
/// pedigrees, largest pedigree's individual count, etc. 
class PedigreeAnalysisSample
{
    typedef std::vector<const FPED::Subpedigree*> SpedVectorType;
  
  public:

  /// \name Object Maintenance
  //@{
    /// Constructor
    ///
    /// \param source_rped     The original source RPED::RefMultipedigree
    /// \param analysis_params The analysis parameters to be used to determine
    ///                        inclusion or exclusion of particular members
    PedigreeAnalysisSample(const RPED::MultiPedigree& source_rped,
                           const AnalysisParameters&  analysis_params);
                           
    /// Destructor
    ///
    ~PedigreeAnalysisSample();
  //@}
    
  /// \name Statistics
  //@{
    
    /// Returns the number of members in the original data set
    ///
    size_t get_original_member_count      () const;
    
    /// Returns the number of subpedigrees in the original data set
    ///
    size_t get_original_subpedigree_count () const;
    
    /// Returns the number of pedigrees in the original data set
    ///
    size_t get_original_pedigree_count    () const;
    
    /// Returns the number of members in the filtered data set
    ///
    size_t get_member_count      () const;
    
    /// Returns the number of subpedigrees in the filtered data set
    ///
    size_t get_subpedigree_count () const;
    
    /// Returns the number of pedigrees in the filtered data set
    ///
    size_t get_pedigree_count    () const;
    
    /// Returns the number of bits of the largest subpedigree.  This
    /// is the minimum number of bits that MLOD will have to support to
    /// make this analysis work.
    size_t get_largest_bit_count () const;
  //@}
  
  /// \name Member Inclusion Information
  //@{
  
    /// The MemberInclusionType describes the reason each member is included or
    /// excluded from the filtered data.
    enum MemberInclusionType
    {
      MIT_UNDEFINED,                 ///< Not known.  Used for initialization purposes. Should never appear in data
      MIT_SINGLETON,                 ///< Member is a singleton (irrelevant for MLOD)
      MIT_DIRECTLY_INFORMATIVE,      ///< Member is informative due to its data.
      MIT_STRUCTURALLY_INFORMATIVE,  ///< Member data is uninformative, but it is structurally informative
      MIT_UNINFORMATIVE,             ///< Member is uninformative both structurally and data based.
      MIT_SINGLETON_AFTER,           ///< Member is a singleton after uninformative members have been removed.
      MIT_SPED_TOO_LARGE,            ///< Member's subpedigree is too large to be processed with the analysis
      MIT_SPED_TOO_SMALL,            ///< Member's subpedigree is too small to be processed with the analysis
      MIT_ERROR                      ///< Something has gone seriously wrong with this member.  Shouldn't ever happen.
    };
  
    /// This boolean function simply returns if a member from the original multipedigree
    /// was included in the filtered data.
    ///
    /// \param mem The member whose status we want.
    bool is_member_included(const RPED::Member& mem);
    
    /// Returns the MemberInclusion type for the specified member.
    ///
    /// \param mem The member whose type we want.
    MemberInclusionType get_member_type(const RPED::Member& mem);
    
    /// Returns the MemberInclusion type for the specified member.
    ///
    /// \param mem The member whose type we want.
    MemberInclusionType get_member_type(const FPED::Member& mem);
  //@}
  
  /// \name Data Sample Access
  //@{
    /// Returns the multipedigree on which the filtering was done.
    ///
    const RPED::MultiPedigree& get_original_multipedigree() const;

    /// Returns the filtered multipedigree
    ///
    const FPED::Multipedigree& get_filtered_multipedigree() const;
  //@}
  
  /// \name Subpedigree Access
  ///
  /// Once a sample is generated, the MLOD analysis is performed on each subpedigree
  /// in the sample.  Rather than requiring the analysis to iterate over the 
  /// pedigrees, then the subpedigrees, the PedigreeAnalysisSample provides a 
  /// set of iterators for iterating over the subpedigrees.
  ///@{

    /// Typedef for an iterator over all subpedigrees.  This is not provided
    /// by the multipedigree.
    ///
    /// \internal
    ///
    /// This should likely be replaced with an actual iterator class, possibly
    /// boost's dereferencing iterator so that subpedigrees are returned, rather
    /// than subpedigree pointers, but boost's custom iterators do not compile under
    /// KCC.
    typedef SpedVectorType::const_iterator SpedIterator;
    
    /// Returns the beginning of our set of subpedigrees
    ///
    SpedIterator get_subpedigree_begin() const;

    /// Returns the end of our set of subpedigrees.
    ///
    SpedIterator get_subpedigree_end()   const;
    
    /// Returns the subpedigree at the index specified by spid.
    const FPED::Subpedigree& get_subpedigree(size_t spid) const;
  ///@}
  
  /// \name Display options
  ///
  /// Common File output options of the results of the processing of the data
  ///@{
    
    /// Print the summary statistics (# members, pedigrees and speds before and after filtering)
    ///
    /// \param o The ostream where the output should go
    void print_summary_statistics_table(ostream& o);
    
    /// Print the table of individuals.  Default output only includes those individuals
    /// who were removed due to problems.  Detailed output prints for all individuals
    ///
    /// \param o        The ostream where the output should go
    /// \param detailed Is detailed output desired?
    void print_member_table(ostream& o, bool detailed = false);
  ///@}
  
  private:
  
  /// \name Disabled object management options declared as private
  //@{
    //lint -e{1704} We don't want these built incorrectly
    /// \internal
    /// Disabled Constructor
    PedigreeAnalysisSample();
    //lint -e{1704} We don't want these built incorrectly
    /// \internal
    /// Disabled copy constructor
    PedigreeAnalysisSample(const PedigreeAnalysisSample&);
    
    /// \internal
    /// Disabled Copy Operator
    PedigreeAnalysisSample& operator=(const PedigreeAnalysisSample&);
  //@}

    void initialize_member_type_vector();
    void classify_members_informativity(const AnalysisParameters&  analysis_params);
    void create_fped                   (const AnalysisParameters&  analysis_params);

    void count_speds();
    
    size_t calculate_sped_bit_count(const FPED::Subpedigree& sped) const;

    /// Reference to the complete data set
    ///
    const RPED::MultiPedigree&  my_original_rped;
    
    /// Multipedigree filtered on analysis
    ///
    FPED::Multipedigree         my_fped;
    
    /// Details for each individual on reasons for inclusion/exclusion
    vector<MemberInclusionType> my_member_types;
    
    /// Subpedigrees in the filtered data set
    ///
    SpedVectorType              my_spedigrees;
    
    /// Stored quantity of the largest number of bits for any pedigree in the
    /// sample
    size_t                      my_largest_bit_count;

    /// Number of subpedigrees in my_original_rped
    ///
    size_t                      my_original_sped_count;
    
    /// Number of subpedigrees in my_fped
    ///
    size_t                      my_filtered_sped_count;
};

}
}

#include "mlod/PedigreeAnalysisSample.ipp"

#endif
