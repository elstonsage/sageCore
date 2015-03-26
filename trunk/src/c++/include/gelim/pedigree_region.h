#ifndef PEDIGREE_REGION_H
#define PEDIGREE_REGION_H

#include <vector>
#include "error/errorstream.h"
#include "mped/mp.h"
#include "mlocus/imodel.h"
#include "rped/genome_description.h"
#include "gelim/ped_imodel_gen.h"

namespace SAGE
{

/// \brief Constructs and stores subpedigree-specific inheritance models.
///
/// The pedigree region is a storage and construction mechanism for 
/// subpedigree-specific inheritance_models.  Given subpedigrees and
/// genome regions, it constructs and then provides access to, inheritance_models
/// which are specific to the region and the subpedigree.  
  
class pedigree_region
{
  private:
  
    typedef std::vector<MLOCUS::inheritance_model> inheritance_model_vector;

  public:

    typedef FPED::Subpedigree                        subpedigree;
    typedef RPED::genome_description::region_type    region;
    typedef inheritance_model_vector::size_type      size_type;

    /// Constructor:
    ///
    /// This constructor does no build, since it has no subpedigree or region.
    ///
    /// \param err       The errorstream to which errors should be reported, if any
    /// \param quiet     Should errors be reported?
    pedigree_region(cerrorstream& err = sage_cerr, bool quiet = false);

    /// Constructor:
    ///
    /// Given a subpedigree and a region, create supedigree specific inheritance
    /// models for that pedigree in that region.
    ///
    /// \param sped      The subpedigree to construct this pedigree_region for
    /// \param r         The region to use
    /// \param err       The errorstream to which errors should be reported, if any
    /// \param quiet     Should errors be reported?
    /// \param eliminate Should genotype elimination be done after creating the
    ///                  inheritance_models?
    pedigree_region(const subpedigree& sped,
                    const region&      r,
                    cerrorstream&      err       = sage_cerr,
                    bool               quiet     = false,
                    bool               eliminate = true);

    /// Constructor:
    ///
    /// Based upon an existing pedigree region, \c pr, create a new one with
    /// \c sped as the subpedigree.  Uses \c pr_ids as the supedigree indices
    /// of the members of \c sped within the subpedigree contained in \c pr.  The
    /// phenotypes of those individuals are copied as the phenotypes of \c sped's
    /// members.  This function allows easy transfer of phenotypes and penetrances 
    /// from one subpedigree to another (say a more restricted subpedigree).
    ///
    /// \param sped      The subpedigree to construct this pedigree_region for
    /// \param pr        The previous pedigree_region
    /// \param pr_ids    The indicies of the members of \c sped within the subpedigree
    ///                  contained in \c pr
    /// \param err       The errorstream to which errors should be reported, if any
    /// \param quiet     Should errors be reported?
    /// \param eliminate Should genotype elimination be done after creating the
    ///                  inheritance_models?
    pedigree_region(const subpedigree&     sped,
                    const pedigree_region& pr,
                    const vector<uint>&    pr_ids,
                    cerrorstream&          err       = sage_cerr,
                    bool                   quiet     = false,
                    bool                   eliminate = true);

    /// Copy Constructor
    ///
    pedigree_region(const pedigree_region&);
    
    /// Destructor
    ///
    ~pedigree_region();
    
    /// Copy operator
    ///
    pedigree_region& operator=(const pedigree_region&);

    /// Error Reporting
    //@{
      
    /// Returns the error stream to which the pedigree_region sends any errors
    /// that might be generated when doing builds (when !quiet)
    cerrorstream& get_errorstream();

    /// Sets the error stream to which the pedigree_region sends any errors
    /// that might be generated when doing builds (when !quiet)
    ///
    /// \param err The error stream to which the pedigree_region should send errors
    void set_errorstream(cerrorstream& err);
    
    /// Returns true if the quiet flag is set, false otherwise.
    ///
    bool is_quiet() const;
    
    /// Sets the quiet flag.
    ///
    /// \param q What the quiet flag should be set to.
    void set_quiet(bool q);
    //@}
    
    /// Build functions
    //@{

    /// Given a subpedigree and a region, create supedigree specific inheritance
    /// models for that pedigree in that region.
    ///
    /// \param sped      The subpedigree to construct this pedigree_region for
    /// \param r         The region to use
    /// \param eliminate Should genotype elimination be done after creating the
    ///                  inheritance_models?
    bool build(const subpedigree& sped,
               const region&      r,
               bool               eliminate = true);

    /// Based upon an existing pedigree region, \c pr, rebuild this one with
    /// \c sped as the subpedigree.  Uses \c pr_ids as the supedigree indices
    /// of the members of \c sped within the subpedigree contained in \c pr.  The
    /// phenotypes of those individuals are copied as the phenotypes of \c sped's
    /// members.  This function allows easy transfer of phenotypes and penetrances 
    /// from one subpedigree to another (say a more restricted subpedigree).
    ///
    /// \param sped      The subpedigree to construct this pedigree_region for
    /// \param pr        The previous pedigree_region
    /// \param pr_ids    The indicies of the members of \c sped within the subpedigree
    ///                  contained in \c pr
    /// \param eliminate Should genotype elimination be done after creating the
    ///                  inheritance_models?
    bool build(const subpedigree&     sped,
               const pedigree_region& pr,
               const vector<uint>&    pr_ids,
               bool                   eliminate = true);

    /// Returns \c true if last build was successful, \c false otherwise 
    ///
    bool is_built() const;
    
    //@}

    /// Current subpedigree and region information
    //@{
    
    /// Returns the current subpedigree.  Calls an internal error if there
    /// isn't one available.
    const subpedigree& get_subpedigree() const;

    /// Returns the current region.  Calls an internal error if there 
    /// isn't one available
    const region& get_region() const;
    //@}

    /// Inheritance Model Retrieval
    //@{
          MLOCUS::inheritance_model& operator[] (size_type t);
    const MLOCUS::inheritance_model& operator[] (size_type t) const;
    //@}

    /// Indicators for informativity and consistency.  Note that if the model
    /// is inconsistent, it is always uninformative
    //@{
    bool model_consistent  (size_type t) const;
    bool model_informative (size_type t) const;
    //@}

    // Misc.
    
    /// Returns the number of inheritance models in the 
    size_type inheritance_model_count() const;

  private:

    void reset();

    void create_model_storage();
    bool check_build_invariants();

    void set_marker_status(size_t i,const pedigree_imodel_generator& generator);

    // Data members

    /// Error Reporting Variables
    //@{
    cerrorstream my_errors;

    bool my_is_quiet;
    //@}

    bool my_is_built;

    const subpedigree*             my_subpedigree;

    region                         my_region;
    inheritance_model_vector       my_markers;

    vector<bool>                   my_model_consistencies;
    vector<bool>                   my_model_informatives;
};

#include "gelim/pedigree_region.ipp"

}

#endif
