#ifndef MEMBER_CALCULATOR_H
#define MEMBER_CALCULATOR_H
//=========================================================================
//  File:	member_calculator.h
//
//  Author:	Stephen Gross
//
//  History:	0.1  sag  Initial implementation	Jul 11 01
//              1.0  gcw  Split into cont. & bin.       Jul    02
//
//  Notes:	Calculates analysis traits, expected means, and expected
//		variates.
//
//  Copyright (c) 2002 R. C. Elston
//=========================================================================

/// \file
///
/// The member calculators take care of calculating various values for each
/// member of each pedigree.  These values vary depending on the type of
/// model being fit (continuous, binary or onset) and sometimes on the
/// analysis type (FPMM, class A vs. class D, etc)


#include "segreg/model.h"
#include "segreg/cov_sub_model.h"
#include "segreg/type_sub_model.h"
#include "segreg/sub_model_base.h"
#include "segreg/segreg_datatypes.h"
#include "fped/fped.h"
#include <math.h>
#include <cmath>
#include <vector>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

namespace SAGE {
namespace SEGREG {

/// \brief Base class common to all the member calculators
///
/// This class provides core functionality of all the member calculators, regardless
/// of model.
class member_calculator_base
{
  public:

    /// Member classes determine the type of each individual, which in turn
    /// determines how the member penetrance and other functions are to be
    /// calculated by calculator functions both internal to and external to
    /// the member calculators.  The list of member classes that are used
    /// for any given class varies, so overlap occurs as indicated below.
    enum member_class
    {
      missing           = 0,        ///< Member is missing            (all)
      actual            = 1,        ///< Member to use actual values  (cont)
      available         = 1,        ///< Member is available          (binary)

      age_of_onset      = 1,        ///< Member uses age of onset     (onset)
      age_at_exam_aff   = 2,        ///< Member is affected at exam   (onset)
      age_at_exam_unaff = 3,        ///< Member is unaffected at exam (onset)
      age_at_onset      = 4,        ///< Member uses age *at* onset   (onset) (ascertained)
      
      gte_thresh        = 2,        ///< Member uses threshold equ.s  (cont)  (ascertained)
      lte_thresh        = 3         ///< Member uses threshold equ.s  (cont)  (ascertained)
    };

    typedef FPED::Member member_type;

    /// Destructor for base class
    virtual ~member_calculator_base() { }

    /// Returns if member is valid (ie, has data, is in ascertained subset, etc)
    bool         is_member_valid    (const member_type& member)  const;

    /// Returns the member_class for a specific individual
    member_class get_member_class   (const member_type& member)  const;

  protected:

    typedef FPED::PedigreeConstPointer        pedigree_const_pointer;
    typedef FPED::Multipedigree::pedinfo_type pedigree_info;

    size_t get_abs_mem_ref(const member_type& member) const;

    /// Returns if member is valid (ie, has data, is in ascertained subset, etc)
    bool         is_member_valid    (size_t mem_ref) const;

    /// Returns the member_class for a specific individual
    member_class get_member_class   (size_t mem_ref) const;

    // Import helpers

    /// Looks up the trait for the individual and puts it into the double
    inline bool import_trait(size_t trait_index, pedigree_const_pointer ped,
                             size_t member_index, size_t mem_ref,
                             double& val);

    /// Looks up the trait for the individual and puts it into the bool vector
    //lint -e{1015,1013,601} <- spurious errors caused by vector<bool>
    inline bool import_trait(size_t trait_index, pedigree_const_pointer ped,
                             size_t member_index, size_t mem_ref,
                             vector<bool>::reference val);

    void import_covariate_data(pedigree_const_pointer   ped, size_t member_index,
                               size_t abs_mem_ref, vector<vector<double> >& cov_data,
                               const CovariateSubmodel& cov_model);

    inline void center_trait(vector<double>& v, double mean) const;

    /// Private Variables common to all MCs
    //@{

    const FPED::Multipedigree & my_ped_data;
    const model               & mod;

    // pedigree_index_map stores quick reference pointer info for looking up
    // covariate trait data based on a pedigree-relative member number

    std::map<pedigree_const_pointer,size_t> pedigree_index_map;

    // The following two variables are for optimizing the lookup into the
    // pedigree_index_map

    mutable pedigree_const_pointer last_pedigree_accessed;
    mutable size_t                 last_pedigree_location;

    vector<member_class>             member_classes;

    /// Is this MCC using ascertainment?
    bool use_ascertainment;

    //@}

    /// Constructor is protected to prevent instantiation of class
    
    member_calculator_base (const FPED::Multipedigree & ped_data,
                            const model               & modp,
                            bool                        use_asc = false);

    //lint --e{1712} <- no default constructor
};

class continuous_member_calculator : public member_calculator_base
{
    friend class segreg_utilities;

  public:

    continuous_member_calculator(const FPED::Multipedigree & multipedigree_data_param,
                                 const model               & modp,
                                 bool                        use_asc = false);

    // Update the values

    int update();

    // Diagnostic functions:

    void dump() const;

    // Get all the values calculated

    double get_composite_trait      (const member_type& member)  const;
    double get_expected_mean        (const member_type& member, genotype_index genotype) const;
    double get_expected_mean        (const member_type& member, genotype_index genotype, size_t polygenotype) const;
    double get_expected_variance    (const member_type& member, genotype_index genotype) const;
    double get_expected_sd          (const member_type& member, genotype_index genotype) const;

    double get_standardization             (const member_type& member,genotype_index genotype) const;
    double get_ascertained_standardization (const member_type& member,genotype_index genotype) const;
    double get_estimated_standardization   (const member_type& member,genotype_index genotype_mother,genotype_index genotype_father) const;
    double get_polygenic_standardization   (const member_type& member,genotype_index genotype,size_t polygenotype) const;

  protected:

    double get_composite_trait         (size_t mem_ref                                            ) const;
    double get_expected_mean           (size_t mem_ref,genotype_index genotype                    ) const;
    double get_expected_mean           (size_t mem_ref,genotype_index genotype,size_t polygenotype) const;      
    double get_expected_variance       (size_t mem_ref,genotype_index genotype                    ) const;
    double get_expected_sd             (size_t mem_ref,genotype_index genotype                    ) const;

    double get_standardization             (size_t mem_ref,genotype_index genotype) const;  
    double get_ascertained_standardization (size_t mem_ref,genotype_index genotype) const;  
    double get_estimated_standardization   (size_t mem_ref,genotype_index genotype_mother,genotype_index genotype_father) const;
    double get_polygenic_standardization   (size_t mem_ref,genotype_index genotype,size_t polygenotype) const;  

    // Initial input routines

    void allocate_memory ();
    void import_data     ();
    void center_data     ();

    // Calculation functions - return 0 if ok, error code otherwise

    int calculate_composite_traits                    ();
    int calculate_expected_means                      ();
    int calculate_expected_polygenic_means            ();
    int calculate_expected_variances                  ();
    int calculate_expected_sds                        ();
    int calculate_expected_polygenic_sds              ();
    int calculate_polygenic_means                     ();
    int calculate_standardizations                    ();
    int calculate_ascertained_standardizations        ();
    int calculate_polygenic_standardizations          ();
    int calculate_estimated_standardizations          ();

    // Calculation routines

    void calc_exp_var    (size_t mem_ref, genotype_index g);
    void calc_exp_sd     (size_t mem_ref, genotype_index g);
    void calc_stan       (size_t mem_ref, genotype_index g);
    void calc_asc_stan   (size_t mem_ref, genotype_index g, double Thigh, double Tlow);
    void calc_pg_stan    (size_t mem_ref, genotype_index g, size_t pg);
    void calc_est_stan   (size_t mem_ref, genotype_index g_m, genotype_index g_f);

    // Fixed quantities:

    vector<double>                   primary_traits;

    vector<vector<double> >          comp_trait_data;
    vector<vector<double> >          mean_cov_data;
    vector<vector<double> >          var_cov_data;

    // Cached calculations:

    vector<double>                   composite_traits;
    vector<vector<double> >          expected_means;
    vector<vector<vector<double> > > expected_polygenic_means;
    vector<vector<double> >          expected_variances;
    vector<vector<double> >          expected_sds;
    vector<vector<double> >          standardizations;
    vector<vector<double> >          ascertained_standardizations;
    vector<vector<vector<double> > > estimated_standardizations;
    vector<vector<vector<double> > > polygenic_standardizations;

    //lint --e{1712}
};

class binary_member_calculator : public member_calculator_base
{
    friend class segreg_utilities;

  public:
  
    binary_member_calculator(const FPED::Multipedigree & ped_data,
                             const model               & modp,
                             bool                        use_asc = false);

    // Update the values

    int update();

    // Diagnostic functions:

    void dump() const;

    /// Accessors to the calculated values
    //@{

    /// Returns affection status
    bool   get_aff_status              (const member_type& member)       const;

    /// Returns expected susceptibility (genotype susc + covariates)
    double get_expected_susc           (const member_type& member,
                                        genotype_index     genotype)     const;

    /// Returns penetrance based on affection status and expected
    /// susceptibility (genotypic)
    double get_penetrance              (const member_type& member,
                                        genotype_index     genotype)     const;

    /// Returns penetrance based on affection status and expected
    /// susceptibility (includes polygenotype for FPMM)
    double get_penetrance              (const member_type& member,
                                        genotype_index     genotype,
                                        size_t             polygenotype) const;
    //@}

    /// calculates unconnected individaul likelihood
    double         unconnected_likelihood(const member_type&, genotype_index) const;

  protected:

    bool   int_get_aff_status    (size_t mem_ref                                            ) const;
    double int_get_expected_susc (size_t mem_ref,genotype_index genotype                    ) const;
    double int_get_penetrance    (size_t mem_ref,genotype_index genotype                    ) const;
    double int_get_penetrance    (size_t mem_ref,genotype_index genotype,size_t polygenotype) const;      

    // Initial input routines

    void allocate_memory ();
    void import_data     ();
    void center_data     ();

    // Calculation functions - return 0 if ok, error code otherwise

    int calculate_susceptibilities();
    int calculate_penetrances ();
    int calculate_polygenic_penetrances ();

    // Calculation routines

    double calculate_penetrance(double mean, double affection) const;

    double calc_exp_susc (size_t mem_ref, genotype_index genotype);

    // Fixed quantities:

    vector<bool>            affections;
    vector<vector<double> > susc_cov_data;

    // Cached calculations:

    vector<vector<double> >          expected_susceptibilities;
    vector<vector<double> >          penetrances;
    vector<vector<vector<double> > > polygenic_penetrances;

    //lint --e{1712}
};

class onset_member_calculator : public member_calculator_base
{
    friend class segreg_utilities;

  public:

    onset_member_calculator(const FPED::Multipedigree & ped_data,
                            const model               & modp,
                            bool                        use_asc = false);

    // Update the values

    int update();

    // Diagnostic functions:

    void dump() const;

    // Get all the values calculated

    /// Returns affection status
    bool   get_aff_status           (const member_type& member)  const;

    /// Returns transformed age of onset.  Returns QNAN if not available
    double get_age_onset            (const member_type& member)  const;

    /// Returns transformed age of exam.  Returns QNAN if not available
    double get_age_exam             (const member_type& member)  const;

    /// Returns transformed expected age of onset.  Note that it only uses
    /// the polygenotype if Multi-Dependent is A.
    double get_expected_age_onset   (const member_type& member, genotype_index genotype, size_t polygenotype) const;

    /// Returns expected susceptibility. Note that it only uses the
    /// polygenotype if Multi-Dependent is S.
    double get_expected_susc        (const member_type& member, genotype_index genotype, size_t polygenotype) const;

    /// Returns \f$\alpha_i\f$, where
    ///
    /// \f[ \alpha_i = \frac {\pi} 
    ///                      { \sqrt { 3 } \sigma_{a_i} }
    /// \f]
    double get_alpha_i              (const member_type& member, genotype_index genotype) const;

  protected:

    bool         get_aff_status         (size_t mem_ref) const;
    double       get_age_onset          (size_t mem_ref) const;
    double       get_age_exam           (size_t mem_ref) const;

    double       get_expected_age_onset (size_t mem_ref,genotype_index genotype,size_t polygenotype) const;      
    double       get_alpha_i            (size_t mem_ref,genotype_index genotype                    ) const;
    double       get_expected_susc      (size_t mem_ref,genotype_index genotype,size_t polygenotype) const;      

    // Initial input routines

    void allocate_memory ();
    void import_data     ();
    void center_data     ();

    // Calculation functions - return 0 if ok, error code otherwise

    int calculate_transf_age_onset ();
    int calculate_transf_age_exam  ();
    int calculate_expected_means   ();
    int calculate_expected_suscs   ();
    int calculate_expected_alphas  ();

    // Calculation routines

    double calc_exp_mean  (size_t mem_ref, genotype_index);
    double calc_exp_susc  (size_t mem_ref, genotype_index);
    double calc_exp_alpha (size_t mem_ref, genotype_index);

    /// Did we fail to calculate the geometric mean of the age of onset?
    bool valid_geom_mean;

    // Fixed quantities:

    vector<bool>                     affections;
    vector<double>                   age_onsets;
    vector<double>                   age_exams;

    vector<vector<double> >          mean_cov_data;
    vector<vector<double> >          var_cov_data;
    vector<vector<double> >          susc_cov_data;

    // Calculated transforms

    vector<double>                   transf_age_onsets;
    vector<double>                   transf_age_exams;

    // Cached calculations:

    vector<vector<double> >          temp_calc_space;

    vector<vector<vector<double> > > expected_means;
    vector<vector<double> >          expected_alphas;
    vector<vector<vector<double> > > expected_susceptibilities;

    //lint --e{1712}
};

}} // End namespace

#include "segreg/member_calculator.ipp"

#endif
