#ifndef REG_PENETRANCE_CONTEXT_H
#define REG_PENETRANCE_CONTEXT_H

#include <math.h>
#include "segreg/member_calculator.h"
#include "segreg/sub_model_base.h"
#include "segreg/resid_sub_model.h"
#include "segreg/segreg_datatypes.h"
#include "segreg/RegPenetranceCommonAdjustments.h"
#include "segreg/model.h"
#include "segreg/types/TypeDescription.h"
#include "numerics/functions.h"

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

namespace SAGE {
namespace SEGREG {

class PenetranceContext;

class PenetranceContextStatistics
{
  public:
    friend class PenetranceContext;
    
    PenetranceContextStatistics
                     (const residual_correlation_sub_model& rsm,
                      const continuous_member_calculator&   cmc,
                      const RegPenetranceCommonAdjustments& corr_adj);
  
    /// \name Calculated Value Accessors
    ///
    /// These routines return values as calculated on the present context.
    //@{
    double get_parental_mean_adjustment() const;
    double get_spousal_mean_adjustment() const;
    double get_sibling_mean_adjustment(size_t                   prior_sib_count,
                                       const PenetranceContext& context) const;

    bool   get_mother_informativity() const;
    bool   get_father_informativity() const;
    //@}
    
  private:

    void calculate_family_statistics   (const PenetranceContext&      context);
    void calculate_spousal_statistics  (const PenetranceContext&      context,
                                        const TypeDescription::State& spouse_state);
    void calculate_maternal_statistics (const PenetranceContext&      context,
                                        const TypeDescription::State& spouse_state);
    void calculate_paternal_statistics (const PenetranceContext&      context,
                                        const TypeDescription::State& spouse_state);
    
    void calculate_sibling_mean_adjustments(const PenetranceContext& context) const;

    /// \name Basic Data Members
    //@{
    const residual_correlation_sub_model& my_corr_sub_model;
    const continuous_member_calculator&   my_member_calc;
    const RegPenetranceCommonAdjustments& my_corr_adjustments;
    //@}

    bool my_mother_informative;    ///< Is mother informative?
    bool my_father_informative;    ///< Is father informative?
    
    double my_alpha_m;    /// \f$\alpha_m\f$ given parental informativity.
    double my_alpha_f;    /// \f$\alpha_f\f$ given parental informativity.

    /// \name Adjustments
    //@{
    double my_spousal_mean_adjustment;
    double my_maternal_mean_adjustment;
    double my_paternal_mean_adjustment;

    /// Storage for sibling mean adjustments.
    /// 
    /// \note This data member is mutable to allow for lazy calculation
    ///       during the get_sibling_mean_adjustment() function
    mutable vector<double> my_sibling_mean_adjustments;
    //@}
  
};
  
/// The PenetranceContext stores a context under which penetrance
/// calculations can be done for the regressive models.
///
/// Witin the regressive models, penetrance calculations are dependant upon
/// the context in which they are calculated.  A context is a set of
/// individuals that affect a particular set of penetrance calculations. The
/// PenetranceContext creates and stores these contexts for use by the 
/// RegPenetranceCalculator.  
///
/// The context within which penetrances are calculated can be one of four
/// sets:
///
/// -# Nobody.  In this case, penetrances can be calculated for any
///    individual as if they were independent (eq. 47)
/// -# A single individual with spouse.  Penetrances can be calculated for
///    the individual only, dependent on spouse (eq. 48)
/// -# A family.  Penetrances can be calculated for the children, dependent
///    on parents and possibly siblings (eqs. 48a, 64)
/// -# A family with a primary sibling and their spouse.  Penetrance can be
///    calculated for all siblings, including the primary one, dependent on
///    parents, siblings and the spouse (Only the primary sib will use the
///    spouse data). (eqs. 49a, 66)
///
/// When a PenetranceContext is first instantiated, it is in the first
/// state (nobody in the context).  A family and/or a individual/spouse can
/// be added to change to the appropriate context.  Once set, these should
/// not be changed without first calling a function to clear the current
/// context (resetting to set 1).
///
/// After creation, parents and the spouse (if included) must be given
/// states which will affect the penetrance calculations.  These
/// states may be changed at will without invalidating the context.
///
/// As context is added, adjustments to the penetrance calculations due to
/// correlations with the context individuals are pre-computed.  These
/// adjustments are then accessed as needed by the RegPenetranceCalculator. 
/// Note that these adjustments are portions of equations in the SEGREG
/// Formula Document, and do not have their own equations in that paper. 
/// The equations which use the particular adjustments are listed in the
/// following table.
///
/// <table>
///   <TR>
///     <TD> \b Adjustment  </TD>
///     <TD> \b Formula     </TD>
///     <TD> \b Description </TD>
///   </TR>
///
///   <TR>
///     <TD> Spousal Mean Adjustment </TD>
///     <TD> \f$-\rho_{fm} d_u(s)\f$ </TD>
///     <TD> Used when spouse is included in penetrance in regressive
///          models. (eqs. 48, 49a, and 66)</TD>
///   </TR>
///
///   <TR>
///     <TD> Maternal Mean Adjustment </TD>
///     <TD> \f$-\alpha_m d_u(m) \f$ </TD>
///     <TD> Used when mother is included in the penetrance in regressive
///          models (eqs. 48a, 49a, 64 and 66)</TD>
///   </TR>
///
///   <TR>
///     <TD> Paternal Mean Adjustment </TD>
///     <TD> \f$-\alpha_f d_u(f) \f$ </TD>
///     <TD> Used when father is included in the penetrance in regressive
///          models (eqs. 48a, 49a, 64 and 66)</TD>
///   </TR>
///
///   <TR>
///     <TD> Sibling Mean Adjustment </TD>
///     <TD> \f[ - \beta_j \sum_{l=1}^{j} \hat{d}_{b_l} \f]
///     </TD> 
///     <TD> Siblings affect the penetrance for class D models.  However,
///          they are dependent upon both parents states.  Therefore, sibling
///          calculations are performed only when a penetrance is called for. 
///          Then all sibling adjustments are calculated and stored for later
///          retrieval. </TD>
///   </TR>
/// </table>
class PenetranceContext
{
  public:
  
    friend class RegPenetranceCalculator;
    friend class PenetranceContextStatistics;

    typedef FPED::Member             member_type;
    typedef FPED::MemberConstPointer member_pointer;
    typedef FPED::Family             family_type;

    PenetranceContext(const residual_correlation_sub_model& rsm,
                      const continuous_member_calculator&   cmc,
                      const RegPenetranceCommonAdjustments& corr_adj);

    /// \name Context initialization routines
    ///
    /// Notes:  Context are only needed if the penetrances need to be
    ///         computed conditionally on the family and/or spouse of the
    ///         individuals.
    //@{
    void clear_family_info();

    void set_nuclear_family(const family_type& family, bool include_sibs = true);
    void set_link_couple(const member_type& ind, const member_type& spouse);

    void set_link_spouse_state(const TypeDescription::State& spouse_state);
    void set_mother_state     (const TypeDescription::State& mother_state);
    void set_father_state     (const TypeDescription::State& father_state);
    //@}

    /// \name Basic Context Information Accessors
    //@{
    bool includes_family  () const;
    bool includes_spouse  () const;
    bool includes_siblings() const;
    
    bool incorporate_spousal_adjustments(const member_type& i) const;
    //@}
    
    /// \name Basic Accessors
    ///
    /// These routines return basic info that has been set on the context
    //@{
    const member_type& get_individual  () const;
    const member_type& get_mother      () const;
    const member_type& get_father      () const;
    const member_type& get_link_spouse () const;
    
    const family_type& get_family() const;

    const TypeDescription::State& get_mother_state() const;
    const TypeDescription::State& get_father_state() const;
    //@}

  private:

    /// Family of the current context
    FPED::FamilyConstPointer my_current_family;

    /// \name Context Individuals
    //@{
    member_pointer my_link_ind;
    member_pointer my_link_spouse;
    //@}

    /// \name Parental Information
    //@{
    TypeDescription::State my_mother_state;    ///< Current maternal state
    TypeDescription::State my_father_state;    ///< Current paternal state
    //@}

    /// Tells whether this context include sibs in class D models.  Class A
    /// models ignore this variable.
    bool my_include_siblings;
    
    PenetranceContextStatistics my_stats;
};

}
}

#include "segreg/RegPenetranceContext.ipp"

#endif
