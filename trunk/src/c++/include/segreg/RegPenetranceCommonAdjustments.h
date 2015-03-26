#ifndef REG_PENETRANCE_COMMON_ADJUSTMENTS_H
#define REG_PENETRANCE_COMMON_ADJUSTMENTS_H

#include <iostream>
#include "boost/array.hpp"
#include "segreg/definitions.h"
#include "segreg/resid_sub_model.h"

namespace SAGE   {
namespace SEGREG {

/// The RegPenetranceCommonAdjustments calculates and stores adjustments to
/// penetrance calculations in the regressive models due to correlations.
///
/// Within the regressive models, the penetrance calculations include
/// adjustments due to correlations.  Many of these adjustments are shared
/// among many individuals.  The RegPenetranceCommonAdjustments object
/// calculates and stores these adjustments so that they may be more
/// efficiently calculated only once per function evaluation.
///
/// Note that these adjustments are portions of equations in the SEGREG
/// Formula Document, and do not have their own equations in that paper. 
/// The equations which use the particular adjustments are listed in the
/// following table.
///
/// The following table lists the adjustments that are calculated in this fashion.
///
/// <table>
///   <TR>
///     <TD> \b Adjustment  </TD>
///     <TD> \b Formula     </TD>
///     <TD> \b Description </TD>
///   </TR>
///
///   <TR>
///     <TD> Spousal Variance Adjustment </TD>
///     <TD> \f$ -\rho^2_{fm}\f$ </TD>
///     <TD> Used when spouse is included in penetrance in regressive
///          models. (eqs. 48, 49a, and 66)</TD>
///   </TR>
///
///   <TR>
///     <TD> Parental Variance Adjustment </TD>
///     <TD> \f$ -\alpha_m \rho_{ms} -\alpha_f \rho_{fs} \f$ </TD>
///     <TD> Used when parents are included in the penetrance in regressive
///          models.  There are 4 possible values, dependent upon the
///          informativity of the mother and the father (see eqs. 52 and
///          53).  Used by eqs. 48a, 49a, 64 and 66.</TD>
///   </TR>
///
///   <TR>
///     <TD> Prior Sibling Coefficient (eq. 60) </TD>
///     <TD> \f$ \beta_j =
///                   \left\{
///                     \begin{array}{ll}
///                         \frac{\rho_{bb} - \delta}
///                              {1 - \rho_{bb} + j (\rho_{bb} - \delta)} &
///                         \mbox{if $j>0$}
///                     \\  0                                             &
///                         \mbox{if $j=0$}
///                     \end{array}
///                   \right. \f$
///     </TD> 
///     <TD> Used for siblings prior to the present sibling in class D
///          models.  \em j is the sibling's index within the sibling order. 
///          There are 4 possible values for each \em j, dependent upon the
///          validity of the mother and the father. (eqs. 64 and 66) </TD>
///   </TR>
///   <TR>
///     <TD> Sibling Variance Adjustment </TD>
///     <TD> \f$ -j \beta_j \rho_{bb} \f$ </TD>
///     <TD> Used for siblings prior to the present sibling in class D
///          models.  j is the sibling's index within the sibling order.
///          There are 4 possible values for each j, dependent upon the
///          validity of the mother and the father. (eqs. 64 and 66) </TD>
///   </TR>
/// </table>
class RegPenetranceCommonAdjustments
{
public:
    /// \name Constructors & operators
    //@{

    /// Default Constructor
    RegPenetranceCommonAdjustments();

    /// Constructor
    ///
    /// \param m_class          The model class of the model we're computing.
    ///                         Should be either model_A or model_D.
    ///
    /// \param corr_sub_model   The correlation submodel from which
    ///                         correlations (rho's and deltas) will be
    ///                         taken.
    ///
    /// \param max_sibship_size The maximum size of sibship in the sample.
    ///                         You can use the getMaxSibshipSize() function
    ///                         in the RPED::RefMultiPedigree to get this value or
    ///                         calculate it directly from the data sample.
    RegPenetranceCommonAdjustments(model_class                           m_class,
                                   const residual_correlation_sub_model& corr_sub_model,
                                   unsigned int                          max_sibship_size = 0);

    /// Copy constructor
    ///
    /// \param other The object which will be copied.
    RegPenetranceCommonAdjustments(const RegPenetranceCommonAdjustments & other);

    /// Assignment operator
    /// \param other The object which will be copied.
    RegPenetranceCommonAdjustments& operator= (const RegPenetranceCommonAdjustments & other);

    //@}

    /// \brief Updates the internal data based upon the status of the
    ///        residual_correlation_sub_model passed at construction time.
    ///
    /// update() must be invoked once per Maxfun function evaluation. Make
    /// sure you invoke update() \b after updating the
    /// residual_correlation_sub_model.
    ///
    /// \retval 0 The calculator updated its internal values succesfully.
    /// \retval 1 The calculator did not update its internal values succesfully.
    int update();

    /// \name Data Accessors
    //@{

    /// \brief Returns current Spousal Variance Adjustment.
    ///
    /// Calcualted as: \f$\-\rho^2_{fm}\f$.  Used for eqs. 48, 49a and 66.
    double get_spousal_var_adjustment() const;

    /// \brief Returns current Parental Variance Adjustment based on
    ///        parental informativity.
    /// 
    /// Calculated as: \f$-\alpha_m\rho_{ms}-\alpha_f\rho_{fs}\f$.  Used for
    /// eqs. 48a, 49a, 64 and 66.
    ///
    /// \param mother_informative \c true if the mother is informative,
    ///                           \c false if mother is noninformative.
    /// \param father_informative \c true if the father is informative,
    ///                           \c false if father is noninformative.
    double get_parental_var_adjustment(bool mother_informative,
                                       bool father_informative) const;

    /// \brief Returns current Prior Sibling Coefficient based on parental
    ///        informativity and prior sibling count.
    ///
    /// Prior Sibling Coefficient is calculated by eq. 60 as:
    ///
    /// \f$ \beta_j =
    ///               \left\{
    ///                   \begin{array}{ll}
    ///                       \frac{\rho_{bb} - \delta}
    ///                            {1 - \rho_{bb} + j (\rho_{bb} - \delta)} &
    ///                       \mbox{if $j>0$}
    ///                   \\  0                                             &
    ///                       \mbox{if $j=0$}
    ///                   \end{array}
    ///               \right.
    /// \f$
    ///
    /// Used by eqs. 64 and 66 in class D models.
    ///
    /// \param mother_informative \c true if the mother is informative,
    ///                           \c false if mother is noninformative.
    /// \param father_informative \c true if the father is informative,
    ///                           \c false if father is noninformative.
    /// \param prior_sib_count    The number of prior sibs for this individual.
    double get_prior_sib_coefficient(bool          mother_informative,
                                     bool          father_informative,
                                     unsigned int  prior_sib_count    ) const;

    /// \brief Returns current Sibling Variance Adjustment based on parental
    ///        informativity and prior sibling count.
    ///
    /// Calculated as: \f$ -j \beta_j \rho_{bb} \f$.  Used for eqs. 64 and
    /// 66 in class D models.
    ///
    /// \param mother_informative \c true if the mother is informative,
    ///                           \c false if mother is noninformative.
    /// \param father_informative \c true if the father is informative,
    ///                           \c false if father is noninformative.
    /// \param prior_sib_count    The number of prior sibs for this individual.
    double get_sibling_var_adjustment(bool          mother_informative,
                                      bool          father_informative,
                                      unsigned int  prior_sib_count    ) const;

    //@}

private:

    /// Copies from one object to another (used by copy constructor and operator=).
    void copy(const RegPenetranceCommonAdjustments & other);

    /// \name Calculation Functions
    //@{

    void calculate_spousal_var_adjustment   ();
    void calculate_parental_var_adjustments ();

    /// Calcualtes the prior sibling coefficients and the sibling variance
    /// adjustments
    void calculate_sib_values ();

    //@}

    /// \name Necessary Construction Information
    //@{

    model_class                            my_model_class;
    const residual_correlation_sub_model*  my_corr_sub_model;
    unsigned int                           my_max_sibs;

    //@}
    
    /// \name Data Storage
    //@{

    /// The InformativeValueArray is a set of four values (doubles), indexed by the
    /// informativity of the parents.
    class InformativeValueArray
    {
        /// Array storage for the doubles
        typedef boost::array<double, 4> ValueArray;

    public:

        InformativeValueArray();

        /// \brief Set the value specified by the parental informativity
        ///
        /// \param mother_informative \c true if the mother is informative,
        ///                           \c false if mother is noninformative.
        /// \param father_informative \c true if the father is informative,
        ///                           \c false if father is noninformative.
        void set_value(bool mother_informative, bool father_informative, double value);

        /// \brief Get the value specified by the parental informativity
        ///
        /// \param mother_informative \c true if the mother is informative,
        ///                           \c false if mother is noninformative.
        /// \param father_informative \c true if the father is informative,
        ///                           \c false if father is noninformative.
        double  get_value(bool mother_informative, bool father_informative) const;

    private:

        /// \brief Helper function used to calculate the index into
        ///        my_double_array based on parental informativity.
        ///
        /// \param mother_informative \c true if the mother is informative,
        ///                           \c false if mother is noninformative.
        /// \param father_informative \c true if the father is informative,
        ///                           \c false if father is noninformative.
        ValueArray::size_type calculate_index(bool mother_informative,
                                              bool father_informative) const;

        ValueArray my_values;
    };
    
    /// \brief Spousal Variance Adjustment
    double my_spousal_var_adj;

    /// \brief Parental Variance Adjustments
    ///
    /// \b Note: that when maternal and paternal informativity is \c false, the
    ///          adjustment is always 0.0.  We only calculate non-zero
    ///          cases.
    InformativeValueArray my_parental_var_adjs;

    /// \brief Prior Sibling Coefficients
    vector<InformativeValueArray> my_prior_sib_coeffs;

    /// \brief Sibling Variance Adjustments
    vector<InformativeValueArray> my_sib_var_adjs;

    //@}
};

/// \brief Given a ostream, a model_class and a residual_correlation_sub_model,
///        creates a RegPenetranceCommonAdjustments object and updates it, then dumps
///        the calculated values to the ostream.
///
///  If the model is class D, then it assigns a max_sibship_size of 5 for the
///  test.
void test_common_reg_corr_adjustments(std::ostream& o, model_class m,
                                      const residual_correlation_sub_model& r);

} // End namespace SEGREG
} // End namespace SAGE

#include "segreg/RegPenetranceCommonAdjustments.ipp"

#endif
