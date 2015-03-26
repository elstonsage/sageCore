#ifndef CONTINUOUS_PENETRANCE_FUNCTION_H
#define CONTINUOUS_PENETRANCE_FUNCTION_H

#include "numerics/functions.h"
#include "numerics/normal_pdf.h"
#include "segreg/member_calculator.h"

namespace SAGE {
namespace SEGREG {

/// This function calculates the penetrance for continuous data given the z
/// and w values.  Primarily, this means calculating equ. 46 from the SEGREG
/// document.
///
/// However, if the individual is missing or a member of the proband
/// sampling frame under certain kinds of ascertainment, the equation
/// changes.
///
/// The equations are:
///
/// <table>
///   <TR>
///     <TD> \b Formula  </TD>
///     <TD> <B> If individual is ... </B> </TD>
///   </TR>
///   <TR>
///     <TD> \f$ \frac{1} {\sqrt{ 2 \pi w_i } } exp \left( - \frac{ z_i^2 } { 2 w_i } \right) \f$ </TD>
///     <TD> normal.  This is the default case (eq. 46)                                </TD>
///   </TR>
///   <TR>
///     <TD> \f$ 1.0 \f$                                                               </TD>
///     <TD> missing or invalid.                                                       </TD>
///   </TR>
///   <TR>
///     <TD> \f$ \Phi \left( - \frac {z_{i (high)} } { \sqrt{w_i} } \right) \f$        </TD>
///     <TD> a member of the PSF to be computed using the Greater Than Threshold 
///          formula (see Appendix C).                                                 </TD>
///   </TR>
///   <TR>
///     <TD> \f$ \Phi \left( \frac {z_{i (high)} } { \sqrt{w_i} } \right) \f$          </TD>
///     <TD> a member of the PSF to be computed using the Greater Than Threshold 
///          formula (see Appendix C).                                                 </TD>
///   </TR>
// 
double calculate_continuous_penetrance
    (double                               z,
     double                               w,
     member_calculator_base::member_class mt);

}
}

#include "segreg/ContinuousPenetranceFunction.ipp"

#endif
