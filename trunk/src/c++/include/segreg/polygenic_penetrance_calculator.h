#ifndef POLYGENIC_PENETRANCE_CALCULATOR_H
#define POLYGENIC_PENETRANCE_CALCULATOR_H
//============================================================================
// File:      	polygenic_penetrance_calculator.h
//                                                                          
// Author:    	Stephen Gross
//                                                                          
// History:   	djb  Created.			  Jul  9 01
//            	sag  Added probability functions  Jul 20 01
//		sag  Added documentation	  Aug 02 01
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================



#include "segreg/model.h"
#include "segreg/member_calculator.h"
#include "segreg/sub_model_base.h"
#include "segreg/resid_sub_model.h"
#include "segreg/segreg_datatypes.h"
#include "segreg/ContinuousPenetranceFunction.h"
#include <math.h>

namespace SAGE {
namespace SEGREG {


//==================================================================
//
//  Class:    polygenic_penetrance_calculator
//                                                                          
//==================================================================
/**             
 *              This class provides all sorts of penetrance functions.  A
 *		penetrance functions is usually written in the form:
 *
 *			P(t|...) or,
 *
 *                      P(y|...)
 *
 *		There are two polygenic penetrance functions.  Both use the
 *		individual information, while the second also includes the
 *		spouse information.
 *
 *		For continuous traits, the internal functioning of the
 * 		penetrance calculator works like this:
 *
 *		(1) The client calls the appropriate get_penetrance(...)
 *		    function, based on what info will be incorporated
 *		(2) The function calculates a "z" and a "w" value, which
 *		    are then passed on to the internal penetrance_value(...)
 *		    function, which processes all z's and w's through the
 *		    same algorithm, and shoots back the result, which is
 *		    then returned to the client.
 *
 *              For binary traits, it uses values collected from the
 *              binary_member_calculator to perform a logit function.
 */
class polygenic_penetrance_calculator
{
  public:
  
    typedef FPED::Member             member_type;
    typedef FPED::MemberConstPointer member_pointer;

    //==============================================================
    // Constructor/destructor/operators:
    //==============================================================

    polygenic_penetrance_calculator
          (const FPED::Multipedigree& ped_data,
           const model &              modp,
           bool                       use_ascertainmentp = false);

    int update();

    //==============================================================
    // Structs:
    //==============================================================

    struct penetrance_info 
    { 
      penetrance_info(const member_type& index_param,
                      genotype_index     genotype_param     = index_AA,
                      size_t             polygenotype_param = NO_POLYGENIC_DATA)
      { member = &index_param; genotype = genotype_param; polygenotype = polygenotype_param; }
      penetrance_info()
      { member = NULL; genotype = index_AA; polygenotype = NO_POLYGENIC_DATA; }

      penetrance_info(const penetrance_info& other)
      {
        member       = other.member;
        genotype     = other.genotype;
        polygenotype = other.polygenotype;
      }

      member_pointer member; 
      genotype_index genotype; 
      size_t         polygenotype; 
    };

    //================================================================
    // Penetrance functions:
    //================================================================

    double get_polygenic_penetrance (penetrance_info i)    const;

    //================================================================
    // Private data members:
    //================================================================

  private:
    const model & mod;

    continuous_member_calculator my_cont_member_calc;
    binary_member_calculator     my_bin_member_calc;
    onset_member_calculator      my_onset_member_calc;

    //================================================================
    // Private member functions:
    //================================================================

    double prior_sibship(size_t prior_sib_count, double delta) const;

    double continuous_penetrance (penetrance_info i)    const;
    double binary_penetrance     (penetrance_info i)    const;
    double onset_penetrance      (penetrance_info i)    const;
};

// End namespace
} 
} 

#include "segreg/polygenic_penetrance_calculator.ipp"

#endif
