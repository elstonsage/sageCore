#ifndef AO_MEMBER_COVARIATE_CALCULATOR_H
#define AO_MEMBER_COVARIATE_CALCULATOR_H
//======================================================================
//
//  File:	MemberCovariateCalculator.h
//
//  Author:	Stephen Gross
//
//  Copyright 2002 R. C. Elston
//
//  The MemberCovariateCalculator is responsible for updating the
//  composite susceptibilities, means, and variances of every individual
//  in the sample. The update() function is called for every iteration
//  of maxfun, which cause the MemberCovariateCalculator to
//  recalculate all the aforementioned values (as based on the newly 
//  updated Model).
//
//======================================================================


#include <string>
#include <vector>
#include "mped/mp.h"
#include "mped/sp.h"
#include "rped/rped.h"
#include "sampling/sampling.h"
#include "ageon/Datatypes.h"
#include "ageon/Model.h"
#include "ageon/ParamFieldCache.h"

namespace SAGE {
namespace AO   {

//======================================================================
//
// class MemberCovariateCalculator
//
//======================================================================
class MemberCovariateCalculator
{
  public:
    //==================================================================
    // Constructors & operators:
    //==================================================================

    MemberCovariateCalculator(const Model      &, 
                                         const SAMPLING::PartitionedMemberDataSample & sample,
                                         int);
    MemberCovariateCalculator(const MemberCovariateCalculator &);

    //==================================================================
    // Public utility functions:
    //==================================================================

    void reset     ();
    int  update    ();
    void copy      (const MemberCovariateCalculator &);
    int  transform (double &);

    //==================================================================
    // Public accessors:  
    //==================================================================

    const SAMPLING::PartitionedMemberDataSample & getSample() const { return my_sample; }


    double get_genetic_suscept (size_t id) const;
    double get_AO_mean         (size_t id) const;
    double get_AO_var          (size_t id) const;
    double get_AO_stdev        (size_t id) const;
    double get_AO_transf       (size_t id) const;
    double get_AE_transf       (size_t id) const;
    double get_zero_transf     ()          const;

    // Debugging:

    void dumpContents();

  private:
    //==================================================================
    // Private utility functions:
    //==================================================================

    int calculate_genetic_suscepts    ();
    int calculate_AO_means            ();
    int calculate_AO_vars             ();
    int calculate_AO_stdevs           ();
    int precalculate_transformed_vals ();
    int transform_AOs                 ();
    int transform_AEs                 ();
    int calculate_zero_transf         ();

    //==================================================================
    // Data members:
    //==================================================================

    const Model       & my_model;
    const SAMPLING::PartitionedMemberDataSample & my_sample;
    int              my_analysis_type;

    ParamFieldCache my_cache;

          size_t           my_num_of_inds;
          vector<double>   my_genetic_suscept;
          vector<double>   my_AO_mean;
          vector<double>   my_AO_var;
          vector<double>   my_AO_stdev;
	  vector<double>   my_transformed_AO;
	  vector<double>   my_transformed_AE;
          double           my_zero_transf;

	// Performance enhancers:

          vector<double>   my_suscept_intercepts;

	  vector<double>   my_transformed_vals;
};


inline double MemberCovariateCalculator::get_genetic_suscept (size_t id) const { return my_genetic_suscept [id]; }
inline double MemberCovariateCalculator::get_AO_mean         (size_t id) const { return my_AO_mean         [id]; }
inline double MemberCovariateCalculator::get_AO_var          (size_t id) const { return my_AO_var          [id]; }
inline double MemberCovariateCalculator::get_AO_stdev        (size_t id) const { return my_AO_stdev        [id]; }
inline double MemberCovariateCalculator::get_AO_transf       (size_t id) const { return my_transformed_AO  [id]; }
inline double MemberCovariateCalculator::get_AE_transf       (size_t id) const { return my_transformed_AE  [id]; }
inline double MemberCovariateCalculator::get_zero_transf     ()          const { return my_zero_transf;          }


//======================================================================
//
//  transform(...)
//
//======================================================================
inline int
MemberCovariateCalculator::transform(double & x)
{ 
  // 1. Set up local variables:
  
        double lambda1 = my_model.GetParameterMgr()(my_model.lambda1_mxid),
               lambda2 = my_model.GetParameterMgr()(my_model.lambda2_mxid);

  // 2. Adjust x and check it:

	double x_ = x + lambda2;

	if(x_ <= 0.0)
	  return 1;

  // 3. Transform value:
  
	if(lambda1 != 0.0)
	  x = (pow(x_, lambda1) - 1) / lambda1;
	else
	  x = log(x_);

  // 4. Return success:

	return 0;
} 

}} // End namespace

#endif
