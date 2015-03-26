#ifndef MEMBER_COVARIATE_DATA_H
#define MEMBER_COVARIATE_DATA_H
//=======================================================================
//  File:	member_covariate_data.h
//
//  Purpose:	Provide data structures and functionality to calculate
//		likelihoods of variates and covariates.
//
//  Author:	Stephen Gross
//
//  History:	0.1  sag  Initial Implementation	Jul 11 01
//
//  Copyright (c) 2001 R.C. Elston
//  All Rights Reserved
//=======================================================================


#include <vector>
#include <cmath>
#include "segreg/segreg_datatypes.h"
#include "segreg/ascertainment_sub_model.h"

using namespace std;

namespace SAGE   {
namespace SEGREG {

//=======================================================================
//  Class:	member_covariate_data
//
//  Purpose:    This class stores information about variate and covariate
//		traits.
//=======================================================================
class member_covariate_data
{
  public:

    /** the member_type (from ascertainment) determines what type of calculation
     * is to be performed on an individual under ascertainment.  By default,
     * this type is set to 'actual', which means that nothing special has to
     * be done if ascertainment isn't being performed. */
    typedef ascertainment_sub_model::v_sm_option member_type;

    member_covariate_data();
    void set_number_of_covariates(size_t Number_of_Composite_Covariates,
                                  size_t Number_of_Susceptibility_Covariates,
                                  size_t Number_of_Covariate_Covariates,
                                  size_t Number_of_Variance_Covaraites);

    size_t get_composite_cov_count      () const;
    size_t get_susceptibility_cov_count () const;
    size_t get_mean_cov_count           () const;
    size_t get_variance_cov_count       () const;

    double get_analysis_trait     () const;
    double set_analysis_trait     (double value      = numeric_limits<double>::quiet_NaN());
    double get_composite_cov      (size_t TraitIndex) const; 
    bool   set_composite_cov      (size_t TraitIndex,
                                   double value      = numeric_limits<double>::quiet_NaN());
    double get_susceptibility_cov (size_t TraitIndex) const; 
    bool   set_susceptibility_cov (size_t TraitIndex,
                                   double value      = numeric_limits<double>::quiet_NaN());
    double get_mean_cov           (size_t TraitIndex) const;
    bool   set_mean_cov           (size_t TraitIndex,
                                   double value      = numeric_limits<double>::quiet_NaN());
    double get_variance_cov       (size_t TraitIndex) const;
    bool   set_variance_cov       (size_t TraitIndex,
                                   double value      = numeric_limits<double>::quiet_NaN());

    void        set_member_type(member_type t);
    member_type get_member_type() const;

    // The valid options are really boolean options on the member type, returning
    // !missing.
    //@{
    bool   validate();       //< Verifies that there are no missing trait/cov values
    bool   is_valid() const; //< Returns !missing.
    //@}

  private:
    vector<double>     composite_covariates;
    vector<double>     susceptibility_covariates;
    vector<double>     mean_covariates;
    vector<double>     variance_covariates;
    double             Analysis_Trait;
    member_type        my_member_type;
};

}} // End namespace

#include "member_covariate_data.ipp"

#endif
