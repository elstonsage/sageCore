#ifndef SEGREG_FPMM_SUB_MODEL_H
#define SEGREG_FPMM_SUB_MODEL_H
//============================================================================
// File:      fpmm_sub_model.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   8/22/01  - created.                          djb
//                                                                          
// Notes:     defines the finite polygenic mixed model sub-model for SEGREG.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include "numerics/binomial_dist.h"
#include "sub_model_base.h"

namespace SAGE
{

namespace SEGREG
{

/// @name fpmm sub-model constants
//@{
extern const std::string  FPMM_NAME;
extern const double       FPMM_DEFAULT_FREQ;
extern const size_t       FPMM_DEFAULT_LOCI;
extern const size_t       FPMM_MAX_LOCI;
extern const double       FPMM_DEFAULT_VALUE;              // 1;
extern const double       FPMM_EPSILON;
extern const double       FPMM_LB;
extern const double       FPMM_UB;
extern const bool         FPMM_DEFAULT_FIXED;
//@}

//----------------------------------------------------------------------------
//  Class:    finite_polygenic_mixed_model_sub_model
//                                                                          
//----------------------------------------------------------------------------
//
class finite_polygenic_mixed_model_sub_model : public SegregSubmodel
{
  public:
    
    // Constructor/destructor.  
    finite_polygenic_mixed_model_sub_model(cerrorstream& errors = sage_cerr);
    finite_polygenic_mixed_model_sub_model(const finite_polygenic_mixed_model_sub_model& other);
    finite_polygenic_mixed_model_sub_model&  operator=(const finite_polygenic_mixed_model_sub_model& other);
    virtual ~finite_polygenic_mixed_model_sub_model();
    
    // Gets.
    virtual string  name() const;
    string          option_description() const;
    double          variance() const;
    double          frequency() const;
    size_t          loci() const;
    size_t          max_pgt() const;

    double          mean(size_t polygeno) const;
    double          pop_freq(size_t polygeno) const;
    double          pop_freq(size_t polygeno, size_t spouse_pg, double corr) const;

    bool variance_fixed(); // due to JA
    
    void  dump(std::ostream& out) const;

    // Sets.
    bool  set(model_input var, double freq, size_t loci); 

    bool is_complete() const;
  
  protected:
    virtual int update();
    
    virtual int finalizeConfiguration();

    void  calculate_means();
    
    // Ancillary functions.
    bool  input_meets_constraints(model_input& input);
    
    // Data members. 
    double   my_variance;       // Variance of polygenic variable.
    bool     my_variance_fixed; ///< true if the variance is fixed, false otherwise.
    double   my_frequency;      // Allele frequency of polygenic loci.
    size_t   my_max_pgt;        // Max polygenotypes  2 * loci + 1

    vector<double> my_means; // Polygenic means calculated using D11
};


} 
} 
#include "segreg/fpmm_sub_model.ipp"

#endif


