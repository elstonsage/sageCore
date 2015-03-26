#ifndef MEMBER_COVARIATE_CALCULATOR_H
#define MEMBER_COVARIATE_CALCULATOR_H
//=======================================================================
//
//  File:	MemberCovariateCalculator.h
//
//  Author:	Stephen Gross
//
//  Copyright 2002 R. C. Elston
//
//=======================================================================


#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include "mped/mp.h"
#include "mped/sp.h"
#include "rped/rped.h"
#include "assoc/Datatypes.h"
#include "assoc/Configuration.h"
#include "assoc/Residuals.h"

namespace SAGE  {
namespace ASSOC {

/// \class MemberCovariateCalculator
///
/// The MemberCovariateCalculator (MCC) calculates the following equations
/// (for every individual) from the internal documentation:
///
/// -# D2 \f$ y_i \f$ (main phenotype)
/// -# D9 \f$ (x_{i1}, x_{i2}, ..., x_{ip}) \f$ (vector of centered covariates)
/// -# D10 \f$ \mu_i \f$ (transformed phenotypic mean)
///  Note \f$ \mu_i = h(\beta_0 + \beta_1x_{i1} + ... + \beta_px_{ip}) \f$
///  and \f$ \beta \f$ 's can be obtained from the regression model
///  \f$ E(y_i) = \beta_0 + \beta_1x_{i1} + ... + \beta_px_{ip}, i = 1,2,...,n \f$
/// -# D11 \f$ t_i \f$ (transformed phenotype)
/// -# D24 \f$ z_i = t_i - \mu_i \f$ (mean adjusted phenotype)
///
/// In addition, the MCC invokes the transformation equations. The implemenentation
/// of the transformation, however, takes place inside the transformation_sub_model.

class MemberCovariateCalculator
{
  public:
    MemberCovariateCalculator(const Configuration& config, 
                              const MAXFUN::ParameterMgr& mgr, const Sampledata& sample, 
                              MFSUBMODELS::Transformation::Facade& facade, double divisor = 1.0);
                              
    virtual ~MemberCovariateCalculator();

    virtual  int update(MAXFUN::ParameterMgr&) = 0;    
    double  getZScore(size_t id) const;
    const vector<double>&  getDiffs() const;
    virtual Residuals  getResiduals() const = 0;
    virtual double  getMeanResidualsAffected() const;
    virtual double  getMeanResidualsUnaffected() const;
    
  protected:
    void calculatePhenotypicMeans();
    virtual bool  doTransformations() = 0;
    virtual void  calculateZScore(int id) = 0;
    void dumpValues() const;    

    // Data members.
    const Configuration&  my_config;
    const MAXFUN::ParameterMgr&  my_mgr;
    MAXFUN::ParameterConstIterator  my_cov_param_begin;
    MAXFUN::ParameterConstIterator  my_cov_param_end;
    int  my_intercept_idx;    
    
    const Sampledata&  my_sampledata;
    SAMPLING::FieldConstIterator  my_cov_field_begin;
    const SAMPLING::Field&  my_main_phenotype_field;        
    
    MFSUBMODELS::Transformation::Facade&  my_transf_facade;

    vector<double>  my_phenotypes;
    vector<double>  my_phenotypic_means;
    vector<double>  my_diffs;       
    vector<double>  my_z_scores;
    double  my_divisor;          // Used to adjust residuals during transformation.
};


// - For a continuous primary trait where the primary trait and the
//   expected value of the primary trait are to be transformed separately.
//
  class MccContinuousBoth : public MemberCovariateCalculator
  {
    public:
      MccContinuousBoth(const Configuration& config,
                        const MAXFUN::ParameterMgr& mgr, const Sampledata& sample,
                        MFSUBMODELS::Transformation::Facade& facade, double divisor = 1.0);
                       
      virtual ~MccContinuousBoth();
      
      int update(MAXFUN::ParameterMgr&);      
      
      Residuals  getResiduals() const;
  
    protected:
      void  calculateDiffs();
      bool  doTransformations();
      void  calculateZScore(int id);     
    
      // Data Members
      vector<double>  my_transformed_phenotypes;
      vector<double>  my_transformed_phenotypic_means;
    
    private:
  };


// - For a continuous primary trait where the difference between the primary
//   trait and the expected value of the primary trait are to be transformed.
//
class MccContinuousDiff : public MemberCovariateCalculator
{
  public:
    MccContinuousDiff(const Configuration& config, 
                      const MAXFUN::ParameterMgr& mgr, const Sampledata& sample, 
                      MFSUBMODELS::Transformation::Facade& facade, double divisor = 1.0);
                      
    virtual ~MccContinuousDiff();
    
    int update(MAXFUN::ParameterMgr&);
    
    Residuals  getResiduals() const;
  
  protected:
    virtual void  calculateDiffs();
    bool  doTransformations();
    void  calculateZScore(int id);
    
    // Data Members

    vector<double>  my_transformed_diffs;  
};

class MccBinaryDiff : public MccContinuousDiff
{
  public:
    MccBinaryDiff(const Configuration& config, 
                  const MAXFUN::ParameterMgr& mgr, const Sampledata& sample, 
                  MFSUBMODELS::Transformation::Facade& facade, double divisor = 1.0);
                      
    virtual ~MccBinaryDiff();
    
    double getMeanResidualsAffected() const;
    double getMeanResidualsUnaffected() const;
    
  protected:
    virtual void  calculateDiffs();
      virtual void  calculateNewPhenotypes();       // (y + 1 / 2n) 
      void  calculateBinaryPhenotypicMeans();       // e ** Mu / (1 + e ** Mu)
      
    // Data Members
    vector<double>  my_binary_phenotypic_means;     // e ** Mu / (1 + e ** Mu)
    vector<double>  my_new_phenotypes;
    vector<double>  my_affected_diffs;
    vector<double>  my_unaffected_diffs;
};


#include "assoc/MemberCovariateCalculator.ipp"

} 
} 

#endif
