//=======================================================================
//
//  File:  MemberCovariateCalculator.cpp
//
//  Author:  Stephen Gross
//
//  Copyright 2002 R. C. Elston
//=======================================================================

#include "assoc/MemberCovariateCalculator.h"

namespace SAGE  {
namespace ASSOC {

//=======================================================================
//  MemberCovariateCalculator 
//=======================================================================
MemberCovariateCalculator::MemberCovariateCalculator(const Configuration& config,
                                                     const MAXFUN::ParameterMgr& mgr,
                                                     const Sampledata& sample,
                                                     MFSUBMODELS::Transformation::Facade& facade,
                                                     double divisor)
    : my_config(config), my_mgr(mgr), 
      my_cov_param_begin(my_mgr.getParamBegin("Standardized covariates")),
      my_cov_param_end(my_mgr.getParamEnd("Standardized covariates")), 
      my_intercept_idx(mgr.getParameter("Intercept", "Intercept").getIndex()), my_sampledata(sample), 
      my_cov_field_begin(sample.getFieldBegin("Covariates")),
      my_main_phenotype_field(sample.getField("Core traits", "Main phenotype")),
      my_transf_facade(facade),
      my_divisor(divisor)
{
  my_phenotypic_means.resize((size_t)my_sampledata.getTotalIndividualCount(), QNAN);
  my_z_scores.resize((size_t)my_sampledata.getTotalIndividualCount(), QNAN);

  //cout << "\nmain trait" << endl;
  my_phenotypes.resize((size_t)my_sampledata.getTotalIndividualCount());
  
  //double  max = my_main_phenotype_field.getSummaryInfo().max;
  for(size_t i = 0; i < my_phenotypes.size(); ++i)
  {
    my_phenotypes[i] = my_main_phenotype_field.getAdjValue(i);
    //cout << my_phenotypes[i] << endl;
  }
}


// - Applicable only to binary case.
//
double  
MemberCovariateCalculator::getMeanResidualsAffected() const
{
  return  QNAN;
}

// - Applicable only to binary case.
//
double  
MemberCovariateCalculator::getMeanResidualsUnaffected() const
{
  return  QNAN;
}
        
void
MemberCovariateCalculator::calculatePhenotypicMeans()
{
  double  intercept = my_mgr(my_intercept_idx);
  for(size_t id = 0; id < (size_t)my_sampledata.getTotalIndividualCount(); ++id)
  {
    my_phenotypic_means[id] = intercept;     
  }      
  
  //cout << "intercept  " << intercept << endl;

  MAXFUN::ParameterConstIterator coefficient = my_cov_param_begin;
  SAMPLING::FieldConstIterator   field       = my_cov_field_begin;

  for(; coefficient != my_cov_param_end; ++coefficient, ++field)
  {
    //cout << "coeficient  " << coefficient->getCurrentEstimate() << endl;
    size_t  individual_count = (size_t)my_sampledata.getTotalIndividualCount();
    for(size_t id = 0; id < individual_count; ++id)
    {
      if(my_sampledata.isValid(id))
      {
        my_phenotypic_means[id] += coefficient->getCurrentEstimate() * field->getAdjValue(id);
        //if(id < 10)
          //cout << "covariate " << field->getAdjValue(id) << "  phenotypic means  " << my_phenotypic_means[id] << endl;
      }
    }
  }
}

  
void 
MemberCovariateCalculator::dumpValues() const
{
  OUTPUT::Table t("member_covariate_calculator");
  
  t << OUTPUT::TableColumn("Pedigree ID")
    << OUTPUT::TableColumn("Individual ID")
    << OUTPUT::TableColumn("Residual");

  for(size_t id = 0; id < (size_t)my_sampledata.getTotalIndividualCount(); ++id)
  {
    if(my_sampledata.isValid(id))
    {
      const FPED::Member&  member = my_sampledata.getIndividual(id);
    
      t << (OUTPUT::TableRow() << member.pedigree()->name() << member.name() << my_z_scores[id]);
    }
  }
  
  cout << t;
}


//======================================================================
// IMPLEMENTATION:  MccContinuousBoth       
//======================================================================
//
MccContinuousBoth::MccContinuousBoth(const Configuration& config,
                                     const MAXFUN::ParameterMgr& mgr,
                                     const Sampledata& sample,
                                     MFSUBMODELS::Transformation::Facade& facade, double divisor)
      : MemberCovariateCalculator::MemberCovariateCalculator(config, mgr, sample, facade, divisor)
{
  my_transformed_phenotypes.resize((size_t)my_sampledata.getTotalIndividualCount(), QNAN);
  my_transformed_phenotypic_means.resize((size_t)my_sampledata.getTotalIndividualCount(), QNAN);
  my_diffs.resize((size_t)my_sampledata.getTotalIndividualCount(), QNAN);
}
 
int
MccContinuousBoth::update(MAXFUN::ParameterMgr&)
{
  calculatePhenotypicMeans();
  calculateDiffs();
  
  // Transform phenotypes and phenotypic means and makes sure the transformation works:
  if(! doTransformations())
    return 1;
 
  for(size_t id = 0; id < my_sampledata.getTotalIndividualCount(); id++)
    calculateZScore(id);
 
  return  0;
}

void
MccContinuousBoth::calculateDiffs()
{
  for(size_t id = 0; id < (size_t)my_sampledata.getTotalIndividualCount(); ++id)
  {
    if(my_sampledata.isValid(id))
    {
      my_diffs[id] = my_phenotypes[id] - my_phenotypic_means[id];
    }
  }
}
 
bool
MccContinuousBoth::doTransformations()
{
  my_transformed_phenotypes = my_phenotypes;
  my_transformed_phenotypic_means = my_phenotypic_means;
 
  if(my_config.getTransConfig().get_type() != MFSUBMODELS::Transformation::NONE)
  {
    if(! my_transf_facade.calculate_geometric_mean(my_phenotypes))
      return  false;

    for(size_t i = 0; i < my_transformed_phenotypes.size(); ++i)
    {
      if(my_sampledata.isValid(i))
      {
        if(! my_transf_facade.transform(my_transformed_phenotypes[i]))
          return false;

        if(! my_transf_facade.transform(my_transformed_phenotypic_means[i]))
          return false;
      }
    }
  }

  return true;
}
 
void
MccContinuousBoth::calculateZScore(int id)
{
  // tr_mean = h(beta_0 + beta_1*x_1 + beta_2*x_2 + ...) --> h = transformation function
  double  tr_mean = my_transformed_phenotypic_means[(size_t)id];
      
  // tr_y = h(Y) --> h = transformation function
  double tr_y = my_transformed_phenotypes[(size_t)id];

  my_z_scores[(size_t)id] = tr_y - tr_mean;
}

Residuals
MccContinuousBoth::getResiduals() const
{
  
  return  Residuals(my_diffs, my_divisor);
}

//======================================================================
// IMPLEMENTATION:  MccContinuousDiff        
//======================================================================
//
MccContinuousDiff::MccContinuousDiff(const Configuration& config,
                                     const MAXFUN::ParameterMgr& mgr,
                                     const Sampledata& sample,
                                     MFSUBMODELS::Transformation::Facade& facade, double divisor)
      : MemberCovariateCalculator::MemberCovariateCalculator(config, mgr, sample, facade, divisor)
{
  my_diffs.resize((size_t)my_sampledata.getTotalIndividualCount(), QNAN);
  my_transformed_diffs.resize((size_t)my_sampledata.getTotalIndividualCount(), QNAN);
}

int
MccContinuousDiff::update(MAXFUN::ParameterMgr&)
{
  calculatePhenotypicMeans();
  calculateDiffs();
  
  if(! doTransformations())
  {
    return  1;
  }

  for(size_t id = 0; id < my_sampledata.getTotalIndividualCount(); id++)
  {
    calculateZScore(id);
    //cout << "z score  " << my_z_scores[id] << endl;
  }

  return  0;
}

void
MccContinuousDiff::calculateDiffs()
{
  for(size_t i = 0; i < my_phenotypes.size(); ++i)
  {
    if(my_sampledata.isValid(i))
    {
      my_diffs[i] = (my_phenotypes[i] - my_phenotypic_means[i]) / my_divisor;  
    }
    
    //cout << my_diffs[i] << endl;
  }
}

bool
MccContinuousDiff::doTransformations()
{
  my_transformed_diffs = my_diffs;

  if(my_config.getTransConfig().get_type() != MFSUBMODELS::Transformation::NONE)
  {
    /* Don't exclude geometric mean from transformation - RCE, March 2012.
    if(dynamic_cast<MccBinaryDiff*>(this))
    {
      my_transf_facade.calculate_geometric_mean(my_diffs, false);
      
      cout << "\n\nGEOMETRIC MEAN IS  " << my_transf_facade.geometric_mean() << endl;
    }
    */
    //else
    //{
      if(! my_transf_facade.calculate_geometric_mean(my_diffs))
      {
        return  false;
      }
    //}
    
    for(size_t i = 0; i < my_transformed_diffs.size(); ++i)
    {
      if(my_sampledata.isValid(i))
      {
        if(! my_transf_facade.transform(my_transformed_diffs[i]))
          return  false;
        //cout << "transformed diff " << my_transformed_diffs[i] << endl;
      }
    }    
  }

  return true;
}

void
MccContinuousDiff::calculateZScore(int id)
{
  my_z_scores[(size_t)id] = my_transformed_diffs[(size_t)id];
}

Residuals
MccContinuousDiff::getResiduals() const
{
  return  Residuals(my_diffs, my_divisor);
}

//======================================================================
// IMPLEMENTATION:  MccBinaryDiff        
//======================================================================
//
MccBinaryDiff::MccBinaryDiff(const Configuration& config,
                             const MAXFUN::ParameterMgr& mgr,
                             const Sampledata& sample,
                             MFSUBMODELS::Transformation::Facade& facade, double divisor)
      : MccContinuousDiff::MccContinuousDiff(config, mgr, sample, facade, divisor)  //, my_evaluation_count(0)
{
  my_binary_phenotypic_means.resize((size_t)my_sampledata.getTotalIndividualCount(), QNAN);
  my_new_phenotypes.resize((size_t)my_sampledata.getTotalIndividualCount(), QNAN);
  calculateNewPhenotypes();
}

void
MccBinaryDiff::calculateNewPhenotypes()
{
  double  one_over_two_n = .5 / my_sampledata.getValidIndividualCount(); 
  size_t  total_individual_count = (size_t)my_sampledata.getTotalIndividualCount();
  for(size_t id = 0; id < total_individual_count; ++id)
  {
    
    //if(id < 10)
      //cout << endl;
    if(my_sampledata.isValid(id))
    {
      my_new_phenotypes[id] = (my_phenotypes[id] + one_over_two_n);
      
      //if(id < 10)
        //cout << "phenotype  " << my_phenotypes[id] << "  new phenotype   " << my_new_phenotypes[id] << endl;
    }
  }
}

void
MccBinaryDiff::calculateDiffs()
{
  calculateBinaryPhenotypicMeans();
  
  //  cout << "\n my divisor " << my_divisor << endl;

  size_t  individual_count = (size_t)my_sampledata.getTotalIndividualCount();
  for(size_t id = 0; id < individual_count; ++id)
  {
    if(my_sampledata.isValid(id))
    {
      
      //double  divisor = 1.0; 
            double  divisor = sqrt(my_binary_phenotypic_means[id] * (1 - my_binary_phenotypic_means[id]));
            /*
            double  intercept = my_mgr(my_intercept_idx);
            double  p = exp(intercept) / (1 + exp(intercept));
            double  divisor = sqrt(p * (1-p));    
            */
      
            /*
            double  trait_mean = my_sampledata.getField("Core traits", "Main phenotype").getMean();
            
            trait_mean += .5 / my_sampledata.getValidIndividualCount();
            double  divisor = sqrt(trait_mean * (1 - trait_mean));
            */
      
      my_diffs[id] = (my_new_phenotypes[id] - my_binary_phenotypic_means[id]) / divisor / my_divisor;
      
      if(my_phenotypes[id] == 1.0)
      {
        my_affected_diffs.push_back(my_diffs[id] * divisor * my_divisor);
      }
      
      if(my_phenotypes[id] == 0.0)
      {
        my_unaffected_diffs.push_back(my_diffs[id] * divisor * my_divisor);
      }
    }
  }
}


void
MccBinaryDiff::calculateBinaryPhenotypicMeans()
{
  size_t  individual_count = (size_t)my_sampledata.getTotalIndividualCount();
  for(size_t id = 0; id < individual_count; ++id)
  {
    if(my_sampledata.isValid(id))
    {
      my_binary_phenotypic_means[id] = exp(my_phenotypic_means[id]) / (1 + exp(my_phenotypic_means[id]));
      
     //if(id < 10)
        //cout << "phenotypic mean   " << my_phenotypic_means[id] << "  " << "binary phenotypic mean   " << my_binary_phenotypic_means[id] << endl;      
    }
  }
}

double  
MccBinaryDiff::getMeanResidualsAffected() const
{
  return  Residuals(my_affected_diffs, 1.0).mean();
}

double  
MccBinaryDiff::getMeanResidualsUnaffected() const
{
  return  Residuals(my_unaffected_diffs, 1.0).mean();
}


} // End namespace ASSOC
} // End namespace SAGE
