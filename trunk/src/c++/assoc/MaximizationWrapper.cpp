//=============================================================================
//
//  File:  MaximizationWrapper.cpp
//
//  Author:  Stephen Gross
//
//  Copyright 2002 R. C. Elston
//=============================================================================

#include "boost/lambda/bind.hpp"
#include "assoc/MaximizationWrapper.h"

namespace SAGE  {
namespace ASSOC {

const double  TOTAL_BETA_CHANGE_THRESHOLD = .00001;
const double  LH_CHANGE_THRESHOLD         = .000001;
const double  INITIAL_PREVALENCE = .9;


//==========================================================================
// IMPLEMENTATION:  vc_order
// Note:  Trac #1780.  Defines order for removal of variance components
//        fixed at lower bound by MAXFUN.
//==========================================================================
//
vc_order::vc_order()
{
  my_order["Sibling"]   = 0;
  my_order["Marital"]   = 1;
  my_order["Family"]    = 2;
  my_order["Polygenic"] = 3;
}

bool
vc_order::operator()(const string& left, const string& right)
{
  assert(my_order.count(left) == 1 && my_order.count(right) == 1);
  
  return  my_order[left] < my_order[right];
}


//==========================================================================
//  maximize()
//==========================================================================
void 
MaximizationWrapper::maximize(const string& model_name, Configuration& config,
                              const FPED::Multipedigree& fped, const Sampledata& sampledata,
                              AnalysisResults::ModelResults& results, ostream& messages,
                              bool base_initialized, cerrorstream errors)
{
  MaximizationWrapper(model_name, config, fped, sampledata, results, messages, base_initialized, errors);
}

//==========================================================================
//  MaximizationWrapper() 
//==========================================================================
MaximizationWrapper::MaximizationWrapper(const string& model_name, Configuration& config,
                                         const FPED::Multipedigree& fped, const Sampledata& sampledata,
                                         AnalysisResults::ModelResults& results, 
                                         ostream& messages, bool base_initialized, cerrorstream errors)
      : my_model_name(model_name), my_config(config), my_messages(messages), null_results_available(base_initialized),
        my_sampledata(sampledata), my_fped(fped), my_errors(errors), alt_covariates_included(false), score_test(false),
        my_divisor(1.0)
{
  const TransConfig&  trans_config = my_config.getTransConfig();

  if(! base_initialized)
  {
    cout << "    Maximizing model without test covariates..." << endl;  
  
    // - Do model maximization in two passes -- first w/o estimating
    //   transformation parameters, then with.
    //
    if(trans_config.either_lambda_estimated())
    {
      TransConfig   first_pass_trans_config(trans_config);
      
      fixToNoTrans(first_pass_trans_config);

      my_phase = ONE_OF_TWO;
      
      doMaximization(results, first_pass_trans_config);
      
      null_results_available = results.maxfun_results_null.getConverged();
      if(! null_results_available)
      {
        throw(NonConvergenceWithTransformation(my_model_name));
      }
      
      my_phase = TWO_OF_TWO;
      
      doMaximization(results, trans_config);
      
      null_results_available = results.maxfun_results_null.getConverged();
    }
    else
    {
      my_phase = ONE_OF_ONE;
      doMaximization(results, trans_config);
      null_results_available = results.maxfun_results_null.getConverged();
    }
  }

  if(my_config.getAltCovCount(my_model_name))
  {
    alt_covariates_included = true;
    
    cout << "    Maximizing model with test covariates..." << endl;  
    my_phase = ONE_OF_ONE;
    doMaximization(results, trans_config);   
  }
  
  results.residuals_scaling_factor = my_divisor;
}

//==========================================================================
//  doMaximization()
//==========================================================================
void
MaximizationWrapper::doMaximization(AnalysisResults::ModelResults& results, 
                                    const TransConfig& trans_config)
{
  string  extra_variance_component = "dummy";
  while(! extra_variance_component.empty())
  {     
    MAXFUN::ParameterMgr mgr;
    MAXFUN::Function func(mgr);
     
    addParams(results, func);
    
    // Create a Transformation::Facade bound to other MAXFUN::Function.
    MFSUBMODELS::Transformation::Facade  trans_facade(trans_config, func);
    MemberCovariateCalculator*  mcc;
    
    if(my_config.getDependentTraitType() == BINARY)
    {
      mcc = new MccBinaryDiff(my_config, mgr, my_sampledata, trans_facade, my_divisor);
    }
    else
    {
      if(my_config.getTransformBothSides())
      {
        mcc = new MccContinuousBoth(my_config, mgr, my_sampledata, trans_facade);
      }
      else
      {
        mcc = new MccContinuousDiff(my_config, mgr, my_sampledata, trans_facade, my_divisor);
      }
    }
    
    Calculator  calc(my_messages, my_config, mgr, my_fped, 
                    my_sampledata, *mcc, trans_facade, my_model_name, alt_covariates_included);
    
    establishInitialEstimates(mgr, INITIAL_PREVALENCE, results);
    
    // Add the necessary steps to the function from the calculator
    func.addStep(boost::bind(&MemberCovariateCalculator::update, boost::ref(calc.getMcc()), _1), "");
    func.setEvaluator(boost::bind(&Calculator::calculateLh, boost::ref(calc), _1));

    MAXFUN::DebugCfg  dbg(my_config.getMaxfunDebug() ? MAXFUN::DebugCfg::COMPLETE : MAXFUN::DebugCfg::NO_DEBUG_INFO);

    if(dbg.hasAnyReportedOutput())
      dbg.setDebugOutput(my_config.getTitle() + ".max", true);

    MAXFUN::SequenceCfg  cfg(MAXFUN::SequenceCfg::USER_DEFINED, 
      my_model_name + (alt_covariates_included ? " with" : " without") + " test covariates", 
      "ln likelihood");
      
    setMaximizationSequence(cfg);
    
    // Maximize
    *(alt_covariates_included ? &results.maxfun_results_alt : &results.maxfun_results_null) = 
                                                                    MAXFUN::Maximizer::maximize(func, cfg, dbg);
    checkResults(*(alt_covariates_included ? &results.maxfun_results_alt : &results.maxfun_results_null));                                                                     
  
    if(my_phase != ONE_OF_TWO)
    {
      const Configuration::Model& current_model = my_config.getModel(my_model_name);
      if(! alt_covariates_included && current_model.null_residuals)
      {
        results.residuals_null = calc.getMcc().getResiduals();
        results.ind_residuals_null = calculateIndResiduals(mgr, calc);
      }
      else if(alt_covariates_included && current_model.test_residuals)
      {
        results.residuals_alt = calc.getMcc().getResiduals();
        results.ind_residuals_alt = calculateIndResiduals(mgr, calc);
      }                                                                      
                                                                      
      if(my_config.getDependentTraitType() == BINARY) 
      {
        *(alt_covariates_included ? &results.affected_mean_residuals_alt : &results.affected_mean_residuals_null) = 
                                                                      calc.getMcc().getMeanResidualsAffected();    
        *(alt_covariates_included ? &results.unaffected_mean_residuals_alt : &results.unaffected_mean_residuals_null) = 
                                                                      calc.getMcc().getMeanResidualsUnaffected(); 
      }  
    }    
    if(my_phase == ONE_OF_TWO) 
    {
      if(my_config.getDependentTraitType() == QUANTITATIVE && ! my_config.getTransformBothSides())
      {
        my_divisor = calc.getMcc().getResiduals().std_dev();
      }
      else if(my_config.getDependentTraitType() == BINARY)
      {
        //my_divisor = calc.getMcc().getResiduals().std_dev();   
        my_divisor = BINARY_DIVISOR;
      }
      
      //cout << "My divisor " << my_divisor << endl; 
    }

    delete  mcc;    
    
    // - Trac #1780
    //
    extra_variance_component = checkVarianceComponents(mgr);
    if(! extra_variance_component.empty())
    {
      removeVarianceComponent(extra_variance_component);
    }
    else
    {
      checkRandomVariance(mgr);
    }
  }
}

Residuals
MaximizationWrapper::calculateIndResiduals(const MAXFUN::ParameterMgr& mgr, const Calculator& calc) const
{
  const map<string, FortranMatrix<double> >&  shared_effects = calc.getSharedEffects();
  const vector<double>&  diffs = calc.getMcc().getDiffs();
  const vector<FPED::MemberConstPointer>&  member_lookup = calc.getMemberLookup();
  const map<FPED::MemberConstPointer, size_t>&  index_lookup = calc.getIndexLookup();
  size_t  member_count = member_lookup.size();
  
  FortranMatrix<double>  var_cov;
  populateVarCovMatrix(var_cov, shared_effects, mgr, member_count);
  
  FortranMatrix<double>  var_cov_inverse;
  SVD  matrix_inverter;
  matrix_inverter.inverse_of(var_cov, var_cov_inverse);
  
  FortranMatrix<double>  exp_values;
  populateExpValues(exp_values, diffs, index_lookup, member_count);
  
  FortranMatrix<double>  phenotypes;
  populatePhenotypes(phenotypes, index_lookup, member_count);
  
  FortranMatrix<double>  ind_residuals;
  ind_residuals.resize_fill(member_count, 1, 0.0);
  
  FortranMatrix<double>  accumulator;
  accumulator.resize_fill(member_count, 1, 0.0);
  
  map<string, FortranMatrix<double> >::const_iterator  eff_iter = shared_effects.begin();
  map<string, FortranMatrix<double> >::const_iterator  eff_end_iter = shared_effects.end();
  for(; eff_iter != eff_end_iter; ++eff_iter)
  {
    string  effect_name = eff_iter->first;
    double  effect_value = mgr.getParameter("Variance components", effect_name).getFinalEstimate();
    
    accumulator += eff_iter->second * effect_value * var_cov_inverse * exp_values;
  }  
  
  ind_residuals = phenotypes - accumulator;
  
  //print_matrix(phenotypes, cout);
  //print_matrix(accumulator, cout);
  //print_matrix(ind_residuals, cout); 
  
  vector<double>  ind_residuals_vector;
  
  reindexResiduals(ind_residuals_vector, ind_residuals, member_lookup);
  
  /*
  for(size_t i = 0; i < ind_residuals_vector.size(); ++i)
    cout << "mped index  " << i << "residual  " << ind_residuals_vector[i] << endl;
  */
  
  return  Residuals(ind_residuals_vector, 1.0);
}


// - Extract residuals from a matrix and store them in a vector indexed by mped
//   index.
//
void
MaximizationWrapper::reindexResiduals(vector<double>& ind_residuals_vector,
                                      const FortranMatrix<double>& ind_residuals,
                                      const vector<FPED::MemberConstPointer>& member_lookup) const
{
  ind_residuals_vector.resize(my_fped.member_count(), QNAN);
  
  size_t  valid_member_count = member_lookup.size();
  for(size_t i = 0; i < valid_member_count; ++i)
  {
    double  residual = ind_residuals(i, 0);
    size_t  mp_index = member_lookup[i]->mpindex();
    
    ind_residuals_vector[mp_index] = residual;
  }
}


void 
MaximizationWrapper::populatePhenotypes(FortranMatrix<double>& phenotypes, 
                                   const map<FPED::MemberConstPointer, size_t>& index_lookup, size_t member_count) const
{
  phenotypes.resize_fill(member_count, 1, 0.0);

  const SAMPLING::Field&  phenotype_field = my_sampledata.getField("Core traits", "Main phenotype");

  size_t  total_count = my_sampledata.getTotalIndividualCount();
  for(size_t i = 0; i < total_count; ++i)
  {
    if(my_sampledata.isValid(i))
    {
      FPED::MemberConstPointer  member_ptr = &(my_fped.member_index(i));
      map<FPED::MemberConstPointer, size_t>::const_iterator result  = index_lookup.find(member_ptr);
      
      assert(result != index_lookup.end());
      
      size_t  sample_index = result->second;
      
      phenotypes(sample_index, 0) = phenotype_field.getAdjValue(i);
    }
  }
  
  //print_matrix(phenotypes, cout);  
}                                   
                                   


void 
MaximizationWrapper::populateExpValues(FortranMatrix<double>& exp_values, const vector<double>& diffs, 
                              const  map<FPED::MemberConstPointer, size_t>& index_lookup, size_t member_count) const
{
  exp_values.resize_fill(member_count, 1, 0.0);

  size_t  diff_count = diffs.size();
  for(size_t i = 0; i < diff_count; ++i)
  {
    if(my_sampledata.isValid(i))
    {
      FPED::MemberConstPointer  member_ptr = &(my_fped.member_index(i));
      map<FPED::MemberConstPointer, size_t>::const_iterator result  = index_lookup.find(member_ptr);
      
      assert(result != index_lookup.end());
      
      size_t  sample_index = result->second;
      
      exp_values(sample_index, 0) = diffs[i] * my_divisor;
    }
  }
  
  //print_matrix(exp_values, cout);
}                              
                              


void  
MaximizationWrapper::populateVarCovMatrix(FortranMatrix<double>& var_cov, 
                                          const map<string, FortranMatrix<double> >& shared_effects,
                                          const MAXFUN::ParameterMgr& mgr, size_t member_count      ) const
{
  var_cov.resize_fill(member_count, member_count, 0.0);
  
  map<string, FortranMatrix<double> >::const_iterator  eff_iter = shared_effects.begin();
  map<string, FortranMatrix<double> >::const_iterator  eff_end_iter = shared_effects.end();
  for(; eff_iter != eff_end_iter; ++eff_iter)
  {
    string  effect_name = eff_iter->first;
    double  effect_value = mgr.getParameter("Variance components", effect_name).getFinalEstimate();
    var_cov += eff_iter->second * effect_value;
  }
  
  // - Set every diagonal element to the total variance.
  //
  double  total_variance = mgr.getParameter("Other parameters", "Total variance").getFinalEstimate();
  for(size_t i = 0; i < member_count; ++i)
    for(size_t j = 0; j < member_count; ++j)
      if(i == j)
        var_cov(i, j) = total_variance;
  
  //print_matrix(var_cov, cout);
}
                                   


void
MaximizationWrapper::doScoreTest(AnalysisResults::ModelResults& results, 
                                          const TransConfig& trans_config)
{
  results.score_test_parameter_mgr = results.maxfun_results_null.getParameterMgr();
  MAXFUN::Function  func(results.score_test_parameter_mgr);
   
  addTestCovCoefficients(func);
  
  // Create a Transformation::Facade bound to other MAXFUN::Function.
  MFSUBMODELS::Transformation::Facade  trans_facade(trans_config, func);
  MemberCovariateCalculator*  mcc;
  
  if(my_config.getDependentTraitType() == BINARY)
  {
    mcc = new MccBinaryDiff(my_config, results.score_test_parameter_mgr, my_sampledata, trans_facade);
  }
  else
  {
    if(my_config.getTransformBothSides())
    {
      mcc = new MccContinuousBoth(my_config, results.score_test_parameter_mgr, my_sampledata, trans_facade);
    }
    else
    {
      mcc = new MccContinuousDiff(my_config, results.score_test_parameter_mgr, my_sampledata, trans_facade);
    }
  }
  
  Calculator  calc(my_messages, my_config, results.score_test_parameter_mgr, my_fped, 
                  my_sampledata, *mcc, trans_facade, my_model_name, alt_covariates_included);
  
  // Add the necessary steps to the function from the calculator
  func.addStep(boost::bind(&MemberCovariateCalculator::update, boost::ref(calc.getMcc()), _1), "");
  func.setEvaluator(boost::bind(&Calculator::calculateLh, boost::ref(calc), _1));

  // Calculate score test.
  MAXFUN::DebugCfg  dbg;
  
  cout << "\ncalling calculateScoreTestStatistic() ..." << endl;
  
  results.score_statistic = MAXFUN::Maximizer::calculateScoreTestStatistic(func, dbg);
  
  delete  mcc;
}

void  
MaximizationWrapper::setMaximizationSequence(MAXFUN::SequenceCfg& cfg) const
{
  cfg.addRunCfg(MAXFUN::RunCfg::DIRECT_WITHOUT, 1);
  cfg.getLatestRunCfg().epsilon1 = 1e-2;   
  cfg.getLatestRunCfg().epsilon2 = 1e-7;  
                                                           
  cfg.addRunCfg(MAXFUN::RunCfg::VAR_METRIC_IDENTITY, 20);
  cfg.getLatestRunCfg().epsilon1       = 1e-3;   
  cfg.getLatestRunCfg().epsilon2       = 1e-7;  
  cfg.getLatestRunCfg().var_cov        = MAXFUN::RunCfg::FINAL;    
  cfg.getLatestRunCfg().control_option = MAXFUN::RunCfg::PREVIOUS_NONCONVERGENCE;
                                                                                                                                                                                                       
  cfg.addRunCfg(MAXFUN::RunCfg::DIRECT_WITHOUT, 50);
  cfg.getLatestRunCfg().epsilon1       = 1e-4;
  cfg.getLatestRunCfg().epsilon2       = 1e-7;
  cfg.getLatestRunCfg().control_option = MAXFUN::RunCfg::PREVIOUS_NONCONVERGENCE;
                                                                                                                                                                                                                     
  cfg.addRunCfg(MAXFUN::RunCfg::VAR_METRIC_ESTIMATE, 20);
  cfg.getLatestRunCfg().epsilon1       = 1e-5;  
  cfg.getLatestRunCfg().epsilon2       = 1e-8;  
  cfg.getLatestRunCfg().var_cov        = MAXFUN::RunCfg::FINAL;
                                
  cfg.addRunCfg(MAXFUN::RunCfg::VAR_METRIC_ESTIMATE, 20);
  cfg.getLatestRunCfg().epsilon1       = 1e-6;   
  cfg.getLatestRunCfg().epsilon2       = 1e-8;  
  cfg.getLatestRunCfg().var_cov        = MAXFUN::RunCfg::FINAL;
  cfg.getLatestRunCfg().control_option = MAXFUN::RunCfg::PREVIOUS_NONCONVERGENCE;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
  cfg.addRunCfg(MAXFUN::RunCfg::VAR_METRIC_ESTIMATE, 20);
  cfg.getLatestRunCfg().epsilon1       = 1e-6;   
  cfg.getLatestRunCfg().epsilon2       = 1e-10;  
  cfg.getLatestRunCfg().var_cov        = MAXFUN::RunCfg::FINAL;
  cfg.getLatestRunCfg().control_option = MAXFUN::RunCfg::PREVIOUS_NONCONVERGENCE;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
  cfg.addRunCfg(MAXFUN::RunCfg::VAR_METRIC_ESTIMATE, 20);
  cfg.getLatestRunCfg().epsilon1       = 1e-4;   
  cfg.getLatestRunCfg().epsilon2       = 1e-8;  
  cfg.getLatestRunCfg().var_cov        = MAXFUN::RunCfg::FINAL;
  cfg.getLatestRunCfg().control_option = MAXFUN::RunCfg::PREVIOUS_NONCONVERGENCE;
}


// - Fix unfixed transformation parameters to values which correspond to
//   no transformation.
//
void
MaximizationWrapper::fixToNoTrans(TransConfig& trans_config) const
{
  if(trans_config.both_lambdas_estimated())
  {
    trans_config.set_type(MFSUBMODELS::Transformation::NONE);
  }
  else
  {
    switch(trans_config.get_type())
    {
      case MFSUBMODELS::Transformation::GEORGE_ELSTON:
        if(trans_config.get_lambda1().initial_type != MAXFUN::Parameter::FIXED)
        {
          trans_config.set_lambda1_init_est(1.0);
          trans_config.set_lambda1_fixed(true);
        }
        
        if(trans_config.get_lambda2().initial_type != MAXFUN::Parameter::FIXED)
        {
          trans_config.set_lambda2_init_est(0.0);
          trans_config.set_lambda2_fixed(true);
        }
        
        break;
        
      default:
        assert(false);
    }
  }
}


void
MaximizationWrapper::removeVarianceComponent(const string& component)
{
  string  maximization_desc = my_model_name == "Baseline" ? "model 'Baseline'" :
                                  string(alt_covariates_included ? "test" : "non-test") + " maximization " +
                                  "of model '" + my_model_name + "'"; 
    
  my_errors << priority(information) << "Variance of " << component << " effect converged to "
            << VARIANCE_COMPONENT_LB << " for " << maximization_desc << ".  Remaximizing without " 
            << component << " effect ..." << endl;
            
  my_config.removeEffectFromModel(component, my_model_name, alt_covariates_included);
  
  // - Remove variance component from test model also if one exists Trak #1845.
  //
  if(! alt_covariates_included)
  {
    if(my_config.getAltCovCount(my_model_name))
    {
      my_config.removeEffectFromModel(component, my_model_name, true);
    }
    
    if(my_config.getAllowAveraging())    // Remove the variance component from all models Trak #1846.
    {
      const vector<string>&  model_names = my_config.getModelList();
      
      vector<string>::const_iterator  model_iter     = model_names.begin();
      vector<string>::const_iterator  model_end_iter = model_names.end();
      for(; model_iter != model_end_iter; ++model_iter)
      {
        if(my_config.modelHasFixedEffect(*model_iter, false, string2effect(component)))
        {
          my_config.removeEffectFromModel(component, *model_iter, false);
        }
        
        if(my_config.modelHasFixedEffect(*model_iter, true, string2effect(component)))
        {
          my_config.removeEffectFromModel(component, *model_iter, true);
        }
      }
    }
  }
}


// - Previously, before MAXFUN::Maximizer maximized likelihood it evaluated the likelihood
//   function once.  If the resultant value was not finite, it silently gave up
//   and returned the invalid results as they were.  This function checks for that
//   condition.  Added 10/25/7 -djb
//
void
MaximizationWrapper::checkResults(const MAXFUN::Results& results) const
{
  if(! results.getValidInitialValue())
  {
    throw BadLikelihood(my_model_name);
  }
}

// - Trac #1780.  Are any variance components fixed to a lower bound?  If so, issue
//   a warning and return one. 
//
string
MaximizationWrapper::checkVarianceComponents(const MAXFUN::ParameterMgr& mgr)
{
  string  special_component = checkSpecialVarianceComponents(mgr);
  
  if(! special_component.empty())
  {
    return  special_component;
  }
  
  set<string, vc_order>  fixed_vcs;

  MAXFUN::ParameterConstIterator  vc_iter     = mgr.getParamBegin("Variance components");
  MAXFUN::ParameterConstIterator  vc_end_iter = mgr.getParamEnd("Variance components");
  for(; vc_iter != vc_end_iter; ++vc_iter)
  {
    MAXFUN::Parameter::ParamTypeEnum  final_status = vc_iter->getFinalType();
    //cout << "\n" << vc_iter->getName() << "  " 
    //     << MAXFUN::ParamTypeEnum2str(final_status) << endl;
    if(final_status == MAXFUN::Parameter::IND_FUNC_FIXED_AT_BOUND ||
       final_status == MAXFUN::Parameter::IND_FIXED_AT_BOUND        )
    {
      string  name = vc_iter->getName();
      if(name != "Random")
      {
        fixed_vcs.insert(name);
      }
    }
  }
  
  return  fixed_vcs.empty() ? "" : *(fixed_vcs.begin());  
}  


string
MaximizationWrapper::checkSpecialVarianceComponents(const MAXFUN::ParameterMgr& mgr)
{
  MAXFUN::ParameterConstIterator  vc_iter     = mgr.getParamBegin("Variance components");
  MAXFUN::ParameterConstIterator  vc_end_iter = mgr.getParamEnd("Variance components");
  for(; vc_iter != vc_end_iter; ++vc_iter)
  {
    string  param_name = vc_iter->getName();
    if(! (fixed_effect(param_name) || param_name == "0.5 * Polygenic"))
    {
      MAXFUN::Parameter::ParamTypeEnum  final_status = vc_iter->getFinalType();
     
      if(final_status == MAXFUN::Parameter::IND_FUNC_FIXED_AT_BOUND ||
         final_status == MAXFUN::Parameter::IND_FIXED_AT_BOUND        )
      {
        return  vc_iter->getName();
      }
    }
  }
  
  return  "";
}

// - Trac #1780.  If the random variance component is fixed to a lower bound, issue a
//   special warning and leave it in the model.
//
void
MaximizationWrapper::checkRandomVariance(const MAXFUN::ParameterMgr& mgr)
{
  MAXFUN::Parameter::ParamTypeEnum  final_status = 
              mgr.getParameter("Variance components", "Random").getFinalType();
              
  if(final_status == MAXFUN::Parameter::IND_FUNC_FIXED_AT_BOUND ||
     final_status == MAXFUN::Parameter::IND_FIXED_AT_BOUND        )
  {
     my_errors << priority(warning) << "Variance of the individual environmental "
               << "random effect has converged to " << VARIANCE_COMPONENT_LB << ".  This should never "
               << "happen.  Fix at least one of the other random effects being "
               << "estimated to zero." << endl;
  }
}
                                                  

//==========================================================================
//  addParams()
//==========================================================================
void
MaximizationWrapper::addParams(AnalysisResults::ModelResults& results, MAXFUN::Function& func)
{
  vector<size_t>  fixed_effect_indices(5, (size_t)(-1));
  vector<bool>  fixed_effect_usage;
  size_t  tv_idx = (size_t)(-1);
  
  fixed_effect_usage.push_back(true);    // random effect
  fixed_effect_usage.push_back(my_config.modelHasFixedEffect(my_model_name, 
                                                             alt_covariates_included, Configuration::POLYGENIC));
  fixed_effect_usage.push_back(my_config.modelHasFixedEffect(my_model_name, 
                                                             alt_covariates_included, Configuration::FAMILY));
  fixed_effect_usage.push_back(my_config.modelHasFixedEffect(my_model_name, 
                                                             alt_covariates_included, Configuration::SIBLING));
  fixed_effect_usage.push_back(my_config.modelHasFixedEffect(my_model_name, 
                                                             alt_covariates_included, Configuration::MARITAL));                                                             

  func.getMgr().reset();
  
  addFixedEffects(func, fixed_effect_indices, fixed_effect_usage);
  addUserEffects(func);
  addTotalVariance(results, func, fixed_effect_indices, fixed_effect_usage, tv_idx);
  addHeritability(func, fixed_effect_indices, fixed_effect_usage, tv_idx); 
  addVarianceProportions(func, fixed_effect_indices, fixed_effect_usage, tv_idx);
  addFamilialCorrelations(func, fixed_effect_indices, fixed_effect_usage, tv_idx);
  addEnvironmentalCorrelations(func, fixed_effect_indices, fixed_effect_usage, tv_idx);
  addIntercept(func);
  addNullCovCoefficients(func);

  if(alt_covariates_included)
  {
    addTestCovCoefficients(func);
  }
}


void
MaximizationWrapper::addVarianceProportions(MAXFUN::Function& func,
                                            const vector<size_t>& fixed_effect_indices,
                                            const vector<bool>& fixed_effect_usage, size_t tv_idx) const
{
  MAXFUN::ParameterMgr&  mgr = func.getMgr();

  mgr.addGroup("Variance proportions");
  
  MAXFUN::Parameter&  var_prop_random = 
                    mgr.addParameterAlt("Variance proportions", "Random", MAXFUN::Parameter::DEPENDENT, 0.5, 0.0, 1.0);
  var_prop_random.setIncludeInOutput(my_config.getDependentTraitType() == BINARY);

  func.addDependentParamStep(var_prop_random.getIndex(), 
                             MF_PARAM(fixed_effect_indices[Configuration::RANDOM]) / MF_PARAM(tv_idx));

  if(fixed_effect_usage[Configuration::POLYGENIC])
  {
    MAXFUN::Parameter&  var_prop_poly = 
                    mgr.addParameterAlt("Variance proportions", "Polygenic", MAXFUN::Parameter::DEPENDENT);
    var_prop_poly.setPValueInf(MAXFUN::Parameter::ONE_SIDED);
    var_prop_poly.setIncludeInOutput(my_config.getDependentTraitType() == BINARY);

    func.addDependentParamStep(var_prop_poly.getIndex(), 
                               MF_PARAM(fixed_effect_indices[Configuration::POLYGENIC]) / MF_PARAM(tv_idx));
  }

  if(fixed_effect_usage[Configuration::FAMILY])
  {
    MAXFUN::Parameter&  var_prop_family = 
                    mgr.addParameterAlt("Variance proportions", "Family", MAXFUN::Parameter::DEPENDENT);
    var_prop_family.setPValueInf(MAXFUN::Parameter::ONE_SIDED);
    var_prop_family.setIncludeInOutput(my_config.getDependentTraitType() == BINARY);

    func.addDependentParamStep(var_prop_family.getIndex(), 
                               MF_PARAM(fixed_effect_indices[Configuration::FAMILY]) / MF_PARAM(tv_idx));
  }

  if(fixed_effect_usage[Configuration::MARITAL])
  {
    MAXFUN::Parameter&  var_prop_marital = 
                    mgr.addParameterAlt("Variance proportions", "Marital", MAXFUN::Parameter::DEPENDENT);
    var_prop_marital.setPValueInf(MAXFUN::Parameter::ONE_SIDED);
    var_prop_marital.setIncludeInOutput(my_config.getDependentTraitType() == BINARY);

    func.addDependentParamStep(var_prop_marital.getIndex(), 
                               MF_PARAM(fixed_effect_indices[Configuration::MARITAL]) / MF_PARAM(tv_idx));
  }

  if(fixed_effect_usage[Configuration::SIBLING])
  {
    MAXFUN::Parameter&  var_prop_sib = 
                    mgr.addParameterAlt("Variance proportions", "Sibling", MAXFUN::Parameter::DEPENDENT);
    var_prop_sib.setPValueInf(MAXFUN::Parameter::ONE_SIDED);
    var_prop_sib.setIncludeInOutput(my_config.getDependentTraitType() == BINARY);

    func.addDependentParamStep(var_prop_sib.getIndex(), MF_PARAM(fixed_effect_indices[Configuration::SIBLING]) / MF_PARAM(tv_idx));
  }
  
  MAXFUN::ParameterConstIterator  p_iter     = mgr.getParamBegin("Variance components");
  MAXFUN::ParameterConstIterator  p_end_iter = mgr.getParamEnd("Variance components");
  for(; p_iter != p_end_iter; ++p_iter)
  {
    string  param_name = p_iter->getName();
    
    if(! (fixed_effect(param_name) || param_name == "0.5 * Polygenic"))
    {
      if(my_config.modelHasEffect(my_model_name, alt_covariates_included, param_name))
      {
        MAXFUN::Parameter&  var_prop_user = 
                        mgr.addParameterAlt("Variance proportions", param_name, MAXFUN::Parameter::DEPENDENT);
        var_prop_user.setPValueInf(MAXFUN::Parameter::ONE_SIDED);
        var_prop_user.setIncludeInOutput(my_config.getDependentTraitType() == BINARY);

        size_t  param_index = p_iter->getIndex();
        func.addDependentParamStep(var_prop_user.getIndex(), 
                                   MF_PARAM(param_index) / MF_PARAM(tv_idx));
      }
    }
  }
}                                            


void
MaximizationWrapper::addHeritability(MAXFUN::Function& func, 
                                     const vector<size_t>& fixed_effect_indices, 
                                     const vector<bool>& fixed_effect_usage, size_t tv_idx) const
{
  MAXFUN::ParameterMgr&  mgr = func.getMgr();

  if(fixed_effect_usage[Configuration::POLYGENIC])
  {
    MAXFUN::Parameter&  heritability = 
                  mgr.addParameterAlt("Other parameters", "Heritability", MAXFUN::Parameter::DEPENDENT, 0.0, 0.0, 1.0);
    heritability.setPValueInf(MAXFUN::Parameter::ONE_SIDED);

    func.addDependentParamStep(heritability.getIndex(), 
                               MF_PARAM(fixed_effect_indices[Configuration::POLYGENIC]) / MF_PARAM(tv_idx));
  }
}


// Ok, so you might ask what the heck is going on here. Here's the deal: if you were to express the calculation for
// total variance as an equation, it would look like this:
//
// total var = random + polygenic + sibling + marital + (only one generation ? 1 : 2) * family
//
// The problem is that, because of the way functors work, you can't simply specify this all as one composite functor.
// So instead, I did it with one logical step to separate it out:
//
// If family_eff enabled and more than one generation in sample:
//
//  unadjusted total var = random + polygenic + sibling + marital + family
//  total var            = unadjusted total var + family
//
// Otherwise:
//
//  total var = random + polygenic + sibling + marital
//
void
MaximizationWrapper::addTotalVariance(AnalysisResults::ModelResults& results, MAXFUN::Function& func,
                                      const vector<size_t>& fixed_effect_indices,
                                      const vector<bool>& fixed_effect_usage, size_t& tv_idx) const
{
  MAXFUN::ParameterMgr&  mgr = func.getMgr();
  
  mgr.addGroup("Other parameters");
  
    //if(my_phase != TWO_OF_TWO)
    //{ 
    mgr.addGroup("Other parameters");
    tv_idx = mgr.addParameter("Other parameters", "Total variance",
                                   MAXFUN::Parameter::DEPENDENT, 0.0, 0.0, MAXFUN::MF_INFINITY);
    
    func.addDependentParamStep(tv_idx, MF_PARAM(fixed_effect_indices[Configuration::RANDOM]));
    
    if(fixed_effect_usage[Configuration::SIBLING])
    {
      func.addDependentParamStep(tv_idx, MF_PARAM(tv_idx) + MF_PARAM(fixed_effect_indices[Configuration::SIBLING]));
    }
    
    if(fixed_effect_usage[Configuration::MARITAL])
    {
      func.addDependentParamStep(tv_idx, MF_PARAM(tv_idx) + MF_PARAM(fixed_effect_indices[Configuration::MARITAL]));
    }
    
    if(fixed_effect_usage[Configuration::POLYGENIC])
    { 
      func.addDependentParamStep(tv_idx, MF_PARAM(tv_idx) + MF_PARAM(fixed_effect_indices[Configuration::POLYGENIC]));
    }
    
    if(fixed_effect_usage[Configuration::FAMILY]) 
    {
      func.addDependentParamStep(tv_idx, MF_PARAM(tv_idx) + MF_PARAM(fixed_effect_indices[Configuration::FAMILY]));
    }
                                       
    // - Add user defined effects to total variance.
    //
    MAXFUN::ParameterConstIterator  p_iter     = mgr.getParamBegin("Variance components");
    MAXFUN::ParameterConstIterator  p_end_iter = mgr.getParamEnd("Variance components");
    for(; p_iter != p_end_iter; ++p_iter)
    {
      string  param_name = p_iter->getName();
      if(! (fixed_effect(param_name) || param_name == "0.5 * Polygenic"))
      {
        if(my_config.modelHasEffect(my_model_name, alt_covariates_included, param_name))
        {
          func.addDependentParamStep(tv_idx, MF_PARAM(tv_idx) + MF_PARAM(p_iter->getIndex()));
        }
      }
    }
/*  }
  else   //Second maximazation of transformation procedure.
         //Random variance is a constant minus all other variance components.
  {
    const MAXFUN::ParameterMgr&  first_phase_mgr = results.maxfun_results_null.getParameterMgr();
    double  first_phase_total_variance = first_phase_mgr.getParameter("Other parameters", "Total variance").getFinalEstimate();  
  
    tv_idx = mgr.addParameter("Other parameters", "Total variance",
                                   MAXFUN::Parameter::FIXED, first_phase_total_variance);
    func.addDependentParamStep(tv_idx, MF_CONSTANT(first_phase_total_variance));

    size_t  rand_idx = fixed_effect_indices[Configuration::RANDOM];
    
    mgr.getParameter(rand_idx).setInitialType(MAXFUN::Parameter::DEPENDENT);
    func.addDependentParamStep(rand_idx, MF_CONSTANT(first_phase_total_variance));
    
    if(fixed_effect_usage[Configuration::SIBLING])
    {
      func.addDependentParamStep(rand_idx, MF_PARAM(rand_idx) - MF_PARAM(fixed_effect_indices[Configuration::SIBLING]));
    }
    
    if(fixed_effect_usage[Configuration::MARITAL])
    {
      func.addDependentParamStep(rand_idx, MF_PARAM(rand_idx) - MF_PARAM(fixed_effect_indices[Configuration::MARITAL]));
    }
    
    if(fixed_effect_usage[Configuration::POLYGENIC])
    { 
      func.addDependentParamStep(rand_idx, MF_PARAM(rand_idx) - MF_PARAM(fixed_effect_indices[Configuration::POLYGENIC]));
    }
    
    if(fixed_effect_usage[Configuration::FAMILY]) 
    {
      func.addDependentParamStep(rand_idx, MF_PARAM(rand_idx) - MF_PARAM(fixed_effect_indices[Configuration::FAMILY]));
    }
                                       
    // - Add user defined effects to total variance.
    //
    MAXFUN::ParameterConstIterator  p_iter     = mgr.getParamBegin("Variance components");
    MAXFUN::ParameterConstIterator  p_end_iter = mgr.getParamEnd("Variance components");
    for(; p_iter != p_end_iter; ++p_iter)
    {
      string  param_name = p_iter->getName();
      if(! (fixed_effect(param_name) || param_name == "0.5 * Polygenic"))
      {
        if(my_config.modelHasEffect(my_model_name, alt_covariates_included, param_name))
        {
          func.addDependentParamStep(rand_idx, MF_PARAM(rand_idx) - MF_PARAM(p_iter->getIndex()));
        }
      }
    } 
  } */
}                                      

void
MaximizationWrapper::addFixedEffects(MAXFUN::Function& func,
                                     vector<size_t>& fixed_effect_indices,
                                     const vector<bool>& fixed_effect_usage)
{
  MAXFUN::ParameterMgr&  mgr = func.getMgr();

  mgr.addGroup("Variance components");
  
  fixed_effect_indices[Configuration::RANDOM] = mgr.addParameter(my_config.getFixedEffect(Configuration::RANDOM));
  
  if(fixed_effect_usage[Configuration::POLYGENIC])
  {
    fixed_effect_indices[Configuration::POLYGENIC] = mgr.addParameter(my_config.getFixedEffect(Configuration::POLYGENIC));
    mgr.getParameter(fixed_effect_indices[Configuration::POLYGENIC]).setPValueInf(MAXFUN::Parameter::ONE_SIDED);

    MAXFUN::Parameter&  half_polygenic = mgr.addParameterAlt("Variance components", "0.5 * Polygenic", MAXFUN::Parameter::DEPENDENT);
    half_polygenic.setIncludeInOutput(false);
    
    func.addDependentParamStep(half_polygenic.getIndex(), 
                               (MF_PARAM(fixed_effect_indices[Configuration::POLYGENIC]) * MF_CONSTANT(0.5)));
  }

  if(fixed_effect_usage[Configuration::FAMILY])
  {
    fixed_effect_indices[Configuration::FAMILY] = mgr.addParameter(my_config.getFixedEffect(Configuration::FAMILY));
    mgr.getParameter(fixed_effect_indices[Configuration::FAMILY]).setPValueInf(MAXFUN::Parameter::ONE_SIDED);
  }

  if(fixed_effect_usage[Configuration::MARITAL])
  {
    fixed_effect_indices[Configuration::MARITAL] = mgr.addParameter(my_config.getFixedEffect(Configuration::MARITAL));
    mgr.getParameter(fixed_effect_indices[Configuration::MARITAL]).setPValueInf(MAXFUN::Parameter::ONE_SIDED);  
  }

  if(fixed_effect_usage[Configuration::SIBLING])
  {
    fixed_effect_indices[Configuration::SIBLING] = mgr.addParameter(my_config.getFixedEffect(Configuration::SIBLING));
    mgr.getParameter(fixed_effect_indices[Configuration::SIBLING]).setPValueInf(MAXFUN::Parameter::ONE_SIDED);
  }
}


void
MaximizationWrapper::addUserEffects(MAXFUN::Function& func)
{
  MAXFUN::ParameterMgr&  mgr = func.getMgr();
  
  Configuration::effect_iterator  ue_iter     = my_config.userEffectBegin();
  Configuration::effect_iterator  ue_end_iter = my_config.userEffectEnd();
  for(; ue_iter != ue_end_iter; ++ue_iter)
  {
    if(my_config.modelHasEffect(my_model_name, alt_covariates_included, ue_iter->param_name))
    {
      size_t  effect_idx = mgr.addParameter(*ue_iter);
      mgr.getParameter(effect_idx).setPValueInf(MAXFUN::Parameter::ONE_SIDED);  
    }
  }
}                                     

// Add covariates / standardized covariates:
// Covariate parameters are corrected to put them back on the unstandardized scale.
// The standardized covariates are those parameters estimated using standardized data.
// values.  -djb
//
void
MaximizationWrapper::addNullCovCoefficients(MAXFUN::Function& func) const
{
  MAXFUN::ParameterMgr&  mgr = func.getMgr();
  
  mgr.addGroup("Covariates");
  mgr.addGroup("Standardized covariates");

  Configuration::NullCovConstIter  cov_info_iter     = my_config.nullCovariateBegin(my_model_name);
  Configuration::NullCovConstIter  cov_info_end_iter = my_config.nullCovariateEnd(my_model_name);
  for(; cov_info_iter != cov_info_end_iter; ++cov_info_iter)
  {
    MAXFUN::Parameter::ParamTypeEnum  ptype = cov_info_iter->cfg.initial_type == MAXFUN::Parameter::FIXED ? 
                                              MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL;
                                              
    // - Dependent parameter must be added first or strange results are produced.  -djb 10/3/7
    //
    size_t  cov_idx = mgr.addParameter("Covariates", cov_info_iter->cfg.param_name, MAXFUN::Parameter::DEPENDENT);
    size_t  std_cov_idx = mgr.addParameter("Standardized covariates", cov_info_iter->cfg.param_name, ptype);
    mgr.getParameter(std_cov_idx).setIncludeInOutput(false);
        
    double stdev = my_sampledata.getField("Covariates", cov_info_iter->cfg.param_name).getStdev();
    func.addDependentParamStep(cov_idx, MF_PARAM(std_cov_idx) / MF_CONSTANT(stdev));
  }
}

void
MaximizationWrapper::addTestCovCoefficients(MAXFUN::Function& func) const
{
  MAXFUN::ParameterMgr&  mgr = func.getMgr();

  Configuration::TestCovConstIter  cov_info_iter     = my_config.testCovariateBegin(my_model_name);
  Configuration::TestCovConstIter  cov_info_end_iter = my_config.testCovariateEnd(my_model_name);
  for(; cov_info_iter != cov_info_end_iter; ++cov_info_iter)
  {
    size_t  cov_idx = mgr.addParameter("Covariates", cov_info_iter->cfg.param_name, MAXFUN::Parameter::DEPENDENT);
  
    MAXFUN::Parameter::ParamTypeEnum ptype = cov_info_iter->cfg.initial_type == MAXFUN::Parameter::FIXED ? 
                                        MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL;

    size_t  std_cov_idx   = mgr.addParameter("Standardized covariates", cov_info_iter->cfg.param_name, ptype);
    
    
    if(score_test)
    {
      mgr.getParameter(std_cov_idx).scoretest = true;
      mgr.getParameter(std_cov_idx).setInitialEstimate(0.0);
    }
    
    
    double stdev = my_sampledata.getField("Covariates", cov_info_iter->cfg.param_name).getStdev();

    mgr.getParameter(std_cov_idx).setIncludeInOutput(false);

    func.addDependentParamStep(cov_idx, MF_PARAM(std_cov_idx) / MF_CONSTANT(stdev));
    
  }
}

// - Removed step parent - offspring correlation per rce.  See trac ticket #1779.  -djb
//   Extensively modified 10-11-7 per rce.  See trac ticket #1797.  -djb
//
void
MaximizationWrapper::addFamilialCorrelations(MAXFUN::Function& func, 
                                             const vector<size_t>& fixed_effect_indices,
                                             const vector<bool>& fixed_effect_usage, size_t tv_idx) const
{
  MAXFUN::ParameterMgr&  mgr = func.getMgr();

  if(fixed_effect_usage[Configuration::FAMILY]  || 
     fixed_effect_usage[Configuration::SIBLING] || 
     fixed_effect_usage[Configuration::POLYGENIC] )
  {
    size_t  full_sib_idx = mgr.addParameter("Residual familial correlations", "Res full sibs corr", MAXFUN::Parameter::DEPENDENT);
    func.addDependentParamStep(full_sib_idx, MF_CONSTANT(0));
    
    if(fixed_effect_usage[Configuration::FAMILY])
    {
      func.addDependentParamStep(full_sib_idx, MF_PARAM(full_sib_idx) + MF_PARAM(fixed_effect_indices[Configuration::FAMILY])); 
    }
    
    if(fixed_effect_usage[Configuration::SIBLING])
    {
      func.addDependentParamStep(full_sib_idx, MF_PARAM(full_sib_idx) + MF_PARAM(fixed_effect_indices[Configuration::SIBLING])); 
    }
    
    if(fixed_effect_usage[Configuration::POLYGENIC])
    {
      func.addDependentParamStep(full_sib_idx, MF_PARAM(full_sib_idx) + 
                                 MF_PARAM(fixed_effect_indices[Configuration::POLYGENIC]) * MF_CONSTANT(.5)); 
    }
    
    func.addDependentParamStep(full_sib_idx, MF_PARAM(full_sib_idx) / MF_PARAM(tv_idx));
  }

  if(fixed_effect_usage[Configuration::FAMILY]   || 
     fixed_effect_usage[Configuration::POLYGENIC]  )
  {
    // - Half sib
    //
    size_t  half_sib_idx = mgr.addParameter("Residual familial correlations", "Res half sibs corr", MAXFUN::Parameter::DEPENDENT);
    func.addDependentParamStep(half_sib_idx, MF_CONSTANT(0));
    
    if(fixed_effect_usage[Configuration::FAMILY])
    {
      func.addDependentParamStep(half_sib_idx, MF_PARAM(half_sib_idx) + MF_PARAM(fixed_effect_indices[Configuration::FAMILY])); 
    }
    
    if(fixed_effect_usage[Configuration::POLYGENIC])
    {
      func.addDependentParamStep(half_sib_idx, MF_PARAM(half_sib_idx) + 
                                 MF_PARAM(fixed_effect_indices[Configuration::POLYGENIC]) * MF_CONSTANT(.25)); 
    }
    
    func.addDependentParamStep(half_sib_idx, MF_PARAM(half_sib_idx) / MF_PARAM(tv_idx));
    
    // - Parent - offspring
    //
    size_t  parent_offspring_idx = mgr.addParameter("Residual familial correlations", "Res parent offspring corr", MAXFUN::Parameter::DEPENDENT);
    func.addDependentParamStep(parent_offspring_idx, MF_CONSTANT(0));
    
    if(fixed_effect_usage[Configuration::FAMILY])
    {
      func.addDependentParamStep(parent_offspring_idx, MF_PARAM(parent_offspring_idx) + 
                                 MF_PARAM(fixed_effect_indices[Configuration::FAMILY])); 
    }
    
    if(fixed_effect_usage[Configuration::POLYGENIC])
    {
      func.addDependentParamStep(parent_offspring_idx, MF_PARAM(parent_offspring_idx) + 
                                 MF_PARAM(fixed_effect_indices[Configuration::POLYGENIC]) * MF_CONSTANT(.5)); 
    }
    
    func.addDependentParamStep(parent_offspring_idx, MF_PARAM(parent_offspring_idx) / MF_PARAM(tv_idx));    
  }
  
  if(fixed_effect_usage[Configuration::FAMILY] || 
     fixed_effect_usage[Configuration::MARITAL]  )
  {
    size_t  spouses_idx = mgr.addParameter("Residual familial correlations", "Res marital(spouse) corr", MAXFUN::Parameter::DEPENDENT);
    func.addDependentParamStep(spouses_idx, MF_CONSTANT(0));
    
    if(fixed_effect_usage[Configuration::FAMILY])
    {
      func.addDependentParamStep(spouses_idx, MF_PARAM(spouses_idx) + MF_PARAM(fixed_effect_indices[Configuration::FAMILY])); 
    }
    
    if(fixed_effect_usage[Configuration::MARITAL])
    {
      func.addDependentParamStep(spouses_idx, MF_PARAM(spouses_idx) + MF_PARAM(fixed_effect_indices[Configuration::MARITAL])); 
    }
    
    func.addDependentParamStep(spouses_idx, MF_PARAM(spouses_idx) / MF_PARAM(tv_idx));
  }
}

// - Extensively modified 10-11-7 per rce.  See trac ticket #1797.  -djb
//
void
MaximizationWrapper::addEnvironmentalCorrelations(MAXFUN::Function& func, 
                                                  const vector<size_t>& fixed_effect_indices,
                                                  const vector<bool>& fixed_effect_usage, size_t tv_idx) const
{
  MAXFUN::ParameterMgr&  mgr = func.getMgr();

  if(fixed_effect_usage[Configuration::FAMILY])
  {
    size_t  nuc_family_idx = mgr.addParameter("Environmental intraclass correlations", "Env nuclear family corr", MAXFUN::Parameter::DEPENDENT);
    
    if(fixed_effect_usage[Configuration::POLYGENIC])
    {
      func.addDependentParamStep(nuc_family_idx, MF_PARAM(fixed_effect_indices[Configuration::FAMILY]) / 
                                                (MF_PARAM(tv_idx)- MF_PARAM(fixed_effect_indices[Configuration::POLYGENIC]))); 
    }
    else
    {
      func.addDependentParamStep(nuc_family_idx, MF_PARAM(fixed_effect_indices[Configuration::FAMILY]) / MF_PARAM(tv_idx));     
    }
  }

  if(fixed_effect_usage[Configuration::FAMILY] || fixed_effect_usage[Configuration::MARITAL])
  {
    size_t  marital_corr_idx = mgr.addParameter("Environmental intraclass correlations", "Env marital(spouse) corr", MAXFUN::Parameter::DEPENDENT);
    func.addDependentParamStep(marital_corr_idx, MF_CONSTANT(0));
    
    if(fixed_effect_usage[Configuration::FAMILY])
    {
      func.addDependentParamStep(marital_corr_idx, MF_PARAM(marital_corr_idx) + 
                                                   MF_PARAM(fixed_effect_indices[Configuration::FAMILY])); 
    }
    
    if(fixed_effect_usage[Configuration::MARITAL])
    {
      func.addDependentParamStep(marital_corr_idx, MF_PARAM(marital_corr_idx) + 
                                                   MF_PARAM(fixed_effect_indices[Configuration::MARITAL])); 
    }
    
    if(fixed_effect_usage[Configuration::POLYGENIC])
    {
      func.addDependentParamStep(marital_corr_idx, MF_PARAM(marital_corr_idx) / 
                                                  (MF_PARAM(tv_idx) - MF_PARAM(fixed_effect_indices[Configuration::POLYGENIC])));
    }
    else
    {
      func.addDependentParamStep(marital_corr_idx, MF_PARAM(marital_corr_idx) / MF_PARAM(tv_idx));    
    }
  }
    
  if(fixed_effect_usage[Configuration::FAMILY] || fixed_effect_usage[Configuration::SIBLING])
  {
    size_t  sibship_idx = mgr.addParameter("Environmental intraclass correlations", "Env full sibs corr", MAXFUN::Parameter::DEPENDENT);
    func.addDependentParamStep(sibship_idx, MF_CONSTANT(0));
    
    if(fixed_effect_usage[Configuration::FAMILY])
    {
      func.addDependentParamStep(sibship_idx, MF_PARAM(sibship_idx) + MF_PARAM(fixed_effect_indices[Configuration::FAMILY])); 
    }
    
    if(fixed_effect_usage[Configuration::SIBLING])
    {
      func.addDependentParamStep(sibship_idx, MF_PARAM(sibship_idx) + MF_PARAM(fixed_effect_indices[Configuration::SIBLING])); 
    }
    
    if(fixed_effect_usage[Configuration::POLYGENIC])
    {
      func.addDependentParamStep(sibship_idx, MF_PARAM(sibship_idx) / 
                                             (MF_PARAM(tv_idx) - MF_PARAM(fixed_effect_indices[Configuration::POLYGENIC])));
    }
    else
    {
      func.addDependentParamStep(sibship_idx, MF_PARAM(sibship_idx) / MF_PARAM(tv_idx));    
    }
  }    
}


void
MaximizationWrapper::addIntercept(MAXFUN::Function& func) const
{
  MAXFUN::ParameterMgr&  mgr = func.getMgr();
    
  mgr.addGroup("Intercept");
  mgr.addParameter("Intercept", "Intercept", 
                    MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL, 0.0, -MAXFUN::MF_INFINITY, MAXFUN::MF_INFINITY);
}


// - Override calculated values with  user specified values for covariate coefficients.
//   Added 10/3/7  -djb
//   
void
MaximizationWrapper::setUsersCovVals(MAXFUN::ParameterMgr& mgr) const
{
  MAXFUN::ParameterIterator  p_iter     = mgr.getParamBegin("Standardized covariates");
  MAXFUN::ParameterIterator  p_end_iter = mgr.getParamEnd("Standardized covariates");
  for(; p_iter != p_end_iter; ++p_iter)
  {
    double  init_value = my_config.getCovariateInfo(p_iter->getName()).cfg.initial_estimate;
    
    if(p_iter->getInitialType() == MAXFUN::Parameter::FIXED)
    {
      assert(! isnan(init_value));
    }
      
    if(! isnan(init_value))
    {
      p_iter->setInitialEstimate(init_value);
    }
  }
}

void
MaximizationWrapper::establishInitialEstimates(MAXFUN::ParameterMgr& mgr, double prevalence, 
                                                const AnalysisResults::ModelResults& results)
{
  if(! null_results_available)
  {
    calculateInitialEsts(mgr, prevalence);
    setUsersCovVals(mgr);
  }
  else   // Use estimates from null model as initial values for test model or
         // use estimates from first of two maximizations as initial values for the second maximization.
  {
    usePreviousResults(mgr, results.maxfun_results_null.getParameterMgr());
  }
}


void
MaximizationWrapper::usePreviousResults(MAXFUN::ParameterMgr& mgr, const MAXFUN::ParameterMgr& previous_mgr)
{
  // - Variance components
  //
  MAXFUN::ParameterIterator  p_iter     = mgr.getParamBegin("Variance components");
  MAXFUN::ParameterIterator  p_end_iter = mgr.getParamEnd("Variance components");
  for(; p_iter != p_end_iter; ++p_iter)
  {
    string  param_name = p_iter->getName();
    if(param_name != "0.5 * Polygenic" && p_iter->getInitialType() != MAXFUN::Parameter::FIXED)
    {
      p_iter->setInitialEstimate(previous_mgr.getParameter("Variance components", param_name).getFinalEstimate() / (my_divisor * my_divisor));
    }
  }
  
  // - Covariates 
  //
  p_iter     = mgr.getParamBegin("Standardized covariates");
  p_end_iter = mgr.getParamEnd("Standardized covariates");
  for(; p_iter != p_end_iter; ++p_iter)
  {
    string  param_name = p_iter->getName();
    if(previous_mgr.doesParamExist("Standardized covariates", param_name))    // Null covariate
    {  
      double  previous_value = previous_mgr.getParameter("Standardized covariates", param_name).getFinalEstimate();
      p_iter->setInitialEstimate(previous_value);
      
      
      if(my_config.getDependentTraitType() == BINARY && my_phase == TWO_OF_TWO)
      {
        p_iter->setInitialType(MAXFUN::Parameter::FIXED);
      }
      
    }
    else                                                                   // Test covariate
    {
      const Configuration::CovariateInfo&  user_settings = my_config.getCovariateInfo(param_name);
      if(! isnan(user_settings.cfg.initial_estimate))
      {
        p_iter->setInitialEstimate(user_settings.cfg.initial_estimate);
      }
      else
      {
        p_iter->setInitialEstimate(0.0);
      }
    }
  }
  
  // - Intercept
  //
  double  previous_value = previous_mgr.getParameter("Intercept", "Intercept").getFinalEstimate();
  mgr.getParameter("Intercept", "Intercept").setInitialEstimate(previous_value);
  
  /* Don't fix intercept - RCE, March 2012.
  if(my_config.getDependentTraitType() == BINARY && my_phase == TWO_OF_TWO)
  {
    mgr.getParameter("Intercept", "Intercept").setInitialType(MAXFUN::Parameter::FIXED);
  }
  */
  
  // - Transformation parameters
  // 
  if(alt_covariates_included)
  {
    usePreviousLambdas(mgr, previous_mgr);
  }
}

void  
MaximizationWrapper::usePreviousLambdas(MAXFUN::ParameterMgr& mgr, const MAXFUN::ParameterMgr& previous_mgr) const
{
  if(previous_mgr.doesParamExist("Transformation", "Lambda1")         && 
     mgr.doesParamExist("Transformation", "Lambda1")               &&
     previous_mgr.getParameter("Transformation", "Lambda1").isEstimated() &&
     mgr.getParameter("Transformation", "Lambda1").isEstimated()        )
  {
    mgr.getParameter("Transformation", "Lambda1").
      setInitialEstimate(previous_mgr.getParameter("Transformation", "Lambda1").getFinalEstimate());
  }
  
  if(previous_mgr.doesParamExist("Transformation", "Lambda2")         && 
     mgr.doesParamExist("Transformation", "Lambda2")               &&
     previous_mgr.getParameter("Transformation", "Lambda2").isEstimated() &&
     mgr.getParameter("Transformation", "Lambda2").isEstimated()        )
  {
    mgr.getParameter("Transformation", "Lambda2").
      setInitialEstimate(previous_mgr.getParameter("Transformation", "Lambda2").getFinalEstimate());
  }
}


//=======================================================================
//  calculateInitialEsts()
//=======================================================================
// Use linear regression to calculate initial values for both the beta vector and
// the total variance.
void
MaximizationWrapper::calculateInitialEsts(MAXFUN::ParameterMgr& mgr, double prevalence)
{
  size_t  num_of_valid_inds = (size_t)my_sampledata.getValidIndividualCount();
  checkDataSufficiency(num_of_valid_inds, mgr);

  size_t  num_of_covs = (size_t)mgr.getParamCount("Covariates");

  FortranMatrix<double>  all_covariates_mod(0, 0);
  FortranMatrix<double>  cov(0, 0);
  FortranMatrix<double>  A(0, 0);
  FortranMatrix<double>  residual(0, 0);
  FortranMatrix<double>  main_phenotypes(1, num_of_valid_inds);
  FortranMatrix<double>  all_covariates(1 + num_of_covs, num_of_valid_inds);
  
  populateMatrices(main_phenotypes, all_covariates, mgr, prevalence);

  // Now make all the magical math happen:
  main_phenotypes.transpose();
  all_covariates.transpose();
  XTX(all_covariates, all_covariates_mod);
  cov = SVDinverse(all_covariates_mod);
  A = cov * (transpose(all_covariates) * main_phenotypes);

  // Calculate the total variance:
  residual = main_phenotypes - all_covariates * A;

  double  est_total_variance = 0.0;
  for(size_t i = 0; i < num_of_valid_inds; ++i)
    est_total_variance += residual(i, 0) * residual(i, 0);

  est_total_variance /= (num_of_valid_inds - (1 + num_of_covs));
  setVarianceComponents(est_total_variance, mgr);
  
  if(my_config.getDependentTraitType() != BINARY)
  {
    mgr.getParameter("Intercept", "Intercept").setInitialEstimate(A(0, 0));
  }
  else
  {
    double  mean = my_sampledata.getField("Core traits", "Main phenotype").getMean();
    double  init_intercept = log(mean / (1 - mean));
    mgr.getParameter("Intercept", "Intercept").setInitialEstimate(init_intercept);    
  }

  // Set initial estimates for standardized covariates
  vector<double>  initial_estimates(num_of_covs);  
  for(size_t cov_id = 0; cov_id < num_of_covs; ++cov_id)
    initial_estimates[cov_id] = my_config.getDependentTraitType() == BINARY ? 0.0 : A(1 + cov_id, 0);
    
  MAXFUN::ParameterIterator param     = mgr.getParamBegin("Standardized covariates");
  MAXFUN::ParameterIterator param_end = mgr.getParamEnd("Standardized covariates");
  for(; param != param_end; ++param)
  {
    param->setInitialEstimate(initial_estimates[param->getGroupIndex()]);
  }
}

void
MaximizationWrapper::checkDataSufficiency(size_t num_of_valid_inds, const MAXFUN::ParameterMgr& mgr)
{
  if(num_of_valid_inds + 1 <= (size_t)mgr.getEstimatedParamCount())
  {
    throw(InsufficientData(my_model_name));
  }
}

void
MaximizationWrapper::populateMatrices(FortranMatrix<double>& main_phenotypes, FortranMatrix<double>& all_covariates,  
                                      const MAXFUN::ParameterMgr& mgr, double prevalence)
{
  size_t  num_of_all_inds   = (size_t)my_sampledata.getTotalIndividualCount();    
  for(size_t i = 0, id = 0; i < num_of_all_inds; ++i)
  {
    if(my_sampledata.isValid(i))
    {
      if(my_config.getDependentTraitType() == BINARY)
      {
        main_phenotypes(0, id) =  calculateBinaryPhenotype(prevalence, i, mgr);  
      }
      else
      {
        main_phenotypes(0, id) = my_sampledata.getAdjValue(i, "Core traits", "Main phenotype");
      }
        
      all_covariates(0, id)  = 1.0;

      MAXFUN::ParameterConstIterator param     = mgr.getParamBegin("Covariates");
      MAXFUN::ParameterConstIterator param_end = mgr.getParamEnd("Covariates");
      for(size_t cov_id = 0; param != param_end; ++param, ++cov_id)
      {
        all_covariates(1 + cov_id, id) = my_sampledata.getField("Covariates", param->getName()).getAdjValue(i);
      }

      ++id;
    }
  }
}

double
MaximizationWrapper::calculateBinaryPhenotype(double prevalence, size_t ind_idx, const MAXFUN::ParameterMgr& mgr) const
{
  double  prev = prevalence;
  
  if(isnan(prev))
  {
    prev = calculatePrevalence(ind_idx, mgr);
  }
  
  double  p_present        = prev;
  double  p_absent         = 1.0 - prev;
  double  binary_z_present = (1.0 - p_present) / sqrt(p_present * (1.0 - p_present));
         
  // Altered per discussion w. Courtney 6-26-7.  djb 
  //
  double  binary_z_absent  = (0.0 - p_absent) / sqrt(p_absent * (1.0 - p_absent));
  double  phenotype = my_sampledata.getAdjValue(ind_idx, "Core traits", "Main phenotype");
  
  return  phenotype ? binary_z_present : binary_z_absent;
}

// - Equation 3.8 in Courtney Gray-Mcguire dissertation.
//
double
MaximizationWrapper::calculatePrevalence(size_t ind_idx, const MAXFUN::ParameterMgr& mgr) const
{
  double  exponent = mgr.getParameter("Intercept", "Intercept").getFinalEstimate();
  MAXFUN::ParameterConstIterator param     = mgr.getParamBegin("Standardized Covariates");
  MAXFUN::ParameterConstIterator param_end = mgr.getParamEnd("Standardized Covariates");
  for(; param != param_end; ++param)
  {
    exponent += param->getFinalEstimate() *
                my_sampledata.getField("Covariates", param->getName()).getAdjValue(ind_idx);
  }

  return  exp(exponent) / (1 + exp(exponent));  
}


// Split the total variance into sub-variances (the distribution is given by RCE)
void  
MaximizationWrapper::setVarianceComponents(double est_total_variance, MAXFUN::ParameterMgr& mgr) const
{
  bool  use_P = my_config.modelHasFixedEffect(my_model_name, alt_covariates_included, Configuration::POLYGENIC);
  int  var_component_count = mgr.getParamCount("Variance components");
  
  if(use_P)
  {
    --var_component_count;        // Don't count "0.5 * Polygenic".
  }
  
  double  est_component  = est_total_variance / var_component_count;
  MAXFUN::ParameterIterator  p_iter     = mgr.getParamBegin("Variance components");
  MAXFUN::ParameterIterator  p_end_iter = mgr.getParamEnd("Variance components");
  for(; p_iter != p_end_iter; ++p_iter)
  {
    string  param_name = p_iter->getName();
    if(param_name != "0.5 * Polygenic")
    {
      if(SAGE::isnan(p_iter->getInitialEstimate()))
      {
        p_iter->setInitialEstimate(est_component);
      }
    }
  }
  
  checkVarianceEstimates(mgr);
}
 
// - Make sure the initial estimates are in-bounds.
//
void
MaximizationWrapper::checkVarianceEstimates(MAXFUN::ParameterMgr& mgr) const
{
  MAXFUN::ParameterIterator  p_iter     = mgr.getParamBegin("Variance components");
  MAXFUN::ParameterIterator  p_end_iter = mgr.getParamEnd("Variance components");
  for(; p_iter != p_end_iter; ++p_iter)
  {
    string  param_name = p_iter->getName();
    if(param_name != "0.5 * Polygenic")
    {
      double  param_est         = p_iter->getInitialEstimate();
      double  param_lower_bound = p_iter->getLowerBound();
      if(param_est < param_lower_bound)
      {
        p_iter->setInitialEstimate(2.0 * param_lower_bound);
      } 
    }
  }
}

void
MaximizationWrapper::rememberEstimates(vector<double>& estimates, const MAXFUN::ParameterMgr& mgr) const
{
  estimates.clear();

  estimates.push_back(mgr.getParameter("Intercept", "Intercept").getCurrentEstimate());
  
  MAXFUN::ParameterConstIterator  param_iter     = mgr.getParamBegin("Standardized covariates");
  MAXFUN::ParameterConstIterator  param_end_iter = mgr.getParamEnd("Standardized covariates");
  for(; param_iter != param_end_iter; ++param_iter)
  {
    estimates.push_back(param_iter->getCurrentEstimate());
  }
}    

void
MaximizationWrapper::setNewInitialValues(MAXFUN::ParameterMgr& mgr)
{
  calculateInitialEsts(mgr, QNAN);
}

double
MaximizationWrapper::calculateTotalBetaChange(const vector<double>& previous, const vector<double>& current) const
{
  double  total_beta_change = 0.0;
  for(size_t i = 0; i < previous.size(); ++i)
  {
    total_beta_change += abs(previous[i] - current[i]);
  }
  
  return  total_beta_change;
}                                           
                                                                            

} // End namespace ASSOC
} // End namespace SAGE
