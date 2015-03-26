//=============================================================================
//
//  File:	AnalysisResults.cpp
//
//  Author:	Stephen Gross
//
//  Copyright 2002 R. C. Elston
//=============================================================================

#include "assoc/AnalysisResults.h"

namespace SAGE  {
namespace ASSOC {


//=============================================================================
//  IMPLEMENTATION:  AnalysisResults::ModelResults 
//=============================================================================
//
AnalysisResults::ModelResults::ModelResults()
{
  score_statistic = QNAN;
  residuals_scaling_factor = 1.0;
  affected_mean_residuals_null = QNAN;
  unaffected_mean_residuals_null = QNAN;
  affected_mean_residuals_alt = QNAN;
  unaffected_mean_residuals_alt = QNAN;  
}


// - Do the Null and Alt results have the same components?
//
bool
AnalysisResults::ModelResults::sameVarianceComponents() const
{
  MAXFUN::ParameterMgr  null_parameters = maxfun_results_null.getParameterMgr();
  MAXFUN::ParameterMgr  alt_parameters  = maxfun_results_alt.getParameterMgr();
  
  set<string>  null_components;
  set<string>  alt_components;
  
  MAXFUN::ParameterConstIterator  p_null_iter     = null_parameters.getParamBegin("Variance components");
  MAXFUN::ParameterConstIterator  p_null_end_iter = null_parameters.getParamEnd("Variance components");
  for(; p_null_iter != p_null_end_iter; ++p_null_iter)
  {
    null_components.insert(p_null_iter->getName());
  }
  
  MAXFUN::ParameterConstIterator  p_alt_iter     = alt_parameters.getParamBegin("Variance components");
  MAXFUN::ParameterConstIterator  p_alt_end_iter = alt_parameters.getParamEnd("Variance components");
  for(; p_alt_iter != p_alt_end_iter; ++p_alt_iter)
  {
    alt_components.insert(p_alt_iter->getName());
  }  
  
  return  null_components == alt_components; 
}

//=============================================================================
//  IMPLEMENTATION:  AnalysisResults::ComparisonInfo 
//=============================================================================
//
AnalysisResults::ComparisonInfo::ComparisonInfo()
    : converged(false), baseline_converged(false), lrt_pvalue(QNAN), 
      score_pvalue(QNAN), same_variance_components(false)
{}

double
AnalysisResults::ComparisonInfo::waldPvalue() const
{
  const MAXFUN::Parameter&  param = covs[covs.size() - 1];
  
  double  wald_pvalue = param.isPValueAvailable() ? param.getPValue() : QNAN;
  
  if(! isnan(wald_pvalue))
  {
    // - Don't allow Wald pvalue to be 0.
    //
    wald_pvalue = wald_pvalue == 0 ? numeric_limits<double>::min() : wald_pvalue;  
  }
  
  return  wald_pvalue;
}

// - Return larger of LRT and Wald pvalues accounting for possibility that
//   one or both are missing.
//
double
AnalysisResults::ComparisonInfo::largerPvalue() const
{
  double  wald_pvalue = waldPvalue();

  if(! (isnan(wald_pvalue) || isnan(lrt_pvalue)))
  {
    return  wald_pvalue > lrt_pvalue ? wald_pvalue : lrt_pvalue;
  }
  else if(isnan(wald_pvalue) && ! isnan(lrt_pvalue))
  {
    return  lrt_pvalue;
  }
  else if(! isnan(wald_pvalue) && isnan(lrt_pvalue))
  {
    return  wald_pvalue;
  }
  else              // Both are unknown
  {
    return  QNAN;
  }
}

double
AnalysisResults::ComparisonInfo::pvalueRatio() const
{
  double  wald_pvalue = waldPvalue();

  if(! (isnan(wald_pvalue) || isnan(lrt_pvalue)))
  {
    return  min(-log10(lrt_pvalue), -log10(wald_pvalue)) -
            max(-log10(lrt_pvalue), -log10(wald_pvalue))  ;
  }
  else              
  {
    return  QNAN;
  }
}

//=============================================================================
//  IMPLEMENTATION:  AnalysisResults::LRTLess 
//=============================================================================
//
bool
AnalysisResults::LRTLess::operator()(const ComparisonInfo& lhs,
                                     const ComparisonInfo& rhs) const
{
  return  AnalysisResults::less(lhs.lrt_pvalue, rhs.lrt_pvalue);
}

//=============================================================================
//  IMPLEMENTATION:  AnalysisResults::WaldLess 
//=============================================================================
//
bool
AnalysisResults::WaldLess::operator()(const ComparisonInfo& lhs,
                                      const ComparisonInfo& rhs) const
{
  return  AnalysisResults::less(lhs.waldPvalue(), rhs.waldPvalue());
}

//=============================================================================
//  IMPLEMENTATION:  AnalysisResults::LargerLess 
//=============================================================================
//
// - Use the larger of the LRT and Wald pvalues for this comparison.
//
bool
AnalysisResults::LargerLess::operator()(const ComparisonInfo& lhs,
                                        const ComparisonInfo& rhs) const
{
  return AnalysisResults::less(lhs.largerPvalue(), rhs.largerPvalue());
}

//=============================================================================
//  IMPLEMENTATION:  AnalysisResults::RatioLess 
//=============================================================================
//
// - Use the ratio of the LRT and Wald pvalues for this comparison.
//
bool
AnalysisResults::RatioLess::operator()(const ComparisonInfo& lhs,
                                        const ComparisonInfo& rhs) const
{
  return  AnalysisResults::greater(lhs.pvalueRatio(), rhs.pvalueRatio());
}    

//=============================================================================
//  IMPLEMENTATION:  AnalysisResults::LRTFilter 
//=============================================================================
//
AnalysisResults::LRTFilter::LRTFilter(double cutpoint)
    : my_cutpoint(cutpoint)
{}

bool
AnalysisResults::LRTFilter::operator()(const ComparisonInfo& c) const
{
  return  ! less(c.lrt_pvalue, my_cutpoint);
}


//=============================================================================
//  IMPLEMENTATION:  AnalysisResults::WaldFilter 
//=============================================================================
//
AnalysisResults::WaldFilter::WaldFilter(double cutpoint)
    : my_cutpoint(cutpoint)
{}

bool
AnalysisResults::WaldFilter::operator()(const ComparisonInfo& c) const
{
  return  ! less(c.waldPvalue(), my_cutpoint);
}

//=============================================================================
//  IMPLEMENTATION:  AnalysisResults 
//=============================================================================
//
// - Will put undfined values last if used to put values in ascending order.
//
bool
AnalysisResults::less(double lhs, double rhs)
{
  if(! (isnan(lhs) || isnan(rhs)))
  {
    return  lhs < rhs;
  }
  else if(isnan(lhs) && ! isnan(rhs))
  {
    return  false;
  }
  else if(! isnan(lhs) && isnan(rhs))
  {
    return  true;
  }
  else              // Both are unknown
  {
    return  false;
  }
}

// - Will put undfined values last if used to put values in descending order.
//
bool
AnalysisResults::greater(double lhs, double rhs)
{
  if(! (isnan(lhs) || isnan(rhs)))
  {
    return  lhs > rhs;
  }
  else if(isnan(lhs) && ! isnan(rhs))
  {
    return  false;
  }
  else if(! isnan(lhs) && isnan(rhs))
  {
    return  true;
  }
  else              // Both are unknown
  {
    return  false;
  }
}

//=============================================================================
//  AnalysisResults() 
//=============================================================================
AnalysisResults::AnalysisResults(const Configuration& config)
    : my_config(config)
{ 
  my_cinfos.clear();
  my_detailed_sections.clear();
  my_null_residual_outputs.clear();
  my_alt_residual_outputs.clear();  
}

//=============================================================================
// addModelResults()
//=============================================================================
void 
AnalysisResults::addModelResults(const ModelResults& r, const Sampledata& sampledata) 
{ 
  addModelDetailedSection(r);
  addModelResidualOutput(r, sampledata);
  //addIndependentResiduals(r, sampledata);
  addPValues(r);
}

const Configuration&
AnalysisResults::getConfig() const
{
  return  my_config;
}


const map<string, string>&
AnalysisResults::getNullResidualOutputs() const
{
  return  my_null_residual_outputs;
}


const map<string, string>&
AnalysisResults::getAltResidualOutputs() const
{
  return  my_alt_residual_outputs;
}


void
AnalysisResults::addPValues(const ModelResults& results)
{
  if(my_config.getAltCovCount(results.name) == 1)
  {
    ComparisonInfo cinfo;
      
    cinfo.model_name = results.name;
    cinfo.baseline_converged = results.maxfun_results_null.getConverged();
  
    if(! my_config.getScoreTestOnly())
    {                
      cinfo.converged  = results.maxfun_results_alt.getConverged() && results.maxfun_results_null.getConverged();
      cinfo.intercept  = results.maxfun_results_alt.getParameterMgr().getParameter("Intercept", "Intercept");

      cinfo.covs.clear();
      
      MAXFUN::ParameterConstIterator p     = results.maxfun_results_alt.getParameterMgr().getParamBegin("Covariates");
      MAXFUN::ParameterConstIterator p_end = results.maxfun_results_alt.getParameterMgr().getParamEnd("Covariates");
      for(; p != p_end; ++p)
      {
        cinfo.covs.push_back(*p);
      }
      
      cinfo.same_variance_components = results.sameVarianceComponents();
      if(cinfo.same_variance_components)
      {
        cinfo.lrt_pvalue = MAXFUN::JointTest(results.maxfun_results_null, results.maxfun_results_alt).getPValue();
      }
    }
    else
    {
      cinfo.intercept  = results.score_test_parameter_mgr.getParameter("Intercept", "Intercept");

      cinfo.covs.clear();
      
      MAXFUN::ParameterConstIterator p     = results.score_test_parameter_mgr.getParamBegin("Covariates");
      MAXFUN::ParameterConstIterator p_end = results.score_test_parameter_mgr.getParamEnd("Covariates");
      for(; p != p_end; ++p)
      {
        cinfo.covs.push_back(*p);
      }
    }
    
    cinfo.score_pvalue =  chdtrc(1, results.score_statistic);
    
    my_cinfos.push_back(cinfo);
  }
}

    
// Returns TRUE is r1 should be listed BEFORE r2, for the purpose of generating
// the summary comparison table (generateSummaryOutput).
struct SortCompTableRow
{
  bool  operator()(const OUTPUT::TableRow& r1, const OUTPUT::TableRow& r2) const
  {
    double  r1_lrt  = r1.getVector().getAbs<OUTPUT::Double>(r1.getVector().size() - 1).toVal();
    double  r2_lrt  = r2.getVector().getAbs<OUTPUT::Double>(r2.getVector().size() - 1).toVal();
                
    return  r1_lrt < r2_lrt;
  }
};


void
AnalysisResults::generateDetailedOutput(ofstream& out)
{
  OUTPUT::Section s("Results");
  
  string  pvalue_note = "P-values for variance components and proportions are one-sided.  All others are two-sided.\n";
  pvalue_note        += "The null hypothesis for transformation parameters is 'no transformation.'\n\n";
      
  s << OUTPUT::NamedString("Note", pvalue_note);
  
  s << summarizeConfig()
    << my_config.getTransConfig().summarize_as_table();
    
  if(my_config.getTransConfig().get_type() != MFSUBMODELS::Transformation::NONE)
  {
    string  transform_diff = "  Apply transformation to ";

    transform_diff += (my_config.getTransformBothSides() ? "each side of the regression equation separately." :       
                                                           "the difference between the two sides of the regression equation.");      
    transform_diff += "\n";
                                                               
    s << OUTPUT::NamedString("", transform_diff);
  }
    
  for(size_t i = 0; i < my_detailed_sections.size(); ++i)
  {
    s << my_detailed_sections[i];
  }
  
  out << s;
}


void
AnalysisResults::generateSummaryOutput(ofstream& out)
{
  OUTPUT::Section s("Results");
  
  s << summarizeConfig()
    << my_config.getTransConfig().summarize_as_table();
    
  vector<string>  no_lrts;      // Names of models w/o likelihood ratio test results.   
  s << generateSummaryTable(no_lrts, false);
    
  if(! no_lrts.empty())
  {
    string  note;
       
    note += "No likelihood ratio test was done for the following models because different\n";
    note += "numbers of variance components were estimated under the null and the alternative.\n";
    note += "You may wish to rerun the analysis with the extra variance components eliminated.\n";

    s << OUTPUT::NamedString("Note", note);
    
    size_t  model_count = no_lrts.size();
    for(size_t m = 0; m < model_count; ++m)
    {
      s << OUTPUT::NamedString("Model", no_lrts[m] + "\n");
    }
  }
    
  if(! my_config.getOmitCompleteSummary())
  {
    s << generateSummaryTable(no_lrts, true);
  }
  
  string  note = "See .det file for more information about each model\n";
  
  s << OUTPUT::NamedString("Note", note);
  
  out << s;
}


OUTPUT::Table 
AnalysisResults::generateSummaryTable(vector<string>& no_lrts, bool complete) const
{
  list<ComparisonInfo>  results(my_cinfos.begin(), my_cinfos.end());
  
  no_lrts.clear();
  removeInvalidResults(results, no_lrts);
  
  if(! complete)
  {
    orderResults(results);
    filterResults(results);
  }
  
  OUTPUT::Table t("Summary table");
  
  if(complete)
  {
    t.setTitle("Complete summary table");
  }
  
  if(results.empty())
  {
    return  OUTPUT::Table();
  }

  t << OUTPUT::NamedString("Note", "Models with an asterisk (*) may not have converged.");
  t << OUTPUT::TableColumn("Model")
    << OUTPUT::TableColumn("Intercept");
  
  Configuration::AllNullCovConstIter  cov_iter     = my_config.allNullCovariateBegin();
  Configuration::AllNullCovConstIter  cov_end_iter = my_config.allNullCovariateEnd();
  for(; cov_iter != cov_end_iter; ++cov_iter)
  {
    t << OUTPUT::TableColumn(cov_iter->cfg.param_name);
  }

  // Add columns for the other stuff:
  t << OUTPUT::TableColumn("Test cov.") 
    << OUTPUT::TableColumn("Estimate") 
    << OUTPUT::TableColumn("Std. err.")
    << OUTPUT::TableColumn("P-value (Wald)")
    << OUTPUT::TableColumn("P-value (LRT)");
  //  << OUTPUT::TableColumn("P-value (Score)");

  // Add a row for each test.
  list<ComparisonInfo>::const_iterator comparison_iter     = results.begin();
  list<ComparisonInfo>::const_iterator comparison_end_iter = results.end();
  for(; comparison_iter != comparison_end_iter; ++comparison_iter)
  {
    // Get test covariate
    const MAXFUN::Parameter&  param = comparison_iter->covs[comparison_iter->covs.size() - 1];

    OUTPUT::TableRow r;
      
    // Add model name:
    r << (comparison_iter->model_name + string(comparison_iter->converged ? "" : "*"));

    // Add intercept estimate:
    r << comparison_iter->intercept.getFinalEstimate();

    // Add baseline covariates:
    for(size_t i = 0; i < comparison_iter->covs.size() - 1; ++i)
    {
      r << comparison_iter->covs[i].getFinalEstimate();
    }

    // Add test covariate:
    r << param.getName() << param.getFinalEstimate();
    
    if(param.isStdErrorAvailable())
      r << param.getStdError();
    else
      r << "Unavailable";
      
    double  wald_pvalue = comparison_iter->waldPvalue();
    if(! isnan(wald_pvalue))
    {
      r << wald_pvalue;
    }
    else
    {
      r << "Unavailable";
    }
      
    // Add lrt pvalue:
    r << comparison_iter->lrt_pvalue;
    
    /*
    double  score_pvalue = comparison_iter->score_pvalue;
    if(! isnan(score_pvalue))
    {
      r << score_pvalue;
    }
    else
    {
      r << "Unavailable";
    }
    */
    
    t << r;
  } 
  
  return  t;    
}

void 
AnalysisResults::generateTsvOutput(ofstream& tsv_file)
{
  list<ComparisonInfo>  results(my_cinfos.begin(), my_cinfos.end());

  tsv_file << "Model\tP-value (Wald)\tP-value (LRT)" << endl;
  
  list<ComparisonInfo>::const_iterator comparison_iter     = results.begin();
  list<ComparisonInfo>::const_iterator comparison_end_iter = results.end();
  for(; comparison_iter != comparison_end_iter; ++comparison_iter)
  {
    tsv_file << (comparison_iter->model_name + string(comparison_iter->converged ? "" : "*") + "\t");
    tsv_file << comparison_iter->waldPvalue() << "\t";
    tsv_file << comparison_iter->lrt_pvalue << endl;
  }  
}


// - Comparison is not valid if null and alternate cases use different variance components.
//
void
AnalysisResults::removeInvalidResults(list<ComparisonInfo>& results, vector<string>& no_lrts) const
{
  list<ComparisonInfo>::iterator comparison_iter     = results.begin();
  list<ComparisonInfo>::iterator comparison_end_iter = results.end();
  for(; comparison_iter != comparison_end_iter; ++comparison_iter)
  {
    if(! comparison_iter->same_variance_components)
    {
      no_lrts.push_back(comparison_iter->model_name);
      comparison_iter = results.erase(comparison_iter);
    }
  }
}

void
AnalysisResults::orderResults(list<ComparisonInfo>& results) const
{
  switch(my_config.getDisplayOrder())
  {
    case Configuration::AS_INPUT:
      break;
      
    case Configuration::BY_LRT:
      results.sort(LRTLess());
      break;
      
    case Configuration::BY_WALD:
      results.sort(WaldLess());
      break;
      
    case Configuration::BY_LARGER_PVALUE:
      results.sort(LargerLess());
      break;
      
    case Configuration::BY_PVALUE_RATIO:
      results.sort(RatioLess());
      break;      
      
    default:
      assert(false);
  }
}

void
AnalysisResults::filterResults(list<ComparisonInfo>& results) const
{
  if(my_config.getDisplayAll())
  {
    return;                           // No filtering required.
  }
  
  if(my_config.getFilterByLRT())
  {
    results.remove_if(LRTFilter(my_config.getLRTFilterThreshold()));
  }
  
  if(my_config.getFilterByWald())
  {
    results.remove_if(WaldFilter(my_config.getWaldFilterThreshold()));
  }
  
  if(my_config.getLimitNumberDisplayed())
  {
    size_t  display_number = my_config.getDisplayNumber();
    if(results.size() > display_number)
    {
      results.resize(display_number);
    }
  }
}
    
//=============================================================================
//  addModelDetailedSection()
//=============================================================================
void 
AnalysisResults::addModelDetailedSection(const ModelResults& results)
{
  OUTPUT::Section s("Model '" + results.name + "'");
    
  s << results.sample_summary;

  if(results.covariate_infos.size())
  {
    OUTPUT::Table cov_table("Covariates");
   
    cov_table << OUTPUT::TableColumn("Name")
              << OUTPUT::TableColumn("Mean")
              << OUTPUT::TableColumn("Std. dev.")
              << OUTPUT::TableColumn("Min")
              << OUTPUT::TableColumn("Max");
              
    for(std::vector<SAMPLING::Field::SummaryInfo>::const_iterator i = results.covariate_infos.begin(); i != results.covariate_infos.end(); ++i)
    {
      cov_table << (OUTPUT::TableRow() << i->name << i->mean << i->stdev << i->min << i->max);
    }
    
    s << cov_table;
  }

  s << OutputFormatter::convertEstimates(results.maxfun_results_null, results.residuals_scaling_factor, true);
  s << OutputFormatter::convertMatrix(results.maxfun_results_null, results.residuals_scaling_factor);

  if(my_config.getModel(results.name).null_residuals)
  {
    ostringstream  residual_misc;
    
    residual_misc << "Mean of the residuals:   " << results.residuals_null.mean() << "\n"
                  << "Standard deviation of the residuals:   " << results.residuals_null.std_dev() << "\n" << endl;
                  
    if(my_config.getDependentTraitType() == BINARY)
    {
      residual_misc << "Mean of the affected residuals:   " << results.affected_mean_residuals_null << "\n"
                    << "Mean of the unaffected residuals:  " << results.unaffected_mean_residuals_null << endl;
    }
                  
    s << OUTPUT::NamedString("", residual_misc.str());
  }
    
  if(my_config.getAltCovCount(results.name) != 0 &&
               (! my_config.getScoreTestOnly()))
  {
    s << OutputFormatter::convertEstimates(results.maxfun_results_alt, results.residuals_scaling_factor, true);
    s << OutputFormatter::convertMatrix(results.maxfun_results_alt, results.residuals_scaling_factor);    
    
    if(my_config.getModel(results.name).test_residuals)
    {
      ostringstream  residual_misc;

      residual_misc << "Mean of the residuals:   " << results.residuals_alt.mean() << "\n"
                    << "Standard deviation of the residuals:   " << results.residuals_alt.std_dev() << "\n" << endl;
                    
    if(my_config.getDependentTraitType() == BINARY)
    {
      residual_misc << "Mean of the affected residuals:   " << results.affected_mean_residuals_alt << "\n"
                    << "Mean of the unaffected residuals:  " << results.unaffected_mean_residuals_alt << endl;
    }                    
                    
      s << OUTPUT::NamedString("", residual_misc.str());
    }
  
    s << MAXFUN::JointTest(results.maxfun_results_null, results.maxfun_results_alt).summarizeAsTable();
  }

  my_detailed_sections.push_back(s);
}


const string  
AnalysisResults::determine_sex(const FPED::Member& member) const
{
  MPED::SexCode  sex = member.get_effective_sex();
  
  switch(sex)
  {
    case MPED::SEX_MALE:
      return  member.multipedigree()->info().sex_code_male();
    case MPED::SEX_FEMALE:
      return  member.multipedigree()->info().sex_code_female();
    case MPED::SEX_MISSING:
      return  member.multipedigree()->info().sex_code_unknown();
    default:
      assert(false);
  }
}

void  
AnalysisResults::write_residual(ostringstream& out, double residual, double ind_residual, const Sampledata& sampledata, size_t member_id) const
{
   if(sampledata.isValid(member_id))
   {
     const FPED::Member&  member = sampledata.getIndividual(member_id);
     const string  sex = determine_sex(member);
     const string&  p1 = member.parent1() ? member.parent1()->name() : member.multipedigree()->info().individual_missing_code();
     const string&  p2 = member.parent2() ? member.parent2()->name() : member.multipedigree()->info().individual_missing_code(); 
                                   
     out << member.pedigree()->name() << "\t" << member.name() << "\t" 
         << p1 << "\t" << p2 << "\t" << sex << "\t" << residual << "\t" << ind_residual << "\n";
   }
}

void 
AnalysisResults::addModelResidualOutput(const ModelResults& results, const Sampledata& sampledata)
{
  // - Alternative residuals.
  //
  ostringstream  alt_output;
    
  if(results.residuals_alt.size() && 
     (! my_config.getScoreTestOnly()) )
  {
    assert(results.residuals_alt.size() == results.ind_residuals_alt.size());
    alt_output << "Pedigree_ID\tIndividual_ID\tParent1\tParent2\tSex\tResidual\tIndependent_Residuale\n";    
                
    
    size_t  alt_residual_count = results.residuals_alt.size();
    for(size_t id = 0; id < alt_residual_count; ++id)
    {
      double  residual = results.residuals_alt.standardized_residual(id);
      double  ind_residual = results.ind_residuals_alt.standardized_residual(id);
      write_residual(alt_output, residual, ind_residual, sampledata, id);
    }
  }
  
  my_alt_residual_outputs.insert(make_pair(results.name, alt_output.str()));
  
  // - Null residuals.
  //
  ostringstream  null_output;
  
  if(results.residuals_null.size())
  {
    assert(results.residuals_null.size() == results.ind_residuals_null.size());  
    null_output << "Pedigree_ID\tIndividual_ID\tParent1\tParent2\tSex\tResidual\tIndependent_Residual\n";    
                
    size_t  null_residual_count = results.residuals_null.size();
    for(size_t id = 0; id < null_residual_count; ++id)
    {
      double  residual = results.residuals_null.standardized_residual(id);
      double  ind_residual = results.ind_residuals_null.standardized_residual(id);
      write_residual(null_output, residual, ind_residual, sampledata, id);
    }
  }
  
  my_null_residual_outputs.insert(make_pair(results.name, null_output.str()));
}


//=============================================================================
//  summarizeConfig()
//=============================================================================
OUTPUT::Table 
AnalysisResults::summarizeConfig() const
{
  OUTPUT::Table t("Analysis description");

  t << (OUTPUT::TableRow() << "Title"           << my_config.getTitle())
    << (OUTPUT::TableRow() << "Primary Trait"   << (my_config.getPrimaryTraitName() +
                              std::string(my_config.getDependentTraitType() == BINARY ? " (Binary)" : " (Quantitative)")))
    << (OUTPUT::TableRow() << "Allow averaging" << (my_config.getAllowAveraging() ? "Enabled" : "Disabled"))
  /*  << (OUTPUT::TableRow() << "Score test only" << (my_config.getScoreTestOnly() ? "Enabled" : "Disabled")) */
    << (OUTPUT::TableRow() << "Omit complete summary table" << (my_config.getOmitCompleteSummary() ? "Enabled" : "Disabled"));

  if(my_config.getDependentTraitType() == BINARY)
  {
    string  binary_warning  = "This output contains results for a binary trait ";
    binary_warning         += "using an algorithm still in Beta release form. Please verify results before reporting.";
    t << OUTPUT::NamedString("Note", binary_warning);
  }
 
  summarizeDisplayOptions(t);
  
  return t;
}

void
AnalysisResults::summarizeDisplayOptions(OUTPUT::Table& options) const
{
  options << (OUTPUT::TableRow() << "Summary table display order" << displayOrder2String(my_config.getDisplayOrder()));
  
  if(my_config.getDisplayAll())
  {
    options << (OUTPUT::TableRow() << "Summary filters" << "None");
    return;
  }
  
  if(my_config.getFilterByLRT())
  {
    ostringstream  filter_desc;
    
    filter_desc << "Show only results where LRT p-value < " << my_config.getLRTFilterThreshold();
    options << (OUTPUT::TableRow() << "Summary filter" << filter_desc.str());
  }
  
  if(my_config.getFilterByWald())
  {
    ostringstream  filter_desc;
    
    filter_desc << "Show only results where Wald p-value < " << my_config.getWaldFilterThreshold();
    options << (OUTPUT::TableRow() << "Summary filter" << filter_desc.str());
  }
  
  if(my_config.getLimitNumberDisplayed())
  {
    ostringstream  filter_desc;
    
    filter_desc << "Show no more than " << my_config.getDisplayNumber() << " results";
    options << (OUTPUT::TableRow() << "Summary filter" << filter_desc.str());
  }
}

} 
} 
