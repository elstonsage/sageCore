//* File:      analysis_out.cpp                                              *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0 Initial Implementation.               yjs May. 01 *
//*                                                                          *
//* Notes:     Implements analysis_out classes for SEGREG.                   *
//*                                                                          *
//* Copyright (c) 2001 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include <iostream>
#include <iomanip>
#include <map>
#include <cmath>
#include <cassert>
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "output/Output.h"
#include "LSF/parse_ops.h"
#include "numerics/functions.h"
#include "numerics/print_util.h"
#include "rped/rped.h"
#include "segreg/analysis_out.h"

namespace SAGE
{
namespace SEGREG
{

bool analysis_output::skip_poly_locus = false; // due to JA
double analysis_output::deriv_sum_sq =  0.0; // due to JA

double calculate_akaike(const MAXFUN::Results& mf)
{
  double akaike = -2.0 * mf.getFinalFunctionValue() + 2 * mf.getNumOfVaryingIndependentParams();
  
  return akaike;
}   

void analysis_output::output_likelihood_table_header
(ostream& out, model& test_model) //const primary_analysis_results& test)
{
  std::string title;

  if(test_model.get_type_missing())
    title = "Commingling Analysis";
  else
    title = "Segregation Analysis";

  out<<"=================================================================================="<<endl;
  out<<"                    Likelihoods for "<<title                                       <<endl;
  out<<endl;
  out<<"Option                       LN(Likelihood)  -2 LN(Likelihood)  Akaike's AIC score"<<endl;
  out<<"----------------------------------------------------------------------------------"<<endl;
}

void analysis_output::output_likelihood_table_results
(ostream& out, const primary_analysis_results& test)
{
  std::string description;
  
  const model& md = test.get_final_model();
 
  // Determine if we're doing a segregation or a commingling analysis to pick the
  // right option description.
  if( md.get_type_missing() )
    description = test.get_final_model().type_dependent_sub_model().option_description();
  else
    description = test.get_final_model().transm_sub_model.option_description();
  
  const MAXFUN::Results& results = test.get_maxfun_results();
     
  out << left << setfill(' ') << setw(29)
      << description;
      
  // Determine if the results are valid.  If they are, print them.  If not, don't
  if(test.is_valid())
  {
    out << setfill(' ')
        << fp(results.getFinalFunctionValue(),13,4)
        << fp(-2*results.getFinalFunctionValue(),19,4)
        << fp(calculate_akaike(results),20,4)<<endl;
  }
  else
  {
    out << " ----- No Valid Maximization -----" << endl;
  }
}
void analysis_output::output_likelihood_table_footer
(ostream& out)
{
  out<<"=================================================================================="<<endl;
}

/////////////////////////////////////////////////////////////////////////////////////////////////

void analysis_output::output_model            
    (ostream& out, const model& mod)
{
  out << "=====================================================================================================" << endl;

  out << "  SEGREG Analysis for Trait : ";
  out << mod.get_primary_trait() << endl;

  out << "=====================================================================================================" << endl;

  out << "  # Model Specification  " << endl << endl
      << "    Model Class ";
  switch( mod.get_model_class() )
  {
    case model_A:
         out << "A"; break;
    case model_D:
         out << "D"; break;
    case model_MLM:
         out << "MLM"; break;
    case model_FPMM:
         out << "FPMM"; break;
    case model_INVALID:
         out << "INVALID"; break;
  }
  out << endl;

// back to original code
  primary_type pt = mod.get_primary_trait_type();

  if(pt == pt_CONTINUOUS || pt == pt_ONSET)
  {
    output_model_continuous(out, mod);
  }

  if(pt == pt_BINARY || pt == pt_ONSET)
  {
    output_model_binary(out, mod);
  }

  // FPMM or Residuals Option

  if(mod.get_model_class() == model_FPMM)
  {
       if ( ! analysis_output::skip_poly_locus) // number of polygenic loci set > 0 by user
           {
            out << "    " << setw(28) << mod.fpmm_sub_model.name()
                 << " : "    << mod.fpmm_sub_model.option_description()  << endl;
           } // zero polygenic loci set by user
          if ( (analysis_output::skip_poly_locus))
           {
             out << "    " << setw(28) << mod.fpmm_sub_model.name()
                 << " : "    <<"0 polygenic loci" <<endl;
           }
  }
  else
  {
    cerrorstream resid(out);

    resid.prefix("    Residual correlations        : ");

    resid << mod.resid_sub_model.option_description();

    resid << endl;
  }

  if(pt == pt_CONTINUOUS || pt == pt_ONSET)
  {
    // Transformation Option

    out << "    " << setw(28) << mod.transf_sub_model.name()
        << " : "    << mod.transf_sub_model.option_description()  << endl;
  }

  // Freq Option

  out << "    " << setw(28) << mod.freq_sub_model.name()
      << " : "    << mod.freq_sub_model.option_description()  << endl;
   

  // Transmission Option

  out << "    " << setw(28) << mod.transm_sub_model.name()
      << " : "    << mod.transm_sub_model.option_description()  << endl;

  // Ascertainment option

  if(mod.ascer_sub_model.s_option() != ascertainment_sub_model::none)
  {
    const ascertainment_sub_model& ascer = mod.ascer_sub_model;

    out << "    " << setw(28) << ascer.name() << " : ";

    if(ascer.option_description().size() <= 30)
    {
      out << ascer.option_description() << endl;
    }
    else
    {
      out                    << ascer.s_option_description() << ","  << endl;
      out << setw(36) << " " << ascer.v_option_description() << endl;
    }
  }

  // Output Prevalence Option

  if(pt == pt_BINARY || pt == pt_ONSET)
  {
    const prevalence_sub_model& psm = mod.prev_sub_model;

    if(psm.get_constraint_count() > 0 || psm.get_estimate_count() > 0)
    {
      out << "    " << setw(28) << "Prevalence options" << " : " << endl;

      for(size_t i = 0; i < psm.get_constraint_count(); ++i)
      {
        out << "        " << setw(22) << "Constraint " << setw(2) << i+1 << " : ";

        out <<                      "Sample Size     : " << psm.get_constraint_sample_size(i) << endl;
        out << setw(35) << " "   << "Number Affected : " << psm.get_constraint_number_affected(i) << endl;

        if(mod.get_primary_trait_type() == pt_ONSET)
          out << setw(35) << ' ' << "Age:              " << psm.get_constraint_age(i) << endl;

        for(size_t cov = 0; cov < psm.get_constraint_susc_covariate_count(i); ++cov)
          out << setw(35) << " " << setw(15)
              << psm.get_constraint_susc_covariate_name(i, cov).substr(0,15) << " : "
              << psm.get_constraint_susc_covariate_value(i, cov) << endl;
        for(size_t cov = 0; cov < psm.get_constraint_mean_covariate_count(i); ++cov)
          out << setw(35) << " " << setw(15)
              << psm.get_constraint_mean_covariate_name(i, cov).substr(0,15) << " : "
              << psm.get_constraint_mean_covariate_value(i, cov) << endl;
        for(size_t cov = 0; cov < psm.get_constraint_var_covariate_count(i); ++cov)
          out << setw(35) << " " << setw(15)
              << psm.get_constraint_var_covariate_name(i, cov).substr(0,15) << " : "
              << psm.get_constraint_var_covariate_value(i, cov) << endl;
      }
      for(size_t i = 0; i < psm.get_estimate_count(); ++i)
      {
        out << "        " << setw(22) << "Estimate " << setw(2) << i+1 << " :" << endl;

        if(mod.get_primary_trait_type() == pt_ONSET)
          out << setw(35) << ' ' << "Age:              " << psm.get_estimate_age(i) << endl;

        for(size_t cov = 0; cov < psm.get_estimate_susc_covariate_count(i); ++cov)
          out << setw(35) << " " << setw(15)
              << psm.get_estimate_susc_covariate_name(i, cov).substr(0,15) << " : "
              << psm.get_estimate_susc_covariate_value(i, cov) << endl;
        for(size_t cov = 0; cov < psm.get_estimate_mean_covariate_count(i); ++cov)
          out << setw(35) << " " << setw(15)
              << psm.get_estimate_mean_covariate_name(i, cov).substr(0,15) << " : "
              << psm.get_estimate_mean_covariate_value(i, cov) << endl;
        for(size_t cov = 0; cov < psm.get_estimate_var_covariate_count(i); ++cov)
          out << setw(35) << " " << setw(15)
              << psm.get_estimate_var_covariate_name(i, cov).substr(0,15) << " : "
              << psm.get_estimate_var_covariate_value(i, cov) << endl;
      }
    }
  }

  out << endl;

  if(mod.comp_trait_sub_model.get_covariate_count() ||
     mod.mean_cov_sub_model  .get_covariate_count() ||
     mod.var_cov_sub_model   .get_covariate_count() ||
     mod.susc_cov_sub_model  .get_covariate_count())
  {
    out << endl;

    out << "    NOTE:  All covariates are centered at mean 0.0 for analysis.  The centering\n"
        << "           means listed are the means of all non-missing values of the covariate,\n"
        << "           regardless of the missingness of other variables used in the analysis.\n"
        << "           It follows that, unless the missingness is the same for the trait and\n"
        << "           all covariates used, these means may not be the appropriate centering\n"
        << "           means for the results below.\n" << endl;
  }
}


void analysis_output::output_model_continuous
    (ostream& out, const model& mod)
{
  if(mod.comp_trait_sub_model.get_covariate_count())
  {
    out << "    " << setw(28) << mod.comp_trait_sub_model.name()
        << " :" << endl;

    mod.comp_trait_sub_model.print_covariate_table(out, 10);
  }

  // Mean Option

  out << "    " << setw(28) << mod.mean_sub_model.name()
      << " : "    << mod.mean_sub_model.option_description()  << endl;

  // Mean Covariate Options.

  if(mod.mean_cov_sub_model.get_covariate_count())
  {
    out << "    " << setw(28) << mod.mean_cov_sub_model.name()
        << " :" << endl;

    mod.mean_cov_sub_model.print_covariate_table(out, 10);
  }

  // Var Option

  out << "    " << setw(28) << mod.var_sub_model.name()
      << " : "    << mod.var_sub_model.option_description()  << endl;

  // Var Covariate Options.

  if(mod.var_cov_sub_model.get_covariate_count())
  {
    out << "    " << setw(28) << mod.var_cov_sub_model.name()
        << " :" << endl;

    mod.var_cov_sub_model.print_covariate_table(out, 10);
  }
}


void analysis_output::output_model_binary
    (ostream& out, const model& mod)
{
  // Susc Option

  out << "    " << setw(28) << mod.susc_sub_model.name()
      << " : "    << mod.susc_sub_model.option_description()  << endl;

  // Susc Covariates

  if(mod.susc_cov_sub_model.get_covariate_count())
  {
    out << "    " << setw(28) << mod.susc_cov_sub_model.name()
        << " :" << endl;

    mod.susc_cov_sub_model.print_covariate_table(out, 10);
  }

}


void analysis_output::output_initial_estimates
    (ostream& out, const primary_analysis_results& test)
{
  const MAXFUN::ParameterMgr& fin_params = test.get_maxfun_results().getParameterMgr();
     
  out << "  # Initial Estimates : " << endl;
  out << endl;
  out << "  Parameter           Lower Bound.    Initial Est.    Upper Bound.   Status" << endl;
  out << "-----------------------------------------------------------------------------------------------------" << endl;
  for( int i = 0; i < fin_params.getParamCount(); ++i )
  {
    const MAXFUN::Parameter& param = fin_params.getParameter(i);
    
    out << left    << "  " << setw(19) << param.getName().substr(0,19);
    out << right   << "  " << fp(param.getLowerBound(), 11, 8);
    out << setw(5) << "  " << fp(param.getInitialEstimate(), 11, 8);
    out << setw(5) << "  " << fp(param.getUpperBound(), 11, 8);
    out << setw(5) << "  " << MAXFUN::ParamTypeEnum2str(param.getInitialType());
    out << endl;
  }
  out << "-----------------------------------------------------------------------------------------------------" << endl;
  out << endl;
}


void analysis_output::output_final_estimates  
    (ostream& out, const primary_analysis_results& test)
{
  const MAXFUN::ParameterMgr& fin_params = test.get_maxfun_results().getParameterMgr();

  analysis_output::deriv_sum_sq = 0;
  // Determine if std errors are available
  bool stde = false;

  for( int i = 0; !stde && i < fin_params.getParamCount(); ++i )
    if(std::abs(fin_params.getParameter(i).getStdError()) > 1e-12)
      stde = true;

  out << "  # Final Estimates : " << endl;
  out << endl;
  out << "  Parameter           Parameter Est.  Standard Err.   First Deriv.   Status" << endl;
  out << "-----------------------------------------------------------------------------------------------------" << endl;

  bool has_infinite_std_err = false;

  for( int i = 0; i < fin_params.getParamCount(); ++i )
  {
    const MAXFUN::Parameter& param = fin_params.getParameter(i);
    
   if ( (param.getName().substr(0,5) == "polyg") && (analysis_output::skip_poly_locus) )continue; // addition by JA

    out << left  << setw(2) << " " << setw(15) << param.getName().substr(0,15);
    out << right << setw(5) << " " << fp(param.getFinalEstimate(), 11, 8);

    out << setw(5) << " ";

    if(stde)
    {
      double std_err = param.getStdError();
      if( std_err < 0. )
      {
        out << " INF.      ";

        has_infinite_std_err = true;
      }
      else if(std_err < 1e-12)          // Should only happen if the
      {                                 // parameter is fixed or at a bound.
        out << "           ";
      }
      else
      {
        out << fp(param.getStdError(), 11, 8);
      }
    }
    else
      out << "-----------";

    // Print out the first derivative for all independent variables
    
    if(param.getInitialType() < 3)
    {
      double this_deriv = param.getDeriv();
      analysis_output::deriv_sum_sq += this_deriv*this_deriv;
      out << setw(5) << " " << fp(param.getDeriv(), 11, 8);
    }
    else
    {
      out << setw(16) << " ";
    }
    
    // Print out the final type
    
    out << setw(5) << " " << MAXFUN::ParamTypeEnum2str(param.getFinalType());;
    out << endl;
  }
  out << "-----------------------------------------------------------------------------------------------------" << endl;

  out << endl;


  if(has_infinite_std_err)
  {
    cerrorstream err(out);
    
    err.prefix("%%SEGREG-%p: ");

    err << priority(information) << "Infinite standard errors indicate the "
        << "detection of a flat or near-flat likelihood for certain parameters.  "
        << "The standard errors of such parameters cannot be computed.  "
        << "If two or more standard errors are equal, the corresponding variables "
        << "may be highly correlated, in which case neither the estimates nor "
        << "their standard errors are individually meaningful."
        << endl;
  }
  if(!stde)
  {
    cerrorstream err(out);
    
    err.prefix("%%SEGREG-%p: ");

    err << priority(information) << "The likelihood surface is flat and standard "
        << "errors cannot be computed." << endl;
  }
}

void analysis_output::output_likelihoods      
    (ostream& out, const primary_analysis_results& test)
{
  const MAXFUN::Results& results = test.get_maxfun_results();

  out << endl;
  out << "        LN(Likelihood) : " << results.getFinalFunctionValue()    << endl;
  out << "     -2 LN(Likelihood) : " << -2*results.getFinalFunctionValue() << endl;
  out << "    Akaike's AIC score : " << calculate_akaike(results) << endl;
  out << "=====================================================================================================" << endl;

  out << endl;

  output_final_error(out, test);
}

void analysis_output::output_final_error
    (ostream& out, const primary_analysis_results& test)
{
  const MAXFUN::Results& results = test.get_maxfun_results();

  size_t code = results.getExitFlag();

//  if(code < 6) return; // original code
  if(code < 4) return;

  cerrorstream err(out);
  
  err.prefix("%%SEGREG-%p: ");

  switch(code)
  {
//    case 6 : //original code 
    case 5: // new addition 
      err << priority(warning) << "Accumulation of round off errors or boundary problem " << endl
          << "Check results carefully " << endl;
      break;
    case 6: // new addition
      err << priority(warning) << " Problems with variable metric methods " << endl
          << "Results may not be fully maximized " << endl;
           break;
    case 7 :
      err << priority(warning) << "Maximization procedure did not complete cleanly. "
          << "Results may not be totally maximized.  Check results carefully." << endl;
      break;
    case 8 :
      if (analysis_output::deriv_sum_sq > 1.E-03) { // activate warning only if derivatives are not small (due to JA)
      err << priority(error) << "Likelihood may not be maximized because of model constraints on parameters. "
          << endl << "Estimating fewer parameters may help." << endl;
      }
      break;
    case 9 :
      err << priority(information) << "Variance-covariance matrix could not be computed."
          << endl;
      break;
    case 10:
      err << priority(warning) << "All independent parameters converged to bounds."
          << endl;
      break;
    case 11:
      err << priority(warning) << "Too many parameters to estimate using derivatives. "
          << "Reduce the number of parameters." << endl;
          break; // new addition
     case 12: // new addition
       err << priority(warning) << " Initial parameter estimates not in domain of function to be " 
           << " maximized " << endl;
  }

  out << endl << endl;
}


void analysis_output::output_vc_matrix        
    (ostream& out, const primary_analysis_results& test)
{
  const MAXFUN::Results& results = test.get_maxfun_results();

  if(results.getExitFlag() > 7)
  {
    out << "    Variance-covariance matrix not available due to computational problems." << endl 
        << "    See below." << endl;
    return;
  }

  // If VC Matrix Unavailable, let them know.
  if(results.getCovMatrixStatus() > 2)
  {
    out << "    Variance-covariance matrix not available due to computational problems." << endl;
    return;
  }

  const MAXFUN::ParameterMgr params = results.getParameterMgr();
  
  out << endl;
  out << "  # Variance-Covariance Matrix : " << endl;
  out << endl;

  out << setw(13) << " " << "|" << left;
  for( int i = 0; i < params.getParamCount(); ++i )
    if(params.getParameter(i).getInitialType() != 4)
      out << setw(5) << " " << setw(11) << params.getParameter(i).getName().substr(0,11);
  out << endl;

  out << "--------------";
  for( int i = 0; i < params.getParamCount(); ++i )
    if(params.getParameter(i).getInitialType() != 4)
      out << "----------------";
  out << endl;

  
  size_t ip = 0;
  for( int i = 0; i < params.getParamCount(); ++i )
  {
    if(params.getParameter(i).getInitialType() != 4)
    {
      out << setw(2) << " " << left << setw(11) << params.getParameter(i).getName().substr(0,11) << "|";

      size_t jp = 0;      

      for( int j = 0; j < params.getParamCount(); ++j )
        if(params.getParameter(j).getInitialType() != 4)
        {
          out << setw(5) << " " << right << fp(results.getCovarianceMatrix().getCovariance(ip,jp), 11, 8);

          ++jp;
        }

      ++ip;

      out << endl;
    }
  }

  out << "--------------";
  for( int i = 0; i < params.getParamCount(); ++i )
    if(params.getParameter(i).getInitialType() != 4)
      out << "----------------";
  out << endl << endl;

  if(results.getCovMatrixStatus() == 1){
   out << " Likelihood Surface Flat: variance matrix near singular and " << endl;
   out << " may be affected by rounding error " << endl;
   out << " Reducing the number of parameters to be estimated is recommended " << endl;
  }
}


////////////////////////////////////////////////////////////////////////////
//    Implementation of analysis_result_file (non-Inline)
////////////////////////////////////////////////////////////////////////////

void
analysis_viewer::print_header(const primary_analysis_results& test)
{
  ostream& out = output_stream();

  analysis_output::output_model(out, test.get_final_model());

  // Only print this stuff if the model was actually maximized
  if(test.is_valid())
  {
    out << "  # Number of constituent pedigrees   : " << test.get_subpedigree_count() << endl;
    out << "  # Number of unconnected individuals : " << test.get_unconnected_count() << endl;
    out << endl;
  }
}

void
analysis_result_file::print_results(const primary_analysis_results& test)
{
  if(test.is_valid())
  {
    analysis_output::output_final_estimates(output_stream(), test);
    if(test.get_final_model().get_model_class() == model_MLM )
      analysis_output::output_mlm_resid_corr_results(output_stream(), test);
  }
  else
  {
    output_stream() << "This model was unable to be maximized due to "
                    << "complications with the initial" << endl
                    << "estimates.  See the .inf file progress table "
                    << "for more information." << endl << endl;
  }
}

void
analysis_result_file::print_footer(const primary_analysis_results& test)
{
  if(test.is_valid())
  {
    analysis_output::output_likelihoods(output_stream(), test);
  }
}

void
analysis_detailed_file::print_results(const primary_analysis_results& test)
{
  if(test.is_valid())
  {
    analysis_output::output_final_estimates   (output_stream(), test);
    analysis_output::output_vc_matrix         (output_stream(), test);
    if(test.get_final_model().get_model_class() == model_MLM )
      analysis_output::output_mlm_resid_corr_results(output_stream(), test);
  
    if(my_debug)
    {
      output_stream() << "  MAXFUN Num. Function Eval: "
                      << test.get_maxfun_results().getIterations() << endl;
      output_stream() << "  MAXFUN Final return type:  "
                      << test.get_maxfun_results().getExitFlag() << endl;
    }
  }
  else
  {
    output_stream() << "This model was unable to be maximized due to "
                    << "complications with the initial" << endl
                    << "estimates.  See the .inf file progress table "
                    << "for more information." << endl << endl;
  }
}

void
analysis_detailed_file::print_footer(const primary_analysis_results& test)
{
  if(test.is_valid())
  {
    analysis_output::output_likelihoods(output_stream(), test);
  }
}

////////////////////////////////////////////////////////////////////////////
//    Implementation of analysis_intermediate_file (non-Inline)
////////////////////////////////////////////////////////////////////////////

void
analysis_intermediate_file::print_results(const primary_analysis_results& test)
{
  analysis_output::output_final_estimates(output_stream(), test);
  if(test.get_final_model().get_model_class() == model_MLM )
    analysis_output::output_mlm_resid_corr_results(output_stream(), test);  
}

void
analysis_intermediate_file::print_footer(const primary_analysis_results& test)
{
  ostream& out = output_stream();

  out << endl;
  out << "  Final Results    : " << endl;
  out << endl;
  out << "  Number of Function Evaluations : " << test.get_maxfun_results().getIterations() << endl;
  out << "  Last Return Code               : " << test.get_maxfun_results().getExitFlag()           << endl;

  analysis_output::output_likelihoods(output_stream(), test);
}

std::string spellout(string st)
{
  if(st == "FM")        return "father_mother";
  if(st == "FS")        return "father_offspring";
  if(st == "MS")        return "mother_offspring";
  if(st == "SS")        return "sib_sib";
  else                  return "no this relation";
}

void output_mlm_correlation
    (ostream& out,
     string corr_name,
     MlmResidCorrelationCalculator& mlm_corrs,
     residual_correlation_sub_model::corr ctype)
{
  out << left << setfill(' ') << setw(20) << corr_name;

  for(size_t j = 0; j < 3; ++j)
  {
    out << setfill(' ') << setw(20) 
        << fp(mlm_corrs.get_corr_ij_u(ctype,j),9,4);
  }
  out << setfill(' ') << setw(30) << mlm_corrs.get_corr_ij(ctype)
      << endl;
}

void analysis_output::output_mlm_resid_corr_results
(ostream& out, const primary_analysis_results& test)
{
  MlmResidCorrelationCalculator  mrc = test.get_mlm_resid_corr();

  //start at 10, end at 110
  out<<endl;
  out<<"-----------------------------------------------------------------------------------------------------"<<endl;
  out<<"                                        Residual Correlations"<<endl;
  out<<endl;
  //out<<setfill(' ')<<setw(20)
  out<<"     "
     <<setfill(' ')<<setw(20)
     <<"AA"
     <<setfill(' ')<<setw(20)
     <<"AB"
     <<setfill(' ')<<setw(20)
     <<"BB"
     <<setfill(' ')<<setw(25)
     <<"Overall"<<endl;
  out<<"-----------------------------------------------------------------------------------------------------"<<endl;

  output_mlm_correlation(out, "father_mother",    mrc, residual_correlation_sub_model::fm);
  
  if(test.get_final_model().resid_sub_model.option() == residual_correlation_sub_model::arb)
  {
    output_mlm_correlation(out, "parent_offspring", mrc, residual_correlation_sub_model::fs);
  }
  else
  {
    output_mlm_correlation(out, "father_offspring", mrc, residual_correlation_sub_model::fs);
    output_mlm_correlation(out, "mother_offspring", mrc, residual_correlation_sub_model::ms);
  }

  output_mlm_correlation(out, "sib_sib",          mrc, residual_correlation_sub_model::ss);

  out<<"-----------------------------------------------------------------------------------------------------"<<endl;
  out<<endl;
}

inline 
string get_transmission_type_short_name(transmission_sub_model::sm_option sm_opt)
{
  switch(sm_opt)
  {
    case transmission_sub_model::homog_no_trans  : return "homo_no_trans";
    case transmission_sub_model::homog_mendelian : return "homo_mendelian";
    case transmission_sub_model::homog_general   : return "homo_general";
    case transmission_sub_model::general         : return "general";
    case transmission_sub_model::tau_ab_free     : return "tau_AB_free";
    case transmission_sub_model::no_trans        : return "no_trans";
    case transmission_sub_model::mitochondrial   : return "mitochondrial";
  }
    
  return "";
}

/// Calculates and returns the two entries for the pair of results objects.
/// The first string is for the upper triangle part of the matrix, the
/// second is for the lower triangular.
///
/// By definition, this function should only be called if rv1 and rv2
/// have different models.
std::pair<string,string> get_asymptotic_p_value_table_entries
  (const primary_analysis_results& rv1,
   const primary_analysis_results& rv2)
{
  // Get the two options
  transmission_sub_model::sm_option sm_opt1 = rv1.get_final_model().transm_sub_model.option();
  transmission_sub_model::sm_option sm_opt2 = rv2.get_final_model().transm_sub_model.option();

  // Function invariants
  assert(sm_opt1 != sm_opt2);
  
  // Since we're doing this in a triangular fashion, we sort so that the 
  // first value is always the smaller.  
  if(sm_opt1 > sm_opt2) std::swap(sm_opt1, sm_opt2);

  double value1 = rv1.get_maxfun_results().getFinalFunctionValue();
  double value2 = rv2.get_maxfun_results().getFinalFunctionValue();

  // Calculate the absolute difference
  double abs_diff = abs(value1 - value2);

  // Number of degrees of freedom in the comparison.  0 is used as an invalid
  // state
  size_t deg_freedom = 0;
  
  double p_value = 0;
        
  switch(sm_opt1)
  {
    case transmission_sub_model::homog_no_trans  :
      switch(sm_opt2)
      {
        case transmission_sub_model::homog_mendelian :
            deg_freedom = 0;
            break;
        case transmission_sub_model::homog_general   : 
            deg_freedom = 2;
            p_value = 1.0 - chi_square_cdf(2, 2*abs_diff);
            break;
        case transmission_sub_model::general         :
            deg_freedom = 3;
            p_value = 1.0 - chi_square_cdf(3, 2*abs_diff);
            break;
        case transmission_sub_model::tau_ab_free     :
            deg_freedom = 0;
            break;
        case transmission_sub_model::no_trans :
        case transmission_sub_model::homog_no_trans : 
            // Not used
            break;
        default :
            // This is not a model we should be using in this table!
            SAGE_internal_error();
      }
      break;
    case transmission_sub_model::homog_mendelian :
      switch(sm_opt2)
      {
        case transmission_sub_model::homog_general   :
        {
            deg_freedom = 2;
            double linear = 0.25 + 
                            0.5  * chi_square_cdf(1, 2*abs_diff) +
                            0.25 * chi_square_cdf(2, 2*abs_diff);
            p_value = 1.0 - linear;
            break;
        }
        case transmission_sub_model::general         :
        {
            deg_freedom = 3;
            double linear = 0.25 * chi_square_cdf(1, 2*abs_diff) +
                            0.5  * chi_square_cdf(2, 2*abs_diff) +
                            0.25 * chi_square_cdf(3, 2*abs_diff);
            p_value = 1.0 - linear;
            break;
        }
        case transmission_sub_model::tau_ab_free     :
        {
            deg_freedom = 1;
            double linear = chi_square_cdf(1, 2*abs_diff);
            p_value = 1.0 - linear;
            break;
        }
        case transmission_sub_model::no_trans :
        case transmission_sub_model::homog_no_trans : 
        case transmission_sub_model::homog_mendelian : 
            // Not used
            break;
        default :
            // This is not a model we should be using in this table!
            SAGE_internal_error();
      }
      break;
    case transmission_sub_model::homog_general   :
      switch(sm_opt2)
      {
        case transmission_sub_model::general         :
            deg_freedom = 1;
            p_value = 1.0 - chi_square_cdf(1, 2*abs_diff);
            break;
        case transmission_sub_model::tau_ab_free     :
            deg_freedom = 0;
            break;
        case transmission_sub_model::no_trans :
        case transmission_sub_model::homog_no_trans : 
        case transmission_sub_model::homog_mendelian : 
        case transmission_sub_model::homog_general : 
            // Not used
            break;
        default :
            // This is not a model we should be using in this table!
            SAGE_internal_error();
      }
      break;
    case transmission_sub_model::general         :
        // sm_opt2 must == tau_ab_free
        deg_freedom = 2;
        p_value = 1.0 - chi_square_cdf(2, 2*abs_diff);
        break;
    case transmission_sub_model::no_trans :
    case transmission_sub_model::tau_ab_free :
        // Not used
        break;
    default :
        // This is not a model we should be using in this table!
        SAGE_internal_error();
  }
  
  if(!deg_freedom)
  {
    return make_pair(string("-----"),string("-----"));
  }
  else
  {
    std::stringstream s;

    s << fp(2*abs_diff,5,3) << '[' << deg_freedom << ']';
    
    return make_pair(s.str(), fp(p_value,5,3));
  }
}

void create_asymptotic_p_value_table
  (vector<string>& table,
   const vector<primary_analysis_results>& rv)
{
  table.resize(rv.size() * rv.size());
  
  for(size_t it1 = 0; it1 != rv.size(); ++it1)
  {
    table[it1*rv.size() + it1] = "-----";
    
    for(size_t it2 = it1+1; it2 != rv.size(); ++it2)
    {
      // Get the table entries
      std::pair<string,string> entries = 
        get_asymptotic_p_value_table_entries(rv[it1], rv[it2]);
        
      table[it1*rv.size() + it2] = entries.first;
      table[it2*rv.size() + it1] = entries.second;
    }
  }
}

void analysis_output::output_asymptotic_p_value_results
(ostream& out, const vector<primary_analysis_results>& rv)
{
  typedef vector<primary_analysis_results> ResultVector;

  vector<string> p_val_table;
  
  create_asymptotic_p_value_table(p_val_table, rv);
  
  // Store the old fill char and set the fill to ' ';
  char fill_char = out.fill();
  
  out << setfill(' ');
  
  // Output Table Header
  out<<endl<<endl;
  out<<"========================================================================================="<<endl;
  out<<"    Likelihood Ratio Criteria(above diagonal) and Asymptotic P-values(below diagonal)"<<endl;
  out<<"               Rows: null hypothesis   Columns: model (alternate hypothesis)"<<endl;

  // Create the column headers.

  out<<"               ";

  // For each element, look up the type
  for(ResultVector::const_iterator it = rv.begin(); it != rv.end(); ++it)
  {
    string name = get_transmission_type_short_name
        (it->get_final_model().transm_sub_model.option());

    out << setw(16) << name;
  }

  out<<endl;
  
  // 
  out<<"-----------------------------------------------------------------------------------------"<<endl;

  // Output entries
  for(size_t it1 = 0; it1 != rv.size(); ++it1)
  {
    string name = get_transmission_type_short_name
        (rv[it1].get_final_model().transm_sub_model.option());
    
    out << setw(15) << name;

    for(size_t it2 = 0; it2 != rv.size(); ++it2)
    {
      out << setw(16) << p_val_table[it1*rv.size() + it2];
    }
    out<<endl;
  }
  
  out<<"========================================================================================="<<endl;
  out<<"Note: The quoted P-values assume large samples and that the difference in the number of"<<endl
     <<"functionally independent parameters estimated in the two models is as indicated in "    <<endl
     <<"bracket []. Because bounds placed on other parameters in the model used may result in a"<<endl
     <<"number that is different, you are cautioned to check this before quoting the corresponding"<<endl
     <<"large-sample P-value."<<endl;

  out << setfill(fill_char);
}



} // End SEGREG namespace
} // End SAGE namespace

