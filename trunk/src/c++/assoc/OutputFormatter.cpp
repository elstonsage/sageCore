#include "LSF/parse_ops.h"
#include "assoc/OutputFormatter.h"

namespace SAGE   {
namespace ASSOC {

//======================================================================
// getWasSkippedMessage()
//======================================================================
string OutputFormatter::getWasSkippedMessage()
{
  ostringstream os;

  os << "This maximization step was skipped, therefore no results are available." << endl;

  return os.str();
}

//======================================================================
//  convertMatrix()
//======================================================================
OUTPUT::Table
OutputFormatter::convertMatrix(const MAXFUN::Results& results, double residuals_scaling_factor)
{
  // Set up local variables:

  const MAXFUN::ParameterMgr&  param_mgr = results.getParameterMgr();
  const MAXFUN::CovarianceMatrix&  matrix = results.getCovarianceMatrix();

  OUTPUT::Table t("VARIANCE-COVARIANCE MATRIX " + results.getSequenceName());

  // Check WasSkipped:

  if(results.getWasSkipped())
  {
    t << OUTPUT::NamedString("Error", OutputFormatter::getWasSkippedMessage());
  }

  // If matrix is unavailable...
  else if(results.getCovMatrixStatus() >= 3) 
  {
    t << OUTPUT::NamedString("Error", "Matrix is not available.");
  }

  // If matrix is available...
  else
  {
    // Columns:

    t << OUTPUT::TableColumn(""); // Empty top-left corner

    for(vector<pair<string, string> >::const_iterator name  = matrix.getNames().begin ();
                                                      name != matrix.getNames().end   (); ++name)
    {
      const MAXFUN::Parameter& param = param_mgr.getParameter(name->first, name->second);

      if(param.getIncludeInOutput() == true)
        t << OUTPUT::TableColumn(param.getName());
    }

    // Matrix rows:

    for(size_t i = 0; i < matrix.getNames().size(); i++)
    {
      string  group_name1 = matrix.getNames()[i].first;
      string  param_name1 = matrix.getNames()[i].second;
      const MAXFUN::Parameter& param = param_mgr.getParameter(group_name1, param_name1);

      if(param.getIncludeInOutput() == true)
      {
        OUTPUT::TableRow r;
        
        r << OUTPUT::String(param.getName());

        for(size_t j = 0; j < matrix.getNames().size(); j++)
        {
          string  group_name2 = matrix.getNames()[j].first;
          string  param_name2 = matrix.getNames()[j].second;
          const MAXFUN::Parameter& param2 = param_mgr.getParameter(group_name2, param_name2);

          if(param2.getIncludeInOutput() == true)
          {
            double  covariance = matrix.getCovariance(i, j);
            if(residuals_scaling_factor != 1.0)
            {
              covariance *= covariance_adjustment(group_name1, param_name1, group_name2, param_name2, residuals_scaling_factor);
            }
            
            r << OUTPUT::Double(covariance);
          }
        }
        
        t << r;
      }
    }
    
    if((results.getCovMatrixStatus()) == 1) 
    { 
      t.insertBlankRow();
      t << (OUTPUT::TableRow() << OUTPUT::String("Likelihood Surface flat "))
        << (OUTPUT::TableRow() << OUTPUT::String("variance matrix near singular and"))
        << (OUTPUT::TableRow() << OUTPUT::String("may be affected by rounding error.")) 
        << (OUTPUT::TableRow() << OUTPUT::String("Reducing the number of parameters "))
        << (OUTPUT::TableRow() << OUTPUT::String("to be estimated is recommended"));
    }
  }
  
  return t;
}

double  
OutputFormatter::covariance_adjustment(const string& group_name1, const string& param_name1,
                                       const string& group_name2, const string& param_name2,
                                                             double residuals_scaling_factor)
{
  double  adjustment = 1.0;
  
  // - At least one is a variance component
  if(group_name1 == "Variance components" || param_name1 == "Total variance" ||
     group_name2 == "Variance components" || param_name2 == "Total variance"   )
  {
    adjustment *= residuals_scaling_factor * residuals_scaling_factor;
  }
  
  // - Both are variance components
  if((group_name1 == "Variance components" || param_name1 == "Total variance") &&
     (group_name2 == "Variance components" || param_name2 == "Total variance")    )
  {
    adjustment *= residuals_scaling_factor * residuals_scaling_factor;
  }
  
  return  adjustment;
}                                                             
                                                                                      

//=======================================================
//  convertEstimates()
//=======================================================
OUTPUT::Table
OutputFormatter::convertEstimates(const MAXFUN::Results& results, double residuals_scaling_factor, bool detailed)
{
  const MAXFUN::ParameterMgr& info = results.getParameterMgr();

  OUTPUT::Table t("MAXIMIZATION RESULTS " + results.getSequenceName());

  // Check WasSkipped:

  if(results.getWasSkipped())
  {
    t << OUTPUT::NamedString("Error", "This maximization step was skipped, therefore no results are available.");
  }

  // Check for valid initial value:

  else if(! results.getValidInitialValue())
  {
    t << OUTPUT::NamedString("Error", "No valid initial function evaluation was available with the given parameter values.");
  }

  // Check for convergence:

  else
  {
    if(!results.getConverged())
      t << OUTPUT::NamedString("Error", "Please note that this maximization did not fully converge. The final estimates may be incorrect - check that the first derivatives are close to zero.");

  // Print parameters:

    t << OUTPUT::TableColumn("Parameter", OUTPUT::TableColumn::RIGHT)
      << OUTPUT::TableColumn("Estimate")
      << OUTPUT::TableColumn("S.E.");

    OUTPUT::RenderingRules pvalue_rules;
    pvalue_rules.addLowerThreshold(0.0000001);
    pvalue_rules.addLowerThreshold(0.000001);
    pvalue_rules.addLowerThreshold(0.00001);
    t << (OUTPUT::TableColumn("P-value") << pvalue_rules);

    if(detailed)
      t << (OUTPUT::TableColumn("Deriv") << OUTPUT::RenderingRules(OUTPUT::HasNumberFormat::FIXED, 10));

  // Get param group names:

    vector<string>  names  = info.getOrderedGroupNames();
    string  group_name = "";
    string  separator  = " ";

  // Loop through groups:

    for(vector<string>::iterator group = names.begin(); group != names.end(); group++)
    {
      group_name = *group;

      if((group_name == "ALL") || (group_name == "__INTERNALS__") || (group_name == "") || (info.getParamCount(group_name) == 0))
        continue;

      // If there is only one parameter in the group, and param_name == group_name, then don't print
      // the group name header.

      if(info.getParamBegin(group_name)->getName() != group_name)
      {
        t.beginRowGroup(group_name);
      }
      else
      {
        t.endRowGroup();
      }
        
  // Loop through parameters:

      for(MAXFUN::ParameterConstIterator p = info.getParamBegin(group_name); p != info.getParamEnd(group_name); p++)
      {
        if(! p->getIncludeInOutput())
          continue;

  // Create TableRow and add parameter name:

        OUTPUT::TableRow r;
        
        // - Correct variance components for residual scaling.
        //
        double  estimate = p->getFinalEstimate();
        string  parameter_name = p->getName();
        
        if(group_name == "Variance components" ||  parameter_name == "Total variance")
        {
          estimate *= residuals_scaling_factor * residuals_scaling_factor;
        }

        r << OUTPUT::String(parameter_name);          // Parameter name
        r << OUTPUT::Double(estimate); // Final estimate

  // 3.3.1.3. Now output the parameter's info:

        if(p->getFinalType() == MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL ||
           p->getFinalType() == MAXFUN::Parameter::INDEPENDENT            ||
           p->getFinalType() == MAXFUN::Parameter::DEPENDENT)
        {

  // Standard error:

          if(! p->getIncludeStdError())
          {
            r << OUTPUT::UnavailableCell(false);
          }
          else if(p->isStdErrorAvailable())
          {
            double  std_error = p->getStdError();
            
            if(group_name == "Variance components" || parameter_name == "Total variance")
            {
              std_error *= residuals_scaling_factor * residuals_scaling_factor;
            }            
          
            r << OUTPUT::Double(std_error);
          }
          else
          {
            r << OUTPUT::UnavailableCell();
          }

  // P-value:

          if(p->getIncludePValue())
          {
            if(p->isPValueAvailable())
            {
              r << OUTPUT::Double(p->getPValue());
            }
            else
            {
              r << OUTPUT::UnavailableCell();
            }
          }
          else
          {
            r << OUTPUT::UnavailableCell(false);
          }

  // Derivative:

          if(detailed)
          {
            if(!p->getIncludeDeriv())
            {
              r << OUTPUT::UnavailableCell(false);
            }
            else if(p->isDerivAvailable())
            {
              r << OUTPUT::Double(p->getDeriv());
            }
            else
            {
              r << OUTPUT::UnavailableCell();
            }
          }
        }

  // Information unavailable:

        else
        {
          r << OUTPUT::String(MAXFUN::ParamTypeEnum2str(p->getFinalType()));
          r.spanLatestCell(detailed ? 3 : 2);
        }

        t << r;

      } // End parameter loop

    } // End group loop

    t << OUTPUT::NamedDouble("Final " + results.getFunctionName(), results.getFinalFunctionValue());

  } // End if-converged section

  return t;
}

} 
} 

