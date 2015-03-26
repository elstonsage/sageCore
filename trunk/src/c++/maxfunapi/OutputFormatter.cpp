#include "LSF/parse_ops.h"
#include "maxfunapi/OutputFormatter.h"

namespace SAGE   {
namespace MAXFUN {

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
OutputFormatter::convertMatrix(const Results& results)
{
  // Set up local variables:

  const ParameterMgr     & param_mgr = results. getParameterMgr     ();
  const CovarianceMatrix & matrix    = results. getCovarianceMatrix ();

  OUTPUT::Table t("VARIANCE-COVARIANCE MATRIX " + results.getSequenceName());

  // Check WasSkipped:

  if(results.getWasSkipped())
  {
    t << OUTPUT::NamedString("Error", OutputFormatter::getWasSkippedMessage());
  }

  // If matrix is unavailable...

//  else if(results.getCovMatrixStatus() != 0) 
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
      const Parameter & param = param_mgr.getParameter(name->first, name->second);

      if(param.getIncludeInOutput() == true)
        t << OUTPUT::TableColumn(param.getName());
    }

    // Matrix rows:

    for(size_t i = 0; i < matrix.getNames().size(); i++)
    {
      const Parameter & param = param_mgr.getParameter(matrix.getNames()[i].first, matrix.getNames()[i].second);

      if(param.getIncludeInOutput() == true)
      {
        OUTPUT::TableRow r;
        
        r << OUTPUT::String(param.getName());

        for(size_t j = 0; j < matrix.getNames().size(); j++)
        {
          const Parameter & param2 = param_mgr.getParameter(matrix.getNames()[j].first, matrix.getNames()[j].second);

          if(param2.getIncludeInOutput() == true)
            r << OUTPUT::Double(matrix.getCovariance(i, j));
        }
        
        t << r;
      }
    }
   if((results.getCovMatrixStatus()) == 1) { 
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

//=======================================================
//  convertEstimates()
//=======================================================
OUTPUT::Table
OutputFormatter::convertEstimates(const Results& results, bool detailed, bool p_val)
{
  const ParameterMgr& info = results.getParameterMgr();

  OUTPUT::Table t("MAXIMIZATION RESULTS " + results.getSequenceName());

  // Check WasSkipped:

  if(results.getWasSkipped())
  {
    t << OUTPUT::NamedString("Error", "This maximization step was skipped, therefore no results are available.");
  }

  // Check for valid initial value:

  else if(!results.getValidInitialValue())
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

    if( p_val )
    {
      OUTPUT::RenderingRules pvalue_rules;
      pvalue_rules.addLowerThreshold(0.0000001);
      pvalue_rules.addLowerThreshold(0.000001);
      pvalue_rules.addLowerThreshold(0.00001);

      t << (OUTPUT::TableColumn("P-value") << pvalue_rules);
    }

    if(detailed)
      t << (OUTPUT::TableColumn("Deriv") << OUTPUT::RenderingRules(OUTPUT::HasNumberFormat::FIXED, 10));

  // Fetch param group names:

    vector<string> names      = info.getOrderedGroupNames();
    string         group_name = "",
                   separator  = " ";

  // Loop through groups:

    for(vector<string>::iterator group = names.begin(); group != names.end(); group++)
    {
      group_name = *group;

      if((group_name == "ALL") || (group_name == "__INTERNALS__") || (group_name == "") || (info.getParamCount(group_name) == 0))
        continue;

// Make sure that there is at least one InUse parameter in the group:
//      int N = 0;
//      for(ParameterConstIterator p  = info.getParamBegin (group_name); p != info.getParamEnd   (group_name); p++)
//        if(p->isInUse() && p->getIncludeInOutput()) N++;
//      if(!N) continue;

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

      for(ParameterConstIterator p = info.getParamBegin(group_name); p != info.getParamEnd(group_name); p++)
      {

  // If param not in use, skip it:

//        if(!p->isInUse() || !p->getIncludeInOutput())
        if(! p->getIncludeInOutput())
          continue;

  // Create TableRow and add parameter name:

        OUTPUT::TableRow r;

        r << OUTPUT::String(p->getName());          // Parameter name
        r << OUTPUT::Double(p->getFinalEstimate()); // Final estimate

  // 3.3.1.3. Now output the parameter's info:

        if(p->getFinalType() == Parameter::INDEPENDENT_FUNCTIONAL ||
           p->getFinalType() == Parameter::INDEPENDENT            ||
           p->getFinalType() == Parameter::DEPENDENT)
        {

  // Standard error:

          if(! p->getIncludeStdError())
          {
            r << OUTPUT::UnavailableCell(false);
          }
          else if(p->isStdErrorAvailable())
          {
            r << OUTPUT::Double(p->getStdError());
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
          else if( p_val )
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
          r << OUTPUT::String(ParamTypeEnum2str(p->getFinalType()));
          r.spanLatestCell(detailed ? 3 : 2);
        }

        t << r;

      } // End parameter loop

    } // End group loop

    t << OUTPUT::NamedDouble("Final " + results.getFunctionName(), results.getFinalFunctionValue());

  } // End if-converged section

  return t;
}


JointTest::JointTest(const Results& results1, const Results& results2)
{
  my_H0_name  = "H0 " + results1.getFunctionName() + " " + results1.getSequenceName();
  my_H1_name  = "H1 " + results2.getFunctionName() + " " + results2.getSequenceName();

  my_H0_valid = results1.getConverged() && !results1.getWasSkipped();
  my_H1_valid = results2.getConverged() && !results2.getWasSkipped();

  my_H0_val   = results1.getFinalFunctionValue();
  my_H1_val   = results2.getFinalFunctionValue();

  my_degrees_of_freedom = abs(results1.getNumOfIndependentParams() - results2.getNumOfIndependentParams());
  my_comp_val           = 2 * (abs(my_H0_val - my_H1_val)),
  my_p_value            = chdtrc(my_degrees_of_freedom, my_comp_val);

  if(SAGE::isnan(my_p_value)) 
    my_p_value = 1.0;
}
  
JointTest::JointTest(const JointTest& other) 
    : my_H0_name            (other.my_H0_name),
      my_H1_name            (other.my_H1_name),
      my_H0_valid           (other.my_H0_valid),
      my_H1_valid           (other.my_H1_valid),
      my_H0_val             (other.my_H0_val),
      my_H1_val             (other.my_H1_val),
      my_degrees_of_freedom (other.my_degrees_of_freedom),
      my_comp_val           (other.my_comp_val),
      my_p_value            (other.my_p_value)
{}
    
JointTest&
JointTest::operator= (const JointTest& other)
{
  if(this != &other)
  {
    my_H0_name            = my_H0_name;
    my_H1_name            = my_H1_name;
    my_H0_valid           = my_H0_valid;
    my_H1_valid           = my_H1_valid;
    my_H0_val             = my_H0_val;
    my_H1_val             = my_H1_val;
    my_degrees_of_freedom = my_degrees_of_freedom;
    my_comp_val           = my_comp_val;
    my_p_value            = my_p_value;
  }
  
  return *this;
}
  
OUTPUT::Table 
JointTest::summarizeAsTable(string title) const
{
  if( title == "" )
    title = "Likelihood Ratio Test";

  OUTPUT::Table t(title);

  OUTPUT::TableRow row0 = OUTPUT::TableRow() << my_H0_name << my_H0_val;

  if(!my_H0_valid)
    row0 << OUTPUT::String("(possible non-convergence)");

  t << row0;

  OUTPUT::TableRow row1 = OUTPUT::TableRow() << my_H1_name << my_H1_val;

  if(!my_H1_valid)
    row1 << OUTPUT::String("(possible non-convergence)");

  t << row1;

  t.insertBlankRow();

  t << (OUTPUT::TableRow() << OUTPUT::String("2 * |H0 - H1|")      <<      my_comp_val)
    << (OUTPUT::TableRow() << OUTPUT::String("Degrees of freedom") << (int)my_degrees_of_freedom)
    << (OUTPUT::TableRow() << OUTPUT::String("P-value")            <<      my_p_value);

  return t;
}

} 
} 

