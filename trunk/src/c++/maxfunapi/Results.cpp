#include <iostream>
#include <sstream>
#include "LSF/parse_ops.h"
#include "maxfunapi/Results.h"

namespace SAGE   {
namespace MAXFUN {

//======================================================================
//  CONSTRUCTOR #1
//======================================================================
Results::Results()
{
  my_sequence_name     = "";
  my_function_name     = "function value";
  my_WasSkipped        = false;
  my_ValidInitialValue = true;

  score                = numeric_limits<double>::quiet_NaN();
  DebugCfg debug_cfg(DebugCfg::NO_DEBUG_INFO);
  Maxfun_Data maxfun_data;

  flat_dir.resize(100,""); // vector of strings with names of flat parameters
  inputMaxfunData(maxfun_data, debug_cfg);
}

//======================================================================
//  CONSTRUCTOR #2
//======================================================================
Results::Results(const Maxfun_Data& data, const DebugCfg & debug)
{
  my_sequence_name     = "";
  my_function_name     = "function_value";
  my_WasSkipped        = false;
  my_ValidInitialValue = true;
  score                = numeric_limits<double>::quiet_NaN();

  flat_dir.resize(100,"");
  inputMaxfunData(data, debug);
}

//======================================================================
//  COPY CONSTRUCTOR
//======================================================================
Results::Results(const Results& other)
{
  copy(other);
}

//======================================================================
//  DESTRUCTOR
//======================================================================
Results::~Results() 
{}

//======================================================================
//  operator=()
//======================================================================
Results& 
Results::operator= (const Results& other)
{
  if(&other != this)
  {
    copy(other);
  }

  return *this;
}

//======================================================================
//  copy()
//======================================================================
void 
Results::copy (const Results & other)
{
  my_AgeOfCovMatrix                       = other.my_AgeOfCovMatrix;                  
  my_AgeOfGradientVector                  = other.my_AgeOfGradientVector;             
  my_ChangeInFunctionValue                = other.my_ChangeInFunctionValue;           
  my_CovMatrixStatus                      = other.my_CovMatrixStatus;                 
  my_CovarianceMatrix                     = other.my_CovarianceMatrix;                 
  my_ERM                                  = other.my_ERM;                             
  my_ExitFlag                             = other.my_ExitFlag;                        
  my_FinalFunctionValue                   = other.my_FinalFunctionValue;              
  my_FirstDerivApproxIndic                = other.my_FirstDerivApproxIndic;           
  my_function_name                        = other.my_function_name;
  my_GradientNorm                         = other.my_GradientNorm;                    
  my_GradientStatusFlag                   = other.my_GradientStatusFlag;              
  my_Iterations                           = other.my_Iterations;                      
  my_MaximumDifference                    = other.my_MaximumDifference;               
  my_NewtonGradientNorm                   = other.my_NewtonGradientNorm;              
  my_NumOfBoundConvergedIndependentParams = other.my_NumOfBoundConvergedIndependentParams; 
  my_NumOfDependentParams                 = other.my_NumOfDependentParams;            
  my_NumOfEstimatedParams                 = other.my_NumOfEstimatedParams;            
  my_NumOfIndependentParams               = other.my_NumOfIndependentParams;          
  my_NumOfTerms                           = other.my_NumOfTerms;
  my_NumOfTrialSearch                     = other.my_NumOfTrialSearch;                
  my_NumOfVaryingIndependentParams        = other.my_NumOfVaryingIndependentParams;   
  my_ParameterMgr                         = other.my_ParameterMgr;
  my_PreviousFunctionValue                = other.my_PreviousFunctionValue;           
  my_sequence_name                        = other.my_sequence_name;
  my_StepSizeChange                       = other.my_StepSizeChange;                  
  my_TerminationAtConstraint              = other.my_TerminationAtConstraint;         
  my_ValidInitialValue                    = other.my_ValidInitialValue;
  my_WasSkipped                           = other.my_WasSkipped;
  score                                   = other.score;
  flat_dir                                = other.flat_dir;
}

//======================================================================
//  InputMaxfunData()
//======================================================================
void 
Results::inputMaxfunData(const Maxfun_Data& data, const DebugCfg& debug)
{
  // 1. Read in Maxfun data:

        my_ExitFlag                             = data.last_error();
        my_FinalFunctionValue                   = data.value     ();
        my_MaximumDifference                    = data.difmax    ();
        my_ERM                                  = data.erm       ();
        my_ChangeInFunctionValue                = data.fch       ();
        my_PreviousFunctionValue                = data.fpr       ();
        my_GradientNorm                         = 0; // data.gtg       (); ACK!
        my_NewtonGradientNorm                   = 0; // data.ptg       (); ACK!
        my_FirstDerivApproxIndic                = data.idif      ();
        my_AgeOfGradientVector                  = data.igage     ();
        my_GradientStatusFlag                   = data.igfl      ();
        my_TerminationAtConstraint              = data.impbnd    ();
        my_Iterations                           = data.it        ();
        my_AgeOfCovMatrix                       = data.ivage     ();
        my_CovMatrixStatus                      = data.ivfl      ();
        my_NumOfTerms                           = data.nt        ();
        my_NumOfBoundConvergedIndependentParams = data.nb        ();
        my_NumOfDependentParams                 = data.nd        ();
        my_NumOfEstimatedParams                 = data.ne        ();
        my_NumOfIndependentParams               = data.ni        ();
        my_NumOfVaryingIndependentParams        = data.nv        ();
        flat_dir                                = data.flat_dir;

  // 2. Read in covariance matrix:

  if(getCovMatrixStatus() < 3)
    getCovarianceMatrix().inputData(data);

  // 3. Self-check if requested:

  if(debug.reportMaxfunOutput())
  {
    if(my_ExitFlag != 0)      // May result in a segfault otherwise.  -djb
    {
      checkSelf(debug, data);
    }
  }
}

//======================================================================
//  inputParameterMgr()
//======================================================================
void 
Results::inputParameterMgr(const ParameterMgr& info)
{
  // 1. Copy over the ParameterMgr object:

  my_ParameterMgr.copyInUseParameterMgr(info);
}

//======================================================================
//  checkSelf()
//======================================================================
int
Results::checkSelf(const DebugCfg& debug, const Maxfun_Data& data) const
{
  // 0. Set up local variables:

  int err_code = 0;

  // 1. Header:

  OUTPUT::Table table("MAXFUN OUTPUT");

  // 2. Exit flag:

  OUTPUT::TableRow exit_flag_row;
  
  exit_flag_row << "Exit flag" << my_ExitFlag;

  switch(my_ExitFlag)
  {
    case 0:  exit_flag_row << "No errors found, but zero iterations requested"; break;
    case 1:  exit_flag_row << "Convergence by criterion 1"; break;
    case 2:  exit_flag_row << "Convergence by criterion 2"; break;
    case 3:  exit_flag_row << "Convergence by criterion 3"; break;
    case 4:  exit_flag_row << "Reached maximum number of iterations without convergence"; break;
    case 5:  exit_flag_row << "Accumulation of round-off errors or boundary problem prevents further progress (in variable metric methods)"; break;
    case 6:  exit_flag_row << "Computed search direction not upward (in variable metric methods)"; break;
    case 7:  exit_flag_row << "All significant digits lost in optimal conditioning during matrix update (in variable metric methods)"; break;
    case 8:  exit_flag_row << "Gradient could not be computed (because too close to a bound or other constraint)"; break;
    case 9:  exit_flag_row << "Variance-covariance matrix could not be computed (in Newton-Raphson methods); either second partial derivatives could not be computed (because too close to a bound or other constraint) or matrix of second partials could not be inverted."; break;
    case 10: exit_flag_row << "All independent parameters converged to bounds."; break;
    case 11: exit_flag_row << "Too many parameters to estimate using derivatives (N(i) > NPV)"; break;
    case 12: exit_flag_row << "Initial parameter estimates not in domain of the function to be maximized."; break;
    case 13: exit_flag_row << "Control input invalid."; break;
  }

  table << exit_flag_row;

  // 3. Final value, iterations, break:

  table << (OUTPUT::TableRow() << "Final value"                                 << my_FinalFunctionValue)
        << (OUTPUT::TableRow() << "Iterations"                                  << my_Iterations)
        << OUTPUT::Table::INSERT_BLANK_ROW()
        << (OUTPUT::TableRow() << "Num of parameters"                           << my_NumOfTerms)
        << (OUTPUT::TableRow() << "Num of independent params converged @ bound" << my_NumOfBoundConvergedIndependentParams)
        << (OUTPUT::TableRow() << "Num of dependent params"                     << my_NumOfDependentParams)
        << (OUTPUT::TableRow() << "Num of estimated params"                     << my_NumOfEstimatedParams)
        << (OUTPUT::TableRow() << "Num of independent params"                   << my_NumOfIndependentParams)
        << (OUTPUT::TableRow() << "Num of varying independent params"           << my_NumOfVaryingIndependentParams);

  if(debug.reportMaxfunOutput())
    debug.getOutputStream() << table;

  // 6. Param info:

  OUTPUT::Table            ptable;
  Parameter::ParamTypeEnum ptype           =  Parameter::NO_PARAMTYPE;
  int                      deriv_idx       = -1,
                           var_idx         = -1;
  bool                     deriv_available =  false;

  vector<vector<double> > var_covs;

  ptable << OUTPUT::TableColumn("Parameter")
         << OUTPUT::TableColumn("Final type")
         << OUTPUT::TableColumn("Final estimate")
         << OUTPUT::TableColumn("Std. Err.")
         << OUTPUT::TableColumn("Derivative")
         << OUTPUT::TableColumn("Matrix index");
  
  for(size_t i = 0; i < (size_t)my_NumOfTerms; i++)
  {
    OUTPUT::TableRow r;

    deriv_available = false;

    r << data.label(i);

    ptype = (Parameter::ParamTypeEnum)data.ist(i);

    r << ParamTypeEnum2str(ptype);

    r << data.param(i);

    r << data.stde(i);

    if((ptype == Parameter::INDEPENDENT_FUNCTIONAL) || 
       (ptype == Parameter::DEPENDENT)              || 
       (ptype == Parameter::INDEPENDENT))
    {
      deriv_idx++;
      deriv_available = true;
    }
  
    if(deriv_available)
      r << data.g(deriv_idx);
    else
      r << OUTPUT::UnavailableCell();

    if(ptype != (int)Parameter::FIXED)
    {
      var_idx++;

      r << var_idx;

      if(data.ivfl() == 0)
      {
        vector<double> covs(0);

        for(size_t j = 0; j < (size_t)data.ne(); j++)
          covs.push_back(data.av(var_idx, j));

        var_covs.push_back(covs);
      }
    }
    
    ptable << r;
  }

  if(debug.reportMaxfunOutput())
    debug.getOutputStream() << ptable;

  // 6. Var-Cov Matrix status:

  OUTPUT::Table matrix("Variance-Covariance Matrix");

  matrix << OUTPUT::NamedInt("Var-cov matrix status code", my_CovMatrixStatus);

  switch(my_CovMatrixStatus)
  {
    case 0: matrix << OUTPUT::NamedString("Var-cov matrix status", "No problem"); break;
    case 1: matrix << OUTPUT::NamedString("Var-cov matrix status", "Round-off error in Hessian matrix"); break;
    case 2: matrix << OUTPUT::NamedString("Var-cov matrix status", "Variance-covariance matrix has negative values on the diagonal (round-off error may also be a$ problem)"); break;
    case 3: matrix << OUTPUT::NamedString("Var-cov matrix status", "Hessian matrix cannot be inverted (round-off error may also be a problem)"); break;
    case 4: matrix << OUTPUT::NamedString("Var-cov matrix status", "Hessian matrix could not be computed (too close to a bounse or other constraint)"); break;
    case 5: matrix << OUTPUT::NamedString("Var-cov matrix status", "Not attempted"); break;
  }

  // 7. Var-Cov matrix:

  if(my_CovMatrixStatus == 0)
  {
    matrix << OUTPUT::TableColumn("");

    for(size_t j = 0; j < var_covs.size(); j++)
      matrix << OUTPUT::TableColumn(doub2str(j));

    for(size_t j = 0; j < var_covs.size(); j++)
    {
      OUTPUT::TableRow r;

      r << j;

      for(size_t k = 0; k < var_covs[j].size(); k++)
        r << var_covs[j][k];

      matrix << r;
    }
  }

  if(debug.reportMaxfunOutput())
    debug.getOutputStream() << matrix;

  // X-2. Check for errors:

  if(debug.reportMaxfunOutput())
    if(err_code)
      debug.getOutputStream() << OUTPUT::NamedString("Warning", "Errors detected in this maximization.");
    else
      debug.getOutputStream() << OUTPUT::NamedString("Note", "No errors detected in this maximization.");

  // X. Return error code:

  return err_code;
}

//======================================================================
//  getParameterMgr()
//======================================================================
ParameterMgr&
Results::getParameterMgr()
{
  return my_ParameterMgr;
}

//======================================================================
//  getConverged()
//======================================================================
bool 
Results::getConverged() const
{
  if(getExitFlag() > 0 && getExitFlag() < 4 && getValidInitialValue())
   return true;
  else
   return false;
}

//======================================================================
//  setFinalFunctionValue()
//======================================================================
int 
Results::setFinalFunctionValue(double val)
{
  my_FinalFunctionValue = val;

  return 0;
}

//======================================================================
//  getCovarianceMatrix()
//======================================================================
CovarianceMatrix&
Results::getCovarianceMatrix()
{
  return my_CovarianceMatrix;
}

//======================================================================
//  getCovarianceMatrix() 
//======================================================================
const CovarianceMatrix&
Results::getCovarianceMatrix() const
{
  return my_CovarianceMatrix;
}

bool Results::compscorestat()
{
    bool needscore = false;

    for(ParameterIterator p = my_ParameterMgr.getParamBegin();
           p != my_ParameterMgr.getParamEnd(); p++) 
       {
         if (p->scoretest == true) needscore = true;
         if (p->scoretest == true) break;

       }

     return needscore;
}
}} // End namespace
