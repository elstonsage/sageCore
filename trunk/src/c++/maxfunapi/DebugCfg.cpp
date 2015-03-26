#include "maxfunapi/DebugCfg.h"

namespace SAGE   {
namespace MAXFUN {

// iter_end introduced by JA for limited output dumps

//======================================================================
//  CONSTRUCTOR #1
//======================================================================
DebugCfg::DebugCfg(DebugLevelEnum debug_level) 
      : my_debug_level(debug_level)
{
  my_ostream_ptr    = NULL;
  my_ofstream_shptr = boost::shared_ptr<ofstream>();

  setType(debug_level);
}

//======================================================================
//  COPY CONSTRUCTOR
//======================================================================
DebugCfg::DebugCfg(const DebugCfg & other)
{
  copy(other);
}

//======================================================================
//  DESTRUCTOR
//======================================================================
DebugCfg::~DebugCfg() 
{}

//======================================================================
//  operator=()
//======================================================================
DebugCfg& 
DebugCfg::operator= (const DebugCfg& other)
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
DebugCfg::copy(const DebugCfg& other)
{
  my_debug_level               = other.my_debug_level;
  my_ostream_ptr               = other.my_ostream_ptr;
  my_ofstream_shptr            = other.my_ofstream_shptr;
  my_ReportBasicConfig         = other.my_ReportBasicConfig;
  my_ReportNonParametricInput  = other.my_ReportNonParametricInput;
  my_ReportParametricInput     = other.my_ReportParametricInput;
  my_ReportNonParametricOutput = other.my_ReportNonParametricOutput; 
  my_ReportParametricOutput    = other.my_ReportParametricOutput;
  my_ReportEvaluations         = other.my_ReportEvaluations;
  my_ReportFinalResults        = other.my_ReportFinalResults;
  iter_end                     = other.iter_end;
}

//======================================================================
//  reset()
//======================================================================
void
DebugCfg::reset()
{
  my_ReportBasicConfig         = false;
  my_ReportNonParametricInput  = false;
  my_ReportParametricInput     = false;
  my_ReportNonParametricOutput = false;
  my_ReportParametricOutput    = false;
  my_ReportFinalResults        = false;
  my_ReportEvaluations         = 0;
  iter_end                     = false;
}

//======================================================================
//  setType()
//======================================================================
void
DebugCfg::setType(DebugLevelEnum debug_level)
{
  reset();

  my_debug_level = debug_level;

  switch(my_debug_level)
  { 
    case NO_DEBUG_INFO:		setReportBasicConfig         (false);
				setReportNonParametricInput  (false);
				setReportParametricInput     (false);
				setReportNonParametricOutput (false);
				setReportParametricOutput    (false);
				break;
    case BASIC:			setReportBasicConfig         (true);
				setReportFinalResults        (true);				
				break;
    case PER_RUN:		setReportBasicConfig         (true);
				setReportNonParametricInput  (true);
				setReportParametricInput     (true);
				setReportNonParametricOutput (true);
				setReportParametricOutput    (true);
				setReportFinalResults        (true);				
				break;
    case COMPLETE:		setReportBasicConfig         (true);
				setReportNonParametricInput  (true);
				setReportParametricInput     (true);
				setReportNonParametricOutput (true);
				setReportParametricOutput    (true);
				setReportFinalResults        (true);
				setReportEvaluations         (1);
				break;
  }
}
 
//======================================================================
//  getOutputType()
//======================================================================
int 
DebugCfg::getOutputType()
{
  return my_ofstream_shptr.get() != NULL; // 0 = ostream; 1 = file
}

//======================================================================
//  setDebugOutput() #1
//======================================================================
void
DebugCfg::setDebugOutput(ostream & output_stream)
{
  my_ofstream_shptr = boost::shared_ptr<ofstream>();
  my_ostream_ptr    = &output_stream;
}

//======================================================================
//  setDebugOutput() #2
//======================================================================
void
DebugCfg::setDebugOutput(const string & ofilename, bool append)
{
  my_ofstream_shptr = boost::shared_ptr<ofstream>(new ofstream(ofilename.c_str(), append ? (ios::out | ios::app) : ios::out));
  my_ostream_ptr    = my_ofstream_shptr.get();
}

//======================================================================
//  hasAnyReportedOutput()
//======================================================================
bool 
DebugCfg::hasAnyReportedOutput() const
{
  return getReportBasicConfig         () ||
         getReportNonParametricInput  () ||
         getReportParametricInput     () ||
         getReportNonParametricOutput () ||
         getReportParametricOutput    () ||
         getReportFinalResults        () ||
         getReportEvaluations         ();
}

bool DebugCfg::getReportBasicConfig         () const { return my_ReportBasicConfig;         }
bool DebugCfg::getReportNonParametricInput  () const { return my_ReportNonParametricInput;  }
bool DebugCfg::getReportParametricInput     () const { return my_ReportParametricInput;     }
bool DebugCfg::getReportNonParametricOutput () const { return my_ReportNonParametricOutput; }
bool DebugCfg::getReportParametricOutput    () const { return my_ReportParametricOutput;    }
bool DebugCfg::getReportFinalResults        () const { return my_ReportFinalResults;        }
int  DebugCfg::getReportEvaluations         () const { return my_ReportEvaluations;         }

int  DebugCfg::setReportBasicConfig         (bool b) { my_ReportBasicConfig         = b; return 0; }
int  DebugCfg::setReportNonParametricInput  (bool b) { my_ReportNonParametricInput  = b; return 0; }
int  DebugCfg::setReportParametricInput     (bool b) { my_ReportParametricInput     = b; return 0; }
int  DebugCfg::setReportNonParametricOutput (bool b) { my_ReportNonParametricOutput = b; return 0; }
int  DebugCfg::setReportParametricOutput    (bool b) { my_ReportParametricOutput    = b; return 0; }
int  DebugCfg::setReportFinalResults        (bool b) { my_ReportFinalResults        = b; return 0; }
int  DebugCfg::setReportEvaluations         (int  i) { my_ReportEvaluations         = i; return 0; }

bool DebugCfg::reportMaxfunOutput() const
{
  return getReportNonParametricOutput() || getReportParametricOutput();
}

bool DebugCfg::reportEachRunCfg() const
{
  return getReportNonParametricInput  () || getReportParametricInput  () ||
         getReportNonParametricOutput () || getReportParametricOutput () ||
         getReportEvaluations               ();
}

//=============================================================================
//  parseDebugParams()
//=============================================================================
bool 
parseDebugParams(
         DebugCfg     & debug_cfg,
   const LSFBase      * param,
         cerrorstream & errors)
{
  if(param->List())
  { 
    for(LSFList::const_iterator iter = param->List()->begin(); iter != param->List()->end(); ++iter)
    {
      string param_name = toUpper((*iter)->name());
     
      if(param_name == "LEVEL")
      {
        AttrVal a = attr_value(*iter, 0);

        if(a.has_value())
        {
          string option_string = toUpper(a.String());

               if(option_string == "NO_DEBUG_INFO") debug_cfg.setType(DebugCfg::NO_DEBUG_INFO);
          else if(option_string == "BASIC")         debug_cfg.setType(DebugCfg::BASIC);
          else if(option_string == "PER_RUN")       debug_cfg.setType(DebugCfg::PER_RUN);
          else if(option_string == "COMPLETE")      {debug_cfg.setType(DebugCfg::COMPLETE);debug_cfg.iter_end = false;}
          else if(option_string == "ITERATION_END") {debug_cfg.setType(DebugCfg::COMPLETE);debug_cfg.iter_end = "true";}

          else // Unrecognized parameter
          {
            errors << priority(warning) << "Unrecognized parameter value " << option_string << endl;
            return false;
          }
        }

      } // End param_name == "LEVEL"

      else // Unrecognizied parameter
      {
        errors << priority(warning) << "Unrecognized parameter value " << param_name << endl;
        return false;
      }

    } // End LSFList iteration

  } // End if(param->List())

  return true;
}

}
}
