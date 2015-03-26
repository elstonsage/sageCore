#ifndef AO_ANALYSIS_OUTPUT_H
#define AO_ANALYSIS_OUTPUT_H
//======================================================================
//
//  File:	AnalysisOutput.h
//
//  Author:	Stephen Gross
//
//  Copyright 2002 R. C. Elston
//======================================================================


#include <string>
#include <fstream>
#include "maxfunapi/maxfunapi.h"
#include "sampling/sampling.h"
#include "output/Output.h"
#include "ageon/Datatypes.h"
#include "ageon/Model.h"
#include "ageon/MemberCovariateCalculator.h"
#include "ageon/ExtraOutput.h"

namespace SAGE {
namespace AO   {

//======================================================================
//
//  AnalysisOutput
//
//======================================================================
class AnalysisOutput
{
  public:
    //==================================================================
    // Constructors & operators:
    //==================================================================

    AnalysisOutput (const SAMPLING::PartitionedMemberDataSample &);
    AnalysisOutput (const AnalysisOutput &);

    //==================================================================
    // Public accessors:
    //==================================================================

    //==================================================================
    // Public utility functions:
    //==================================================================

    OUTPUT::Section generate_output(bool detailed = true);

    const vector<Model>                         & get_models () const;
    const Model                                 & get_model         (int t) const;
    const Maxfun_Data                           & get_maxfun_data   (int t) const;
    const MAXFUN::Results                       & get_Results       (int t) const;
    const SAMPLING::PartitionedMemberDataSample & get_sample        ()      const;
    const ModelTraitsVector                     & get_model_traits  ()      const { return my_model_traits; }

    void input(const Model                     &, int t);
    void input(const Maxfun_Data               &, int t);
    void input(const MAXFUN::Results           &, int t);
    void input(const MemberCovariateCalculator &, int t);

  private:
    //==================================================================
    // Private utility functions:
    //==================================================================

    AnalysisOutput& operator= (const AnalysisOutput &);

    string bool_out          (bool);
    
    OUTPUT::Table generateSampleHeader();
    OUTPUT::Table generateHeader();
    
    OUTPUT::Table generateClassSystem();
    OUTPUT::Section generateClasses();

    //==================================================================
    // Data members:
    //==================================================================

    const SAMPLING::PartitionedMemberDataSample & my_sample;

    vector<Model>                     my_models;
    vector<Maxfun_Data>               my_maxfun_datas;
    vector<MAXFUN::Results *>         my_Results;
    ModelTraitsVector                 my_model_traits;
};

inline string
AnalysisOutput::bool_out(bool x)
{
  if(x) return "YES";
  else  return "NO";
}

} // End namespace AGEON
} // End namespace SAGE

#endif
