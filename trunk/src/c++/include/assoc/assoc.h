#ifndef ASSOC_ASSOC_H
#define ASSOC_ASSOC_H
//=======================================================================
//
//  File:	assoc.h
//
//  Author:	Stephen Gross
//
//  Copyright 2002 R. C. Elston
//=======================================================================


#include <string>
#include <fstream>
#include "LSF/LSFinit.h"
#include "LSF/LSFfile.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "app/SAGEapp.h"
#include "fped/fped.h"
#include "mped/mp.h"
#include "mped/mp_utilities.h"
#include "mped/sp.h"
#include "sampling/sampling.h"
#include "output/Output.h"
#include "sampling/IndividualValidator.h"
#include "assoc/MatrixDefs.h"
#include "assoc/AnalysisResults.h"
#include "assoc/MaximizationWrapper.h"
#include "assoc/Calculator.h"
#include "assoc/Datatypes.h"
#include "assoc/AppData.h"
#include "assoc/Configuration.h"

namespace SAGE  {
namespace ASSOC {

const double EFFECTIVELY_ZERO = 1.0E-7;

class AssocIndividualValidator : public SAMPLING::IndividualValidator
{
  public:
    bool  isValid(size_t i, const SAMPLING::IndividualTraitData& trait_data) const;
};


class assoc : public APP::SAGEapp
{
  public:
    assoc(int argc=0, char **argv=NULL);

    virtual int  main();
    
    enum  residual_type  { NULL_RESIDUALS, TEST_RESIDUALS };

  private:
    void  createSampleFields(SAMPLING::MemberDataSample& sample, const Configuration& config, const string& model_name) const;
    bool  checkSample(const SAMPLING::MemberDataSample& sample, const FPED::Multipedigree& fped, 
                      const Configuration& config, const string& model_name, cerrorstream& errors) const;
      bool  checkCovariateVariances(const SAMPLING::MemberDataSample& sample, 
                                    const string& model_name, cerrorstream& errors) const;
      bool  checkPrimaryTraitVariance(const SAMPLING::MemberDataSample& sample, const string& model_name,
                                      const string& trait_name, cerrorstream& errors) const;
      bool  checkUserCreatedFields(const SAMPLING::MemberDataSample& sample, const Configuration& config,
                                   const string& model_name, cerrorstream& errors) const;
      bool  checkForFamilies(const SAMPLING::MemberDataSample& sample, const FPED::Multipedigree& fped,
                             const Configuration& config, const string& model_name, cerrorstream& errors) const;
      bool  checkForSibships(const SAMPLING::MemberDataSample& sample, const FPED::Multipedigree& fped,
                             const Configuration& config, const string& model_name, cerrorstream& errors) const;
      bool  checkForMatePairs(const SAMPLING::MemberDataSample& sample, const FPED::Multipedigree& fped,
                              const Configuration& config, const string& model_name, cerrorstream& errors) const;
                     
    void perform_analyses(AppData& data, const FPED::FilteredMultipedigree& fmp);
    
    // - Added 6-12-7. djb
    //
    OUTPUT::Table createSampleSummary(const Sampledata& sd) const;
      size_t getSubpedigreeCount(const Sampledata& sd) const;
    void  generateResidualOutput(const Configuration& config, const map<string, string>& residuals, residual_type r_type);
};

}
} 

#endif
