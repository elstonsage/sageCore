#ifndef SEGREG_ANALYSIS_CONTROLLER_H
#define SEGREG_ANALYSIS_CONTROLLER_H

#include "segreg/analysis.h"
#include "segreg/PedigreeDataSet.h"
#include "segreg/model.h"
#include "rped/rped.h"
#include "app/output_streams.h"
#include <fstream>

namespace SAGE
{
namespace SEGREG
{

/// \brief Decides what needs to be done to perform a particular analysis
///
/// The AnalysisController is responsible for processing analyses given a
/// raw dataset and a model.  It determines which analyses to perform, performs
/// them and produces the output related to them.
class AnalysisController
{
  public:

    AnalysisController(APP::Output_Streams&);

    void do_analysis(const RPED::MultiPedigree& mp, const model& model) const;
    unsigned  fpmm_and_no_poly_loci; // due to JA trial for zero polygenic loci)
    
    vector<std::pair<string,double> >
    do_single_evaluation \
    (const PedigreeDataSet& ped_data, const model& m, double& func_val, vector<string>& par_type) const; 
    // for single likelihood evaluation
   void produce_allfixed_output \
   (vector<std::pair<string,double> > name_val,vector<string> par_type, double value, string file) const;

     double like_cutoff; // due to JA, see comments in parser.h 
  protected:

    void process_analysis        (const PedigreeDataSet& ped_data, const model& m) const;

    void do_commingling_analysis (const PedigreeDataSet& ped_data, model& m) const;
    void do_segregation_analysis (const PedigreeDataSet& ped_data, model& m) const;

    void do_reg_analysis         (const PedigreeDataSet& ped_data, model& m) const;

    void print_analysis_header   (const std::string&     title) const;
    void print_analysis_footer   ( )                            const;
    
    // for single likelihood evaluation

    APP::Output_Streams&     my_out;
    mutable primary_analysis my_analysis;

    bool check_all_fixed(const PedigreeDataSet& ped_data, const model& test_model) const; 
    // for single likelihood evaluation
   void produce_allfixed_output(vector<std::pair<string,double> > name_val,\
        vector<string> par_type,double& value,string file);
};

}
}

#include "segreg/AnalysisController.ipp"

#endif 
