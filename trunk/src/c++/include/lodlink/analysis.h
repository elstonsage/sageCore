#ifndef LODLINK_ANALYSIS_H
#define LODLINK_ANALYSIS_H
//============================================================================
// File:      analysis.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   12/20/2 - created.                                   djb
//                                                                          
// Notes:     declaration of a analysis class. 
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include <vector>
#include <ostream>
#include "boost/smart_ptr.hpp"
#include "rped/rped.h"
#include "fped/fped.h"
#include "app/SAGEapp_version_bank.h"
#include "app/SAGEapp.h"
#include "lodlink/instructions.h"
#include "lodlink/linkage_tests.h"
#include "lodlink/homogeneity_tests.h"
#include "lodlink/lods.h"
#include "lodlink/genotypes.h"

using std::vector;
using std::ostream;

namespace SAGE
{

namespace LODLINK
{

void  build_headers();

//----------------------------------------------------------------------------
//  Class:    analysis
//                                                                          
//  Purpose:  execute and write results of a set of user instructions 
//            corresponding to a parameter file lodlink analysis block.
//                                                                          
//----------------------------------------------------------------------------
//
class analysis
{
  public:
    analysis(cerrorstream& errors, const RPED::RefMultiPedigree& mped, const instructions& instr);
    ~analysis();
    
    void  build();
    void  analyze();
    void  write(ostream& summary, ostream& detail, APP::SAGEapp& app);
    
    static void  clear();
    
  private:
    static void  write_results_label(ostream& out);
    void  build_filtered_mped();
    void  write_file_header(ostream& out, const string& file_type, APP::SAGEapp& app);
    bool  missing_sex_parents() const;
    void  assign_arbitrary_sexes();
  
    // Data members.
    cerrorstream&  my_errors;
    const RPED::RefMultiPedigree&  my_original_mped;
    FPED::FilteredMultipedigree  my_mped;
    const instructions&  my_instructions;
    
    typedef boost::shared_ptr<task>  task_ptr;
    vector<task_ptr>  my_tasks;
    
    bool  aborted;
};

#include "lodlink/analysis.ipp"
}
}

#endif

