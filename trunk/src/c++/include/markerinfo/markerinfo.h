#ifndef __MARKERINFO_NEW_H
#define __MARKERINFO_NEW_H

//============================================================================
//  File:     markerinfo.h
//
//  Author:   Yeunjoo Song
//
//  History:  Initial Interface design                               Mar. 2002
//
//  Notes:    MARKERINFO app definition
//
//  Copyright (c) 2002 R. C. Elston
//    All Rights Reserved
//============================================================================

#include <iostream>
#include "app/SAGEapp.h"
#include "markerinfo/input.h"

namespace SAGE
{

class MARKERINFO : public APP::SAGEapp
{
  public:

    MARKERINFO(int argc=0, char **argv=NULL);

    ~MARKERINFO();

    // Run the application
    virtual int main();

    void run_analyses(const markerinfo_data& m);

    void parse_analyses(LSFBase* params, const RPED::RefMultiPedigree& mp, cerrorstream& errors);
    
  private:

    void print_resetted_pedigree(const inconsistency_handler& h, const FPED::Multipedigree& fped);
    void print_marker_pedigree_table(const inconsistency_handler& h, const FPED::Multipedigree& fped);

    string get_pheno_string(size_t p, const MLOCUS::inheritance_model& model) const;

    string         my_outfile_name;
    vector<size_t> my_sample_ids;
    set<size_t>    my_pedigree_skips;

    bool           my_consistent_out;
    bool           my_pedigree_out;

    ostream*       my_summary_file;
    ostream*       my_pedigree_file;
    ostream*       my_parameter_file;
    ostream*       my_table_file;
};

} // end of namespace

#endif
