#ifndef __TEST_IVG_NEW_H
#define __TEST_IVG_NEW_H

//
//  TEST_IVG app definition
//
//  Author: Yeunjoo Song
//
//  History:  Initial Interface design  Mar. 2002
//
//  Copyright (c) 2002 R. C. Elston

#include <iostream>
#include "app/SAGEapp.h"
#include "lvec/test_ivg_input.h"

namespace SAGE
{

class TEST_IVG : public APP::SAGEapp
{
  public:

    TEST_IVG(int argc=0, char **argv=NULL);

    ~TEST_IVG();

    // Print program information
    virtual void print_help(std::ostream &);

    // Run the application
    virtual int main();

    void run_analyses(const test_ivg_data& m);

    void parse_analyses(LSFBase* params, const RPED::RefMultiPedigree& mp, cerrorstream& errors);
    
  private:

    string         my_outfile_name;
    vector<size_t> my_sample_ids;
    set<size_t>    my_pedigree_skips;

    bool           my_consistent_out;

    ostream*       my_summary_file;
};

} // end of namespace

#endif
