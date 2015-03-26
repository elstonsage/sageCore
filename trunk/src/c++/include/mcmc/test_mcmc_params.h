#ifndef TEST_MCMC_PARAMS_H
#define TEST_MCMC_PARAMS_H

//========================================================================== 
//  File:       test_mcmc parameters
//
//  Author:     Yeunjoo Song
//
//  History:    Initial implementation.                              May. 04
//
//  Copyright (c) 2004 R. C. Elston
//  All Rights Reserved
//========================================================================== 

#include "mcmc/definitions.h"

namespace SAGE
{

namespace MCMC
{

struct test_mcmc_region_type
{
  test_mcmc_region_type() : name(), output() { }
  test_mcmc_region_type(const string& n, const string& o) : name(n), output(o) { }

  string name;
  string output;
};

typedef std::list<test_mcmc_region_type>             test_mcmc_region_list;
typedef test_mcmc_region_list::const_iterator        test_mcmc_region_iterator;

class test_mcmc_parameters
{
  public:

    test_mcmc_parameters();

    test_mcmc_parameters(const test_mcmc_parameters&);
    
    test_mcmc_parameters& operator=(const test_mcmc_parameters&);

    ~test_mcmc_parameters();

    bool                       is_multipoint()             const; 

    test_mcmc_region_iterator  region_begin()              const;
    test_mcmc_region_iterator  region_end()                const;

    size_t                     region_count()              const;

    // Modification functions

    void set_multipoint(bool);

    void build_analysis_region(RPED::genome_description*, const test_mcmc_region_list&, cerrorstream&  err = sage_cerr);

    void add_region(const string&, const string& = string());
    void remove_region(const string&);

    // Printout functions
    
    void dump_parameters(ostream& out = cout) const;

  protected:

    void dump_multipoint(SAGE::cerrorstream& out) const;
    void dump_singlepoint(SAGE::cerrorstream& out) const;

    void dump_title(SAGE::cerrorstream& out) const;
    void dump_regions(SAGE::cerrorstream& out) const;

    bool                 my_multipoint;

    test_mcmc_region_list   my_regions;
};

#include "mcmc/test_mcmc_params.ipp"

} // end of namespace MCMC

} // end of namespace SAGE

#endif
