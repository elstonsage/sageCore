#ifndef  MCMC_PARSER_H
#define  MCMC_PARSER_H

//==========================================================================
//  File:       test_mcmc_parser.h
//
//  Author:     Yeunjoo Song
//
//  History:    Initial implementation.                              May. 04
//
//  Notes:      This file defines a parser for test_mcmc analysis.
//
//  Copyright (c) 2004 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "app/aparser.h"
#include "mcmc/test_mcmc_params.h"

namespace SAGE
{

namespace MCMC
{

class test_mcmc_parser : public SAGE::APP::BasicParser
{
  public:

    test_mcmc_parser(cerrorstream&  err = sage_cerr);
    test_mcmc_parser(const test_mcmc_parser& f);

    ~test_mcmc_parser();

    void parse_symbols(const SymbolTable* syms) {}
    void parse_parameter(const LSFBase* param) {}
    void parse_test_parameter_section(const LSFBase* params);
    void parse_test_parameter(const LSFBase* params);

    const test_mcmc_parameters& get_parameters() const {return my_parameters;}
          test_mcmc_parameters& get_parameters()       {return my_parameters;}

    const test_mcmc_region_list& get_regions() const {return my_regions;}

  private:

    void                     parse_region(const LSFBase* param);
    void                     parse_mode(const LSFBase* param);
//    void                     parse_scan_type(const LSFBase* param);

    test_mcmc_parameters     my_parameters;

    test_mcmc_region_list    my_regions;

    cerrorstream          errors;
};

} // end of namespace MCMC

} // end of namespace SAGE
                          
#endif
