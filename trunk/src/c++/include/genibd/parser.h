#ifndef  GENIBD_PARSER_H
#define  GENIBD_PARSER_H

//==========================================================================
//  File:       parser.h
//
//  Author:     Yeunjoo Song
//
//  History:    Initial implementation.                              Nov. 03
//
//  Notes:      This file defines a parser for genibd analysis.
//
//  Copyright (c) 2003 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "app/aparser.h"
#include "genibd/params.h"

namespace SAGE
{

namespace GENIBD
{

class genibd_parser : public APP::BasicParser
{
  public:

    genibd_parser(cerrorstream&  err = sage_cerr);
//    genibd_parser(const genibd_parser& f);

    ~genibd_parser();

    void parse_symbols(const SymbolTable* syms) {}
    void parse_parameter(const LSFBase* param) {}
    void parse_test_parameter_section(const LSFBase* params);
    void parse_test_parameter(const LSFBase* params);

    const genibd_parameters& get_parameters() const {return my_parameters;}
          genibd_parameters& get_parameters()       {return my_parameters;}

    const genibd_region_list& get_regions() const {return my_regions;}

  private:

    void parse_region(const LSFBase* param);
    void parse_mode(const LSFBase* param);
    void parse_scan_type(const LSFBase* param);
    void parse_loops(const LSFBase* param);
    void parse_simulation(const LSFBase* param);
    void parse_output_ibd_state(const LSFBase* param);
    void parse_family(const LSFBase* param);
    void parse_pair_types(const LSFBase* param);

    void parse_choice(const LSFBase* param, choice_type& c, const string& name);

    genibd_parameters     my_parameters;

    genibd_region_list    my_regions;

    mcmc_parser           my_mcmc_parser;

    cerrorstream          errors;
};

} // end of namespace GENIBD

} // end of namespace SAGE
                          
#endif
