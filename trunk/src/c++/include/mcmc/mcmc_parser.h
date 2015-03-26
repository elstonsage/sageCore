#ifndef  MCMC_PARSER_H
#define  MCMC_PARSER_H

//==========================================================================
//  File:       mcmc_parser.h
//
//  Author:     Qing Sun
//              Geoff Wedig
//
//  History:    Version 0.90                                   Aug. 28, 1998
//              1.0 Completely rewritten so makes sense        Feb   4, 2000
//              2.0 Updated to new libraries.                  yjs Aug. 2004
//
//  Notes:      This header gives basic parameters that are used in both
//              mcmc_ibd simulation program and mcmc haplotyping program.
//
//  Copyright (c) 1998 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "app/aparser.h"
#include "mcmc/mcmc_params.h"

namespace SAGE
{

namespace MCMC
{

class mcmc_parser : public SAGE::APP::BasicParser
{
  public:

    mcmc_parser(mcmc_parameters& p, SAGE::cerrorstream& err = SAGE::sage_cerr);
    mcmc_parser(mcmc_parameters& p, const SymbolTable* syms,
                                    SAGE::cerrorstream& err = SAGE::sage_cerr);
    mcmc_parser(mcmc_parameters& p, const LSFBase* params,
                                    SAGE::cerrorstream& err = SAGE::sage_cerr);
    mcmc_parser(mcmc_parameters& p, const SymbolTable* syms, const LSFBase* params,
                                    SAGE::cerrorstream& err = SAGE::sage_cerr);

    void parse_symbols(const SymbolTable* syms);
    void parse_parameter(const LSFBase* param);
    void parse_test_parameter_section(const LSFBase* params);
    void parse_test_parameter(const LSFBase* param);

    const mcmc_parameters& parameters() const { return my_parameters; }
          mcmc_parameters& parameters()       { return my_parameters; }

  protected:

    void parse_mode           (const LSFBase* param);
    void parse_dememorization (const LSFBase* param);
    void parse_sim_steps      (const LSFBase* param);
    void parse_batches        (const LSFBase* param);
    void parse_factor         (const LSFBase* param);
    void parse_base_factor    (const LSFBase* param);
    void parse_demem_factor   (const LSFBase* param);
    void parse_sim_factor     (const LSFBase* param);
    void parse_batch_factor   (const LSFBase* param);
    void parse_local_marker   (const LSFBase* param);
    void parse_local_ind      (const LSFBase* param);
    void parse_weights        (int i, const LSFBase* param, const string&);
    void parse_tunnel         (const LSFBase* param);
    void parse_random_seed    (const LSFBase* param);

  private:
    
    mcmc_parameters&      my_parameters;

    cerrorstream         errors;
};

} // end of namespace MCMC

} // end of namespace SAGE
                          
#endif
