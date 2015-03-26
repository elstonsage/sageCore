//==========================================================================
//  File:     test_mcmc_input.cpp
//
//  Author:   Yeunjoo Song
//
//  History:  Initial implementation.                                May. 04
//
//  Copyright (c) 2004 R. C. Elston
//  All Rights Reserved
//==========================================================================

#include "mcmc/test_mcmc_input.h"

namespace SAGE
{

namespace MCMC
{

test_mcmc_data::test_mcmc_data(const string& program_name, bool debug)
              : APP::SAGE_Data(program_name, debug)
{}

test_mcmc_data::~test_mcmc_data()
{}

bool
test_mcmc_data::input(int argc, char** argv)
{
  //read parameter file
  //
  read_parameter_file(argv[1]);

  // Read in allele delimiter & marker locus file.
  read_locus_description_file(argv[3]);

  //read pedigree data file
  //
  read_family_data_file(argv[2], false, true, true, false, false);

  //read genome file
  //
  if( argc == 5 )
    read_genome_description_file(argv[4]);

  //read test_mcmc analysis
  //
  read_analysis();

  return true;
}

bool
test_mcmc_data::read_analysis()
{
  my_analysis.resize(0);

  LSFList::const_iterator i;
  for( i = my_params->List()->begin(); i != my_params->List()->end(); ++i )
  {
    if( !*i ) continue;

    if( toUpper((*i)->name() ) == "MCMC_ANALYSIS" ||
        toUpper((*i)->name() ) == "MCMC" ||
        toUpper((*i)->name() ) == "IBD_ANALYSIS" )
      my_analysis.push_back(*i);
  }

  return true;
}

} // end of namespace MCMC

} // end of namespace SAGE
