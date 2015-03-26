//==========================================================================
//  File:       test_mcmc_params.cpp
//
//  Author:     Yeunjoo Song
//
//  History:    Initial implementation.                              May. 04
//
//  Notes:      This file implements a param for test_mcmc analysis.
//
//  Copyright (c) 2004 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "mcmc/test_mcmc_params.h"

namespace SAGE
{

namespace MCMC
{

test_mcmc_parameters::test_mcmc_parameters()
                    : 
//                      my_title("Analysis 1"),
                      my_multipoint(true),
//                      my_scan_interval(false),
//                      my_interval_distance(2.0),
                      my_regions()
{}

test_mcmc_parameters::test_mcmc_parameters(const test_mcmc_parameters& p)
                    :
//                      my_title(p.my_title),
//                      my_output(p.my_output),
                      my_multipoint(p.my_multipoint),
//                      my_scan_interval(p.my_scan_interval),
//                      my_interval_distance(p.my_interval_distance),
                      my_regions(p.my_regions)
{}

test_mcmc_parameters&
test_mcmc_parameters::operator=(const test_mcmc_parameters& p)
{
//  my_title             = p.my_title;
//  my_output            = p.my_output;
  my_multipoint        = p.my_multipoint;
//  my_scan_interval     = p.my_scan_interval;
//  my_interval_distance = p.my_interval_distance;
  my_regions           = p.my_regions;

  return *this;
}

test_mcmc_parameters::~test_mcmc_parameters()
{}

void
test_mcmc_parameters::build_analysis_region(RPED::genome_description*    g,
                                            const test_mcmc_region_list& regions,
                                            cerrorstream&             errors)
{
  if( !g )
    return;

  test_mcmc_region_iterator r = regions.begin();
  for( ; r != regions.end(); ++r )
  {
    string s = r->name;

    if(    g->region(s).valid()
        && g->region(s).locus_count() > 0 )
      add_region(s, r->output);

    else
      errors << priority(warning) << "Region name '" << s
             << "' is invalid.  It will be ignored."  << endl;
  }

  if( !region_count() )
  {
    for( size_t r = 0; r < (size_t)g->region_count(); ++r )
    {
      if(  g->region(r).valid()
        && g->region(r).locus_count() > 0 )
      add_region(g->region(r).name());
    }
  }
}

void
test_mcmc_parameters::dump_parameters(ostream& out) const
{
  SAGE::cerrorstream o(out);
  
//  o.prefix(get_title() + " : ");

  dump_title(o);

  o.prefix("        * - ");

  if( is_multipoint() )
    dump_multipoint(o);
  else
    dump_singlepoint(o);

  out << endl;
}

void
test_mcmc_parameters::dump_title(SAGE::cerrorstream& o) const
{
  if( is_multipoint() )
    o << "Multipoint";
  else
    o << "Singlepoint";

  o << " analysis on ";

  dump_regions(o);

  o << ':' << endl;
   
} 
void
test_mcmc_parameters::dump_multipoint(SAGE::cerrorstream& o) const
{
//  if( scan_interval() )
//    o << "Including intervals." << endl;

  o << "Pedigrees will be simulated." << endl;
}

void
test_mcmc_parameters::dump_singlepoint(SAGE::cerrorstream& o) const
{
  o << "Pedigrees will be simulated." << endl;
}

void
test_mcmc_parameters::dump_regions(SAGE::cerrorstream& o) const
{
  if( !region_count() ) return;

  if( region_count() > 1 )
    o << "regions (";
  else
    o << "region (";

  test_mcmc_region_iterator r = region_begin();
  for( ; r != region_end(); ++r )
  {
    if( r != region_begin()) o << ", ";
    
    o << r->name;
  }

  o << ')';
}

} // end of namespace MCMC

} // end of namespace SAGE
