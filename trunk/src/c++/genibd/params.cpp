//==========================================================================
//  File:       params.cpp
//
//  Author:     Yeunjoo Song
//
//  History:    Initial implementation.                              Nov. 03
//
//  Notes:      This file implements a param for genibd analysis.
//
//  Copyright (c) 2003 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "genibd/params.h"

namespace SAGE
{

namespace GENIBD
{

genibd_parameters::genibd_parameters()
                 : my_title(""),
                   my_exact_size(GENIBD_CUTOFF),
                   my_multipoint(true),
                   my_scan_interval(false),
                   my_interval_distance(2.0),
                   my_loops(false),
                   my_output_ibd_state(false),
                   my_simulation(YES),
                   my_family_split(NO),
                   my_pair_type(RELATIVE),
                   my_regions(),
                   my_sim_parameters()
{}

genibd_parameters::genibd_parameters(const genibd_parameters& p)
                 : my_title(p.my_title),
                   my_output(p.my_output),
                   my_exact_size(p.my_exact_size),
                   my_multipoint(p.my_multipoint),
                   my_scan_interval(p.my_scan_interval),
                   my_interval_distance(p.my_interval_distance),
                   my_loops(p.my_loops),
                   my_output_ibd_state(p.my_output_ibd_state),
                   my_simulation(p.my_simulation),
                   my_family_split(p.my_family_split),
                   my_pair_type(p.my_pair_type),
                   my_regions(p.my_regions),
                   my_sim_parameters(p.my_sim_parameters)
{}

genibd_parameters&
genibd_parameters::operator=(const genibd_parameters& p)
{
  my_title             = p.my_title;
  my_output            = p.my_output;
  my_multipoint        = p.my_multipoint;
  my_loops             = p.my_loops;
  my_output_ibd_state  = p.my_output_ibd_state;
  my_simulation        = p.my_simulation;
  my_family_split      = p.my_family_split;
  my_scan_interval     = p.my_scan_interval;
  my_interval_distance = p.my_interval_distance;
  my_exact_size        = p.my_exact_size;
  my_pair_type         = p.my_pair_type;
  my_regions           = p.my_regions;
  my_sim_parameters    = p.my_sim_parameters;

  return *this;
}

genibd_parameters::~genibd_parameters()
{}

void
genibd_parameters::build_analysis_region(genome_description*       g,
                                         const genibd_region_list& regions,
                                         cerrorstream&             errors)
{
  if( !g )
    return;

  genibd_region_iterator r = regions.begin();
  for( ; r != regions.end(); ++r )
  {
    string s = r->name;

    if(    g->region(s).valid()
        && g->region(s).locus_count() > 0 )
    {
      add_region(s, r->output);
      //cout << "adding region from parser region " << s << endl;
    }
    else
      errors << priority(warning) << "Region name '" << s
             << "' is invalid.  It will be ignored."  << endl;
  }

  //cout << "region_count() = " << region_count() << endl;
  //cout << "g->region_count() = " << g->region_count() << endl;

  if( !region_count() )
  {
    for( size_t r = 0; r < (size_t)g->region_count(); ++r )
    {
      if(  g->region(r).valid()
        && g->region(r).locus_count() > 0 )
      {
        add_region(g->region(r).name());
        //cout << "adding region " << g->region(r).name() << endl;
      }
    }
  }
}

void
genibd_parameters::dump_parameters(ostream& out) const
{
  SAGE::cerrorstream o(out);
  
  o.prefix(title() + " : ");

  dump_title(o);

  o.prefix("        * - ");

  if( is_multipoint() )
    dump_multipoint(o);
  else
    dump_singlepoint(o);
}

void
genibd_parameters::dump_title(SAGE::cerrorstream& o) const
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
genibd_parameters::dump_multipoint(SAGE::cerrorstream& o) const
{
  if( allow_simulation() != ALWAYS )
  {
    if( scan_interval() )
      o << "Including intervals." << endl;
  
    o << "Maximum size for exact method is set to " << max_exact_size() << "." << endl;
  }

  if( allow_family_splitting() == YES )
  {
    if( allow_simulation() == NO )
      o << "Pedigrees larger than " << max_exact_size() << " will be split into "
        << "nuclear families." << endl;
  }
  else if( allow_family_splitting() == ALWAYS )
  {
    o << "Pedigrees will be split into nuclear families." << endl;
  }

  if( allow_simulation() == YES )
  {
    o << "Pedigrees larger than " << max_exact_size() << " will be "
      << "simulated." << endl;
  }
  else if( allow_simulation() == ALWAYS )
  {
    o << "Pedigrees will be simulated." << endl;
  }

  if( allow_simulation() != NO )
  {
    if( get_sim_parameters().get_use_factor() )
      o << "Pedigrees will be simulated based on pedigree size and "
        << "the number of markers." << endl;
    else
      o << "Simulated pedigrees will use "
        << get_sim_parameters().get_batch_count()
        << " batch(es) of "
        << get_sim_parameters().get_dememorization_step()
        << " dememorization step(s) and "
        << get_sim_parameters().get_simulation_step()
        << " simulation step(s)." << endl;

//    get_sim_parameters().dump_sim_parameters(o);
  }
}

void
genibd_parameters::dump_singlepoint(SAGE::cerrorstream& o) const
{
  if( allow_simulation() == ALWAYS )
    o << "Pedigrees will be simulated." << endl;
  
  if( allow_loops() )
  {
    if( allow_simulation() == NO )
      o << "Pedigrees with loops will use exact methods.  "
        << "Maximum size for exact method is " << max_exact_size() << '.' << endl;
    else if( allow_simulation() == YES )
      o << "Pedigrees with loops will use exact and simulation methods." << endl
        << "Maximum size for exact method is "
        << max_exact_size() << '.' << endl
        << "Pedigrees larger than " << max_exact_size() << " will be simulated." << endl;

    if( allow_family_splitting() == YES )
    {
      if( allow_simulation() == NO )
        o << "Pedigrees larger than " << max_exact_size() << " will be split into "
          << "nuclear families." << endl;
    }
    else if( allow_family_splitting() == ALWAYS )
    {
      o << "Pedigrees will be split into nuclear families." << endl;
    }
  }
  else
  {
    o << "Pedigrees with loops will not be processed." << endl;

    if( allow_family_splitting() == ALWAYS )
      o << "Pedigrees will be split into nuclear families before processing." << endl;
  }

  if( allow_simulation() != NO && allow_loops() == true )
  {
    if( get_sim_parameters().get_use_factor() )
      o << "Pedigrees will be simulated based on pedigree size and "
        << "the number of markers." << endl;
    else
      o << "Simulated pedigrees will use "
        << get_sim_parameters().get_batch_count()
        << " batch(es) of "
        << get_sim_parameters().get_dememorization_step()
        << " dememorization step(s) and "
        << get_sim_parameters().get_simulation_step()
        << " simulation step(s)." << endl;

//    get_sim_parameters().dump_sim_parameters(o);
  }
}

void
genibd_parameters::dump_regions(SAGE::cerrorstream& o) const
{
  if( !region_count() ) return;

  if( region_count() > 1 )
    o << "regions (";
  else
    o << "region (";

  genibd_region_iterator r = region_begin();
  for( ; r != region_end(); ++r )
  {
    if( r != region_begin()) o << ", ";
    
    o << r->name;
  }

  o << ')';
}

} // end of namespace GENIBD

} // end of namespace SAGE
