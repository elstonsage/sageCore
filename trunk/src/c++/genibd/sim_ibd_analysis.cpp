//==========================================================================
//  File:    sim_ibd_analysis.h
//
//  Author:  Qing Sun
//
//  History: Version 0.90
//           Maternal & paternal bit split for sib pair done.   - yjs Jun 02
//           Updated to new libraries.                          - yjs Sep 04
//
//  Notes:   This file implements detailed ibd implementations related
//           to MCMC ibd simulation program, and provide interfaces to
//           MCMC simulation program.                                     
//
//  Copyright (c) 1998 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "genibd/sim_ibd_analysis.h"

namespace SAGE
{

namespace GENIBD
{

sim_ibd_analysis::sim_ibd_analysis(const mcmc_parameters& p, cerrorstream& e)
                : my_params(p), errors(e)
{
  my_ibds      = NULL;
  my_simulator = NULL;
  
  my_built = false;
  my_valid = false;
  my_verbose = false;
}

sim_ibd_analysis::~sim_ibd_analysis()
{
  if( my_ibds )
    delete my_ibds;

  if( my_simulator )
    delete my_simulator;
}

bool
sim_ibd_analysis::build()
{
  return my_built = true;
}

bool
sim_ibd_analysis::set_pedigree(const meiosis_map& mmap, const pedigree_region& pr)
{
  my_meiosis_map = mmap;

  my_ped_region = pr;

  my_region = my_ped_region.get_region();

  return true;
}

bool
sim_ibd_analysis::build_ibds()
{
  if( my_ibds )
    delete my_ibds;

  my_relative_pairs.resize(0);

  my_ibds = new sim_storage_ibd(my_meiosis_map, my_params, my_relative_pairs);

  gmodel_type m_type = MLOCUS::AUTOSOMAL;
  if( my_region.is_x_linked() )
    m_type = MLOCUS::X_LINKED;

  for( size_t i = 0; i < my_region.locus_count(); ++i )
  {
    string m_name = my_region.locus(i).name();
    double m_dist = my_region.locus(i).point_location(0);

    my_ibds->add_marker(m_name, m_dist, m_type);
  }

  my_ibds->build();

  return true;
}

size_t
sim_ibd_analysis::add_pair(mem_pointer m1, mem_pointer m2, pair_type pt)
{
  if( !my_ibds )
    return false;

  return my_ibds->add_pair(m1, m2, pt);
}

bool
sim_ibd_analysis::compute(ostream& info)
{
  if( !built() )
    return false;

  if( !my_ibds->built() ) 
  {
    my_ibds->build();

    assert(my_ibds->built());
  }

  MCMC::McmcMeiosisMap ped(*my_meiosis_map.get_subpedigree(), my_meiosis_map.is_x_linked());

//  ped.build();

  info << "\nSTARTING PEDIGREE: "     << my_meiosis_map.get_pedigree()->name()
       << " ( " << ped.get_family_count() << " nuclear families, "
       << ped.get_individual_count()      << " individuals. )" << endl;

#if 0
  ped.dump_map(cout);
#endif

  my_simulator = new ibd_mcmc_simulator(ped, my_ped_region, &my_params, errors);

  assert(my_simulator!=NULL);

  my_simulator->set_pairs(&my_relative_pairs);

  if( my_simulator->start(info) ) return my_valid = true;
  else                            return my_valid = false;
}


} // end of namespace GENIBD

} // end of namespace SAGE
