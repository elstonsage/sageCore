#ifndef PALBASE_DEF_H
#define PALBASE_DEF_H

//****************************************************************************
//* File:      definitions.h                                                 *
//*                                                                          *
//* Author:    Kevin Jacobs & Yeunjoo Song                                   *
//*                                                                          *
//* History:   Version 0.0 Initial implementation                            *
//*            1. Added more relationship test functions.      - yjs Apr. 02 *
//*                                                                          *
//* Notes:     This header file defines PALBASE utility structure.           *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "mped/mp.h"
#include "rped/rped.h"
#include "pairs/relpair.h"
#include "ibd/prior_ibd.h"
#include "ibd/ibd.h"
#include "ibd/ibdfile.h"

namespace SAGE    {
namespace PALBASE {

typedef ibd_pair_info  rel_pair_data;

typedef map<sped_pointer, ibd_state_info> subped_ibd_state_map;

struct sibcluster_info
{
  set<mem_pointer> parent_set;

  vector<size_t>   fsib_pairs;
  vector<size_t>   hsib_pairs;

  map< id_pair, vector<size_t> > full_sibship_map;
};

typedef vector<sibcluster_info>         sibship_cluster;
typedef sibship_cluster::iterator       sibship_cluster_iterator;
typedef sibship_cluster::const_iterator sibship_cluster_const_iterator;

typedef set< mem_pointer >                         sib_set;
typedef map< mem_pointer, pair< double, double > > sib_mean_map;

} // end of namespace PALBASE
} // end of namespace SAGE

#endif
