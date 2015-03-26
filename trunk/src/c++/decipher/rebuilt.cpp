//============================================================================
// File:      rebuilt.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   created 2/18/5                               djb
//                                                                          
// Notes:     
//                                                                          
// Copyright (c) 2005 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "decipher/rebuilt.h"

using namespace std;

namespace SAGE
{

namespace DECIPHER
{

//============================================================================
// IMPLEMENTATION:  rebuilt_em_phenotype_map
//============================================================================
//
rebuilt_em_phenotype_map::rebuilt_em_phenotype_map(output_state& ostate, APP::Output_Streams& streams, const locus_group& loci,
                                                   const vector<member>& members,  
                                                   const vector<const member_em_phenotype_map*>& phenotype_maps)
      : base_em_phenotype_map(ostate, streams, loci)
{
  // - Members in each subpopulation.
  //
  vector<set<member, member_order<member> > >  map_members;
  
  size_t  map_count = phenotype_maps.size();
  for(size_t map = 0; map < map_count; ++map)
  {
    assert(phenotype_maps[map]->loci() == my_loci);
    map_members.push_back(phenotype_maps[map]->members());
  }

  size_t  member_count = members.size();
  for(size_t m = 0; m < member_count; ++m)
  {
    // - Combinations in phenotypes are stored as indices into the haplotype_map.  The
    //   haplotype sequences must be retrieved so that the new haplotype map and phenotypes can
    //   be built.
    //
    pair<const em_phenotype*, const em_haplotype_map*>  phenotype_and_map = 
                                  get_phenotype(members[m], map_members, phenotype_maps);
    set<hap_seq_comb>  hs_combs;
    const vector<pair<em_phenotype::combination, double> >&  emp_combs = phenotype_and_map.first->combinations();
    size_t  comb_count = emp_combs.size();
    size_t  hap_count = 0;
    for(size_t c = 0; c < comb_count; ++c)
    {
      hap_seq_comb  hs_comb;
      
      hap_count = emp_combs[c].first.size();
      for(size_t h = 0; h < hap_count; ++h)
      {
        hs_comb.insert(phenotype_and_map.second->index_to_hap_seq(emp_combs[c].first[h]));
      }
 
      hs_combs.insert(hs_comb);
    }
    
    my_total_hap_count += hap_count;
    my_phenotypes.push_back(em_phenotype(&my_haplotypes, hs_combs));
  }
  
  my_haplotypes.update_counts(*this);
}

pair<const em_phenotype*, const em_haplotype_map*>  
rebuilt_em_phenotype_map::get_phenotype(const member ind, 
                                        const vector<set<member, member_order<member> > >& map_members,
                                        const vector<const member_em_phenotype_map*>& phenotype_maps   )
{
  size_t  map_count = phenotype_maps.size();

  size_t  m = 0;
  bool  member_found = false;
  for(; m < map_count; ++m)
  {
    if(map_members[m].find(ind) != map_members[m].end())
    {
      member_found = true;
      break;
    }
  }
  
  assert(member_found);
  
  return  make_pair(&((*phenotype_maps[m])[ind]), &(phenotype_maps[m]->haplotypes()));
}

}
}
