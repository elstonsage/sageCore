#ifndef DECIPHER_REBUILT_H
#define DECIPHER_REBUILT_H

//============================================================================
// File:      rebuilt.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   2/18/5 - created.                                   djb
//                                                                          
// Copyright (c) 2005 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "decipher/em.h"

namespace SAGE
{

namespace DECIPHER
{

//----------------------------------------------------------------------------
//  Class:    rebuilt_em_phenotype_map
//                                                                         
//  Purpose:  em_phenotype_map which is constructed of prebuilt em_phenotypes.
//                                                                          
//----------------------------------------------------------------------------
//
class rebuilt_em_phenotype_map : public base_em_phenotype_map
{
  public:
  
    // Constructor/destructor.
    rebuilt_em_phenotype_map(output_state& ostate, APP::Output_Streams& streams, const locus_group& loci, 
                             const vector<member>& members, 
                             const vector<const member_em_phenotype_map*>& phenotype_maps); 

  private:
    static pair<const em_phenotype*, const em_haplotype_map*>  get_phenotype(const member ind, 
                                              const vector<set<member, member_order<member> > >& map_members,
                                              const vector<const member_em_phenotype_map*>& phenotype_maps);
};

#include "decipher/rebuilt.ipp"

}
} 

#endif

