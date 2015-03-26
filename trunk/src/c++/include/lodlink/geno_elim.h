#ifndef LODLINK_GENO_ELIM_H
#define LODLINK_GENO_ELIM_H
//============================================================================
// File:      geno_elim.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   3/14/3 created         -djb
//                                                                          
// Notes:     Defines a class for creating and storing genotype-eliminated
//            penetrance models. 
//                                                                          
// Copyright (c) 2003 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include <map>
#include <utility>
#include "fped/fped.h"
#include "mlocus/penmodel.h"
#include "gelim/ped_imodel_gen.h"

using std::pair;

namespace SAGE
{

namespace LODLINK
{

//----------------------------------------------------------------------------
//  Class:    ge_models
//                                                                          
//  Purpose:  create and cache genotype-eliminated penetrance models.
//                                                                          
//----------------------------------------------------------------------------
//
class ge_models
{
  public:
    typedef FPED::FilteredMultipedigree::subpedigree_type               subpedigree;
    typedef FPED::FilteredMultipedigree::subpedigree_const_pointer      subped_ptr;
    typedef std::map<pair<subped_ptr, size_t>, MLOCUS::penetrance_model>  model_map;

    static void clear_models();
    
    static const MLOCUS::penetrance_model&  get_model(const subpedigree& subped, size_t locus_idx);

  private:  
    static pedigree_imodel_generator  my_generator;
    static model_map                  my_models;
};

}
}

#endif
