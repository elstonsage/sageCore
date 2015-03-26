#ifndef LODLINK_DEFINITIONS_H
#define LODLINK_DEFINITIONS_H
//============================================================================
// File:      definitions.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   10/18/2 - created.                                   djb
//                                                                          
// Notes:     structs and classes common to more than one LODLINK file.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include <ostream>
#include <string>
#include <vector>
#include <utility>
#include "numerics/log_double.h"
#include "mlocus/penmodel.h"
#include "rped/rped.h"
#include "lodlink/instructions.h"
#include "lodlink/mle_sub_model.h"

using std::ostream;
using std::string;
using std::vector;
using std::pair;

namespace SAGE
{

namespace LODLINK
{

extern const double  NULL_THETA;
extern const double  DERIV_EPS;
extern const double  BOUND_EPS;

struct joint_pen_iter;

struct joint_genotype
{
  joint_genotype(MLOCUS::phased_genotype t, MLOCUS::phased_genotype m);
  joint_genotype(const joint_pen_iter& jpi); 
  
  void  print() const;
  double  frequency() const;
  
  MLOCUS::phased_genotype  tg;         // Trait genotype.
  MLOCUS::phased_genotype  mg;         // Marker genotype.
};

ostream& operator<<(ostream& out, const joint_genotype& jg);

struct haplotype
{
  haplotype(MLOCUS::allele t, MLOCUS::allele m);

  MLOCUS::allele  ta;             // Trait MLOCUS::allele.
  MLOCUS::allele  ma;             // Marker allele.
};

ostream& operator<<(ostream& out, const haplotype& h);

typedef MLOCUS::penetrance_model::phased_penetrance_iterator  pen_iter;

struct joint_pen_iter
{
  joint_pen_iter(pen_iter ti, pen_iter mi);
  log_double  penetrance() const;
  
  pen_iter  trait_iter;
  pen_iter  marker_iter;
};

#include "lodlink/definitions.ipp"
}
}

#endif



