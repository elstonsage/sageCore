#include "tdtex/ScoreCalculator.h"

namespace SAGE  {
namespace TDTEX {

//===================================================================
//
//  score_child(...)
//
//===================================================================
TransmissionList
ScoreCalculator::score_child(
                  size_t        marker,
            const RPED::Member& parent1,
            const RPED::Member& parent2,
            const RPED::Member& child)
{
  // Create TransmissionList:
  TransmissionList xmit;

  // Extract genotypes for parents and child
  const MLOCUS::unphased_genotype parent1_geno = get_genotype(marker, parent1),
                                  parent2_geno = get_genotype(marker, parent2),
                                  child_geno   = get_genotype(marker, child);

  // Are the parents missing genotypes?
  bool parent1_missing = is_missing_genotype(parent1_geno),
       parent2_missing = is_missing_genotype(parent2_geno);

  // No scoring can be done when the child genotype is missing
  // OR when both parents are missing genotype information
  if(is_missing_genotype(child_geno) || (parent1_missing && parent2_missing))
    return xmit;

  // Compute all possible allele sources:
  //   c{x}_in_p{y} means that the child's x allele can possibly come from
  //                parent y's genotype.
  bool c1_in_p1 = is_allele_in_genotype(child_geno.allele1(), parent1_geno),
       c1_in_p2 = is_allele_in_genotype(child_geno.allele1(), parent2_geno),
       c2_in_p1 = is_allele_in_genotype(child_geno.allele2(), parent1_geno),
       c2_in_p2 = is_allele_in_genotype(child_geno.allele2(), parent2_geno);

  // Declare variables that will indicate the source of each child allele
  const MLOCUS::unphased_genotype * geno_source1   = NULL;
  const MLOCUS::unphased_genotype * geno_source2   = NULL;
  const RPED::Member              * allele_source1 = NULL;
  const RPED::Member              * allele_source2 = NULL;

  // Deal with one-missing parent case
  if(parent1_missing || parent2_missing)
  {
    // Set up variables:
    const RPED::Member              * parent  = NULL;
    const MLOCUS::unphased_genotype * geno    = NULL;
          bool                        c1_in_p = false,
                                      c2_in_p = false;

    // Get information from informative parent
    if(parent1_missing)
    {
      parent  = &parent2;
      geno    = &parent2_geno;
      c1_in_p =  c1_in_p2;
      c2_in_p =  c2_in_p2;
    }
    else // parent1 is present
    {
      parent  = &parent1;
      geno    = &parent1_geno;
      c1_in_p =  c1_in_p1;
      c2_in_p =  c2_in_p1;
    }

    // Check to see if the child is a possible descendent of parent

    if(!(c1_in_p || c2_in_p))
    {
      xmit.set_error(true);
      return xmit;
    }

    else if(c1_in_p)
    {
      // Do not score offsping that share both alleles with their parent
      // (to avoid bias, as per Am. J. Hum. Genet. 56:811-812)
      if(c2_in_p)
      {
        return xmit;
      }

      // Score that the child shares only its first allele with the
      // informative parent
      else
      {
        allele_source1 = parent;
        geno_source1   = geno;
      }
    }

    // Else, score that the child shares only its second allele with the
    // informative parent

    else // c2 MUST be in p, c1 NOT in p
    {
      allele_source2 = parent;
      geno_source2   = geno;
    }
  }

  // Both parents have known genotypes:
  else
  {
    // Check to see if the child is a possible descendent of parents

    // If c1 didn't come from either parent, or c2 didn't come from either parent, then this child
    // is not possible:
    if((!c1_in_p1 && !c1_in_p2) || (!c2_in_p1 && !c2_in_p2))
    {
      xmit.set_error(true);
      return xmit;
    }

    // Child allele 1 must come from parent 1 OR
    // child allele 2 must come from parent 2

    else if(c1_in_p1 && !c1_in_p2 || c2_in_p2 && !c2_in_p1)
    {
      allele_source1 = &parent1;
      geno_source1   = &parent1_geno;
      allele_source2 = &parent2;
      geno_source2   = &parent2_geno;
    }

    // Child allele 1 must come from parent 2 OR
    // child allele 2 must come from parent 1

    else if(c1_in_p2 && !c1_in_p1 || c2_in_p1 && !c2_in_p2)
    {
      allele_source1 = &parent2;
      geno_source1   = &parent2_geno;
      allele_source2 = &parent1;
      geno_source2   = &parent1_geno;
    }

    // Otherwise, the source of both child alleles can be assigned arbitrarily

    else
    {
      // Only score parent-of-origin if the transmissions are distinct

      if(child_geno.allele1() != get_other_allele(parent1_geno, child_geno.allele2()) || 
         child_geno.allele2() != get_other_allele(parent2_geno, child_geno.allele1()))
      {
        allele_source1 = &parent1;
        geno_source1   = &parent1_geno;
        allele_source2 = &parent2;
        geno_source2   = &parent2_geno;
      }

      // Otherwise, score alleles with NULL source
      else
      {
        allele_source1 =  NULL;
        geno_source1   = &parent1_geno;
        allele_source2 =  NULL;
        geno_source2   = &parent2_geno;
      }
    }
  }

  // Add the transmissions, if they were in fact detected:
  
  if(geno_source1)
    xmit.get_list().push_back(
      Transmission(
        allele_source1, 
        child_geno.allele1(), 
        get_other_allele(*geno_source1, child_geno.allele1())));

  if(geno_source2)
    xmit.get_list().push_back(
      Transmission(
        allele_source2, 
        child_geno.allele2(), 
        get_other_allele(*geno_source2, child_geno.allele2())));

  /// Return the transmission list:
  return xmit;
}

} // End namespace TDTEX
} // End namespace SAGE
