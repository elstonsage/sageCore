#ifndef TDTEX_SCORE_CALCULATOR
#define TDTEX_SCORE_CALCULATOR

#include <string>
#include "mlocus/genotype.h"
#include "rped/rped.h"
#include "tdtex/Transmission.h"

namespace SAGE  {
namespace TDTEX {

/** \brief Calculates the TDT scores
  *
  */
class ScoreCalculator
{
public:

  /// @name Scoring
  //@{
  
    ///
    /// Looks at the alleles (at the given marker) for parent1, parent2, and the child in question. Figures out
    /// if either (or both!) of the child's alleles are "source-able". That is, figures out if one can conclude
    /// that a given allele did in fact come from a specific parent. Returns a list of 0, 1, or 2 "source-able"
    /// transmissions.
    static TransmissionList score_child(
            size_t        marker,
      const RPED::Member& parent1,
      const RPED::Member& parent2,
      const RPED::Member& child);

  //@}
  
  /// @name Genotype queries
  //@{
      
    ///
    /// Indicates whether or not a given genotype represents a missing genotype.
    /// \param geno The genotype in question
    /// \retval true The genotype is missing
    /// \retval false The genotype is \b not missing
    static bool is_missing_genotype(
      const MLOCUS::unphased_genotype& geno);

    ///
    /// Returns the individual's genotype for the given marker
    /// \param marker The id number of the marker in question
    /// \param The individual in question
    static MLOCUS::unphased_genotype get_genotype(
      size_t marker, 
      const RPED::Member& person);

    ///
    /// Returns the genotype name ("X/X") for the given individual on the
    /// given marker. If missing, returns "?/?".
    /// \param marker The id number of the marker in question
    /// \param person The individual in question
    static string get_genotype_name(
      size_t marker, 
      const RPED::Member& person);

    ///
    /// If the given allele is present in the given genotype, returns true. Otherwise,
    /// returns false.
    /// \param a The allele in question
    /// \param g The genotype in question
    static bool is_allele_in_genotype(
      const MLOCUS::allele& a, 
      const MLOCUS::unphased_genotype& g);

    ///
    /// Given one of the alleles in a diallelic genotype, returns the other allele.
    /// \param geno The genotype in question
    /// \param a The given allele
    static MLOCUS::allele get_other_allele(
      const MLOCUS::unphased_genotype& geno, 
      const MLOCUS::allele& a);
    
  //@}

};

//===========================================================
//  INLINE FUNCTIONS
//===========================================================

//===============================================================
//
//  get_genotype(...)
//
//===============================================================
inline
MLOCUS::unphased_genotype
ScoreCalculator::get_genotype(size_t marker, const RPED::Member& person)
{
  // Fetch marker info and phenotype id for that marker:
  const RPED::RefMarkerInfo & marker_info = person.multipedigree()->info().marker_info(marker);
        uint                  id          = person.pedigree()->info().phenotype(person.index(), marker);

  // If it's out of bounds, or it's the missing id:
  if(id >= marker_info.phenotype_count() || id == marker_info.get_missing_phenotype_id())
  {
    return MLOCUS::unphased_genotype();
  }
  // I guess it's valid; return the unphased genotype for this id:
  else
  {
    return marker_info.unphased_penetrance_begin(id).unphased_geno();
  }
}

//===============================================================
//
//  is_missing_genotype(...)
//
//===============================================================
inline
bool
ScoreCalculator::is_missing_genotype(const MLOCUS::unphased_genotype& geno)
{
  return geno == MLOCUS::unphased_genotype();
}

//===============================================================
//
//  get_genotype_name(...)
//
//===============================================================
inline
string
ScoreCalculator::get_genotype_name(size_t marker, const RPED::Member& person)
{
  const MLOCUS::unphased_genotype & geno = get_genotype(marker, person);

  return is_missing_genotype(geno) ? "?/?" : geno.name();
}

//===============================================================
//
//  is_allele_in_genotype(...)
//
//===============================================================
inline
bool
ScoreCalculator::is_allele_in_genotype(const MLOCUS::allele& a, const MLOCUS::unphased_genotype& g)
{
  return is_missing_genotype(g) ? false : a == g.allele1() || a == g.allele2();
}

//===============================================================
//
//  get_other_allele(...)
//
//===============================================================
inline
MLOCUS::allele
ScoreCalculator::get_other_allele(const MLOCUS::unphased_genotype& geno, const MLOCUS::allele& a)
{
  return geno.allele1() == a ? geno.allele2() : geno.allele1();
}

} // End namespace TDTEX
} // End namespace SAGE

#endif
