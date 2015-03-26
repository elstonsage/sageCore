#ifndef FREQ_PHENOSET_H
#define FREQ_PHENOSET_H

#include <limits>
#include <vector>
#include <iostream>
#include <algorithm>
#include "gelim/ped_imodel_gen.h"
#include "mlocus/penmodel.h"
#include "fped/fped.h"
#include "output/Output.h"
#include "globals/SAGEConstants.h"
#include "freq/AlleleRemapping.h"

namespace SAGE {
namespace FREQ {

/// Represents, for an individual, a possible genotype and the penetrance for that individual (given the genotype).
struct GenotypeInfo
{
  MLOCUS::unphased_genotype genotype;
  double penetrance;
};

/// A vector of genotype ids.
typedef std::vector<GenotypeInfo> GenotypeInfoVector;

/// \brief Stores genotype vectors for all individuals on a particular marker.
///
/// This object is pretty simple: when constructed, it's given a multipedigree and a marker id. It
/// stores, for each individual in the multipedigree, a vector of possible genotypes.
///
/// When constructed, the Sample performs simple genotype elimination (no allele remapping yet) on
/// the given multipedigree. It also creates a vector of VALID subpedigrees (the subpedigrees
/// with genotype inconsistencies are thrown out).
///
class Sample
{
  public:

  /// @name Typedefs and structs
  //@{
  
    struct SpedInfo
    {
      FPED::SubpedigreeConstPointer sped;
    
      MLOCUS::genotype_model remapped_gmodel;
      
      AlleleRemapping allele_remapping;
    };

    typedef std::vector<SpedInfo> SpedInfoVector;

    typedef std::vector<FPED::MemberConstPointer> MemberVector;

  //@}

  /// @name Constructors / operators 
  //@{
  
    ///
    /// Constructor.
    ///
    /// NOTE: Throws a std::exception is there is no one to stick in the sample (uninformative individuals / mendelian inconsistencies).
    /// \param mp The multipedigree that this sample will represent.
    /// \param marker The id number of a marker whose genotypes will be stored.
    Sample(const FPED::FilteredMultipedigree & mp, size_t marker);
    
    ///
    /// Copy constructor.
    Sample(const Sample & other);
    
  //@}

  /// @name Multipedigree & inheritance model
  //@{
  
    ///
    /// Returns a const reference to the multipedigree represented by this Sample.
    const FPED::FilteredMultipedigree & getMultipedigree() const { return my_mp; }

    ///
    /// Returns a const reference to the RefMarkerInfo corresponding to the represented marker.    
    const RPED::RefMarkerInfo & getMarkerInfo() const { return my_mp.info().marker_info(my_marker_id); }
    
    ///
    /// Returns the marker id that this Sample represents.
    size_t getMarkerId() const { return my_marker_id; }
    
  //@}

  /// @name Misc
  //@{
  
    ///
    /// Returns the vector of subpedigrees with NO mendelian inconsistencies.
    const SpedInfoVector & getValidSpedInfos() const { return my_valid_speds; }
    
    const MemberVector & getUnconnecteds() const { return my_unconnecteds; }
    
  //@}

  /// @name Phenotypes / genotypes
  //@{
  
    ///
    /// Indicates whether or not a phenotype is present for this individual.
    bool phenotypePresent(const FPED::Member & ind) const 
    {
      return !(ind.pedigree()->info().phenotype_missing(ind.index(), my_marker_id, my_mp.info().marker_info(my_marker_id)));
    }
  
    ///
    /// Returns the phenotype id for this individual. (Remember to use phenotypePresent() to see
    /// if a phenotype is available!).
    size_t getPhenotype(const FPED::Member & ind) const
    {
      return ind.pedigree()->info().phenotype(ind.index(), my_marker_id);
    }

    ///
    /// Returns a const reference to the genotype vector for the given individual.
    const GenotypeInfoVector & getGenotypeInfoVector(const FPED::Member & ind) const { return my_genotype_infos[ind.mpindex()]; }

  //@}
  
  /// @name Mate itrs
  //@{
  
    const FPED::MateConstIterator & getMateBegin(const FPED::Member & ind) const { return my_mate_begin_itrs[ind.mpindex()]; }

    const FPED::MateConstIterator & getMateEnd(const FPED::Member & ind) const { return my_mate_end_itrs[ind.mpindex()]; }
  
  //@}
  
  /// @name Debugging
  //@{
  
    void dump() const;
  
  //@}

private:

  Sample & operator= (const Sample &); // Disallowed

  /// The multipedigree from which data is drawn.
  const FPED::FilteredMultipedigree & my_mp;
  
  /// The marker id to be analyzed
  size_t my_marker_id;

  /// A vector of genotype vector's (idx is mpindex)
  std::vector<GenotypeInfoVector> my_genotype_infos;
  
  /// List of subpedigrees valid for analysis.
  SpedInfoVector my_valid_speds;
  
  /// List of unconnecteds
  MemberVector my_unconnecteds;

  /// A vector of mate const itr's (idx is mpindex)
  std::vector<FPED::MateConstIterator> my_mate_begin_itrs;

  /// A vector of mate const itr's (idx is mpindex)
  std::vector<FPED::MateConstIterator> my_mate_end_itrs;
};

} // End namespace FREQ
} // End namespace SAGE

#endif
