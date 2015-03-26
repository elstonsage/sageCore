#ifndef LOD_TABLE_H
#define LOD_TABLE_H

//
//  Lod Score table classes for MLOD
//
//  Copyright (C) 2005 R. C. Elston


#include <vector>
#include <limits>
#include "math.h"

namespace SAGE
{
namespace MLOD
{

/// \brief Stores lod scores for a subpedigree over a region for a set of traits.
///
/// Included in the table is a dimension parameter.  This is the size
/// of the likelihood vector in elements that was used to calculate
/// the lod score (ie 2^(2n-f)).  This is stored for later use in calculating
/// summary statistics.
///
/// It is important to note that the SpedLodTable actually knows nothing
/// about the subpedigree to which it belongs.  It is, in this sense, generic
/// data, but is named based upon how it is to be used.
class SpedLodTable
{
public:

  friend class SummaryLodTable;

  /// \brief Stores the lod score/information content pair.
  ///
  struct LodScoreInfoType
  {
    /// \name Object Management
    //@{
    /// Default Constructor
    ///
    LodScoreInfoType();
    
    /// Constructor
    ///
    /// \param lscore The lod score
    /// \param icontent The information content
    LodScoreInfoType(double lscore, double icontent);
    
    /// Copy Constructor
    ///
    LodScoreInfoType(const LodScoreInfoType&);
    
    /// Copy Operator
    ///
    LodScoreInfoType& operator=(const LodScoreInfoType&);
    //@}
    
    double lod_score;     ///< Storage for lod score
    double info_content;  ///< Storage for information content
  };
  
  /// \name Object Management
  //@{
    SpedLodTable();
    
    /// Constructor 
    ///
    /// \param num_traits The number of traits analyzed
    /// \param num_points The number of points in the region
    /// \param dim        The dimension (2^(2n-f)) of the likelihood vector
    SpedLodTable(size_t num_traits, size_t num_points, size_t dim);
    SpedLodTable(const SpedLodTable&);
    
    SpedLodTable& operator= (const SpedLodTable&);
  //@}

  /// \name Basic Info
  //@{
    
    /// Returns the number of traits analyzed
    ///
    size_t get_trait_count() const;

    /// Returns the number of points in the region analyzed
    ///
    size_t get_point_count() const;
    
    /// Returns the dimension of the likelihood vector
    size_t   get_dimension() const;

  //@}
  /// \name Data access
  //@{

    /// Sets a particular element in the matrix (indexed by trait and point)
    /// to the LodScoreInfoType given
    ///
    /// \param trait The index of the trait
    /// \param point The index of the point
    /// \param data  The data to be set at (trait,point)
    void set_lod_score_info (size_t trait, size_t point, const LodScoreInfoType& data);

    /// Returns the information currently stored at (trait,point)
    ///
    /// \param trait The index of the trait
    /// \param point The index of the point
    const LodScoreInfoType& get_lod_score_info(size_t trait, size_t point) const;
  //@}

protected:

  /// Main Data Storage
  ///
  std::vector<LodScoreInfoType> my_lod_scores;

  size_t my_num_traits; ///< Number of traits analyzed
  size_t my_num_points; ///< Number of points in region analyzed
  size_t my_dimension;  ///< Size of likelihood vector
};

/// \brief Stores a summary of lod scores
///
/// Like the SpedLodTable, this stores a set of lod scores at various points along
/// a chromosome for various traits of interest.  It stores the summarized
/// values of a set of SpedLodTables.  For each entry, it stores the number of
/// tables that were used in composing that entry.
class SummaryLodTable
{
public:

  typedef SpedLodTable::LodScoreInfoType LodScoreInfoType;

  /// \name Object Management
  //@{
  /// Constructor
  ///
  /// \param num_traits The number of traits in the analysis
  /// \param num_points The number of points in the region analyzed
  SummaryLodTable(size_t num_traits = 0, size_t num_points = 0);
  
  /// Copy Constructor
  ///
  SummaryLodTable(const SummaryLodTable&);

  /// Copy Operator
  SummaryLodTable& operator= (const SummaryLodTable&);
  //@}
  
  /// Addition operation.  This adds a lod table into the Summary.  Lod
  /// scores are added (on log scale), while information content are a weighted average
  /// based on dimension.  When either value is unknown (quiet NaN) the
  /// previous values are not modified.  It also keeps track of the number of
  /// values that have been accumulated.  If the table is the wrong size
  /// (traits or points) the table is not added and a false is returned. 
  /// Otherwise returns true.
  bool add_table(const SpedLodTable&);

  /// Returns the LodScoreInfoType associated with a particular (trait,point)
  /// location
  ///
  /// \param trait index of the trait
  /// \param point index of the point
  const LodScoreInfoType& get_lod_score_info (size_t trait, size_t point) const;

  /// Returns the number of tables that had valid entries at (trait,point) and
  /// were therefore composed into the entry at that location.
  ///
  /// \param trait index of the trait
  /// \param point index of the point
  size_t table_count(size_t trait, size_t point) const;

  /// Returns number of traits in the analysis
  ///
  size_t get_trait_count() const;
  /// Returns number of points in the region analyzed
  ///
  size_t get_point_count() const;

protected:

  /// \brief Stores data about lod scores and composition statistics
  struct LodScoreData
  {
    LodScoreData();
    LodScoreData(const LodScoreData&);
    LodScoreData& operator=(const LodScoreData&);
    
    LodScoreInfoType info;
    size_t           dimension;
    size_t           table_count;
  };

  /// Primary Storage
  ///
  std::vector<LodScoreData> my_lod_scores;

  size_t my_num_traits; ///< Number of traits in analysis
                        ///
  size_t my_num_points; ///< Number of points in analyzed region
                        ///
};

}}

#include "mlod/lod_table.ipp"

#endif

