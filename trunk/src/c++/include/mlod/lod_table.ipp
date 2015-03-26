//
//  Lod Score table classes for MLOD
//
//  Copyright (C) 2005 R. C. Elston

#ifndef LOD_TABLE_H
#include "mlod/lod_table.h"
#endif

namespace SAGE
{
namespace MLOD
{

inline
SpedLodTable::LodScoreInfoType::LodScoreInfoType()
  : lod_score    (std::numeric_limits<double>::quiet_NaN()),
    info_content (std::numeric_limits<double>::quiet_NaN())
{ }
      
inline   
SpedLodTable::LodScoreInfoType::LodScoreInfoType(double lscore, double icontent)
  : lod_score    (lscore),
    info_content (icontent)
{ }
      
inline   
SpedLodTable::LodScoreInfoType::LodScoreInfoType(const LodScoreInfoType& l)
  : lod_score    (l.lod_score),
    info_content (l.info_content)
{ }
   
inline   
SpedLodTable::LodScoreInfoType& 
    SpedLodTable::LodScoreInfoType::operator=(const LodScoreInfoType& l)
{
  if(this != &l)
  {
    lod_score    = l.lod_score;
    info_content = l.info_content;
  }

  return *this;
}

inline SpedLodTable::SpedLodTable()
  : my_num_traits(0),
    my_num_points(0), 
    my_dimension (0)
{ }

inline
SpedLodTable::SpedLodTable(size_t num_traits, size_t num_points, size_t dim)
  : my_num_traits(num_traits),
    my_num_points(num_points),
    my_dimension (dim)
{
  my_lod_scores.resize(num_traits * num_points);
}

inline
SpedLodTable::SpedLodTable(const SpedLodTable& l)
  : my_lod_scores  (l.my_lod_scores),
    my_num_traits  (l.my_num_traits),
    my_num_points  (l.my_num_points),
    my_dimension   (l.my_dimension)
{ }

inline
SpedLodTable& SpedLodTable::operator= (const SpedLodTable& l)
{
  if(this == &l) return *this;
  
  my_lod_scores   = l.my_lod_scores;
  my_num_traits   = l.my_num_traits;
  my_num_points   = l.my_num_points;
  my_dimension    = l.my_dimension;

  return *this;
}

inline size_t SpedLodTable::get_trait_count() const { return my_num_traits; }
inline size_t SpedLodTable::get_point_count() const { return my_num_points; }
inline size_t SpedLodTable::get_dimension() const   { return my_dimension;  }

inline
void SpedLodTable::set_lod_score_info (size_t trait, size_t point, const LodScoreInfoType& l)
{
  my_lod_scores[trait * my_num_points + point] = l;
}

inline
const SpedLodTable::LodScoreInfoType& 
    SpedLodTable::get_lod_score_info(size_t trait, size_t point) const
{
  return my_lod_scores[trait * my_num_points + point];
}


inline
SummaryLodTable::SummaryLodTable(size_t num_traits, size_t num_points)
  : my_num_traits(num_traits),
    my_num_points(num_points)
{
  my_lod_scores.resize(num_traits * num_points);
}

inline
SummaryLodTable::SummaryLodTable (const SummaryLodTable& l)
  : my_lod_scores(l.my_lod_scores),
    my_num_traits(l.my_num_traits),
    my_num_points(l.my_num_points)
{ }

inline
SummaryLodTable& SummaryLodTable::operator= (const SummaryLodTable& l)
{
  if(this == &l) return *this;
  
  my_lod_scores   = l.my_lod_scores;

  my_num_traits   = l.my_num_traits;
  my_num_points   = l.my_num_points;

  return *this;
}

inline
const SummaryLodTable::LodScoreInfoType&
  SummaryLodTable::get_lod_score_info (size_t trait, size_t point) const
{
  return my_lod_scores[trait * my_num_points + point].info;
}

inline
size_t SummaryLodTable::table_count(size_t trait, size_t point) const
{
  return my_lod_scores[trait * my_num_points + point].table_count;
}

inline size_t SummaryLodTable::get_trait_count() const { return my_num_traits; }
inline size_t SummaryLodTable::get_point_count() const { return my_num_points; }

inline 
SummaryLodTable::LodScoreData::LodScoreData()
  : info(),
    dimension(0),
    table_count(0)
{ }

inline 
SummaryLodTable::LodScoreData::LodScoreData(const LodScoreData& l)
  : info(l.info),
    dimension(l.dimension),
    table_count(l.table_count)
{ }

inline 
SummaryLodTable::LodScoreData& SummaryLodTable::LodScoreData::operator=(const LodScoreData& l)
{
   if(this != &l)
   {
     info = l.info;
     dimension = l.dimension;
     table_count = l.table_count;
   }

   return *this;
}

}}

