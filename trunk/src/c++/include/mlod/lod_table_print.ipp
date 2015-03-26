//
//  Lod Score table classes for MLOD
//
//  Copyright (C) 2005 R. C. Elston

#ifndef LOD_TABLE_H
#include "mlod/lod_table_print.h"
#endif

namespace SAGE
{
namespace MLOD
{

template<class TABLE_TYPE>
OUTPUT::Table LodTableFormatter::formatTable
    (const TABLE_TYPE& table,
     const string&     tname,
     bool              include_intervals) const
{
  OUTPUT::Table out_table("Lod Scores for " + tname);

  // Create the headings
  out_table << (OUTPUT::TableColumn("Position") << OUTPUT::RenderingRules(OUTPUT::RenderingRules::FIXED, 2))
            << OUTPUT::TableColumn("Marker");
                   
  std::for_each(my_params.get_trait_list().begin(),
                my_params.get_trait_list().end(),
                boost::bind(&LodTableFormatter::insert_trait_column_headers,
                            *this, boost::ref(out_table), 
                            boost::bind(get_trait_model_pair_first_element, _1)));

  // Print each row

  // Keep track of row/point
  size_t point_idx = 0;

  // print the first marker.
  insert_table_row(out_table, table, point_idx,
                   my_params.get_region().locus(0).location(),
                   my_params.get_region().locus(0).name());

  ++point_idx;

  size_t marker_idx = 0;  
  for( ; marker_idx < my_params.get_region().locus_count() - 2; ++marker_idx )
  {
    // Print interval locations for either interval or both cases.
    if( my_params.get_scan_type() != AnalysisParameters::ST_MARKER )
    {
      // Determine how many points exist between marker marker_idx and marker_idx+1
      size_t num_pts_in_intval =  my_params.get_region().locus(marker_idx).interval_point_count(1);

      // Only print the interval rows if specified, otherwise, skip them
      if( include_intervals )
      {
        for(size_t pt = 1; pt < num_pts_in_intval; ++pt)
        {
          insert_table_row(out_table, table, point_idx,
                           my_params.get_region().point_location(point_idx));
          ++point_idx;
        }
      }
      else
      {
        point_idx += num_pts_in_intval - 1;
      }
    }

    // Print marker for either marker or both cases.
    if( my_params.get_scan_type() != AnalysisParameters::ST_INTERVAL )
    {
      if(    tname == "Analysis"
          || my_params.get_ped_output_detail_option() != AnalysisParameters::PD_INTERVALS )
      {
        insert_table_row(out_table, table, point_idx,
                         my_params.get_region().locus(marker_idx+1).location(),
                         my_params.get_region().locus(marker_idx+1).name());
      }

      ++point_idx;
    }
    else
    {
      double current_location = my_params.get_region().locus(marker_idx+1).location();

      if( fabs(remainder(current_location, my_params.get_scan_distance())) < 1.0e-13  )
      {
        insert_table_row(out_table, table, point_idx,
                         my_params.get_region().locus(marker_idx+1).location(),
                         my_params.get_region().locus(marker_idx+1).name());
      }

      ++point_idx;
    }
  }

  // print the intervals before last marker.
  if( my_params.get_scan_type() != AnalysisParameters::ST_MARKER )
  {
    // Determine how many points exist between marker marker_idx and marker_idx+1
    size_t num_pts_in_intval =  my_params.get_region().locus(marker_idx).interval_point_count(1);

    // Only print the interval rows if specified, otherwise, skip them
    if( include_intervals )
    {
      for(size_t pt = 1; pt < num_pts_in_intval; ++pt)
      {
        insert_table_row(out_table, table, point_idx,
                         my_params.get_region().point_location(point_idx));
        ++point_idx;
      }
    }
    else
    {
      point_idx += num_pts_in_intval - 1;
    }
  }

  // print the last marker.
  insert_table_row(out_table, table, point_idx,
                   my_params.get_region().locus(marker_idx+1).location(),
                   my_params.get_region().locus(marker_idx+1).name());

  return out_table;
}

template<class TABLE_TYPE>
void LodTableFormatter::insert_table_row
    (OUTPUT::Table&    otable,
     const TABLE_TYPE& info,
     size_t            point_idx,
     double            position,
     string            marker_name) const
{
  OUTPUT::TableRow row;

  row << position;

  if(marker_name == "")
    row << OUTPUT::UnavailableCell(false);
  else
    row << marker_name;

  for(size_t i = 0; i != my_params.get_trait_list().size(); ++i)
  {
    row << info.get_lod_score_info(i, point_idx).lod_score
        << info.get_lod_score_info(i, point_idx).info_content;
  }

  otable << row;
}

inline 
void LodTableFormatter::insert_trait_column_headers
    (OUTPUT::Table& table,
     size_t trait) const
{
  string trait_name = my_trait_info->marker_info(trait).name();

  table.beginColumnGroup(trait_name);

  OUTPUT::RenderingRules rules(OUTPUT::RenderingRules::FIXED, 4);

  table << (OUTPUT::TableColumn("Lod Score")   << rules)
        << (OUTPUT::TableColumn("Information") << rules);
}

}
}

