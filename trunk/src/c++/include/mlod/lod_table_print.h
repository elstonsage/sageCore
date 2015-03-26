#ifndef LOD_TABLE_PRINT_H
#define LOD_TABLE_PRINT_H

//
//  Lod Score table classes for MLOD
//
//  Copyright (C) 2005 R. C. Elston


#include "output/Output.h"
#include "rped/rped.h"
#include "rped/genome_description.h"
#include "mlod/lod_table.h"

namespace SAGE
{
namespace MLOD
{
  
class LodTableFormatter
{
  public:
  
    LodTableFormatter(const RPED::RefMPedInfo&  trait_info,
                      const AnalysisParameters& params)
      : my_trait_info(&trait_info),
        my_params    (params)
      { }

                      
    template<class TABLE_TYPE>
    OUTPUT::Table formatTable(const TABLE_TYPE& table,
                              const string&     tname,
                              bool              include_intervals) const;
  
  protected:
  
    template<class TABLE_TYPE>
    void insert_table_row(OUTPUT::Table&    otable,
                          const TABLE_TYPE& info,     
                          size_t            point_idx,
                          double            position,
                          string            locus      = "") const;
                          
    void insert_trait_column_headers(OUTPUT::Table& table, size_t trait) const;

    const RPED::RefMPedInfo* my_trait_info;
    AnalysisParameters       my_params;
};
  

}}

#include "mlod/lod_table_print.ipp"

#endif

