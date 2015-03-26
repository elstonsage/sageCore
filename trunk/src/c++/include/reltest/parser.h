#ifndef  RELTEST_PARSER_H
#define  RELTEST_PARSER_H

//==========================================================================
//  File:       parser.h
//
//  Author:     Yeunjoo Song
//
//  History:    Initial implementation.                              Jul. 03
//
//  Notes:      This file defines a parser for reltest analysis.
//
//  Copyright (c) 2003 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "reltest/definitions.h"
#include "reltest/input.h"

namespace SAGE
{

namespace RELTEST
{

class reltest_parser
{
  public:

    reltest_parser();
    reltest_parser(const reltest_parser& f);

    ~reltest_parser();
   
    void                   build_default_reltest_analysis(reltest_data& rd);
    
    void                   parse_reltest_analysis(reltest_data& rd,
                                                  const LSFBase*      params,
                                                  cerrorstream&       err = sage_cerr);

    const string&          get_analysis_name()        const;

    const vector<bool>&    get_analysis_pairtypes()   const;
    const vector<size_t>&  get_analysis_regions()     const;
    const vector<double>&  get_preset_cutpoints()     const;

    bool                   calculate_cutpoints()      const;
    bool                   generate_nucfam_output()   const;
    bool                   generate_detailed_output() const;

    void                   view_parameter()           const;

  private:

    void                   build_region_list();

    void                   parse_pairtype(const LSFBase* param);
    void                   parse_region(const LSFBase* param);
    void                   parse_cutpoints(const LSFBase* param);
    void                   parse_nucfam(const LSFBase* param);
    void                   parse_detailed(const LSFBase* param);

    RPED::genome_description*   my_genome;

    string                my_analysis_name;         //name for current analysis.It is set to
                                                    //"default_analysis" unless user specified
                                                    //analysis is provided.
    
    vector<bool>          my_analysis_pairtypes;    //flags for the putative pair type to analyzed:

    vector<size_t>        my_analysis_regions;      //the indice of the regions will be used
                                                    //in analysis.

    vector<double>        my_preset_cutpoints;      //indicating cutpoints pre-setting status.

    bool                  my_calculate_cutpoints;
    bool                  my_nucfam_out;
    bool                  my_detailed_out;

    cerrorstream          errors;
};

#include "reltest/parser.ipp"

} // end of namespace RELTEST

} // end of namespace SAGE
                          
#endif
