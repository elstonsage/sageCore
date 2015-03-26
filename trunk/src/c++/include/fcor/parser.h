#ifndef FCORPARSER_H
#define FCORPARSER_H

//****************************************************************************
//* File:      parser.h                                                      *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0                                                   *
//*                    0.1  yjs  Added parsing for var-cov option  Jan 23 01 *
//*                                                                          *
//* Notes:     This header file defines functions to parse fcor parameters.  *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "fcor/definitions.h" 

namespace SAGE {
namespace FCOR {

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     FcorParser                                                   ~
// ~                                                                         ~
// ~ Purpose:   Defines functions to parse fcor analysis parameters.         ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class FcorParser
{
  public:
  
    typedef std::map<string, string, less<string> > name_map_type;
    
    FcorParser(const RPED::RefMultiPedigree& mp, cerrorstream& err = sage_cerr);
    FcorParser(const FcorParser& f);

    void    build_default_fcor_analysis();

    void    parse_fcor_analysis(const LSFBase* fcor_params);

    void    build_relationship_name_map();
    
    string                         get_weight(bool intraclass) const;

    const RPED::RefMultiPedigree*  get_multi_pedigree()        const;

    const vector<name_index_type>& get_trait_list()            const;  

    const analysis_option_type&    get_analysis_options()      const;

    const name_map_type&           get_name_map()              const;

    size_t                         get_trait_count()           const;
    size_t                         get_pedigree_count()        const;

    void                           view_parameter()            const;

  private:
    
    void           build_trait_list();

    void           parse_class_weight     (const AttrVal& attrv, weight_type& w);
    void           parse_generation_limit (const AttrVal& attrv);
    void           parse_pairset          (const LSFBase* param);
    void           parse_homogeneity_test (const LSFBase* param);
    void           parse_var_cov          (const LSFBase* param);
    void           parse_output_options   (const LSFBase* param);
    void           parse_boolean_value    (const AttrVal& attrv, bool& b);

    const RPED::RefMultiPedigree*        my_multipedigree;

    vector<name_index_type>              my_trait_list;

    analysis_option_type                 my_analysis_options;

    name_map_type                        my_relationship_name_map;

    cerrorstream                         errors;
};    
    
#include "fcor/parser.ipp"
  
} // end of namespace FCOR
} // end of namespace SAGE

#endif
