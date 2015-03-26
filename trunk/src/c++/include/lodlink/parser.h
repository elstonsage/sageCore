#ifndef LODLINK_PARSER_H
#define LODLINK_PARSER_H
//============================================================================
// File:      parser.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   9/5/2 created                     - djb
//                                                                          
// Notes      Defines class, parser, for parsing lodlink block of
//            a SAGE parameter file.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include <vector>
#include <string>
#include <limits>
#include <algorithm>
#include "rped/rped.h"
#include "app/aparser.h"
#include "LSF/LSF.h"
#include "LSF/Attr.h"
#include "mped/sp.h"
#include "mped/mp.h"
#include "maxfun/sub_model.h"
#include "lodlink/mle_sub_model.h"
#include "lodlink/instructions.h"
#include "lodlink/input.h"

namespace SAGE
{

namespace LODLINK
{

//----------------------------------------------------------------------------
//  Class:    parser
//                                                                          
//  Purpose:  extract information from lodlink  block of SAGE parameter
//            file.
//                                                                          
//----------------------------------------------------------------------------
//
class parser : public SAGE::APP::BasicParser
{
  public:
    
    // Constructor/destructor.
    parser(const RPED::RefMultiPedigree* mp, const lodlink_data& dat, ostream& messages = std::cout,
           SAGE::cerrorstream& errors = SAGE::sage_cerr);
  
    virtual void  parse_symbols(const SymbolTable* syms);
    virtual void  parse_parameter(const LSFBase* param);
    virtual void  parse_test_parameter_section(const LSFBase* params);
    virtual void  parse_test_parameter(const LSFBase* param);
    void  parse(const LSFBase* params);
    
    const instructions&  user_instructions() const;  
    size_t  analysis_id() const;

    // - Repository for lod sub-block pointers.
    //
    struct lods_ptrs
    {
      lods_ptrs();

      LSFBase*  option;
      LSFBase*  sex_specific;
      LSFBase*  male_female;
      LSFBase*  average;
    };
    
    
  private:
  
    void  init_parse();
    void  print_header();
    void  print_footer();

    // General ancillary functions.
    void  reset();
    RPED::RefTraitInfo::trait_t  primary_trait_type() const;
    void  check_trait_or_marker();
    void  check_sf_linkage_type_compatibility();
    void  check_trait_read_in();
    void  check_genotype_model_compatibility();
                         
    // Parsing functions.
    void  parse_out_attr(const LSFBase* param);
    void  parse_title(const LSFBase* param);
    
    void  parse_model(const LSFBase* param);
      bool  parse_marker_attr(const LSFBase* param);
      bool  parse_trait_attr(const LSFBase* param);
    
    void  parse_linkage_tests(const LSFBase* param);
    
    void  parse_homog_tests(const LSFBase* param);
      void  parse_smiths_test(const LSFBase* param);
      void  parse_mortons_test(const LSFBase* param);    
        bool  parse_group(const LSFBase* param);
          void  parse_pedigree_ids(const LSFBase* param, group& ids);
          void  verify_pedigrees(group& ids);
          bool  insert_group(const string& name, const group& g);
          bool  pedigrees_duplicated(const group& g) const;
          bool  groups_comprehensive() const;
        void  build_default_groups();
    
    void  parse_lods(const LSFBase* param);
      void  init_lods_ptrs(const LSFBase* param, lods_ptrs& pointers);
      string  parse_lods_option(const LSFBase* param);
      void  parse_lods_none(const lods_ptrs& pointers);
      void  parse_lods_standard(const lods_ptrs& pointers, bool sex_specific);
      void  parse_lods_specified(const lods_ptrs& pointers, bool sex_specific);
        void  parse_lods_thetas(const LSFBase* param);
          bool  theta_valid(double t);
          bool  not_a_duplicate(double t) const;
        void  parse_lods_theta_pairs(const LSFBase* param);
          void  parse_theta_pair(const LSFBase* param, theta_pair& p);
          bool  theta_pair_valid(theta_pair p);
          bool  not_a_duplicate(theta_pair p) const;
    
    void  parse_genotypes(const LSFBase* param);
    
    // Data members.
    instructions  my_instructions;
    
    size_t  my_analysis_id;       
    const RPED::RefMultiPedigree*  my_ref_mped;
    const lodlink_data&  my_data;
    ostream&  my_messages;
};

#include "lodlink/parser.ipp"

}
}

#endif

