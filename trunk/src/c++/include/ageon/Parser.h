#ifndef AO_PARSER_H
#define AO_PARSER_H
//============================================================================
// File:      Parser.h
//                                                                          
// Author:    Stephen Gross
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
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
#include "ageon/Model.h"
#include "ageon/Datatypes.h"

namespace SAGE {
namespace AO   {

//============================================================================
//
//  class Parser
//
//============================================================================
class Parser : public APP::BasicParser
{
  public:
    //========================================================================
    // Constructor:
    //========================================================================

    Parser(const RPED::RefMultiPedigree & RMP, 
                    ostream          & messages = std::cout,
                    cerrorstream     & errors   = sage_cerr);
  
    //========================================================================
    // Public utility functions:
    //========================================================================

    virtual void parse_test_parameter_section (const LSFBase     * params);
    virtual void parse_symbols                (const SymbolTable * syms)  {}
    virtual void parse_parameter              (const LSFBase     * param) {}
    virtual void parse_test_parameter         (const LSFBase     * param) {}
    
    //========================================================================
    // Public accessors:
    //========================================================================

    const Model & get_model       () const;  
          size_t     analysis_number () const;
    
  private:
    //========================================================================
    // Private utility functions:
    //========================================================================
    
    void init_parse   ();
    void print_title  ();
    void print_footer ();
    void reset        (); 

    void process_parameters        (const LSFBase * param);
    void parse_title               (const LSFBase * param);

    void parse_affectedness_trait  (const LSFBase * param);
    void parse_age_of_onset_trait  (const LSFBase * param);
    void parse_age_of_exam_trait   (const LSFBase * param);

    string parse_trait(const LSFBase * param, string a_name);

    void parse_covariate_sub_block (const LSFBase * param, trait_type t);
    void parse_covariate           (const LSFBase * param, trait_type t);

    void parse_transform_sub_block (const LSFBase * param);
    void parse_allow_averaging     (const LSFBase * param);
    void parse_pool                (const LSFBase * param);

    //========================================================================
    // Data members:
    //========================================================================

    const RPED::RefMultiPedigree & my_RMP;
          ostream          & my_messages;
          cerrorstream     & my_errors;
          Model           my_model;
          size_t             idnum;
};

inline const Model & Parser::get_model       () const { return my_model; }
inline       size_t     Parser::analysis_number () const { return idnum;    }

}} // End namespace

#endif
