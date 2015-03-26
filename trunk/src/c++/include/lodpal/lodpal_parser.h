#ifndef LODPAL_PARSER_H
#define LODPAL_PARSER_H

//=============================================================================
// File:    lodpal_parser.h
//
// Author:  Yeunjoo Song
//
// History: Version 0.0 Initial implementation.                   yjs Jan 01
//
// Notes:   This file contains definition for following data structure.
//
//            class LodpalParser : public APP::BasicParser
//
// Copyright (c) 2001 R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "app/aparser.h"
#include "lodpal/lodpal_params.h"

namespace SAGE   {
namespace LODPAL {

class lodpal_parser : public APP::BasicParser
{
  public:

    lodpal_parser(const RelativePairs& relpairs, cerrorstream& err = SAGE::sage_cerr)
            : APP::BasicParser(err), pairs(relpairs)
    {
      clear();
    }                              

    lodpal_parser(const RelativePairs& relpairs, SymbolTable* syms, 
                  cerrorstream& err = SAGE::sage_cerr) 
            : APP::BasicParser(err), pairs(relpairs)
    {
      clear();
      parse_symbols(syms);
    }                              

    lodpal_parser(const RelativePairs& relpairs, LSFBase* params, 
                  cerrorstream& err = SAGE::sage_cerr)
            : APP::BasicParser(err), pairs(relpairs)
    {
      clear();
      parse_test_parameter_section(params);
    }                              

    lodpal_parser(const RelativePairs& relpairs, SymbolTable* syms, 
                  const LSFBase* params, cerrorstream& err = SAGE::sage_cerr)
            : APP::BasicParser(err), pairs(relpairs)
    {
      clear();
      parse_symbols(syms);
      parse_test_parameter_section(params);
    }                              

    void clear()
    {
      APP::BasicParser::clear();

      my_wide_output              = false;
      my_csv_output               = false;
      my_pval_scientific_notation = false;
      
      my_diagnostic               = false;
      my_discordant               = false;

      my_autosomal_model_exist    = false;
      my_one_parameter_model      = true;
      my_alpha                    = 2.634;
    }
     
    void parse_symbols(const SymbolTable* syms);
    void parse_parameter(const LSFBase* param);
    void parse_test_parameter_section(const LSFBase* params);
    void parse_test_parameter(const LSFBase* param);

    const lodpal_parameters&         parameters() const { return my_parameters; }
    const vector<LSF_ptr<LSFBase> >& get_pair_info_file() const { return my_pair_info_file; }

    bool wide_output()                      const { return my_wide_output; }
    bool csv_output()                       const { return my_csv_output; }
    bool get_pval_scientific_notation()     const { return my_pval_scientific_notation; }

    bool multipoint()                       const { return my_parameters.multipoint(); }
    bool diagnostic()                       const { return my_diagnostic; }
    bool pair_info_file()                   const { return my_pair_info_file.size() ? true : false; }

    void set_multipoint(bool b)                   { my_parameters.set_multipoint(b); }

  protected:

    void parse_trait(const LSFBase* param);
    void parse_subset(const LSFBase* param);
    void parse_marker(const LSFBase* param);         // parse marker & model for this marker
    void parse_covariate(const LSFBase* param);      // parse covariate
    void parse_weight(const LSFBase* param);         // parse weight
    void parse_alpha(const LSFBase* param);          // parse alpha value
    void parse_two_model(const LSFBase* param);      // parse global model for all markers
    void parse_unconstrained(const LSFBase* param);  // parse global model for all markers
    void parse_diagnostic(const LSFBase* param);     // parse diagnostic output option
    void parse_x_linkage_model(const LSFBase* param);
    void parse_autosomal_model(const LSFBase* param);

    bool   discordant()                     const { return my_discordant; }
    bool   one_parameter_model()            const { return my_one_parameter_model; }
    double alpha()                          const { return my_alpha; }

    void set_discordant(bool b)                   { my_discordant = b; }
    void set_one_parameter_model(bool b)          { my_one_parameter_model = b; }
    void set_alpha(double d)                      { my_alpha = d; }

    void set_wide_output(bool b)                  { my_wide_output = b; }
    void set_csv_output(bool b)                   { my_csv_output = b; }
    void set_pval_scientific_notation(bool b)     { my_pval_scientific_notation = b; }

    void set_diagnostic(bool b)                   { my_diagnostic = b; }

  private:

    const RelativePairs&   pairs;

    bool   my_wide_output;
    bool   my_csv_output;
    bool   my_pval_scientific_notation;

    bool   my_diagnostic;
    bool   my_discordant;

    bool   my_autosomal_model_exist;
    bool   my_one_parameter_model;
    double my_alpha;
      
    lodpal_parameters my_parameters;

    vector<LSF_ptr<LSFBase> >  my_pair_info_file;
};

} // end of namespace LODPAL
} // end of namespace SAGE

#endif
