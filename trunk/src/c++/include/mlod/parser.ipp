//===================================================
// File:      parser.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   5/28/02 created        -djb
//                                                                          
// Notes:     Inline implementation for class, parser.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

//============================================================================
// IMPLEMENTATION:  parser
//============================================================================
//
inline
Parser::Parser(const Data& data, ostream& messages, cerrorstream& err)
  : BasicParser(err),
    my_data(data),
    my_messages(messages),
    my_analysis_id(-1)
{
  reset_parameters();
}              

inline void
Parser::reset_parameters()
{
  ++my_analysis_id;
  my_parameters = AnalysisParameters((unsigned int)my_analysis_id);
}

inline void
Parser::parse_symbols(const SymbolTable* )
{}

inline void 
Parser::parse_parameter(const LSFBase* )
{}

inline void 
Parser::parse_test_parameter(const LSFBase* )
{}

inline const AnalysisParameters&
Parser::parameters() const
{
  return my_parameters;
}

inline size_t
Parser::analysis_id() const
{
  return (size_t) my_analysis_id;
}

