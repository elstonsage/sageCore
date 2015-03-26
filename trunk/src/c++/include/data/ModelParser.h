#ifndef DATA_MODEL_PARSER_H
#define DATA_MODEL_PARSER_H

#include <vector>
#include "error/internal_error.h"
#include "app/aparser.h"

namespace SAGE {
namespace APP  {

class ModelParser : public APP::BasicParser
{
public:

  ModelParser (cerrorstream & err = SAGE::sage_cerr);

  template<class ANALYSIS_TYPE>
  void operator() (ANALYSIS_TYPE & analysis, const LSFBase * block) const;

  virtual void  parse_symbols (const SymbolTable *syms)              {}
  virtual void  parse_parameter (const LSFBase *param)               {}
  virtual void  parse_test_parameter_section (const LSFBase *params) {}
  virtual void  parse_test_parameter (const LSFBase *param)          {}

private:
};

inline ModelParser::ModelParser(cerrorstream & err) : BasicParser(err)
{}

template<class ANALYSIS_TYPE>
inline void
ModelParser::operator() (ANALYSIS_TYPE & analysis, const LSFBase * block) const
{
  // Note: This non-specialized function should never be invoked. If it is invoked,
  // that means the programmer forgot to implement his own version.
   
  SAGE_internal_error();
}

}} /// End namespace

#endif
