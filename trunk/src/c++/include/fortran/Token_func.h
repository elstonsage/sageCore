#ifndef __TOKEN_FUNC_H_
#define __TOKEN_FUNC_H_

#include "fortran/Tokenizer.h"

#ifdef NAMESPACE_STD
NAMESPACE_STD
#endif

AttrVal& get_token(Tokenizer::token& a, Tokenizer* t);
AttrVal put_token(Tokenizer::token& a, OutputTokenizer* t);

#endif
