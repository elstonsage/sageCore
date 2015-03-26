#ifndef __LSF_FACTORY_H
#define __LSF_FACTORY_H

#include "LSF/LSF.h"

#ifdef NAMESPACE_STD
NAMESPACE_STD
#endif

class LSFFactory : public GenericFactory
{
public:
  LSFFactory(const char *n = NULL) : GenericFactory(n) 
                              { add_factory(NULL, lsf_mappings); }
  
protected:
  virtual LSFBase *build_local( lsf_t t, const char* name = NULL,
                                const char *type = NULL, AList *s = NULL, 
                                LSFBase* p = NULL );
};

extern LSFFactory* Factory; 

#endif
