#ifndef __LSF_INIT_H
#define __LSF_INIT_H

#include <iostream>
#include "LSF.h"
#include "LSFfactory.h"
#include "error/errorstream.h"

using namespace std;

void LSFInit()  
{ 
  Factory = new LSFFactory("Global Factory");

#if NEEDS_ERRORSTREAM_INIT
  SAGE::manual_errorstream_init();
#endif
}

#endif
