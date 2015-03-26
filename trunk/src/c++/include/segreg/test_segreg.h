#ifndef SEGREGTEST_H
#define SEGREGTEST_H

#include <iostream>
#include <map>
#include "app/SAGEapp.h"
#include "error/errorstream.h"
#include "segreg/segreg_input.h"
#include "segreg/sub_model_base.h"

namespace SAGE {

class SEGREGTEST : public APP::SAGEapp
{
public:

  SEGREGTEST(int argc=0, char **argv=NULL);

  virtual int main();

private:

};

} // End namespace

#endif
