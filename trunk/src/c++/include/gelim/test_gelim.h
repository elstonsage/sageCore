#ifndef GELIMTEST_H
#define GELIMTEST_H

#include <iostream>
#include <map>
#include "app/SAGEapp.h"
#include "error/errorstream.h"
#include "gelim/test_gelim_input.h"

namespace SAGE
{

class GELIMTEST : public APP::SAGEapp
{
public:

  GELIMTEST(int argc=0, char **argv=NULL)
    : APP::SAGEapp(APP::APP_SEGREG, false, argc, argv) { }

  ~GELIMTEST();

  virtual void print_help(std::ostream &);

  virtual int main();

private:

  void init();
  void init_output();

  gelim_data data;
};

// ================
// Inline functions
// ================

inline GELIMTEST::~GELIMTEST() { }

} // SAGE_NS_END

#endif
