#ifndef GELIMTEST_H
#define GELIMTEST_H

#include <iostream>
#include <map>
#include "app/SAGEapp.h"
#include "error/errorstream.h"
#include "peeling/test_peeler_input.h"

namespace SAGE
{

class TESTPEELER : public APP::SAGEapp
{
public:

  TESTPEELER(int argc=0, char **argv=NULL)
    : APP::SAGEapp(APP::APP_SEGREG, false, argc, argv) { }

  ~TESTPEELER();

  // Print program information
  virtual void print_help(std::ostream &);

  virtual int main();

private:

  void init();
  void init_output();

  test_peeler_data data;
};

// ================
// Inline functions
// ================

inline TESTPEELER::~TESTPEELER() { }

}

#endif
