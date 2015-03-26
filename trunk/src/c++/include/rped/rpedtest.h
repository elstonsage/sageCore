#ifndef RPEDTEST_H
#define RPEDTEST_H

#include <iostream>
#include <map>
#include "app/SAGEapp.h"
#include "error/errorstream.h"
#include "rpedtest_input.h"

namespace SAGE {
namespace RPED {

class RPEDTEST : public SAGEapp
{
public:

  RPEDTEST(int argc=0, char **argv=NULL)
    : SAGEapp(argc, argv) { }

  ~RPEDTEST();

  // Print program information
  virtual void print_title(std::ostream &);
  virtual void print_help(std::ostream &);

  virtual int main();
  
private:

  void init();
  void init_output();

  rptest_data data;
};

// ================
// Inline functions
// ================

inline RPEDTEST::~RPEDTEST() { }

} // End namespace RPED
} // End namespace SAGE

#endif
