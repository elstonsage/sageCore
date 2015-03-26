#include <memory>
#include "error/errorstream.h"

namespace SAGE {

//lint -e{1502} Warning about no non-static members supressed
std::ios::Init init_std_ios;

//lint -e{1502} Warning about no non-static members supressed
static ErrorNoInit no_init = ErrorNoInit();

cerrorstream sage_cerr( no_init );
cerrorstream sage_cout( no_init );
cerrorstream sage_clog( no_init );

//lint -e{1502} Warning about no non-static members supressed
errorstream_init error_init;

long errorstream_init::count;

errorstream_init::errorstream_init()
{
  if (count++ == 0)
    manual_errorstream_init();
}

errorstream_init::~errorstream_init()
{
  if (--count == 0)
  {
    sage_cerr << std::flush;
    sage_cout << std::flush;
    sage_clog << std::flush;
  }
}

void manual_errorstream_init()
{
    //lint --e{522} Suprssed message about expressions

    new ((void*)&sage_cout) cerrorstream(std::cout);
    new ((void*)&sage_clog) cerrorstream(std::clog);
    new ((void*)&sage_cerr) cerrorstream(std::cerr);
}

}
