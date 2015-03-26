#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <cassert>

#include "error/errorstream.h"
#include "error/errormanip.h"
#include "error/bufferederrorstream.h"

using namespace std;
using namespace SAGE;

void p(ostream &o, const string &s)
{
  o << s << flush;
}

void _print(const char* n, errorstream<char>& s)
{
  cout << n << ": " << s.rdbuf()->str() << ", " << s.str() << endl;
}

#define print(S)	_print(#S, S)

int main()
{
  free( malloc( 10 ) );

  string line;

  cout << "Test of string buffer base class and iostream insertion." << endl;
  errorstream<char> err;

  err << 'x' << "yz";
  print(err);

  err.str("");
//  err.rdbuf()->str("injection");
  err.str("injection");
  print(err);
  assert(err.rdbuf()->str() == "injection" && err.str() == "injection");

  err << "!";
  print(err);
//  assert(err.rdbuf()->str() == "injection!" && err.str() == "injection!");

  err.str("");
  err << "!";
  err << "injection3?";
//  assert(err.rdbuf()->str() == "!injection3?" && err.str() == "!injection3?");
  
  err.str("");

// Removed since they should always be valid anyway, and Windows doesn't
// like it.
//  assert(cout);
//  assert(err);

  cout << "Simple error messages: " << flush;
  err << location(__FILE__,__LINE__);
  err << suffix("]");
  err << prefix("[ERROR[%f(%l)] %d:") << priority(error) 
      << error_location << "This is an error"    << flush;
  err << prefix("[WARN[%f(%l)] %t:")  << priority(warning) 
      << error_location << "This is a warning"   << flush;
  err << prefix("[%P: ")      << priority(notice)  
      << "This is notice #"   << 3 << flush;

  err << suffix("- %t");
  err << prefix("%%[%P] ");

  errorstream<char> err0(err);
  err0 << "test of copy constructors..." << flush;
  basic_errorstream<char> err1(err);
  err1 << "test of copy constructors..." << flush;
  
  cout << endl 
       << "Test of mixing copy constructors of base and derived types"
       << endl;

  err1 << "Before copy.";
  err1.prefix("WARNING: ");

        errorstream<char> err2(err1);
  basic_errorstream<char> err3(err1);

  err2.prefix("Bad1");
  err3.prefix("Bad2");
  err1 << prefix("Original: ") << flush;
  err2 << prefix("Copy1: ") << "  After copy" << flush;
  err3 << prefix("Copy2: ") << "  After copy" << flush;

  
  cout << endl 
       << "Test of mixing assignment operators of base and derived types"
       << endl;

  basic_errorstream<char> err4;

//err1 << cooked_mode;
  err1 << "Before copy.";
  err1.prefix("WARNING: ");
  err4 = err1;
  err4 = err4;
  err3 = err1;
  err3 = err3;

  err4.prefix("Bad1");
  err3.prefix("Bad2");
  err1 << prefix("Original: ") << flush;
  err4 << prefix("Copy1: ") << "  After copy" << flush;
  err3 << prefix("Copy2: ") << "  After copy" << flush;

  cout << endl
       << "Testing error multistreams.  COUT = [none...warning]" << endl
       << "                             CERR = [error..critical]" << endl
       << endl;

  // Test to make sure everything is assignment compatible
  errormultistream<char> errs;

  err4 = errs;
  errs = errs;

  errorstream<char> err5(cout);
  errorstream<char> err6(cerr);

  err5 << prefix("COUT-%P: ") << line_width(79) << suffix(" [%f:%L]");
  err6 << prefix("CERR-%P: ") << line_width(79) << suffix(" [%f:%l]");

  errs.insert(err5);
  errs.restrict(r_lt,error);
  errs.insert(err6);
  errs.restrict(r_ge,error);

  errs << error_location;
  errs << priority(no_priority) << "This is a test 0" << flush;
  errs << priority(debug)    << "This is a test 1" << flush;

  errs << priority(information) << "This is a test 2" << flush;
  errs << priority(notice)   << "This is a test 3" << flush;
  errs << priority(warning)  << "This is a test 4" << flush;
  errs << priority(error)    << "This is a test 5" << flush;
  errs << priority(critical) << "This is a test 6" << flush;
  errs << priority(fatal)    << "This is a test 7" << flush;

  cout << endl
       << "Test of conversion to ostream &" << endl
       << endl;

  p(err5, "This is a test");
  p(err6, "This is a test2");
  p(errs, "This is a test3");
  
#if 1
  cout << endl
       << "Try your own error messages: (EOF to end, '.' on a line to flush)" 
       << endl << endl;

  while( cin )
  {
    getline(cin,line);
    if(!cin) break;
    if (line.size() == 0) continue;

    if( line != "." )
    {  if(line.size()) errs << line << '\n'; }
    else 
      errs << filename(__FILE__) << linenumber(__LINE__) << flush;
  }
#endif

  return 0;
}


