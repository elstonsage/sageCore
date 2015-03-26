#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <cassert>

#pragma hdrstop

#include "error/errorstream.h"
#include "error/errormanip.h"
#include "error/bufferederrorstream.h"

using namespace std;
using namespace SAGE;

void test_stream(cerrorstream& errs)
{
  string message = "This is a %P message.";
  errs << priority(information) << message << flush;
  errs << priority(notice)   << message << flush;
  errs << priority(warning)  << message << flush;
  errs << priority(error)    << message << flush;
  errs << priority(critical) << message << flush;
  errs << priority(fatal)    << message << flush;
}


void buff0()
{
  cout << "Test 0: Testing standard error stream.  All messages should print."
       << endl << endl;

  errorstream<char> buff;
  cout << "Generating Error Messages" << endl;
  test_stream(buff);

  cout << endl << endl;
}

void buff1()
{
  cout << "Buffer Test 1: " << endl << endl;

  bufferederrorstream<char> buff(sage_cout);

  cout << "Generating Error Messages" << endl;

  test_stream(buff);

  cout << "Messages Generated" << endl;

  buff.flush_buffer();

  cout << endl;
}

void buff2()
{
  cout << "Buffer Test 2: " << endl << endl;

  bufferederrorstream<char> buff(sage_cout);

  errormultistream<char> multi;

  multi.insert(buff);

  multi.restrict(r_lt, error);

  multi.insert(sage_cout);

  multi.restrict(r_ge, error);

  cout << "Generating Error Messages" << endl;

  test_stream(multi);

  cout << "Messages Generated" << endl;

  buff.flush_buffer();

  cout << "Now?" << endl;

  cout << endl;
}

int main()
{
  free(malloc(1));
  buff0();
  buff1();
  buff2();
}
