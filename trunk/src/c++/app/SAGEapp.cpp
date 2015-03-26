#include <errno.h>
#include <iostream>
#include <iomanip>

#ifndef MSDOS
#include <unistd.h>
#else
#include "LSF/getopt.h"
#endif

#ifndef WINDOWS
#include <stdlib.h>
#endif

#include "app/SAGEapp.h"
#include "app/SAGEapp_version_bank.h"
#include "app/VersionNumber.h"

using namespace std;

namespace SAGE {
namespace APP  {

#ifndef CXXFLAGS
  const char *SAGEapp::cxxflags = NULL;
#else
  const char *SAGEapp::cxxflags = CXXFLAGS;
#endif

#ifndef BUILD
  const char *SAGEapp::build = NULL;
#else
  const char *SAGEapp::build = BUILD;
#endif

#ifndef LDFLAGS
  const char *SAGEapp::ldflags = NULL;
#else
  const char *SAGEapp::ldflags = LDFLAGS;
#endif

#ifndef LDLIBS
  const char *SAGEapp::ldlibs = NULL;
#else
  const char *SAGEapp::ldlibs = LDLIBS;
#endif

#ifndef BUILD_DATE
  const char *SAGEapp::build_date = NULL;
#else
  const char *SAGEapp::build_date = BUILD_DATE;
#endif


std::string SAGEapp::release_string;

std::string SAGEapp::getReleaseString() { return release_string; }

//======================================
//
//  checkExpirationStatus()
//
//======================================
void 
SAGEapp::checkExpirationStatus() const
{
  if(getYearsSinceBuild() > 4)
  {
    std::cout << "This version of S.A.G.E. has expired. Please contact S.A.G.E. to obtain the latest version." << std::endl;
    exit(0);
  }
}

//======================================
//
//  getYearsSinceBuild()
//
//======================================
int 
SAGEapp::getYearsSinceBuild() const
{
  time_t cur_time = time(NULL);

  struct tm * now = localtime(&cur_time);
    
  int current_year = 1900 + now->tm_year;
  
  std::string build_date_str(build_date);

  int build_year = atoi(build_date_str.substr(build_date_str.length() - 4).c_str());
          
  return current_year - build_year;
}


//===============================================
//
//  calculateReleaseString()
//
//===============================================
void
SAGEapp::calculateReleaseString() const
{
  // Process current time:

  time_t cur_time = time(NULL);

  struct tm * now = localtime(&cur_time);

  std::ostringstream cur_time_str;

  cur_time_str << now->tm_mday        << " ";

       if(now->tm_mon ==  0) cur_time_str << "Jan";
  else if(now->tm_mon ==  1) cur_time_str << "Feb";
  else if(now->tm_mon ==  2) cur_time_str << "Mar";
  else if(now->tm_mon ==  3) cur_time_str << "Apr";
  else if(now->tm_mon ==  4) cur_time_str << "May";
  else if(now->tm_mon ==  5) cur_time_str << "Jun";
  else if(now->tm_mon ==  6) cur_time_str << "Jul";
  else if(now->tm_mon ==  7) cur_time_str << "Aug";
  else if(now->tm_mon ==  8) cur_time_str << "Sep";
  else if(now->tm_mon ==  9) cur_time_str << "Oct";
  else if(now->tm_mon == 10) cur_time_str << "Nov";
  else if(now->tm_mon == 11) cur_time_str << "Dec";

  cur_time_str << " "
               << 1900 + now->tm_year << " "
               << (now->tm_hour < 10 ? "0" : "") << now->tm_hour        << ":"
               << (now->tm_min  < 10 ? "0" : "") << now->tm_min         << ":"
               << (now->tm_sec  < 10 ? "0" : "") << now->tm_sec;

  // Generate release string:

  release_string =
           toUpper(name)           + " -- " +
           cur_time_str.str()      + " -- " +
    "["  + version.toPrettyPrint() + 
    "; " + build_date              + "]\n";
}

//================================================
//
//  CONSTRUCTOR
//
//================================================
SAGEapp::SAGEapp(const app_index_type n, bool suppress_help_screen, int argc, char **argv)
{
  free(malloc(1));

  app_names::construct();

  name             = app_names::get_application_name(n);
  my_micro_version = app_names::get_micro_version_number(n);
  this->argc       = argc;
  this->argv       = argv;
  hlp              = 0;
  dbg              = 0;

  // Check license:
  checkExpirationStatus();

  // Calculate static release string:
  calculateReleaseString();

  // Print the title of the program
  print_title(cout);

  // Otherwise, check the parameters for flags
  size_t first_nonflag_arg = parse_params();
  
  // If it's greater than 1, then there WERE flags:
  if(first_nonflag_arg > 1)
  {
    size_t flag_count = first_nonflag_arg - 1;

    this->argv[first_nonflag_arg - 1]  = this->argv[0];  // Make the argument before the operands be the name of the program
    this->argc                        -= flag_count;     // Subtract the flag arguments
    this->argv                        += flag_count;     // Advance the array pointer past the flag arguments
  }

  // Set arg_count:
  arg_count = this->argc - 1;

  // Display debug, if requested:
  if(debug())
  {
    print_debug(cerr);
  }
  
  // If help print help message and exit
  if((suppress_help_screen == false) && help())
  {
    print_help(cerr);

    exit(EXIT_SUCCESS);
  }
}

//==================================
//
//  print_title(...)
//
//==================================
void SAGEapp::print_title(ostream &o)
{
  o << getReleaseString() << endl;

  o << "Remember you have agreed to add an appropriate statement (including the" << endl
    << "NIH grant number) under \"acknowledgments\" in any publication of results" << endl
    << "obtained by using this program. Suggested wording is:" << endl << endl
    << "\"(Some of)The results of this paper were obtained by using the software" << endl
    << "package S.A.G.E., which was supported by a U.S. Public Health Service" << endl
    << "Resource Grant (RR03655) from the National Center for Research Resources.\"" << endl
    << endl;
}

//===================================
//
//  print_debug(...)
//
//===================================
void SAGEapp::print_debug(ostream &o)
{
  o << "Application Information: ";

  if( !build_date && !cxxflags && !ldflags && !ldlibs && !build )
  {
    cout << "none. " << endl;
    return;
  }

  cout << endl;

  if(my_micro_version.size())          o << "    Version: " << version.toPrettyPrint()                 << endl;
  if(build_date && strlen(build_date)) o << "      Built: " << build_date                              << endl;
  if(build      && strlen(build))      o << "      BUILD: " << build                                   << endl;
  if(cxxflags   && strlen(cxxflags))   o << "   CXXFLAGS: " << cxxflags                                << endl;
  if(ldflags    && strlen(ldflags))    o << "    LDFLAGS: " << ldflags                                 << endl;
  if(ldlibs     && strlen(ldlibs))     o << "     LDLIBS: " << ldlibs                                  << endl;

  o << endl;
}

//=========================================
//
//  parse_params(...)
//
//=========================================
int SAGEapp::parse_params(const char *opts)
{
  // Process the args:
  bool help_found    = false,
       debug_found   = false;

  // Start at index #1, not index #0, because we have to skip the executed filename:
  size_t first_nonflag_index = 1;
  
  for(; first_nonflag_index < (size_t)argc; ++first_nonflag_index)
  {
    std::string arg = argv[first_nonflag_index];
    
         if(arg == "-h" || arg == "-?") { help_found    = true; }
    else if(arg == "-@")                { debug_found   = true; }
    else                                { break;                }
  }

  hlp  += help_found;
  dbg  += debug_found;
  
  // Return the begin index for the operand portion:
  return first_nonflag_index;
}    

//==================================
//
//  print_inf_banner(...)
//
//==================================
void SAGEapp::print_inf_banner(ostream &o)
{
  o << endl
    << "********************************************************************" << endl
    << "*****                                                          *****" << endl
    << "*****    Always check the INF file before viewing results.     *****" << endl
    << "*****                                                          *****" << endl
    << "********************************************************************" << endl
    << endl;

}

} // End namespace APP
} // End namespace SAGE
