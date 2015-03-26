#ifndef __SAGE_APP_H
#define __SAGE_APP_H

#include <iostream>
#include <string>
#include <time.h>
#include "LSF/Attr.h"
#include "app/SAGEapp_version_bank.h"

// One of the application's main duties is to maintain the licensing
// agreement. This is done through FlexLM's license manager.  However,
// licensing only needs to be done for certain build types.  Thus
// it is all controlled by a environment variable called SAGE_DEFINED

#if defined(SAGE_LICENSED)

// We must include some extra stuff if we're doing a windows (MingW) build

#  if defined(MINGW)

#    include <windows.h>

     // Clean up namespace spewage from Microsoft.

#    ifdef max
#      undef max
#    endif

#    ifdef min
#      undef min
#    endif

#    ifdef ERROR
#      undef ERROR
#    endif

#    ifdef ESTIMATE
#      undef ESTIMATE
#    endif

#    ifdef FIXED
#      undef FIXED
#    endif

#  endif

#else
#  define LM_HANDLE unsigned int
#endif

namespace SAGE {
namespace APP  {

/** \class SAGEapp
  * \brief Application framework for all SAGE programs
  *
  * The SAGEapp is the application framework for all SAGE programs. 
  * Individual programs are derived from this as a baseclass by public
  * inheritance.
  *
  * Functions the SAGE app provides:
  *
  * 1. Basic command line parsing for verbose and debug flags
  *
  * 2. Static information about the SAGE program (name, release, copyright, etc)
  *
  * 3. The basic output streams automatically generated for each SAGE program.
  *
  */
class SAGEapp
{
public:

  /// @name Versioning
  //@{
  
    ///
    /// Returns a pretty-print string containing the following information:
    /// version number, program name, build date, current date, and copyright.
    ///
   /// For example:
    /// \verbatim
    /// AGEON Output -- 18 Mar 2005 08:52:59 -- [S.A.G.E. v5.0.3; bld 17 Mar 2005]
    /// COPYRIGHT (C) 2005 CASE WESTERN RESERVE UNIVERSITY
    /// \endverbatim
    ///
    /// Note: This function is generally called by print_title(); if, however, you
    /// need the release string in a static context you can use this function.
    static std::string getReleaseString();
  
  //@}

  /// @name Constructor
  //@{

    ///
    /// Constructor
    /// \param program_index The index number of the program
    /// \param suppress_help_screen If true, the -h or -? will be ignored; set to true if you intend to use the ArgumentRuleset system in the data library.
    /// \param argc Number of arguments passed on the commandline
    /// \param argv Argument list from the commandline
    SAGEapp(const app_index_type program_index, bool suppress_help_screen = false, int argc=0, char **argv=NULL);

  //@}

  /// @name Required virtual interface
  //@{

    ///
    /// Run the application.
    virtual int main() = 0;
  //@}

    ///
    /// Prints the release version (obtained via getReleaseVersion()) to the
    /// given output stream.
    void print_title(std::ostream & os);

    ///
    /// Prints the banner reminding user to check *.INF file to the
    /// given output stream.
    void print_inf_banner(ostream &o);

  /// @name Option virtual interface
  //@{

    ///
    /// Destructor
    virtual inline ~SAGEapp();
  
    ///    
    /// Print program information
    /// \param os Outputstream to which output will be directed.
    virtual void print_debug(std::ostream & os);

    ///    
    /// Print program information
    /// \param os Outputstream to which output will be directed.
    virtual inline void print_help(std::ostream & os) { ; }

  //@}

  /// @name Accessors
  //@{

    ///
    /// Accessor for user defined command line options
    const AttrList & params();

  //@}

  /// @name Utility functions
  //@{

    ///
    /// Parse command line.
    /// The commandline should follow the form:
    ///
    /// [sage_app] [-v | -h | -?] [-@] [input file(s)...]
    ///
    int parse_params(const char *opts = "vh@");

  //@}

  /// @name Flags for basic run mode options
  //@{

    ///
    /// Returns the help status
    inline bool help() const;

    ///
    /// Sets the help status
    /// \param h The new help status
    inline bool help(bool h);

    ///
    /// Returns the debug status
    inline bool debug() const;

    ///
    /// Sets the debug status
    /// \param d The new debug status
    inline bool debug(bool d);

  //@}

  /// @name Basic program information
  //@{

    ///
    /// The date on which the program was built (compiled)
    static const char *build_date;

    ///
    /// Copyright information
    static const char *copyright;

    ///
    /// C++ compiler flags used at compile time
    static const char *cxxflags;

    ///
    /// Linker flags used at compile time
    static const char *ldflags;

    ///
    /// External libraries linked during build
    static const char *ldlibs;

    ///
    /// Type of build ("DEBUG", "PROFILE", "RELEASE", etc.)
    static const char *build;

  //@} 

protected:
  int      argc;
  char**   argv;
  size_t   arg_count;
  AttrList opts;

  string   name;
  string   my_micro_version;
  
private:

  static std::string release_string;

  int hlp;
  int dbg;

  /// \internal
  /// Populates the static string 'release_string' with the app-specific value.
  /// It is important to invoke this function somewhere in SAGEapp's construction,
  /// because subsequent calls to getReleaseString() assume that the string has
  /// actually been populated with content.
  void calculateReleaseString() const;

  /// @name New license management
  //@{

    ///
    /// Finds the installation path for \b this version of SAGE.
    /// This work by the following algorithm:
    ///
    /// First, search for a sagerc file (see SageRcFile for more information).
    ///
    /// If a sagerc file is found, use the SAGE_PATH declared therein as the installation path.
    ///
    /// If a sagerc file is not found, check the environment's SAGE_PATH variable (or registry key,
    /// if this is a windows environment).
    ///
    /// Otherwise, scan through the sage installations listed in the SAGE_PATH variable; match
    /// the correct SAGE installation to this version, and return that string.
    ///
    /// If the SAGE_PATH variable is not found, or there is no entry in it corresponding to this
    /// version, produce an error.
    /// \param check_for_rcfile By default, this function will scan for a .sagerc file. You can set
    /// this parameter to \c false, however, to override this feature (and thus skip checking for
    /// an rc file).
    std::string getInstallationPath(bool check_for_rcfile = true) const;
    
  //@}

  /// @name Expiration stuff
  //@{
  
    ///
    /// Checks that the program has not expired; exits if it has.
    void checkExpirationStatus() const;
  
    ///
    /// Returns the number of years since the build date.
    int getYearsSinceBuild() const;
    
  //@}

};

//==================================================================================
// INLINE FUNCTIONS
//==================================================================================

inline SAGEapp::~SAGEapp() 
{
  argv = NULL;
}
  
//inline void SAGEapp::print_help (std::ostream &) {}

inline const AttrList & SAGEapp::params     ()             { return opts;         }
inline       bool       SAGEapp::help       ()       const { return !!hlp;        }
inline       bool       SAGEapp::help       (bool h)       { return !!(hlp = h);  }
inline       bool       SAGEapp::debug      ()       const { return !!dbg;        }
inline       bool       SAGEapp::debug      (bool d)       {return !!(dbg = d);   }

} // End namespace APP
} // End namespace SAGE

#endif
