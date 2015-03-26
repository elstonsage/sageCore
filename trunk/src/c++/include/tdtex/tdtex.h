#ifndef TDTEX_H
#define TDTEX_H
//============================================================================
// File:      tdtex.h
//
// Author:    Kevin Jacobs (jacobs@theopalgroup.com)
//            Stephen Gross
//
// CVS:       $Id: tdtex.h 7964 2007-07-18 20:48:53Z kcartier $
//
// History:   6/3/2002 Initial version
//
// Notes:
//
// Copyright (c) 2002 R.C. Elston
// All Rights Reserved
//============================================================================


#include <cassert>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

#include "LSF/LSFsymbol.h"
#include "LSF/LSFinit.h"
#include "LSF/LSF.h"
#include "LSF/LSFfile.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "data/SAGEdata.h"
#include "app/SAGEapp.h"
#include "app/SAGEapp.h"
#include "rped/rpfile.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "numerics/print_util.h"
#include "output/Output.h"
#include "tdtex/Sampler.h"
#include "tdtex/Parser.h"
#include "tdtex/Tests.h"

// #define KCC_DEBUG(x) std::cout << "\n" << #x " = " << x << " at " << __FILE__ << "(" << __LINE__ << ")" << "\n"
// #define KCC_PING() std::cout << "\n--> ***** PING ***** at " << __FILE__ << "(" << __LINE__ << ")" << "\n"

namespace SAGE  {
namespace TDTEX {

class Sampler;
class Parser;

/** \class TdtexApp
  * \brief Defines tdtex application derived from SAGEapp
  */
class TdtexApp : public APP::SAGEapp
{
public:

  /// @name Constructor
  //@{

    ///
    /// Constructor.
    /// \param argc The number of arguments (taken from main())
    /// \param argv The argument list (taken from main())
    TdtexApp(int argc=0, char **argv=NULL);

  //@}

  /// @name Required virtual interface
  //@{


    ///
    /// Main function (starts the program)
    virtual int main();

  //@}

  /// @name Analysis execution functions
  //@{

    ///
    /// Analyzes a single marker on a single sex. Please note this functions assumes the marker for analysis
    /// has \b already been set (via the sampler.set_marker() function), as well as the parent sex of interest
    /// (via the sampler.set_parent_sex() function).
    /// \param mped The multipedigree to analyses
    /// \param parser ???
    /// \param sample The sampler that describes the TDT tables
    /// \param out The output stream to which output should be directed
    void run_single (const Configuration & config, const RPED::RefMultiPedigree& mped, MPED::SexCode sex, cerrorstream& err, OUTPUT::Section & section);

    ///
    /// Analyzes a single marker & single trait on each possible sex (unknown, male, and female).
    /// Please note this functions assumes the marker for analysis
    /// has \b already been set (via the sampler.set_marker() function).
    /// \param mped The multipedigree to analyses
    /// \param parser ???
    /// \param sample The sampler that describes the TDT tables
    /// \param out The output stream to which output should be directed
    void runOneMarkerOneTrait(const Configuration & config, const RPED::RefMultiPedigree & mped,  cerrorstream& err, OUTPUT::Section & section);

    void runAllMarkersAllTraits(Configuration & config, const RPED::RefMultiPedigree & mped,  cerrorstream& err, OUTPUT::Section & section);

  //@}
};


/* ************************************************************************** **
** *****************          Class Specification           ***************** **
** ************************************************************************** */
class AppData : public APP::SAGE_Data {
   /* --------------------------------------------------------------------------
      Friend Functions and Classes
   -------------------------------------------------------------------------- */
   // None

   public:
      /* -----------------------------------------------------------------------
         Methods
      ----------------------------------------------------------------------- */
      // Constructor
      AppData(const string & program_name, bool debug);

      // Utility functions
      void process_input(int argc, char** argv);
      virtual bool read_analysis() { return(true); }


      /* -----------------------------------------------------------------------
         Members
      ----------------------------------------------------------------------- */
      // None

   protected:
      /* -----------------------------------------------------------------------
         Methods
      ----------------------------------------------------------------------- */
      // None

      /* -----------------------------------------------------------------------
         Members
      ----------------------------------------------------------------------- */
      // None

   private:
      /* -----------------------------------------------------------------------
         Methods
      ----------------------------------------------------------------------- */
      // None

      /* -----------------------------------------------------------------------
         Members
      ----------------------------------------------------------------------- */
      // None


}; /* class AppData */


} // End namespace TDTEX
} // End namespace SAGE

#endif
