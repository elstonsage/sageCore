#ifndef OUTPUT_STREAMS_H
#define OUTPUT_STREAMS_H

//============================================================================
// File:      output_streams.h
//                                                                          
// Author:    
//                                                                          
// History:   3-27-01 modified to add test of genome_description class.  - djb                                                   
//            May  01 Generalized so all programs use the same interface - gcw
//            Dec  01 Split out the Output options into their own class  - gcw
//                                                                          
// Notes:     
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include <iostream>
#include <fstream>
#include <string>
#include "error/errorstream.h"
#include "error/errorbuf.h"

namespace SAGE {
namespace APP  {

/** \brief Provides basic output streams for SAGE programs
 * 
 * \par INTRODUCTION
 *
 * The Output_Streams class is the basic output streams needed for most SAGE programs.
 * 
 * It provides four basic 'channels' (screen, information, messages, errors)
 * to which messages can be sent.  The messages are multiplexed to several
 * output channels (Screen, Inf file, Debug file, sage_cerr)
 *
 * Of the four channels, only the errors() stream is an error stream.  This
 * means that only the errors stream is formatted.  The other streams should
 * be considered to be raw streams.
 *
 * Note also that the messages stream is simply a multiplex of the screen and information
 * streams for ease of use.  If a message should go to both of these streams it should
 * be sent to messages, which will take care of the issues.
 *
 * \par FUTURE WORK
 *
 * We would like to extend this class to allow for additional streams to be
 * added at runtime.  In other words, it'd be nice to be able to add string
 * identified streams and/or error streams to the class that can be looked
 * up and utilized by the user classes.  This would allow the addition of
 * analysis streams to the mix, or specialized streams that could then be
 * passed around easily.
 *
 * In the more distant future, it will be possible to have the streams
 * create dialog boxes and other such fetures.  This would require making a
 * module to do this when certain kinds of messages are received.
 *
 * \par MORE DETAIL
 *
 * The main outputs from any S.A.G.E. program are/should be:
 *
 * 1 Screen output (cout)
 *
 * 2 Screen error output (sage_cerr)
 *
 * 3 Debug stream
 *
 * 4 Information Output (.inf file)
 *
 * 5 Analysis Output(s) (various)
 *
 * 6 One or more specialized streams for specific elements (ex: genome.inf)
 *
 * What goes to each stream:
 *
 * 1. Screen Output
 *
 * Primarily, the screen output provides current status to the user.  This includes which data files are read in and 
 * what analyses are to be performed, as well as status of each analysis.
 *
 * 2. Screen error output 
 *
 * Any error with a code = error should go to the screen.  These errors are things which may invalidate the analysis.  
 * Errors < error are warning, information and debugging messages.  If Debugging mode is on, this should include debug 
 * messages as well.
 *
 * 3. Debug Stream
 *
 * When Debugging is on, there is an additional stream.  This should be the info file, plus debug messages.  The debug 
 * messages should also be passed directly to the screen.
 *
 * 4. Information Output
 *
 * The information output file includes information on the running of the program.  As such, it should include 
 * everything that goes to the screen, and perhaps some additional details.  Errors with a code = information should go to this 
 * file.
 *
 * 5. Analysis Outputs
 *
 * The results of the analyses.  These files should include the tables of results as well as any errors which might 
 * affect the validity of those results.
 *
 * Diagram of the standard streams and where output goes:
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * Note:  The error stream is split among the various outputs, based on the criteria given.  If debugging is not 
 * active, then the debugging critera are removed and the debug stream is not created.
 */  
class Output_Streams
{
public:

  //lint --e{1712} no default constructor ok
  //lint --e{1536} returning access to the streams is ok

  /// @name Constructor & destructor
  //@{

    ///
    /// Constructor
    /// \param program_name The name of the program
    /// \param debug_on Indicates whether or not to include debugging information
    Output_Streams(const std::string& program_name, bool debug_on = false);

    ///
    /// Destructor
    ~Output_Streams();

  //@}

  /// @name Streams available for use
  //@{

    ///
    /// Returns an outputstream for screen output.
    std::ostream& screen();

    ///
    /// Returns an outputstream for informative output.
    std::ostream& info();

    ///
    /// Returns an outputstream for messages.
    std::ostream& messages();

    ///
    /// Returns an outputstream for debugging output.
    std::ostream& debug();

    ///
    /// Returns an outputstream for error output.
    cerrorstream& errors();

  //@}

  /// @name Debugging info
  //@{

    ///
    /// Returns the debugging status.
    bool get_debug_status() const;
  
  //@}

private:

  // Default initialization of output streams
  
  void init_output_streams();

  std::string  my_program_name;
  bool         my_debug;

  cerrormultistream my_screen_stream;
  cerrormultistream my_information_stream;
  cerrormultistream my_message_stream;
  cerrormultistream my_error_stream;
  cerrormultistream my_debug_stream;

  std::ofstream inf_file;
  std::ofstream dbg_file;

  // Prohibited functions

  Output_Streams(const Output_Streams&);

  Output_Streams& operator=(const Output_Streams&);
};

} // End APP  namespace
} // End SAGE namespace

#include "app/output_streams.ipp"

#endif
