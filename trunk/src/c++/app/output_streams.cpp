#include "app/output_streams.h"
#include "error/errormanip.h"
#include "LSF/parse_ops.h"

namespace SAGE {
namespace APP  {

void Output_Streams::init_output_streams()
{
  // This routine initializes the output streams that all SAGE programs need
  // consistently

  // Create output streams (Screen, Information)

  cerrorstream Screen = cerrorstream(std::cout);

  if(!inf_file)
  {
    sage_cerr << priority(fatal)
              << "Cannot open information file.  Exiting..." << endl;
    exit(EXIT_FAILURE);
  }

  inf_file.open((my_program_name + ".inf").c_str());

  cerrorstream Information(inf_file);

  // Format the multistreams (see programmer documentation /usr/sage/docs/SAGE Output.doc)

  //   screen      -> Screen

  my_screen_stream.insert(Screen);

  //   information -> Information (+ Debug)

  my_information_stream.insert(Information);

  //   messages    -> Screen + Information (+ Debug)

  my_message_stream.insert(Screen);
  my_message_stream.insert(Information);

  // error stream: set the prefix for all error messages produced.
  // Each stream added to the error stream must be prefixed with this code.

  string prefix = "%%" + toUpper(my_program_name) + "-%p: ";

  //   errors      -> sage_cerr      (if priority >= error)

  sage_cerr.prefix(prefix);

  my_error_stream.insert(sage_cerr);
  my_error_stream.restrict(r_ge, error);

  //   errors      -> Information    (if priority >= information)

  Information.prefix(prefix);

  my_error_stream.insert(Information);
  my_error_stream.restrict(r_ge, information);

  // if debug is true we construct the debug stream and add it to the others

  if(my_debug)
  {
    //   create Debug stream

    dbg_file.open((my_program_name + ".dbg").c_str());

    cerrorstream Debug(dbg_file);

    //   add Debug to information, messages and errors (priority >= debug)

    my_information_stream.insert(Debug);
    my_message_stream.insert(Debug);
    my_debug_stream.insert(Debug);

    Debug.prefix(prefix);

    my_error_stream.insert(Debug);
    my_error_stream.restrict(r_ge, SAGE::debug);

    //   add Screen to errors (priority == debug)

    Screen.prefix(prefix);

    my_error_stream.insert(Screen);
    my_error_stream.restrict(r_eq, SAGE::debug);
  }

  // Classify the non-error streams as raw output.

  my_screen_stream.set_raw_mode();
  my_message_stream.set_raw_mode();
  my_information_stream.set_raw_mode();
}

} // End namespace APP
} // End namespace SAGE
