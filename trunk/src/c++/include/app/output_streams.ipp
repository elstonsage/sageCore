//============================================================================
// File:      output_streams.ipp
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

namespace SAGE {
namespace APP  {

// ================
// Inline Functions
// ================

inline
Output_Streams::Output_Streams(const std::string& p, bool debug_on)
  : my_program_name(p), my_debug(debug_on)
{
  init_output_streams();
}

inline
Output_Streams::~Output_Streams()
{ }

//lint 1536 returning access to the streams is ok

inline std::ostream& Output_Streams::screen()
{
  //lint --e{1536} returning access to the streams is ok

  return my_screen_stream;
}

inline std::ostream& Output_Streams::info()
{
  //lint --e{1536} returning access to the streams is ok

  return my_information_stream;
}

inline std::ostream& Output_Streams::messages()
{
  //lint --e{1536} returning access to the streams is ok

  return my_message_stream;
}
inline cerrorstream& Output_Streams::errors()
{
  //lint --e{1536} returning access to the streams is ok

  return my_error_stream;
}

inline std::ostream& Output_Streams::debug()
{
  //lint --e{1536} returning access to the streams is ok

  return my_debug_stream;
}

inline bool Output_Streams::get_debug_status() const
{
  return my_debug;
}

} // End namespace APP
} // End namespace SAGE

