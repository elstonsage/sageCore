#ifndef SIBPAL_MEANTEST_OUTPUT_H
#define SIBPAL_MEANTEST_OUTPUT_H

//=============================================================================
// File:    meantest_out.h
//
// Author:  Kevin Jacobs
//
// History: Version 0.0 Initial implementation
//
// Notes:   This file contains definition for following data structures.
//
//            class mean_test_viewer
//            class mean_test_textfile : public mean_test_viewer
//            class mean_test_csvfile  : public mean_test_viewer
//
// Copyright (c) 2001   R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "sibpal/meantest.h"

namespace SAGE   {
namespace SIBPAL {

class meantest_viewer
{
  public:

    typedef vector<meantest_parameters>     meantest_param_vector;

    meantest_viewer(ostream& output, const relative_pairs& pairs)
      : my_pairs(pairs), my_output(output) 
    { 
      set_wide_output(false);
      set_print_only_pi(false);
    }

    meantest_viewer(ostream& output, const relative_pairs& pairs, bool wide_out, bool only_pi)
      : my_pairs(pairs), my_output(output) 
    { 
      set_wide_output(wide_out);
      set_print_only_pi(only_pi);
    }

    virtual ~meantest_viewer() {}

    virtual void print_results(const meantest_parameters&) = 0;
    virtual void print_results(const meantest_parameters& , size_t i) = 0;
    virtual void print_results(const meantest_param_vector &, const string&) = 0;
    virtual void print_results(const meantest_param_vector &, size_t i, const string&) = 0;

    ostream&               output_stream()   { return my_output; }
    const relative_pairs&  sib_pairs() const { return my_pairs;   }

    void     set_print_only_pi(bool p)  { my_print_only_pi = p; }
    void     set_wide_output(bool w)    { my_wide_output   = w; }

    bool     print_only_pi()  const     { return my_print_only_pi; }
    bool     wide_output()    const     { return my_wide_output; }
    
  protected:

    const relative_pairs&     my_pairs;

    ostream&                  my_output;

    bool                      my_wide_output;
    bool                      my_print_only_pi;
};

class meantest_textfile : public meantest_viewer
{
  public:

    meantest_textfile(ostream& output, const relative_pairs& pairs)
      : meantest_viewer(output, pairs) 
    { }

    meantest_textfile(ostream& output, const relative_pairs& pairs, bool wide_out, bool only_pi)
      : meantest_viewer(output, pairs, wide_out, only_pi) 
    { }

    virtual ~meantest_textfile() {}

    virtual void print_results(const meantest_parameters&);
    virtual void print_results(const meantest_parameters&, size_t i);
    virtual void print_results(const meantest_param_vector &, const string&);
    virtual void print_results(const meantest_param_vector &, size_t i, const string&);

  protected:
    
    void print_result_line(double est, double se, double t, double p,
                           int w, int pr, int pw, bool sci,
                           bool prt_exp = false, double exp_val = 0.5, bool less = false);
};

class meantest_csvfile : public meantest_viewer
{
  public:

    meantest_csvfile(ostream& output, const relative_pairs& pairs)
      : meantest_viewer(output, pairs) 
    { printed = false; }

    meantest_csvfile(ostream& output, const relative_pairs& pairs, bool wide_out, bool only_pi)
      : meantest_viewer(output, pairs, wide_out, only_pi) 
    { printed = false; }

    virtual void print_results(const meantest_parameters&);
    virtual void print_results(const meantest_parameters&, size_t i);
    virtual void print_results(const meantest_param_vector &, const string&);
    virtual void print_results(const meantest_param_vector &, size_t i, const string&);

  private:

    bool printed;
};

} // end of namespace SIBPAL
} // end of namespace SAGE

#endif
