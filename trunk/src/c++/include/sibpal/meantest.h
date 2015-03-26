#ifndef SIBPAL_MEAN_REGRESSION_H
#define SIBPAL_MEAN_REGRESSION_H

//=============================================================================
// File:     meantest.h
//
// Author:   Kevin Jacobs
//
// History:  Version 0.0 Initial implementation
//
// Notes:
//
// Copyright (c) 2001 R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "sibpal/parser.h"

namespace SAGE   {
namespace SIBPAL {

class SibMeanTest
{
  public:

    SibMeanTest(relative_pairs& p, cerrorstream& err = sage_cerr);

    void  regress();

    void  set_parameters(const meantest_parameters& p);

          meantest_parameters& parameters();
    const meantest_parameters& parameters() const;

    void  set_use_pairs(const pair<bool, bool>& up);
    pair<bool, bool> get_use_pairs() const;

  protected:

    void regress_marker(GLS3& gls, size_t m);

    meantest_parameters my_parameters;

    relative_pairs&     pairs;
    cerrorstream&       errors;
};

#include "sibpal/meantest.ipp"

} //end of namespace SIBPAL
} //end of namespace SAGE

#endif
