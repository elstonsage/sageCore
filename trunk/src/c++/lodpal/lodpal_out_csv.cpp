//****************************************************************************
//* File:      lodpal_out_csv.cpp                                            *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0 Initial implementation                kbj         *
//*                    1.0 Modification.                         yjs Aug. 01 *
//*                                                                          *
//* Notes:     This file implemnets lodpal_test_csvfile class.               *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "lodpal/lodpal_out.h"

using namespace std;

#undef SIMULATION_MODE

namespace SAGE   {
namespace LODPAL {

//---------------------------------------------------------------------------

void
lodpal_test_csvfile::print_header(const ARP_base_analysis& test)
{
  ostream& out = output_stream();

  time_t timer;
  time(&timer);

  out << ctime(&timer) << endl;

  if(!printed)
  {
    out << "TRAIT,MARKER,LOD,F.Sibs,H.Sibs,Other,Pairs,Beta1,Beta2,";

    for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
      out << "Delta" << ci+1 << "1,Delta" << ci+1 << "2,";

    out << "NFE,LFL" << endl;
  }
  printed = true;
}

void
lodpal_test_csvfile::print_results(const ARP_base_analysis& test)
{
  ostream& out = output_stream();

  lodpal_parameters::marker_const_iterator pi = test.parameters().marker_begin();
  for( ; pi != test.parameters().marker_end(); ++pi)
  {
    const trait_parameter& trait = test.parameters().trait_parameters(0);
    out << trait.name(test.relative_pairs())             << ","
        << test.relative_pairs().marker_name(pi->marker) << ","
        << test.get_lodpal_result().lod_score()          << ","
        << test.pairs_info().fsib_pair_count()           << ","
        << test.pairs_info().hsib_pair_count()           << ","
        << test.pairs_info().other_pair_count()          << ","
        << test.pairs_info().pair_count()                << ","
        << pi->beta1.value()                             << ","
        << pi->beta2.value()                             << ",";

    lodpal_parameters::covariate_const_iterator ci = test.parameters().covariate_begin();
    for( ; ci != test.parameters().covariate_end(); ++ci )
      out << ci->delta1.value()                          << ","
          << ci->delta2.value()                          << ",";

    out << test.get_lodpal_result().function_evaluations() << ","
        << test.get_lodpal_result().last_error()           << endl;
  }
}

void
lodpal_test_csvfile::print_footer(const ARP_base_analysis& test)
{}

//---------------------------------------------------------------------------

} // end of namespace LODPAL
} // end of namespace SAGE
