//==========================================================================
//  File:       l2_procedure.cpp
//
//  Author:     Yeunjoo Song
//
//  History:    Version 1.0  Initial implementation.             yjs Sep. 01
//                      2.0  Updated to new libraries.           yjs Jul. 03
//
//  Notes:      This class is the main driver performing analysis.
//
//  Copyright (c) 2003 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "reltest/l2_procedure.h"

namespace SAGE
{

namespace RELTEST
{

double
L2_error_procedure::evaluate(vector<double>& theta)
{
  ++nfe;

  double mean     = theta[0];
  double variance = theta[1];

#if 0
  cout << endl << "mean = " << mean << ", variance = " << variance << endl;
#endif

  double pi_variance = 3.1415926535897931 * variance;

  double f_value = 0.;

#if 0
  if( my_Yj )
    cout << "Yj\t\tnormal_density\t\tf_value" << endl;
  else
    cout << "Yjp\t\tnormal_density\t\tf_value" << endl;
#endif

  for( size_t p = 0; p < my_ptt_pairs.size(); ++p )
  {
    double Yj = my_ptt_pairs[p].get_Yj();

    if( !my_Yj )
      Yj = my_ptt_pairs[p].get_Yjp();

    double N_density = exp( -( pow(Yj-mean, 2) / 2.0*variance ) ) / sqrt( 2.0 * pi_variance );

    if( fabs(N_density) < 1.0e-11 )
      N_density = 0.0;

    f_value += N_density;

#if 0
  cout << Yj  << "\t";
  cout << N_density << "\t\t";
  cout << f_value << endl;
#endif
  }
  
  double final_value = 2. * f_value / my_ptt_pairs.size() - 1.0 / (2.0 * sqrt(pi_variance));

  if( fabs(final_value) < 1.0e-11 )
    final_value = 0.0;

  nfe += 1;

#if 0
  cout << "final f_value = " << final_value << endl;
#endif

  return final_value;
}

int
L2_error_procedure::update_bounds(vector<double>& theta)
{
  double variance = theta[1];

  if( variance < 0.5 )//std::numeric_limits<double>::epsilon() )
    return 1;

  return 0;
}

} // end of namespace RELTEST

} // end of naemspace SAGE
