/* MAXEX-EASY:  SAMPLE DRIVER PROGRAM FOR MAXFUN
 *
 * Last Modified:
 *  8-APR-2004 - Steve Gross  - Updated for new Maxfun API
 *  7-JUN-1999 - Kevin Jacobs - Updated for new Fortran runtime
 * 13-APR-1998 - Kevin Jacobs - C++ conversion
 * 24-JAN-1996 - Kevin Jacobs - updated
 *  7-DEC-1995 - Kevin Jacobs - created
*/

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include "maxfun/MaxfunAPI.h"

using namespace std;
using namespace SAGE;

/* nice aliases for the parameters to be estimates */
enum param_t { pA=0, pB=1, pO=2, params=3 };

const int N_AB = 1269;    /* Constants for function   */
const int N_A  = 9123;    
const int N_B  = 2987;
const int N_O  = 7725;

class MyFunction : public MaxFunction
{
  public:
    virtual double evaluate      (parameter_vector& theta);
    virtual int    update_bounds (parameter_vector& theta);

    MAXFUN::MaxfunInfo minfo;
};

int main(int argc, char* argv[])
{
  // 2. Set up MaxfunInfo:

	MyFunction myfun;  

	myfun.minfo.AddParameter("pA", 0.8, MAXFUN::type_independent_functional, 0.0, 1.0);
	myfun.minfo.AddParameter("pB", 0.1, MAXFUN::type_independent_functional, 0.0, 1.0);
	myfun.minfo.AddParameter("pO", 0.1, MAXFUN::type_dependent,              0.0, 1.0);

  // 3. Maximize

	Maximize(MAXFUN::max_default, myfun.minfo, myfun);
//, MAXDebug(opt1, opt2, opt3));

  // 5. Return:

	return 0;    
}

double MyFunction::evaluate(parameter_vector& tr)
{
    double pA = minfo("pA"),
           pB = minfo("pB"),
           pO = minfo("pO");

    double ftr;
    
    /* calculate the function */
    ftr = N_AB * log(        2*pA*pB) 
        + N_A  * log(pA*pA + 2*pA*pO)
        + N_B  * log(pB*pB + 2*pB*pO) 
        + N_O  * log(pO*pO);

    ++nfe;

    return ftr;    /* and return   */
} 

int MyFunction::update_bounds(parameter_vector& tr)
{
  double pA = minfo("pA"),
         pB = minfo("pB");

  minfo("pO") = 1.0 - pA - pB;

  if(!minfo.GetParameter("pO").InBounds())
    return 1;
  else
    return 0; 
}
