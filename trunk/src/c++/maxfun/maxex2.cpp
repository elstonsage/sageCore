/* MAXEX-EASY:  SAMPLE DRIVER PROGRAM FOR MAXFUN
 *
 * Last Modified:
 *  7-JUN-1999 - Kevin Jacobs - Updated for new Fortran runtime
 * 13-APR-1998 - Kevin Jacobs - C++ conversion
 * 24-JAN-1996 - Kevin Jacobs - updated
 *  7-DEC-1995 - Kevin Jacobs - created
*/

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include "maxfun/maxfun.h"

using namespace std;
using namespace SAGE;

/* nice aliases for the parameters to be estimates */
enum index_t { pA=0, pB=1, pO=2, params=3 };

const int N_AB = 1269;    /* Constants for function   */
const int N_A  = 9123;    
const int N_B  = 2987;
const int N_O  = 7725;

class MyFunction : public MaxFunction
{
  virtual double evaluate(parameter_vector& theta);
  virtual int    update_bounds(parameter_vector& theta);
};

int main(int argc, char* argv[])
{
    int method = 0;

    if(argc > 1)
      method = atoi(argv[1]);

    if(method <= 0 || method > 6)
      method = 1;

    MyFunction myfun;  

    // Create MAXFUN instance
    Maxfun maxfun(myfun);

    maxfun.epsc1()   = 1e-7;
    maxfun.ixvc()    = 2;          // Compute var-cov matrix when done
    maxfun.method()  = method;     // direct search with 2^n trials   
    maxfun.nt()      = params;     // 3 parameters total              
    maxfun.thin(pA)  = exp(.8);    // Initial estimate of pA          
    maxfun.thl(pA)   = exp(0.);    //   lower bound                   
    maxfun.thu(pA)   = exp(1.);    //   upper bound                   
    maxfun.istin(pA) = 1;          //   independent, but used in depar
    maxfun.thin(pB)  = exp(.1);    // Initial estimate of pB          
    maxfun.thl(pB)   = exp(0.);    //   lower bound                   
    maxfun.thu(pB)   = exp(1.);    //   upper bound                   
    maxfun.istin(pB) = 1;          //   independent, but used in depar
    maxfun.thin(pO)  = exp(.1);    // Initial estimate of pB          
    maxfun.thl(pO)   = exp(0.);    //   lower bound                   
    maxfun.thu(pO)   = exp(1.);    //   upper bound                   
    maxfun.istin(pO) = 3;          //   depends on pA, pB             
    maxfun.maxit()   = 1000;       // max. interations                
    
    //
    // Call maxfun with the function to be evaluated, the function to 
    // check the bounds of the parameters, the array of initial values,
    // an int for the number of times the function is evaluated,
    // and an error return code.


    cout << "Running MAXFUN with method " << method << endl;

    maxfun.run();

    printf("\nFinal values:\n");
    printf("  pA=%f\tV(pA)=%12.10f\tSE(pA)=%f\n",  log(maxfun.param(pA)), 
                                 maxfun.av(pA,pA), maxfun.stde(pA));
    printf("  pB=%f\tV(pB)=%12.10f\tSE(pB)=%f\n",  log(maxfun.param(pB)), 
                                 maxfun.av(pB,pB), maxfun.stde(pB));
    printf("  pO=%f\tV(pO)=%12.10f\tSE(pO)=%f\n",  log(maxfun.param(pO)), 
                                 maxfun.av(pO,pO), maxfun.stde(pO));

    cout << endl << "MAXFUN completed with " << maxfun.evaluations()
                 << " function evaluations." << endl;

    /* Prepare to run maxfun again to with final estimates as initial   */
    /* estimates, all parameters independent and 0 interations.  Maxfun */
    /* will only recompute the var-cov matrix.  The results had better  */
    /* be the same as above.                                            */
                                      
    maxfun.thin(pA) = maxfun.param(pA);     /* estimate of pA */
    maxfun.thin(pB) = maxfun.param(pB);     /* estimate of pB */
    maxfun.thin(pO) = maxfun.param(pO);     /* estimate of pO */

    maxfun.maxit()  = 0;             /* Set max iterations to 0 so that */
                                     /* only the var-cov matrix is      */
                                     /* computed since initial estimates*/
                                     /* are the values maximized using  */
                                     /* maxfun above                    */

    /* Call maxfun again */
    maxfun.run();

    printf("\nFinal values (again):\n");
    printf("  pA=%f\tV(pA)=%12.10f\tSE(pA)=%f\n",  log(maxfun.param(pA)), 
                                 maxfun.av(pA,pA), maxfun.stde(pA));
    printf("  pB=%f\tV(pB)=%12.10f\tSE(pB)=%f\n",  log(maxfun.param(pB)), 
                                 maxfun.av(pB,pB), maxfun.stde(pB));
    printf("  pO=%f\tV(pO)=%12.10f\tSE(pO)=%f\n",  log(maxfun.param(pO)), 
                                 maxfun.av(pO,pO), maxfun.stde(pO));

    cout << endl << endl;
    return 0;    
}

double MyFunction::evaluate(parameter_vector& tr)
/*
 * The function to be evaluated :
 *  Parameters: tr  - trial parameters
 */
{
    double ftr;

   double tra = log(tr[pA]);
   double trb = log(tr[pB]);
   double tro = log(tr[pO]);

    /* calculate the function */
    ftr = N_AB * log(                2*tra*trb) 
        + N_A  * log(tra*tra + 2*tra*tro)
        + N_B  * log(trb*trb + 2*trb*tro) 
        + N_O  * log(tro*tro);

    ++nfe;

    return ftr;    /* and return   */
} 

int MyFunction::update_bounds(parameter_vector& tr)
{
   tr[pO] = exp(1.0 - log(tr[pA]) - log(tr[pB]));      

   double tra = log(tr[pA]);
   double trb = log(tr[pB]);
   double tro = log(tr[pO]);

   /* Make sure all parameters are within bounds */
   if(tra <= 0 || tra >= 1 || trb <= 0 || trb >= 1 ||
      tro <= 0 || tro >= 1)
     return 1;                /* if not, return an error */

   return 0;                  /* else return no error    */
}
