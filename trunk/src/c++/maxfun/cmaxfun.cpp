// #define KCC_DEBUG 1
// # if KCC_DEBUG
//    #define INSPECT(expr, type) printf( "\n--> %s:["#type"] at %s(%d)", #expr, expr, __FILE__, __LINE__)
//    #define MARK_PROGRESS() printf( "\n>>>>> . <<<<< at %s(%d)", __FILE__, __LINE__)
// #endif


/* MAXFUN VERSION 6.1 (3-OCT-1995)
 *
 * WRITTEN BY ALEXA J. M. SORANT AND ROBERT C. ELSTON,
 * CONVERTED TO C WITH DYNAMIC ARRAY SUPPORY BY KEVIN JACOBS
 * THEN CONVERTED TO C++ BY KEVIN JACOBS
 * WITH ALGORITHMS FOR DIRECT SEARCH AND NEWTON-RAPHSON METHODS ADAPTED
 * FROM MAXLIK (BY E. B. KAPLAN AND R. C. ELSTON) AND ALGORITHM FOR
 * VARIABLE METRIC METHOD ADAPTED FROM GEMINI (BY J. M. LALOUEL)
 *
 * MAXIMIZES A FUNCTION SPECIFIED BY CALLER (IN SUBROUTINE FUN) WITH
 * RESPECT TO PARAMETERS (ARRAY THETA) SUBJECT TO RESTRICTIONS SPECIFIED
 * BY CALLER (IN SUBROUTINE DEPAR AND BOUND ARRAYS THL AND THU)
 *
 *------------------------------------------------------------------------
 *     >>>>> REFER TO FORTRAN SOURCE CODE FOR MORE INFORMATION <<<<<
 *------------------------------------------------------------------------
 *          Copyright(c) R.C. Elston 1995.   All Rights Reserved.
 *------------------------------------------------------------------------
 *
 * RETURNS FLAG LFL:
 *   0:  NO ERRORS FOUND, BUT ZERO ITERATIONS REQUESTED
 *   1:  CONVERGED BY CRITERION 1 (IN BCNVCH OR VCNVCH)
 *   2:  CONVERGED BY CRITERION 2:  NORMALIZED GRADIENT
 *       <= SPECIFIED TOLERANCE EPSC2 (IN NRSTEP OR DIRECT)
 *   3:  CONVERGED BY CRITERION 3:  CHANGE IN FUNCTION VALUE LESS THAN
 *       SPECIFIED TOLERANCE EPSC3 (IN BCNVCH OR VCNVCH)
 *   4:  REACHED MAXIMUM # ITERATIONS
 *   5:  ACCUMULATION OF ROUNDING ERROR OR BOUNDARY PROBLEM PREVENTS
 *       FURTHER PROGRESS (VARIABLE METRIC METHOD) (NO STEP COMPLETED
 *       IN LSRCH)
 *   6:  SEARCH DIRECTION NOT UPWARDS (VARIABLE METRIC METHOD)
 *       (IN DIRECT)
 *   7:  ALL SIGNIFICANT DIGITS LOST IN OPTIMAL CONDITIONING DURING
 *       UPDATE OF B MATRIX (VARIABLE METRIC METHOD) (IN BUPDT)
 *   8:  GRADIENT COULD NOT BE COMPUTED (IN DERIV1)
 *   9:  VARIANCE-COVARIANCE MATRIX COULD NOT BE COMPUTED (IN VCMX)
 *  10:  ALL INDEPENDENT PARAMETERS CONVERGED TO BOUNDS (IN FIXBND)
 *  11:  TOO MANY PARAMETERS TO USE DERIVATIVES (IN PREPD)
 *  12:  UNACCEPTABLE INITIAL ESTIMATES (IN MAXFUN)
 *  13:  INVALID CONTROL INPUT (IN MAXFUN)
 */
/*
 * Last Modified:
 *  2000-01-17 - Kevin Jacobs - Massive conversion to C++
 *  7-JUN-1999 - Kevin Jacobs - Updated for new Fortran runtime
 *  7-DEC-1995 - Kevin Jacobs - created
*/

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <list>
#include <iostream>
#include "maxfun/maxfun.h"
#include "numerics/clapack.h" // for matrix operations
#include "numerics/fmatrix.h"

namespace SAGE {

using std::vector;
using std::min;
using std::max;

#ifndef _MSC_VER
using std::fabs;
using std::sqrt;
#endif

typedef SAGE::FortranMatrix<double> matrix;

int do_SVD_inverse(vector<double>& a, int a_dim, int m, vector<double>& b, int b_dim, int& lex)
{
  #undef  a_ref
  #define a_ref(x, y) a[y*a_dim + x]

  #undef  b_ref
  #define b_ref(x, y) b[y*b_dim + x]

  SAGE::SVD  svd;

  matrix svdA;
  matrix svdA_i;

  svdA.resize((size_t)m, (size_t)m, 0.0);

  for( int i = 0; i < m; i++ )
  {
    for( int j = 0; j < m; j++ )
    {
      svdA(i, j) = a_ref(i, j);
    }
  }

  int svd_return_code = svd.compute(svdA);

  if( svd_return_code == 0 ) // Not invertable.
  {
    return svd_return_code;
  }

  if( svd.positive_diagonal.size() < svdA.cols() )
  {
    svd.general_inverse(svdA_i);
  }
  else
  {
    svd.inverse(svdA_i);
  }

#if 0
  cout << endl;
  cout << "lex = " << lex << endl;
  cout << "svdA" << endl;
  cout << "--------------------------------" << endl;
  for( size_t i = 0; i < svdA.rows(); i++ )
  {
    for( size_t j = 0; j < svdA.cols(); j++ )
    {
      cout << setw(12) << svdA(i,j) << " ";
    }
    cout << endl;
  }

  cout << "svdA inverse" << endl;
  cout << "--------------------------------" << endl;
  for( size_t i = 0; i < svdA_i.rows(); i++ )
  {
    for( size_t j = 0; j < svdA_i.cols(); j++ )
    {
      cout << setw(12) << svdA_i(i,j) << " ";
    }
    cout << endl;
  }
#endif

  for(int r = 0; r < m; r++)
  {
    for(int s = 0; s < m; s++)
    {
      b_ref(r, s) = svdA_i(r, s);
    }
  }

  #undef b_ref
  #undef a_ref

  return svd_return_code;
}

/* ************************************************************************** **
** *************                Maxfun::pseudoInverse            ************ **
** ************************************************************************** **
   ...
----------------------------------------------------------------------------- */
int pseudoInverse(double* A, int m, int n)
{
   /* --------------------------------------------------------------------------
      Syntactic sugar
   -------------------------------------------------------------------------- */
   #undef  A_ref
   #define A_ref(x, y) A[y*n + x]

   #undef  a_ref
   #define a_ref(x, y) a[y*n + x]

   #undef  b_ref
   #define b_ref(x, y) b[y*n + x]

   /* --------------------------------------------------------------------------
      Local data
   -------------------------------------------------------------------------- */
   int M          = m;
   int N          = n;
   int info       = 0;
   int LDWORK     = -1;
   int NRHS       = M;
   int LDA        = M;
   int LDB        = M;
   int IRANK      = -1;
   double* WORK   = new double[1];
   double* SINGV  = new double[min(M,N)];
   double RCOND   = -1.0f;

   /* --------------------------------------------------------------------------
      Init
   -------------------------------------------------------------------------- */
   // MxM temp matrix
   double* a = new double[M*N];
   for(int i = 0; i < M; i++) {
      for(int j = 0; j < N; j++) {
         a_ref(i, j) = A_ref(i, j);
      }
   }

   // MxM identity matrix
   double* b = new double[M*M];
   for(int i = 0; i < m; i++) {
      for(int j = 0; j < m; j++) {
         if(i == j) {
            b_ref(i, j) = 1.0;
         }
         else {
            b_ref(i, j) = 0.0;
         }
      }
   }

   /* --------------------------------------------------------------------------
      Get the inverse
   -------------------------------------------------------------------------- */
   // Query and allocate the optimal workspace

   int res = F77_CALL(dgelss)(&M, &N, &NRHS, &a[0], &LDA, &b[0], &LDB, SINGV, &RCOND, &IRANK, WORK, &LDWORK, &info);

   if(info == 0) {
      // Compute pseudo inverse
      LDWORK = (int)WORK[0];
      delete [] WORK;
      WORK = new double[LDWORK];
      res = F77_CALL(dgelss)( &M, &N, &NRHS, &a[0], &LDA, &b[0], &LDB, SINGV, &RCOND, &IRANK, WORK, &LDWORK, &info);

      if(info == 0) {
         for(int i = 0; i < m; i++) {
            for(int j = 0; j < m; j++) {
               A_ref(i, j) = b_ref(i, j);
            }
         }
      }
   }

   delete [] b;
   delete [] WORK;
   delete [] SINGV;

   return(info);
}

int do_LU_inverse(vector<double>& a, int a_dim, int m, vector<double>& b, int b_dim, int& lex)
{
   // --------------------------------------------------------------------------
   //  We're using a vector to implement a matrix, so the following
   //  "syntactic sugar" helps with the indexing
   // --------------------------------------------------------------------------
   #undef  a_ref
   #define a_ref(x, y) a[y*a_dim + x]

   #undef  b_ref
   #define b_ref(x, y) b[y*b_dim + x]

   #undef  tempA_ref
   #define tempA_ref(x, y) tempA[(y)*m + x]

   // --------------------------------------------------------------------------
   //    For NORMAL RETURN matrix inversion
   // --------------------------------------------------------------------------
   int LWORK         = 10*m;
   int* permutations = new int[2*m];
   double* WORK      = new double[m*m];
   double* tempA     = new double[m*m];
   int INFO          = 0;
   /* --------------------------------------------------------------------------
      Setup for matrix inversion
   -------------------------------------------------------------------------- */
   for(int i = 0; i < m; i++) {
      for(int j = 0; j < m; j++) {
         b_ref(i, j) = a_ref(i, j);
         tempA_ref(i, j) = a_ref(i, j);
      }
   }

   /* --------------------------------------------------------------------------
      Transpose a -> tempA
   -------------------------------------------------------------------------- */
   for(int r = 0; r < m; r++) {
      for(int s = 0; s < m; s++) {
         tempA[r+m*s] = a_ref(s, r);
      }
   }

#if 1
    cout << endl;
    cout << "tempA                        (1)" << endl;
    cout << "--------------------------------" << endl;
    for (size_t i = 0; i < (size_t)m ; i++)
    {
       for (size_t j = 0; j < (size_t)m ; j++)
       {
          cout << setw(12) << tempA_ref(i,j) << " ";
       }
       cout << endl;
    }
#endif

   /* --------------------------------------------------------------------------
      Invert tempA ... use LAPACK for the heavy lifting
   -------------------------------------------------------------------------- */
   // get LU decomposition of a general matrix
   F77_CALL(dgetrf)(&m, &m, tempA , &m, permutations , &INFO );
   if (INFO != 0) {
      lex = 1;
      return 0;
   }

   // compute the matrix inverse using the LU factorization computed by dgetrf
   F77_CALL(dgetri)(&m, tempA , &m, permutations , WORK, &LWORK, &INFO );

   if (INFO != 0)
   {
#if 1
      // Maybe singular?  Try pseudo-inverse
      cout << endl;
      cout << "tempA                        (2)" << endl;
      cout << "--------------------------------" << endl;
      for (size_t i = 0; i < (size_t)m ; i++)
      {
         for (size_t j = 0; j < (size_t)m ; j++)
         {
            cout << setw(12) << tempA_ref(i,j) << " ";
         }
         cout << endl;
      }
#endif

      for(int i = 0; i < m; i++) {
         for(int j = 0; j < m; j++) {
            tempA_ref(i, j) = a_ref(i, j);
         }
      }

      if( pseudoInverse(tempA, m, m) != 0 ) {
          lex = 1;
          return 0;
      }
   }
#if 1
    cout << endl;
    cout << "tempA                        (3)" << endl;
    cout << "--------------------------------" << endl;
    for (size_t i = 0; i < (size_t)m ; i++)
    {
       for (size_t j = 0; j < (size_t)m ; j++)
       {
          cout << setw(12) << tempA_ref(i,j) << " ";
       }
       cout << endl;
    }
#endif
   /* --------------------------------------------------------------------------
      Transpose tempA -> b
   -------------------------------------------------------------------------- */
   for(int r = 0; r < m; r++)
   {
      for(int s = 0; s < m; s++)
      {
         b_ref(r, s) = tempA[r*m+s];
      }
   }

   /* --------------------------------------------------------------------------
      Cleanup
   -------------------------------------------------------------------------- */
   delete [] permutations;
   delete [] WORK;
   delete [] tempA;

  #undef a_ref
  #undef b_ref
  #undef tempA_ref

  return 1;
}

/* Initialized data */

char Maxfun_Data::maxfst_[][80] =
          { "IND-FN, MAY VARY",
            "IND, MAY VARY",
            "DEPENDENT",
            "FIXED EXTERNALLY",
            "IND-FN, FIXED BY MAXFUN AT BOUND",
            "IND, FIXED BY MAXFUN AT BOUND",
            "IND-FN, FIXED BY MAXFUN NEAR BOUND",
            "IND, FIXED BY MAXFUN NEAR BOUND",
            "IND-FN, FIXED BY MAXFUN NOT NEAR BND",
            "IND, FIXED BY MAXFUN NOT NEAR BOUND"    };

/* Table of constant values */

bool Maxfun::iteration_ended = false;

// by JA for djb

double Maxfun::obtain_score()
{
  double return_val = numeric_limits<double>::quiet_NaN();

  bool goodfirstderiv = score_deriv1();

  bool goodsecderiv = false;

  if (goodfirstderiv) goodsecderiv = score_deriv2();

  if ((goodsecderiv) && (goodfirstderiv))
  {
    return_val = calc_score();
  }

  return return_val;
}

int Maxfun::maxfun_()
{
  my_data.maxf2_.nb = 0;
  my_data.maxf2_.nv = 0;

  /* Local variables */
  int iflb;
  int isti;
  int ihess;
  int iupdt;
  int ih;
  double si;
  int ifirst;
  int ifl;
  int lex;
  double thi;
  int ireturn;

  vector<double> dth;
  vector<double> gpr;


/* --INITIALIZE EXIT FLAG */

  dth.resize(my_data.param_NP);
  gpr.resize(my_data.param_NPV);

  /* Function Body */
  lfl = 0;

/* --PRINT INITIAL MESSAGES AND CHECK INPUT */

  if (my_data.maxf1_.method < 1 || my_data.maxf1_.method > 6)
  {
    lfl = 13;
    my_data.maxf2_.igfl = 2;
    return 0;
  }

  if (my_data.maxf1_.nt < 0 || my_data.maxf1_.nt > my_data.param_NP)
  {
    lfl = 13;
    my_data.maxf2_.igfl = 2;
    return 0;
  }

  my_data.maxf2_.ni = 0;
  my_data.maxf2_.ne = 0;
  my_data.maxf2_.nd = 0;
  for (int i = 0; i < my_data.maxf1_.nt; ++i)
  {
    isti = my_data.maxf1_.istin[i];
    switch ((int)isti)
    {
      case 3:  ++my_data.maxf2_.nd;
      case 1:
      case 2:  ++my_data.maxf2_.ne;
              break;
      case 4: break;
      default: lfl = 13;
    }
  }

  if (my_data.maxf2_.ne <= 0)
  {
    lfl = 13;
    my_data.maxf2_.igfl = 2;
    return 0;
  }

  if (my_data.maxf2_.nd >= my_data.maxf2_.ne)
  {
    lfl = 13;
    my_data.maxf2_.igfl = 2;
    return 0;
  }

/* --SUBSTITUTE DEFAULTS IF NECESSARY */

  if (my_data.maxf1_.ixvc < 0 || my_data.maxf1_.ixvc > 2)
    my_data.maxf1_.ixvc = 0;

  if (my_data.maxf1_.ihit < 0 || my_data.maxf1_.ihit > 1)
    my_data.maxf1_.ihit = 0;

  if (my_data.maxf1_.epsc1 <= 0. || my_data.maxf1_.epsc1 >= .5)
    my_data.maxf1_.epsc1 = .001;

  if (my_data.maxf1_.epsc2 <= 0.)
    my_data.maxf1_.epsc2 = 1e-15;

  if (my_data.maxf1_.epsc3 < 0.)
    my_data.maxf1_.epsc3 = 0.;

  if (my_data.maxf1_.epsd <= 1e-9 || my_data.maxf1_.epsd >= .5)
  {
    if (my_data.maxf1_.epsc1 > 1e-9)
      my_data.maxf1_.epsd = my_data.maxf1_.epsc1;
    else
      my_data.maxf1_.epsd = .001;
  }

  if (my_data.maxf1_.yota <= 0. || my_data.maxf1_.yota >= 1.)
    my_data.maxf1_.yota = sqrt(2.8421709430404007e-14) * 10.;

  if (my_data.maxf1_.epst <= 0.)
    my_data.maxf1_.epst = sqrt(my_data.maxf1_.yota) * 10.;

  for (int i = 0; i < my_data.maxf1_.nt; ++i)
  {
    si = my_data.maxf1_.stpin[i];
    if (si <= 0. || si >= 1.)
      my_data.maxf1_.stpin[i] = .1;
  }

  my_data.clear_results();

/* --PRINT OUT CONTROL VALUES */

/* --GENERAL INITIALIZATION */

  my_data.maxf2_.igage = -1;
  my_data.maxf2_.ivage = -1;
  my_data.maxf2_.igfl = 0;
  my_data.maxf2_.ivfl = 5;
  ih = my_data.maxf1_.ihit;
  my_data.maxf2_.nsurf2 = 0;
  my_data.maxf2_.impbnd = 0;
  my_data.maxf2_.it = 0;
  nfe = 0;

  my_data.maxf2_.ni = my_data.maxf2_.ne - my_data.maxf2_.nd;
  my_data.maxf1_.nt = my_data.maxf1_.nt;


  for (int i = 0; i < my_data.maxf1_.nt; ++i)
  {
    isti = my_data.maxf1_.istin[i];
    my_data.maxf2_.ist[i] = isti;
    my_data.maxf2_.stp[i] = my_data.maxf1_.stpin[i];
    thi = my_data.maxf1_.thin[i];

    if (isti != 3)
    {
      if (thi < my_data.maxf1_.thl[i] || thi > my_data.maxf1_.thu[i])
        lfl = 12;
    }

    theta[i] = thi;
  }

  if (lfl == 0)
  {
    evaluate(theta, f, nfe, lex);

    if (lex > 0)
      lfl = 12;
  }
  else
  {   /* --INITIAL ESTIMATES INVALID */

    my_data.maxf2_.igfl = 2;
    return 0;
  }

/* --CHECK FOR POSSIBLE ZERO-ITERATION RUN */

  if (my_data.maxf1_.maxit == 0)
  {
    prepd_(lex);

    if(lex > 0)
    {
      my_data.maxf2_.igfl = 2;
      goto L7080;
    }

    if (my_data.maxf1_.ixvc <= 0)
      goto L7050;

    goto L7000;
  }

/* --SELECT PATH FOR CHOSEN METHOD */

  switch ((int)my_data.maxf1_.method)
  {
    case 1:  goto L1000;
    case 2:  goto L2000;
    case 3:  goto L3000;
    case 4:  goto L4000;
    case 5:  goto L5000;
    case 6:  goto L6000;
  }

/* --METHOD 1:  DIRECT SEARCH METHOD, INCLUDING BASIC SEARCH AND */
/* --2**(NI)-TRIAL SEARCH */

/* --PRINT METHOD IDENTIFICATION AND INITIAL VALUES */

L1000:
/* --START BASIC ITERATION PROCESS */
L1010:
  bsrch_(iflb);

/* --IF REACHED MAXIMUM # ITERATIONS, EXIT */

  if(iflb >= 3)
  {
    lfl = 4;
    goto L7030;
  }

/* --SKIP OTHER TYPES OF SEARCH IF ONLY ONE PARAMETER TO VARY */

  if (my_data.maxf2_.ni <= 1)
    goto L1050;
/* --SEARCH 2(NI) "NEIGHBORING" PARAMETER SETS */

  nsrch_(dth, ifl);

  if (ifl > 0)
    goto L1030;

/* --SIGNIFICANT IMPROVEMENT; IF < MAXIMUM # ITERATIONS, GO BACK TO BASIC */
/* --ITERATION PROCESS */

  if (my_data.maxf2_.it < my_data.maxf1_.maxit)
    goto L1010;

  lfl = 4;
  goto L7030;

/* --CHECK WHETHER HAVE DONE 2**(NI)-TRIAL SEARCH MORE THAN ONCE ALREADY
*/

L1030:
  if (my_data.maxf2_.nsurf2 <= 1)
    goto L1040;

  my_data.maxf2_.nsurf2 = 3;
  goto L1050;

/* --SEARCH 2**(NI) NEIGHBORING PARAMETER SETS; */
/* --IF SIGNIFICANT IMPROVEMENT, GO BACK TO BASIC SEARCH PROCESS */

L1040:
  psrch_(dth, ifl);
  if (ifl <= 0)
    goto L1010;

/* --FINISH UP */

L1050:
  lfl = 1;

  if (iflb > 1)
    lfl = 3;

  prepd_(lex);

  if (lex > 0)
  {
    my_data.maxf2_.igfl = 2;
    goto L7030;
  }

  if (my_data.maxf1_.ixvc <= 0)
    goto L7030;

  goto L7000;

/* --METHOD 2:  DIRECT SEARCH METHOD, SKIPPING 2**(NI)-TRIAL SEARCH */

/* --PRINT METHOD IDENTIFICATION AND INITIAL VALUES */

L2000:
/* --START BASIC ITERATION PROCESS */

L2010:
  bsrch_(iflb);

/* --IF REACHED MAXIMUM # ITERATIONS, EXIT */

  if (iflb < 3)
     goto L2020;

  lfl = 4;
  goto L7030;

/* --SKIP OTHER TYPE OF SEARCH IF ONLY ONE PARAMETER TO VARY */

L2020:
  if (my_data.maxf2_.ni <= 1)
    goto L2030;

/* --SEARCH 2(NI) "NEIGHBORING" PARAMETER SETS */

  nsrch_(dth, ifl);
  if (ifl == 0)
  {

/* --SIGNIFICANT IMPROVEMENT; IF < MAXIMUM # ITERATIONS, GO BACK TO BASIC */
/* --ITERATION PROCESS */

    if (my_data.maxf2_.it < my_data.maxf1_.maxit)
      goto L2010;

    lfl = 4;
    goto L7030;
  }

/* --FINISH UP */

L2030:
  lfl = 1;
  if (iflb > 1)
    lfl = 3;

  prepd_(lex);
  if (lex <= 0)
    goto L2040;

  my_data.maxf2_.igfl = 2;
  goto L7030;

L2040:
  if (my_data.maxf1_.ixvc <= 0)
    goto L7030;

  goto L7000;

/* --METHOD 3:  NEWTON-RAPHSON METHOD WITHOUT RECOMPUTATION OF VARIANCE- */
/* --COVARIANCE MATRIX */

/* --PRINT METHOD IDENTIFICATION AND INITIAL VALUES */

L3000:
/* --PREPARE FOR NEWTON-RAPHSON ITERATION */

  prepd_(lex);
  if (lex <= 0)
    goto L3020;

  lfl = lex + 9;
  my_data.maxf2_.igfl = 2;
  return 0;

/* --USE CENTRAL DIFFERENCE FOR GRADIENT CALCULATIONS */

L3020:
  my_data.maxf2_.idif = 2;

/* --COMPUTE NEW VARIANCE-COVARIANCE MATRIX FIRST TIME OR IF PARAMETERS */
/* --BECOME FIXED */

L3030:
  vcmx_(my_data.maxf1_.ihit);
  if (my_data.maxf2_.ivfl < 3)
    goto L3040;

  lfl = 9;
  goto L7030;

/* --COMPUTE GRADIENT VECTOR */

L3040:
  deriv1_();
  if (my_data.maxf2_.igfl <= 0)
    goto L3050;

  lfl = 8;
  goto L7030;

/* --DO ONE NEWTON-RAPHSON ITERATION */

L3050:
  nrstep_(ifl);

/* --GO BACK TO COMPUTE NEW G UNLESS MET SOME STOPPING CRITERION */

  switch ((int)(ifl + 1))
  {
    case 1:  goto L3060;
    case 2:  goto L3080;
    case 3:  goto L3090;
    case 4:  goto L3080;
    case 5:  goto L3120;
  }

/* --CHECK FOR CONVERGENCE TO BOUNDS BEFORE NEXT ITERATION */

L3060:
  fixbnd_(lex);
  switch ((int)(lex + 1))
  {
    case 1:  goto L3040;
    case 2:  goto L3030;
    case 3:  goto L3070;
  }

L3070:
  lfl = 10;
  goto L7030;

/* --CONVERGED; IF ITERATION DONE, CHECK FOR CONVERGENCE TO BOUNDS BEFORE */
/* --QUITTING */

L3080:
  fixbnd_(lex);

L3090:
  lfl = ifl;
  switch ((int)(my_data.maxf1_.ixvc + 1))
  {
    case 1:  goto L7030;
    case 2:  goto L3110;
    case 3:  goto L3100;
  }

L3100:
  if (my_data.maxf1_.ihit <= 0)
    goto L7000;

L3110:
  if (my_data.maxf1_.ixvc + my_data.maxf2_.ivage <= 2)
    goto L7030;

  goto L7000;

/* --REACHED MAXIMUM # ITERATIONS */

L3120:
  lfl = 4;
  goto L7030;

/* --METHOD 4:  NEWTON-RAPHSON METHOD WITH RECOMPUTATION OF VARIANCE- */
/* --COVARIANCE MATRIX AT EACH ITERATION */

/* --PRINT METHOD IDENTIFICATION AND INITIAL VALUES */

L4000:
/* --PREPARE FOR NEWTON-RAPHSON ITERATION */

  prepd_(lex);
  if (lex <= 0)
    goto L4020;

  lfl = lex + 9;
  my_data.maxf2_.igfl = 2;
  return 0;

/* --USE CENTRAL DIFFERENCE FOR GRADIENT CALCULATIONS */

L4020:
  my_data.maxf2_.idif = 2;

/* --COMPUTE NEW VARIANCE-COVARIANCE MATRIX AND VECTOR OF 1ST PARTIAL */
/* --DERIVATIVES */

L4030:
  vcmx_(my_data.maxf1_.ihit);
  if (my_data.maxf2_.ivfl < 3)
    goto L4040;

  lfl = 9;
  goto L7030;

L4040:
  deriv1_();
  if (my_data.maxf2_.igfl <= 0)
    goto L4050;

  lfl = 8;
  goto L7030;

/* --DO ONE NEWTON-RAPHSON ITERATION */

L4050:
  nrstep_(ifl);

/* --GO BACK TO COMPUTE NEW V AND G UNLESS MET SOME STOPPING CRITERION */

  switch ((int)(ifl + 1))
  {
    case 1:  goto L4060;
    case 2:  goto L4070;
    case 3:  goto L4080;
    case 4:  goto L4070;
    case 5:  goto L4110;
  }

/* --CHECK FOR CONVERGENCE TO BOUNDS BEFORE NEXT ITERATION */

L4060:
  fixbnd_(lex);
  if (lex <= 1)
    goto L4030;

  lfl = 10;
  goto L7030;

/* --CONVERGED; IF ITERATION DONE, CHECK FOR CONVERGENCE TO BOUNDS BEFORE */
/* --QUITTING */

L4070:
  fixbnd_(lex);

L4080:
  lfl = ifl;
  switch ((int)(my_data.maxf1_.ixvc + 1))
  {
    case 1:  goto L7030;
    case 2:  goto L4100;
    case 3:  goto L4090;
  }

L4090:
  if (my_data.maxf1_.ihit <= 0)
    goto L7000;

L4100:
  if (my_data.maxf1_.ixvc + my_data.maxf2_.ivage <= 2)
    goto L7030;

  goto L7000;

/* --REACHED MAXIMUM # ITERATIONS */

L4110:
  lfl = 4;
  goto L7030;

/* --METHOD 5:  VARIABLE METRIC METHOD WITH INITIAL B = IDENTITY */

/* --PRINT METHOD IDENTIFICATION */

L5000:

/* --SET SWITCH FOR INITIAL B:  DON'T COMPUTE HESSIAN */

  ihess = 0;

/* --METHOD 6 JOINS PATH AT THIS POINT */

/* --PRINT INITIAL VALUES */

L5010:

/* --PREPARE FOR DERIVATIVE COMPUTATION */

  prepd_(lex);
  if (lex <= 0)
    goto L5020;

  lfl = lex + 9;
  my_data.maxf2_.igfl = 2;
  return 0;

/* --INITIALIZE B MATRIX (ULT, DIAG) AND FLAGS */

L5020:
  binit_(ihess);
  my_data.maxf2_.idif = 1;
  ifirst = 1;

L5030:
  iupdt = 0;

/* --START ITERATION */

/* --COMPUTE NEW GRADIENT G, FIRST SAVING OLD ONE FOR BUPDT */

L5040:
  gpr = my_data.maxf2_.g;
  deriv1_();
  if (my_data.maxf2_.igfl <= 0)
    goto L5060;

  lfl = 8;
  goto L7030;

/* --IF APPROPRIATE, UPDATE B (ULT, DIAG) */

L5060:
  if (iupdt <= 0)
    goto L5070;


  bupdt_(gpr, lex);
  if (lex <= 0)
    goto L5070;

  lfl = 7;
  goto L5210;

/* --COMPUTE NEW DIRECTION OF SEARCH */

L5070:
  direct_(lex);
  switch ((int)(lex + 1))
  {
    case 1:  goto L5100;
    case 2:  goto L5080;
    case 3:  goto L5090;
  }

L5080:
  lfl = 2;
  goto L5210;

L5090:
  if (my_data.maxf2_.idif == 1)
    goto L5200;

  lfl = 6;
  goto L5210;

/* --DO LINE SEARCH IN CHOSEN DIRECTION TO COMPLETE ONE ITERATION */

L5100:
  lsrch_(ifirst, ifl);
  switch ((int)(ifl + 2))
  {
    case 1:  goto L5110;
    case 2:  goto L5120;
    case 3:  goto L5150;
    case 4:  goto L5160;
    case 5:  goto L5180;
    case 6:  goto L5190;
  }

/* --ITERATION COMPLETED, BUT T < INITIAL T AND < TMIN */

L5110:
  iupdt = 0;
  if (my_data.maxf2_.idif == 2)
    goto L5130;

  my_data.maxf2_.idif = 2;
  goto L5130;

/* --ITERATION SUCCESSFULLY COMPLETED  WITH T >= INITIAL T AND/OR TMIN */

L5120:
  iupdt = 1;

/* --CHECK FOR CONVERGENCE TO BOUNDS BEFORE NEXT ITERATION */

L5130:
  fixbnd_(lex);
  switch ((int)(lex + 1))
  {
    case 1:  goto L5040;
    case 2:  goto L5020;
    case 3:  goto L5140;
  }

L5140:
  lfl = 10;
  goto L7030;

/* --CONVERGED BY STANDARD TEST */

L5150:
  lfl = 1;
  goto L5170;

/* --NEGLIGIBLE FUNCTION CHANGE */

L5160:
  lfl = 3;

/* --CHECK FOR CONVERGENCE TO BOUNDS BEFORE QUITTING */

L5170:
  fixbnd_(lex);
  goto L5210;

/* --REACHED MAXIMUM # ITERATIONS */

L5180:
  lfl = 4;
  goto L7030;

/* --COULD NOT COMPLETE ITERATION */

L5190:
  if (my_data.maxf2_.idif == 1)
    goto L5200;

  lfl = 5;
  goto L5210;

/* --SWITCH TO CENTRAL DIFFERENCE COMPUTATION OF GRADIENT */
/* --AND RETRY THIS ITERATION */

L5200:
  my_data.maxf2_.idif = 2;
  goto L5030;

/* --FINISH UP */

L5210:
  if (my_data.maxf1_.ixvc > 0)
    goto L7000;

  goto L7030;

/* --METHOD 6:  VARIABLE METRIC METHOD WITH INITIAL B = -H */
/* --           COMPUTED AT INITIAL ESTIMATES */

/* --PRINT METHOD IDENTIFICATION */

L6000:

/* --SET SWITCH TO COMPUTE INITIAL HESSIAN */

  ihess = 1;

/* --JOIN METHOD 5 PATH */

  goto L5010;

/* --FINAL OUTPUT AND TERMINATION */

/* --OBTAIN VARIANCE-COVARIANCE MATRIX FOR FINAL ESTIMATES */

L7000:
  if (my_data.maxf1_.ixvc >= 2)
    ih = 1;

  if (ih <= my_data.maxf1_.ihit)
    goto L7020;

  for (int i = 0; i < my_data.maxf1_.nt; ++i)
    my_data.maxf2_.stp[i] = my_data.maxf1_.epsd;

L7020:
//    cout << "About to terminate " << endl;
  vcmx_(ih);

/* --AUGMENT AND PRINT FINAL VARIANCE-COVARIANCE MATRIX (IF ONE IS */
/* --AVAILABLE) AND COMPUTE STANDARD DEVIATIONS */

L7030:
  if (my_data.maxf2_.ivfl >= 3)
    goto L7050;

  ireturn = augv_(ih);

/* --COMPUTE (IF NECESSARY) FINAL GRADIENT */

L7050:
  if (my_data.maxf2_.igage < 0)
    goto L7060;
  else if (my_data.maxf2_.igage == 0)
    goto L7080;
  else
    goto L7070;

L7060:
  if (my_data.maxf2_.igfl > 0)
    goto L7080;

L7070:
  my_data.maxf2_.idif = 2;
  deriv1_();

/* --PRINT FINAL VALUES */

L7080:
  fun->depar(theta);

  return 0;

} /* maxfun_ */

int Maxfun::bsrch_(int& ifl)
{
    /* Local variables */
    double athi, thib, dthi, thic      = 0;
    int ifix                           = 0;
    double thim, thip                  = 0;
    int lret, isti, n                  = 0;
    double dbthi                       = 0;
    double thmin, thmax, umtpn, fc     = 0;
    double fm                          = 0;
    double  fp                         = 0;
    double si                          = 0;
    double  fy                         = 0;
    int ndecrm                         = 0;
    int ndecrp                         = 0;
    double ei8                         = 0;
    double den                         = 0;
    double thi                         = 0;
    int lex                            = 0;
    double tpn                         = 0;

    vector<int>    mex(my_data.param_NP);
    vector<double> equiv_0(my_data.param_NP);
    vector<double> thy(my_data.param_NP);
    vector<double> thp(my_data.param_NP);

#define thm (equiv_0)
#define ths (equiv_0)

/* --BASIC ITERATION PROCESS OF DIRECT SEARCH */

/* --RETURNS FLAG IFL: */
/* --   1:  CONVERGED BY CRITERION 1 */
/* --   2:  CONVERGED BY NEGLIGIBLE FUNCTION CHANGE (CRITERION 3) */
/* --   3:  REACHED MAXIMUM # ITERATIONS */

    /* Function Body */
    for (int i = 0; i < my_data.maxf1_.nt; ++i) {
   mex[i] = 0;
    }

/* --EACH ITERATION BEGINS HERE */

L20:
    my_data.maxf2_.thpr = theta;
    my_data.maxf2_.fpr = f;

/* --INITIALIZE FOR THIS ITERATION */

    my_data.maxf2_.impbnd = 0;
    my_data.maxf2_.difmax = 0.;

/* --LOOP THROUGH INDEPENDENT PARAMETERS THAT ARE STILL ALLOWED TO VARY */

    for (int i = 0; i < my_data.maxf1_.nt; ++i) {
   isti = my_data.maxf2_.ist[i];
   if (isti > 2) {
       goto L430;
   }

/* --INITIALIZE WORKING VALUES FOR THIS PARAMETER */

   ifix = 0;
   fy = f;
   thi = theta[i];
   athi = fabs(thi);
   si = my_data.maxf2_.stp[i];
   thmax = my_data.maxf1_.thu[i];
   thmin = my_data.maxf1_.thl[i];

/* --FIRST SELECT 3 APPROPRIATE TRIAL VALUES */

/* --CHECK FOR CLOSENESS TO A BOUND */

   bndchk_(thi, thmax, lex);
   if (lex <= 0) {
       goto L40;
   }
   if (si > efn_(athi, my_data.maxf1_.epsd)) {
       goto L80;
   }
   thib = thmax;
   goto L50;
L40:
   bndchk_(thi, thmin, lex);
   if (lex <= 0) {
       goto L150;
   }
   if (si > efn_(athi, my_data.maxf1_.epsd)) {
       goto L100;
   }
   thib = thmin;

/* --FIX THIS PARAMETER FOR THE REST OF BASIC ITERATION PROCESS */

L50:
        thy = theta;
   thy[i] = thib;
   evaluate(thy, fy, nfe, lex);
   if (lex > 0 || fy < f) {
       goto L70;
   }
   my_data.maxf2_.ist[i] = isti + 4;
   goto L400;

L70:
   my_data.maxf2_.ist[i] = isti + 6;
   goto L430;

/* --CLOSE TO UPPER BOUND */

L80:
   lret = 3;
   thip = thmax;

/* --COME BACK TO THIS POINT (FOR LRET = 3) IF HAD TO DECREASE SI */

L90:
   dbthi = dfn_(fabs(thmax), si);
   thim = thip - dbthi;
   goto L120;

/* --CLOSE TO LOWER BOUND */

L100:
   lret = 2;
   thim = thmin;

/* --COME BACK TO THIS POINT (FOR LRET = 2) IF HAD TO DECREASE SI */

L110:
   dbthi = dfn_(fabs(thmin), si);
   thip = thim + dbthi;

/* --TWO CLOSE-TO-BOUND CASES JOIN AT THIS POINT */

L120:
   thi = (thip + thim) * .5;
   athi = fabs(thi);
   dthi = dbthi * .5;
   si = dfninv_(athi, dthi);
   ei8 = efn_(athi, my_data.maxf1_.epsd) / 8.;

        thy = theta;
   thy[i] = thi;
   evaluate(thy, fy, nfe, lex);
   if (lex > 0) {
       goto L140;
   }
   my_data.maxf2_.difmax = max(my_data.maxf2_.difmax,fabs(fy - my_data.maxf2_.fpr));
   goto L200;

L140:
   my_data.maxf2_.impbnd = 1;
   si *= .5;
   if (si >= ei8) {
       switch ((int)(lret - 1)) {
      case 1:  goto L110;
      case 2:  goto L90;
       }
   }
   ifix = 1;
   goto L420;

/* --NOT CLOSE TO EITHER BOUND */

L150:
   lret = 1;
   ei8 = efn_(athi, my_data.maxf1_.epsd) / 8.;

/* --INCREASE OR DECREASE SI ACCORDING TO RESULT OF PREVIOUS ITERATION
, */
/* --IF ANY */

   if (mex[i] < 0) {
       goto L160;
   } else if (mex[i] == 0) {
       goto L180;
   } else {
       goto L170;
   }

L160:
   si *= .5;
   if (si >= ei8) {
       goto L180;
   }

/* --SI HAS NOW BECOME SMALL IN NORMAL ITERATION PROCESS (NOT DUE TO */
/* --PROBLEMS), SO GO AHEAD AND USE IT BEFORE FIXING PARAMETER */

   ifix = 1;
   goto L180;

L170:
   si = min(si + si,.5);

L180:
   ndecrp = 0;
   ndecrm = 0;

/* --COME BACK TO THIS POINT (FOR LRET = 1) IF HAD TO DECREASE SI */

L190:
   dthi = dfn_(athi, si);
   thip = thi + dthi;
   thim = thi - dthi;

/* --CLOSE-TO-BOUND CASES JOIN NORMAL CASE AT THIS POINT */

/* --CHECK THAT ALL TRIAL VALUES ARE IN BOUNDS */

L200:
   if (thim < thmin || thip > thmax) {
       goto L260;
   }

/* --GET FUNCTION VALUES FOR UPPER AND LOWER TRIAL VALUES */

        thp = theta;
   thp[i] = thip;
   evaluate(thp, fp, nfe, lex);
   if (lex > 0) {
       goto L220;
   }
   my_data.maxf2_.difmax = max(my_data.maxf2_.difmax,fabs(fp - my_data.maxf2_.fpr));
   if (lret <= 1 && fp > fy) {
       goto L280;
   }
   goto L230;

L220:
   my_data.maxf2_.impbnd = 1;
   if (lret > 1) {
       goto L260;
   }
   ++ndecrp;
   if (ndecrp <= 3) {
       goto L260;
   }
   thmax = thi;
   goto L80;

L230:
        thm = theta;
   thm[i] = thim;
   evaluate(thm, fm, nfe, lex);
   if (lex > 0) {
       goto L250;
   }
   my_data.maxf2_.difmax = max(my_data.maxf2_.difmax,fabs(fm - my_data.maxf2_.fpr));
   goto L270;

L250:
   my_data.maxf2_.impbnd = 1;
   if (lret > 1) {
       goto L260;
   }
   ++ndecrm;
   if (ndecrm <= 3) {
       goto L260;
   }
   thmin = thi;
   goto L100;

/* --REDUCE SI AND TRY AGAIN */

L260:
   si *= .5;
   if (si >= ei8) {
       switch ((int)lret) {
      case 1:  goto L190;
      case 2:  goto L110;
      case 3:  goto L90;
       }
   }
   ifix = 1;
   goto L390;

/* --HAVE THREE TRIAL VALUES AND CORRESPONDING FUNCTION VALUES; */
/* --FIND BEST ESTIMATES */

L270:
   if (fy >= fp && fy >= fm) {
       goto L320;
   }

/* --EITHER THIP OR THIM IS BEST */

   if (fm > fp) {
       goto L290;
   }

/* --THIP IS BEST */

   if (fp >= f) {
       goto L280;
   }
   mex[i] = -1;
   goto L420;
L280:
   mex[i] = 1;
   f = fp;
   goto L370;

/* --THIM IS BEST */

L290:
   if (fm >= f) {
       goto L300;
   }
   mex[i] = -1;
   goto L420;
L300:
   mex[i] = 1;
   f = fm;
        theta = thm;
   goto L420;

/* --THI IS BEST (OR AS GOOD) */

L320:
        // Added 2003-04-16:  If fm or fp are -infinity, we replace
        // Them with something exceedingly small.  This number should be <
        // fy by a large margin (10^10 should do it)

        if(fm == -numeric_limits<double>::infinity())
        {
          fm = 1e10 * fy;

          if(fm > 0) fm = -fm;
        }
        if(fp == -numeric_limits<double>::infinity())
        {
          fp = 1e10 * fy;

          if(fp > 0) fp = -fp;
        }

   mex[i] = -1;
   den = (fm + fp - fy - fy) * 2.;
   if (den == 0.) {
       goto L390;
   }

/* --TRY PARABOLIC APPROXIMATION TO GET A BETTER ESTIMATE */

/* --(USE ARRAY THP ALREADY PARTLY PREPARED) */

   thic = thi - (fp - fm) * dthi / den;
   thp[i] = thic;
   evaluate(thp, fc, nfe, lex);
   if (lex <= 0) {
       goto L330;
   }
   my_data.maxf2_.impbnd = 1;
   goto L350;
L330:
   my_data.maxf2_.difmax = max(my_data.maxf2_.difmax,fabs(fc - my_data.maxf2_.fpr));
   if (fc - fy < 0.) {
       goto L350;
   } else if (fc - fy == 0) {
       goto L340;
   } else {
       goto L360;
   }

/* --THIC EQUALLY GOOD; FIX PARAMETER */

L340:
   ifix = 1;
   goto L390;

/* --THI IS BETTER; RETRY WITH SMALLER STEPSIZE */

L350:
   dthi = min(si * .5,fabs(thi - thic));
   si = dfninv_(athi, dthi);
   if (si < ei8 * 8.) {
       goto L390;
   }
   goto L190;

/* --THIC IS BETTER */

L360:
   if (f > fc) {
       goto L420;
   }
   f = fc;
L370:
        theta = thp;
   goto L420;

/* --CHOOSE BETTER OF THI, THETA(I) (DIFFERENT ONLY IF LRET > 1) */

L390:
   if (lret == 1 || f > fy) {
       goto L420;
   }

/* --SET PARAMETER TO THI */

L400:
   f = fy;
        theta = thy;

/*  SAVE CURRENT STEPSIZE FACTOR */

L420:
   my_data.maxf2_.stp[i] = si;

/* --FIX PARAMETER IF INDICATED */

   if (ifix <= 0) {
       goto L430;
   }
   my_data.maxf2_.ist[i] = isti + 8;

/* --END OF LOOP */

L430:
   ;
    }

/* --FINISH UP THIS ITERATION */

    endit_();
    my_data.maxf2_.difmax = max(my_data.maxf2_.difmax,my_data.maxf2_.fch);

/* --WRITE OUT DETAILS OF THIS ITERATION IF DESIRED */

/* --CHECK FOR CONVERGENCE */

    bcnvch_(lex);
    switch ((int)(lex + 1)) {
   case 1:  goto L440;
   case 2:  goto L500;
   case 3:  goto L500;
   case 4:  goto L540;
    }

/* --DO "SUPPLEMENTARY ITERATION PROCESS" BEFORE DOING ANOTHER */
/* --BASIC ITERATION */

/* --INITIALIZE */

L440:
    ths = theta;
    n = 0;
    tpn = 2.;

/* --GO FARTHER IN SAME DIRECTION */

L460:
    umtpn = 1. - tpn;
    for (int i = 0; i < my_data.maxf1_.nt; ++i) {
   if (my_data.maxf2_.ist[i] > 2) {
       thy[i] = ths[i];
   } else {
       thy[i] = tpn * ths[i] + umtpn * my_data.maxf2_.thpr[i];
       if (thy[i] < my_data.maxf1_.thl[i] || thy[i] >
          my_data.maxf1_.thu[i]) {
      goto L490;
       }
   }
    }
    evaluate(thy, fy, nfe, lex);
    if (lex > 0 || fy <= f) {
   goto L490;
    }

/* --THY CONTAINS IMPROVED ESTIMATES */

    theta = thy;
    f = fy;

    ++n;
    tpn += tpn;
    goto L460;

/* --NO MORE IMPROVEMENT */

L490:
    if (n <= 0) {
   goto L20;
    }

    --n;
    goto L20;

/* --HAVE CONVERGED; */
/* --BEFORE RETURNING, UNFIX ANY TEMPORARILY FIXED PARAMETERS */

L500:
    for (int i = 0; i < my_data.maxf1_.nt; ++i) {
   switch ((int)(my_data.maxf2_.ist[i])) {
       case 1:  goto L530;
       case 2:  goto L530;
       case 3:  goto L530;
       case 4:  goto L530;
       case 5:  goto L510;
       case 6:  goto L520;
       case 7:  goto L510;
       case 8:  goto L520;
       case 9:  goto L510;
       case 10:  goto L520;
   }
L510:
   my_data.maxf2_.ist[i] = 1;
   goto L530;
L520:
   my_data.maxf2_.ist[i] = 2;
L530:
   ;
    }

/* --HAVE CONVERGED OR REACHED MAXIMUM # ITERATIONS; EXIT */

L540:
    ifl = lex;
    return 0;

} /* bsrch_ */

#undef ths
#undef thm


int Maxfun::nsrch_(vector<double>& dth, int& ifl)
{
    /* Local variables */
    double athi;
    double thiy;
    double fs = 0;
    double si = 0;
    double fy = 0;
    double thi;
    int lex;

/* --TRY 2*(NI) "NEIGHBORING" PARAMETER SETS TO FIND A BETTER ONE */

/* --RETURNS FLAG IFL: */
/* --   0:  SIGNIFICANT CHANGE; CONTINUE ITERATING */
/* --   1:  NEGLIGIBLE OR NO CHANGE */


/* --SAVE INCOMING ESTIMATES AND ESTABLISH NEW STEPSIZES */

    vector<double> ths(my_data.param_NP);
    vector<double> thy(my_data.param_NP);

    /* Function Body */
    ths = theta;
    thy = theta;

    for (int i = 0; i < my_data.maxf1_.nt; ++i) {
   athi = fabs(thy[i]);
   si = efn_(athi, my_data.maxf1_.epsd);
   dth[i] = dfn_(athi, si);
   my_data.maxf2_.stp[i] = si;
    }
    fs = f;

/* --INITIALIZE BOUNDARY FLAG */

    my_data.maxf2_.impbnd = 0;

/* --LOOP THROUGH PARAMETERS TO CHECK NEIGHBORING VALUES */

    for (int i = 0; i < my_data.maxf1_.nt; ++i) {
   if (my_data.maxf2_.ist[i] > 2) {
       goto L80;
   }
   thi = thy[i];
   thiy = thi + dth[i];
   if (thiy > my_data.maxf1_.thu[i]) {
       goto L40;
   }
   thy[i] = thiy;
   evaluate(thy, fy, nfe, lex);
   if (lex <= 0) {
       goto L20;
   }
   my_data.maxf2_.impbnd = 1;
   goto L40;
L20:
   if (fy <= f) {
       goto L40;
   }
        theta = thy;
   f = fy;
L40:
   thiy = thi - dth[i];
   if (thiy < my_data.maxf1_.thl[i]) {
       goto L70;
   }
   thy[i] = thiy;
   evaluate(thy, fy, nfe, lex);
   if (lex <= 0) {
       goto L50;
   }
   my_data.maxf2_.impbnd = 1;
   goto L70;
L50:
   if (fy <= f) {
       goto L70;
   }

        theta = thy;
   f = fy;
L70:
   thy[i] = theta[i];
L80:
   ;
    }

/* --IF ANY IMPROVEMENT HAS OCCURRED, CONSIDER THIS ANOTHER ITERATION */

    if (f > fs) {
   goto L90;
    }
    ifl = 1;
    return 0;

L90:
    my_data.maxf2_.thpr = ths;
    my_data.maxf2_.fpr = fs;
    endit_();
    my_data.maxf2_.difmax = my_data.maxf2_.fch;

/* --PRINT ITERATION DETAILS IF DESIRED */

/* --CHECK FOR CONVERGENCE */

    if (my_data.maxf2_.fch > my_data.maxf1_.epsc1 * my_data.maxf1_.epsc1 && my_data.maxf2_.fch > my_data.maxf1_.epsc3) {
   goto L120;
    }
    ifl = 1;
    return 0;

/* --SIGNIFICANT CHANGE */

L120:
    for (int i = 0; i < my_data.maxf1_.nt; ++i) {
   my_data.maxf2_.stp[i] *= 2.;
    }
    ifl = 0;
    return 0;

} /* nsrch_ */

int Maxfun::psrch_(vector<double>& dth, int& ifl)
{
    /* Local variables */
    double thiy;
    int k;
    int nplus;
    double fs = 0;
    double fy = 0;
    int np1, lex;
    int nsw;

/* --SEARCH 2**(NI) NEIGHBORING PARAMETER SETS FOR BETTER ESTIMATES */

/* --RETURNS FLAG IFL: */
/* --   0:  SIGNIFICANT CHANGE; CONTINUE ITERATING */
/* --   1:  NEGLIGIBLE OR NO CHANGE */


/* --SAVE INCOMING VALUES */

    vector<int> list(my_data.param_NPV+1);
    vector<double> ths(my_data.param_NP);
    vector<double> thy(my_data.param_NP);

    /* Function Body */

    ths = theta;
    fs  = f;

/* --INITIALIZE BOUNDARY FLAG */

    my_data.maxf2_.impbnd = 0;

/* --CONSIDER TRIAL PARAMETER SETS WITH NPLUS INDEPENDENT PARAMETERS */
/* --CHANGED IN + DIRECTION (AND OTHERS IN - DIRECTION), */
/* --FOR NPLUS = 0 TO NI */

/* --NPLUS = 0 CASE:  NO PARAMETERS TO CHANGE IN + DIRECTION */

    nplus = 0;
    list[0] = 0;

/* --SET UP TRIAL PARAMETER SET WITH NPLUS PARAMETERS CHANGING IN */
/* -- + DIRECTION (WHICH ONES ARE INDICATED IN LIST ARRAY) */

L20:
    np1 = 0;
    k = 0;
    for (int i = 0; i < my_data.maxf1_.nt; ++i) {
      if (my_data.maxf2_.ist[i] <= 2) {
         if (k == list[np1]) {
             thy[i] = ths[i] + dth[i];
             ++np1;
         } else {
             thy[i] = ths[i] - dth[i];
         }
         ++k;
         if (thy[i] > my_data.maxf1_.thu[i] || thy[i] < my_data.maxf1_.thl[i]) {
             goto L70;
         }
      }
      else
        thy[i] = ths[i];
    }

/* --CHECK OUT THIS TRIAL SET */

    evaluate(thy, fy, nfe, lex);
    if (lex <= 0) {
   goto L50;
    }
    my_data.maxf2_.impbnd = 1;
    goto L70;
L50:
    if (fy <= f) {
   goto L70;
    }

/* --HAVE AN IMPROVEMENT */

    theta = thy;
    f = fy;

/* --FINISHED WITH THIS TRIAL SET; GET NEXT ONE */

L70:
    if (nplus > 0) {
   goto L90;
    }

/* --DONE WITH THIS NPLUS; INCREASE # PARAMETERS TO CHANGE IN + DIRECTION
*/

L80:
    ++nplus;
    if (nplus > my_data.maxf2_.ni) {
   goto L100;
    }
    nsw = 0;

/* --GET LIST OF NPLUS PARAMETERS TO CHANGE IN + DIRECTION */

L90:
    comb_(my_data.maxf2_.ni, nplus, nsw, list);
    if (nsw > 1) {
   goto L80;
    }
    goto L20;

/* --HAVE TRIED EVERY DIRECTION OF CHANGE; ANY IMPROVEMENT? */

L100:
    if (f > fs) {
   goto L110;
    }

    ifl = 1;
    return 0;

/* --INCREMENT # TIMES 2**(NI)-TRIAL SEARCH MADE IMPROVEMENT */

L110:
    ++my_data.maxf2_.nsurf2;

/* --TRY GOING FARTHER IN BEST DIRECTION */

    for (int i = 0; i < my_data.maxf1_.nt; ++i) {
   thiy = theta[i] * 2. - ths[i];
   if (my_data.maxf2_.ist[i] > 2 || thiy < my_data.maxf1_.thl[i] || thiy > my_data.maxf1_.thu[i]) {
       thy[i] = theta[i];
   } else {
       thy[i] = thiy;
   }
    }
    evaluate(thy, fy, nfe, lex);
    if (lex <= 0) {
   goto L130;
    }
    my_data.maxf2_.impbnd = 1;
    goto L150;

L130:
    if (fy <= f) {
   goto L150;
    }

    theta = thy;
    f = fy;

/* --COMPLETE THIS "ITERATION" */

L150:
    my_data.maxf2_.thpr = ths;
    my_data.maxf2_.fpr = fs;
    endit_();
    my_data.maxf2_.difmax = my_data.maxf2_.fch;

/* --PRINT ITERATION DETAILS IF DESIRED */

/* --CHECK FOR CONVERGENCE */

    if (my_data.maxf2_.fch > my_data.maxf1_.epsc1 * my_data.maxf1_.epsc1 && my_data.maxf2_.fch > my_data.maxf1_.epsc3) {
   goto L180;
    }
    ifl = 1;
    return 0;

/* --SIGNIFICANT CHANGE */

L180:
    for (int i = 0; i < my_data.maxf1_.nt; ++i) {
   my_data.maxf2_.stp[i] = my_data.maxf1_.epsd;
    }

    ifl = 0;
    return 0;

} /* psrch_ */

int Maxfun::comb_(int ntotal, int nchoos, int& nsw, vector<int>& list)
{
    /* System generated locals */

    /* Local variables */
    int notch;

/* --CHOOSE NCHOOS ELEMENTS OUT OF NTOTAL; INDICATE WHICH ONES IN LIST */
/* --(GENERATES COMBINATIONS TO BE USED IN THE 2**(NI)-TRIAL SEARCH: */
/* --INDICATE WHICH PARAMETERS TO CHANGE IN + DIRECTION) */

    /* Function Body */
    if (nsw <= 0)
    {
      /* --START WITH A NEW NCHOOS */

      notch = ntotal - nchoos;

      for (int k = 0; k < nchoos; ++k)
   maxlst[k] = notch + k;

      nsw = 1;
      list[nchoos] = 0;
      kv = 0;
      list[0] = 0;
  }
  else
  /* --IN PROGRESS WITH THIS NCHOOS */
  {
    ++list[kv];

    // [kbj] Not sure about <= or < here.
    if (list[kv] <= maxlst[kv])
   return 0;

    while(kv >= 0 && list[kv] > maxlst[kv])
    {
      --kv;
      if(kv >= 0)
        ++list[kv];
    }

    if (kv < 0)
    {
      nsw = 2;
      return 0;
    }
  }

/* --HAVE ESTABLISHED KV'TH CHOSEN ELEMENT; FILL IN REST OF LIST WITH */
/* --LOWEST POSSIBLE VALUES */

    for (int k = kv + 1; k < nchoos; ++k)
      list[k] = list[k - 1] + 1;

    kv = nchoos - 1;
    return 0;

} /* comb_ */

int Maxfun::nrstep_(int& ifl)
{
    /* Local variables */
    int l;
    double fs, pl, tt;
    int lex;

/* --PERFORM ONE NEWTON-RAPHSON ITERATION */

/* --RETURNS FLAG IFL: */
/* --   0:  CONTINUE ITERATING */
/* --   1:  CONVERGED BY CRITERION 1 (IMPLIED TEST) */
/* --   2:  CONVERGED BY CRITERION 2:  PTG <= EPSC2 */
/* --   3:  CONVERGED BY NEGLIGIBLE FUNCTION CHANGE (CRITERION 3) */
/* --   4:  REACHED MAXIMUM # ITERATIONS */


/* --SAVE INCOMING VALUES */

    vector<double> ths = theta;
    fs = f;

/* --COMPUTE DIRECTION OF CHANGE */

    mmult_(my_data.maxf2_.v, my_data.param_NPV, my_data.maxf2_.nv, my_data.maxf2_.nv, my_data.maxf2_.g,
           my_data.param_NPV, 1, my_data.maxf2_.pdir, my_data.param_NPV);

/* --REDUCE AMOUNT OF CHANGE IF NECESSARY TO PRESERVE BOUNDS ON */
/* --INDEPENDENT PARAMETERS */

    my_data.maxf2_.ptg = 0.;
    tt = 1.;
    l = -1;

    for (int i = 0; i < my_data.maxf1_.nt; ++i) {
   if (my_data.maxf2_.ist[i] > 2)
          continue;

   ++l;
   pl = my_data.maxf2_.pdir[l];

   if (pl == 0.)
          continue;

   if (pl > 0.) {
       tt = min(tt,(my_data.maxf1_.thu[i] - ths[i]) / pl);

            while(1)
            {
         if (ths[i] + tt * pl <= my_data.maxf1_.thu[i]) {
                my_data.maxf2_.ptg += pl * my_data.maxf2_.g[l];
                break;
         }
         tt = (1. - my_data.maxf1_.yota) * tt;
            }
   } else {
       tt = min(tt,(my_data.maxf1_.thl[i] - ths[i]) / pl);
            while(1)
            {
         if (ths[i] + tt * pl >= my_data.maxf1_.thl[i]) {
                my_data.maxf2_.ptg += pl * my_data.maxf2_.g[l];
                break;
         }
         tt = (1. - my_data.maxf1_.yota) * tt;
            }
   }
    }

/* --SEE IF PTG SMALL ENOUGH TO STOP */

    if (fabs(my_data.maxf2_.ptg) > my_data.maxf1_.epsc2) {
   goto L30;
    }
    ifl = 2;
    return 0;

/* --INITIALIZE BOUNDARY FLAG */

L30:
    my_data.maxf2_.impbnd = 0;

/* --COMPUTE NEW ESTIMATES */

L40:
    l = -1;
    for (int i = 0; i < my_data.maxf1_.nt; ++i) {
   if (my_data.maxf2_.ist[i] > 2)
          continue;

   ++l;

   theta[i] = ths[i] + tt * my_data.maxf2_.pdir[l];
    }

/* --GET FUNCTIONAL VALUE IF POSSIBLE */

    evaluate(theta, f, nfe, lex);

    if (lex <= 0) {
   goto L60;
    }
    my_data.maxf2_.impbnd = 1;
    goto L70;

/* --CHECK FOR DECREASING VALUE */

L60:
    if (f >= fs) {
   goto L80;
    }

/* --RAN INTO TROUBLE; CUT DOWN STEP SIZE AND TRY AGAIN */

L70:
    tt *= .5;
    goto L40;

/* --OK TO COMPLETE THIS ITERATION */

L80:
    my_data.maxf2_.tstep = tt;
    my_data.maxf2_.thpr = ths;
    my_data.maxf2_.fpr = fs;
    endit_();

/* --INCREMENT NUMBER OF ITERATIONS WITH THIS VARIANCE-COVARIANCE MATRIX
*/
/* --AND THIS GRADIENT VECTOR */

    ++my_data.maxf2_.ivage;
    ++my_data.maxf2_.igage;

/* --WRITE OUT DETAILS OF THIS ITERATION IF DESIRED */
/* --CHECK FOR CONVERGENCE */

    vcnvch_(lex);

/* --SET FLAG AND RETURN */

    if (lex <= 1) {
   ifl = lex;
    } else {
   ifl = lex + 1;
    }
    return 0;

} /* nrstep_ */

int Maxfun::binit_(int& ihess)
{
    /* Local variables */
    double bl, ultl2l;
    int lex;

#define b_ref(a_1,a_2) b[(a_2)*my_data.param_NPV + a_1]
#define h_ref(a_1,a_2) my_data.maxf2_.h[(a_2)*my_data.param_NPV +a_1]
#define ult_ref(a_1,a_2) my_data.maxf2_.ult[(a_2)*my_data.param_NPV +a_1]


/* --OBTAIN INITIAL B MATRIX (APPROXIMATING OPPOSITE OF HESSIAN AT OPTIMAL */
/* --ESTIMATES) IN FACTORIZED FORM B = ULT*DIAG*ULT' */
/* --(NOTE:  OPPOSITES ARE NECESSARY SINCE UPDATING PROCEDURE REQUIRES */
/* --POSITIVE DEFINITE B) */


/* --SELECT TYPE OF INITIAL MATRIX TO BE USED */

    vector<double> b(my_data.param_NPV*my_data.param_NPV);

    /* Function Body */
    if (ihess <= 0) {
   goto L50;
    }

/* --COMPUTE MATRIX OF SECOND PARTIAL DERIVATIVES FOR INITIAL ESTIMATES */
/* --OF PARAMETERS */

    deriv2_(my_data.maxf1_.ihit, lex);
    if (lex > 1) {
   goto L50;
    }

/* --INITIALIZE B TO -H */

    for (int l1 = 0; l1 < my_data.maxf2_.nv; ++l1) {
   for (int l2 = 0; l2 < my_data.maxf2_.nv; ++l2) {
       b_ref(l1, l2) = -h_ref(l1, l2);
   }
    }

/* --DO CHOLESKY FACTORIZATION */

    for (int l = 0; l < my_data.maxf2_.nv; ++l) {
   bl = b_ref(l, l);
   if (bl <= 0.)
     goto L50;

   my_data.maxf2_.diag[l] = bl;
   ult_ref(l, l) = 1.;
   if (l+1 >= my_data.maxf2_.nv)
          return 0;
   for (int l1 = l + 1; l1 < my_data.maxf2_.nv; ++l1) {
       ult_ref(l1, l) = b_ref(l1, l) / bl;
   }
   for (int l2 = l + 1; l2 < my_data.maxf2_.nv; ++l2) {
       ultl2l = ult_ref(l2, l);
       for (int l1 = l2; l1 < my_data.maxf2_.nv; ++l1) {
      b_ref(l1, l2) = b_ref(l1, l2) - ultl2l * ult_ref(l1, l) * bl;
       }
   }
    }
    return 0;

/* --USE IDENTITY MATRIX */

L50:
    for (int l1 = 0; l1 < my_data.maxf2_.nv; ++l1) {
   for (int l2 = 0; l2 < l1; ++l2) {
       ult_ref(l1, l2) = 0.;
   }
   ult_ref(l1, l1) = 1.;
   my_data.maxf2_.diag[l1] = 1.;
    }

    return 0;

} /* binit_ */

#undef ult_ref
#undef h_ref
#undef b_ref


int Maxfun::bupdt_(vector<double>& gpr, int& lex)
{
    /* Local variables */
    double sdbc, sden, rdiv, alph1, alph2;
    int l;
    double  alpha, lambl, gamtl, bettl;
    int iswup;
    double lambl2, aa, bb, cc, al, bl,  sb, sc;
    int ll;
    double sd, sa, nu, sbc, sdb, eta, phl, mul,
             ytp, del1, del2, tau1, tau2, fbcd;

#define ult_ref(a_1,a_2) my_data.maxf2_.ult[(a_2)*my_data.param_NPV + a_1]

/* --UPDATE MATRIX B = ULT*DIAG*ULT' */

/* --RETURNS FLAG LEX: */
/* --   0:  UPDATE SUCCESSFUL */
/* --   1:  UPDATE FAILED:  LOST ALL SIGNIFICANT DIGITS IN OPTIMAL */
/* --       CONDITIONING */

/* --USING PREVIOUS GRADIENT (BEFORE ITERATION IT) IN GPR AND CURRENT */
/* --GRADIENT (AFTER ITERATION IT) IN G */

/* ************************************************************************/
/* * NOTE THAT G, GPR, Y ARE DEFINED OPPOSITE TO QUANTITIES DESCRIBED IN **/
/* * GEMINI DOCUMENTATION, AND FORMULAS ARE CHANGED ACCORDINGLY          **/
/* ************************************************************************/
/* * NOTE ALSO THAT VARIABLES TAU1, TAU2, AND PHL REFER TO QUANTITIES    **/
/* * THET1, THET2, AND THL IN GEMINI; THESE DO NOT REFER TO THE PARAMETER**/
/* * ESTIMATE VECTOR THETA                                               **/
/* ************************************************************************/


/* --CHECK THAT DIRECTION OF CHANGE IS ONE OF DECREASING GRADIENT; */
/* --IF NOT, DON'T UPDATE */

    vector<double> wtil(my_data.param_NPV);
    vector<double> ztil(my_data.param_NPV);
    vector<double> wplp1(my_data.param_NPV);
    vector<double> zplp1(my_data.param_NPV);
    vector<double> s(my_data.param_NPV);
    vector<double> w(my_data.param_NPV);
    vector<double> y(my_data.param_NPV);
    vector<double> z__(my_data.param_NPV);
    vector<double> dp(my_data.param_NPV);
    vector<double> wtplp1(my_data.param_NPV);
    vector<double> ztplp1(my_data.param_NPV);
    vector<double> wpl(my_data.param_NPV);
    vector<double> zpl(my_data.param_NPV);

    /* Function Body */
    ytp = 0.;
    for (int l = 0; l < my_data.maxf2_.nv; ++l) {
   y[l] = my_data.maxf2_.g[l] - gpr[l];
   ytp += y[l] * my_data.maxf2_.pdir[l];
    }

    if (ytp >= 0.) {
   goto L260;
    }

    if (my_data.maxf2_.nv > 1) {
   goto L20;
    }

/* --FOR NV = 1: */

    my_data.maxf2_.diag[0] = -y[0] * my_data.maxf2_.diag[0] / (my_data.maxf2_.tstep * gpr[0]);
    goto L260;

/* --FOR NV > 1 (FORMULA NUMBERS REFER TO GEMINI DOCUMENTATION): */

/* -----------------------------------------------------------------------
 */

/* --GOLDFARB METHOD FOR VARIABLE METRIC B-MATRIX UPDATING; */
/* --UPDATE OF FACTORIZED B-MATRIX BY RANK TWO MODIFICATION */
/* --IN REAL PRODUCT FORM WITH FORMULA (3.25) OR (3.28) OF GEMINI */
/* --DOCUMENTATION, DEPENDING ON WHETHER NONE OR AT LEAST ONE */
/* --LAMBA(L)**2 IS GREATER THAN 10. */
/* --OPTIMAL CONDITIONING OF DAVIDON ALSO INCORPORATED. */

/* -----------------------------------------------------------------------
 */
/* --SOLVE (3.20) AND (3.21) FOR WTIL AND ZTIL */

L20:
    wtil[0] = gpr[0];
    ztil[0] = -y[0];

    for (int l = 1; l < my_data.maxf2_.nv; ++l) {
   wtil[l] = gpr[l];
   ztil[l] = -y[l];
   for (int ll = 0; ll < l; ++ll) {
       wtil[l] -= ult_ref(l, ll) * wtil[ll];
       ztil[l] -= ult_ref(l, ll) * ztil[ll];
   }
    }

/* -----------------------------------------------------------------------
 */
/* --GET SCALARS B, C AND D, DENOTED SB, SC, SD HERE (3.27) */

    sb = 0.;
    sc = 0.;
    sd = 0.;
    for (int l = 0; l < my_data.maxf2_.nv; ++l) {
   sb += my_data.maxf2_.tstep * ztil[l] * wtil[l] / my_data.maxf2_.diag[l];
   sc += my_data.maxf2_.tstep * my_data.maxf2_.tstep * wtil[l] * wtil[l] /
      my_data.maxf2_.diag[l];
   sd += ztil[l] * ztil[l] / my_data.maxf2_.diag[l];
    }

/* -----------------------------------------------------------------------
 */
/* --OPTIMAL CONDITIONING OF UPDATE, ACCORDING TO THEOREM 3 OF DAVIDON */
/* --(1975) -- MAY BE DISABLED BY "GO TO 60" */
/* --------------------- */

/* --PREVENTING UNDERFLOW: */
/* --IF ANY OF SB, SC, SD IS SMALLER THAN 1.D-10, USE SR1-UPDATE FOR B */

    fbcd = min(min(sb,sc),sd);
    if (fbcd <= 1e-10) {
   goto L50;
    }
    fbcd = sc * 2. * (sd / sb) / (sc + sd);
    if (fbcd < 1.) {
   goto L50;
    }

/* --------------------- */
/* --RANK TWO UPDATE */

/* --GET ALPHA */

    aa = sb / sc - sd / sb * 2. + sd / sc;
    bb = sb / sc - 1.;
    cc = 1. - sb / sd;
    del2 = bb * bb - aa * cc;

/* --IF DISCRIMINANT NEGATIVE OR EQUALS ZERO, TAKE ALPHA EQUAL TO ZERO */

    if (del2 <= 1e-8) {
   goto L60;
    }

/* ---------------- */

    del1 = sqrt(del2);
    alph1 = (-bb + del1) / aa;
    alph2 = (-bb - del1) / aa;

/* --FOR NOW, ALWAYS CHOOSE SOLUTION OF SMALLEST MODULUS */

    if (fabs(alph1) <= fabs(alph2)) {
   alpha = alph1;
    } else {
   alpha = alph2;
    }

/* --IF ALPHA VERY SMALL, ALPHA TAKEN EQUAL TO ZERO */

    if (fabs(alpha) <= 1e-5) {
   goto L60;
    }

/* --GET SA */

    sa = (alpha + 1.) * (alpha + 1.) + sc / sb - alpha * alpha * (sc / sb)
                      * (sd / sb) - 1. + alpha * alpha * (sd / sb);
    if (sa <= 0.) {
   sa = 0.;
    } else {
   sa = 1. / (sqrt(sa) * sb);
    }

/* --GET TAU1 AND TAU FOR NON-TRIVIAL ALPHA */

    rdiv = 1. / (alpha * alpha * sd + alpha * 2. * sb + sc);
    tau1 = -(sa * (alpha * (alpha * sd + sb)) + 1.) * rdiv;
    tau2 = sa + (alpha * sa * (sc + alpha * sb) - alpha) * rdiv;
    goto L70;

/* --------------------- */
/* --SR1 (SYMMETRIC RANK-ONE) UPDATE */

L50:
    alpha = -1.;
    sbc = sb - sc;
    if (sbc == 0.) {
   goto L290;
    }
    sdb = sd - sb;
    if (sdb == 0.) {
   goto L290;
    }
    sden = sqrt((fabs(sbc))) * sqrt((fabs(sdb)));
    if (sden == 0.) {
   goto L290;
    }
    sa = 1. / sden;
    sdbc = sd - sb * 2. + sc;
    if (sdbc == 0.) {
   goto L290;
    }
    tau1 = -(sdb * sa + 1.) / sdbc;
    tau2 = sa + (sa * sbc + 1.) / sdbc;
    goto L70;

/* --------------------- */
/* --ALPHA ZERO:  DFP-UPDATE FOR B */

L60:
    alpha = 0.;
    sa = 1. / (sqrt(sb) * sqrt(sc));
    tau1 = -1. / sc;
    tau2 = sa;

/* -----------------------------------------------------------------------
 */
/* --START WITH SWITCH SET TO DO (3.25) TYPE UPDATE */

L70:
    iswup = 1;

/* --SAVE ULT LOWER TRIANGLE IN UPPER TRIANGLE OF ARRAY ULT, IN CASE A */
/* --SWITCH TO (3.28) SHOULD BE REQUIRED */

    for (int l = 1; l < my_data.maxf2_.nv; ++l) {
   for (int ll = 0; ll < l; ++ll) {
       ult_ref(ll, l) = ult_ref(l, ll);
   }
    }

/* --GET W AND Z AS GIVEN IN (3.22), (3.23) */

    for (int l = 0; l < my_data.maxf2_.nv; ++l) {
   w[l] = my_data.maxf2_.tstep * wtil[l] + alpha * ztil[l];
   z__[l] = my_data.maxf2_.tstep * tau1 * wtil[l] + tau2 * ztil[l];
    }

/* -----------------------------------------------------------------------
 */
/* --READY TO APPLY GOLDFARB RECURRENCE 3 TO SOLVE CONCURRENTLY (3.25), */
/* --(3.24) BEING ALSO SOLVED SIMULTANEOUSLY */
/* --(OR (3.28) AND CORRESPONDING ANALOG OF (3.24)) */

/* --GET S FIRST (P.802 (TOP) OF ?) */

    s[my_data.maxf2_.nv-1] = 0.;
    for (int l = my_data.maxf2_.nv-1; l > 0; --l) {
   s[l - 1] = s[l] + w[l] * w[l] / my_data.maxf2_.diag[l];
    }

/* -----------------------------------------------------------------------
 */
/* --INITIALIZE NU, ETA */

L110:
    nu = 1.;
    eta = 0.;

/* -----------------------------------------------------------------------
 */
/* --INITIALIZE WTPLP1, ZTPLP1 OR WPLP1, ZPLP1 */

    if (iswup > 1) {
   goto L130;
    }
    for (int l = 0; l < my_data.maxf2_.nv; ++l) {
   wtplp1[l] = gpr[l];
   ztplp1[l] = -y[l];
    }
    goto L150;
L130:
    for (int l = 0; l < my_data.maxf2_.nv; ++l) {
   wplp1[l] = my_data.maxf2_.tstep * my_data.maxf2_.g[l] - alpha * y[l];
   zplp1[l] = my_data.maxf2_.tstep * tau1 * my_data.maxf2_.g[l] - tau2 * y[l];
    }

/* -----------------------------------------------------------------------
 */

L150:
    for (int l = 0; l+1 < my_data.maxf2_.nv; ++l) {

/* --RECURRENCE FORMULA FROM (3.24) OR ANALOG */

   if (iswup > 1) {
       goto L170;
   }
   for (int ll = l + 1; ll < my_data.maxf2_.nv; ++ll) {
       wtplp1[ll] -= wtil[l] * ult_ref(ll, l);
       ztplp1[ll] -= ztil[l] * ult_ref(ll, l);
   }
   goto L190;
L170:
   for (int ll = l + 1; ll < my_data.maxf2_.nv; ++ll) {
       wpl[ll] = wplp1[ll];
       zpl[ll] = zplp1[ll];
       wplp1[ll] = wpl[ll] - w[l]   * ult_ref(ll, l);
       zplp1[ll] = zpl[ll] - z__[l] * ult_ref(ll, l);
   }

/* --RECURRENCE 3 (?) TO GET AL, BL, ETC... */

L190:
   al = nu * z__[l] - eta * w[l];
   phl = al * w[l] / my_data.maxf2_.diag[l] + 1.;
   lambl2 = phl * phl + al * al * s[l] / my_data.maxf2_.diag[l];

/* --SWITCH TO (3.28) UPDATE IF ANY LAMBL2 GREATER THAN 4 */

   if (iswup <= 1 && lambl2 > 4.) {
       goto L270;
   }

/* --COMPUTE L'TH ELEMENT OF D+ (NEW DIAG) */

   dp[l] = my_data.maxf2_.diag[l] * lambl2;

/* --TAKES SIGN OF LAMBDA OPPOSITE TO THAT OF THETA */

   lambl = sqrt(lambl2);
   if (phl > 0.) {
       lambl = -lambl;
   }

   mul = phl - lambl;
   bl = phl * w[l] + al * s[l];

/* --NOTE:  GAMTL AND BETTL STAND FOR GAMMA TILDE AND BETA TILDE */
/* --RESPECTIVELY IN (3.26) */

   gamtl = bl * nu / (lambl2 * my_data.maxf2_.diag[l]);
   bettl = (al - bl * eta) / (lambl2 * my_data.maxf2_.diag[l]);

/* --UPDATE L'TH COLUMN OF ULT */

   if (iswup > 1) {
       goto L210;
   }
   for (int ll = l + 1; ll < my_data.maxf2_.nv; ++ll) {
       ult_ref(ll, l) = ult_ref(ll, l) + my_data.maxf2_.tstep
                           * (bettl + tau1 * gamtl)
                           * wtplp1[ll]
                           + (alpha * bettl + tau2 * gamtl) * ztplp1[ll];
   }
   goto L230;
L210:
   for (int ll = l + 1; ll < my_data.maxf2_.nv; ++ll) {
       ult_ref(ll, l) = ult_ref(ll, l) / lambl2 + bettl * wpl[ll]
                           + gamtl * zpl[ll];
   }

/* --UPDATE NU, ETA */

L230:
   nu = -nu / lambl;
   eta = -(eta + al * al / (mul * my_data.maxf2_.diag[l])) / lambl;

/* L240: */
    }

/* --GET LAST LAMBDA FOR D+ */

    al = nu * z__[my_data.maxf2_.nv - 1] - eta * w[my_data.maxf2_.nv - 1];
    lambl = al * w[my_data.maxf2_.nv - 1] / my_data.maxf2_.diag[my_data.maxf2_.nv - 1] + 1.;

    dp[my_data.maxf2_.nv - 1] = my_data.maxf2_.diag[my_data.maxf2_.nv - 1] * lambl * lambl;

/* --UPDATE DIAG FROM D+ */

    my_data.maxf2_.diag = dp;

/* --NORMAL EXIT */

L260:
    lex = 0;
    return 0;

/* --SWITCH TO UPDATING PROCEDURE 2 */

L270:
    iswup = 2;
    for (l = 1; l < my_data.maxf2_.nv; ++l) {
   for (ll = 0; ll < l; ++ll) {
       ult_ref(l, ll) = ult_ref(ll, l);
   }
    }
    goto L110;

/* --ERROR EXIT IN CASE IF DIVISION BY ZERO WOULD OCCUR IN OPTIMAL */
/* --CONDITIONING */

L290:
    lex = 1;
    return 0;

} /* bupdt_ */

#undef ult_ref


int Maxfun::direct_(int& lex)
{
#define ult_ref(a_1,a_2) my_data.maxf2_.ult[(a_2)*my_data.param_NPV + a_1]

/* --COMPUTE DIRECTION OF SEARCH (SOLVING (-ULT*DIAG*ULT')*PDIR = -G FOR
*/
/* --PDIR) AND NORMALIZED GRADIENT */

/* --RETURNS FLAG LEX: */
/* --   0:  PDIR SUCCESSFULLY DETERMINED; CONTINUE */
/* --   1:  CONVERGED BY CRITERION 2:  PTG <= EPSC2 */
/* --   2:  PDIR NOT IN ASCENT DIRECTION */

    vector<double> q(my_data.param_NPV);

/* --SEPARATE PATHS FOR ONE OR MORE PARAMETERS */

    if (my_data.maxf2_.nv > 1) {
   goto L10;
    }

/* --FOR NV = 1: */

    my_data.maxf2_.pdir[0] = my_data.maxf2_.g[0] / my_data.maxf2_.diag[0];
    goto L40;

/* --FOR NV > 1: */

/* --FIRST SOLVE   ULT*Q = G   FOR Q */
/* --(USE G RATHER THAN -G TO COMPENSATE FOR THE FACT THAT ULT*DIAG*ULT' */
/* --CORRESPONDS TO OPPOSITE OF ACTUAL HESSIAN MATRIX) */

L10:
    q[0] = my_data.maxf2_.g[0];
    for (int l = 1; l < my_data.maxf2_.nv; ++l) {
   q[l] = my_data.maxf2_.g[l];
   for (int ll = 0; ll < l; ++ll) {
       q[l] -= ult_ref(l, ll) * q[ll];
   }
    }

/* --THEN SOLVE   ULT'*PDIR = DIAG**(-1)*Q   FOR PDIR */

    my_data.maxf2_.pdir[my_data.maxf2_.nv-1] = q[my_data.maxf2_.nv-1] / my_data.maxf2_.diag[my_data.maxf2_.nv-1];
    for (int l = my_data.maxf2_.nv - 1; l >= 0; --l) {
   my_data.maxf2_.pdir[l] = q[l] / my_data.maxf2_.diag[l];
   for (int ll = l + 1; ll < my_data.maxf2_.nv; ++ll) {
       my_data.maxf2_.pdir[l] -= ult_ref(ll, l) * my_data.maxf2_.pdir[ll];
   }
    }

/* --INITIALIZE EXIT FLAG */

L40:
    lex = 0;

/* --COMPUTE NORMALIZED GRADIENT PDIR'G */

    my_data.maxf2_.ptg = 0.;
    for (int l = 0; l < my_data.maxf2_.nv; ++l) {
   my_data.maxf2_.ptg += my_data.maxf2_.pdir[l] * my_data.maxf2_.g[l];
    }

/* --MAKE SURE PTG > 0 (DIRECTION WILL INCREASE FUNCTION) */

/* L60: */
    if (my_data.maxf2_.ptg > 0.) {
   goto L70;
    }
    lex = 2;
    return 0;

/* --SEE IF PTG SMALL ENOUGH TO STOP */

L70:
    if (my_data.maxf2_.ptg > my_data.maxf1_.epsc2) {
   return 0;
    }
    lex = 1;
    return 0;

} /* direct_ */

#undef ult_ref


int Maxfun::lsrch_(int& ifirst, int& ifl)
{
    /* Local variables */
    double scal, tbnd, tmin, thpm;
    double diff1, diff2;
    int l;
    double fs = 0;
    double ft = 0;
    double tt = 0;
    double ft2 = 0;
    double tt2 = 0;
    double tin;
    int lex;

    bool numeric_difference = false;

/* --LINE SEARCH IN CHOSEN DIRECTION (PDIR) */

/* --RETURN FLAG IFL: */
/* --  -1:  ITERATION COMPLETE, BUT TSTEP < TIN AND < TMIN */
/* --   0:  ITERATION COMPLETE, WITH TSTEP >= TIN AND/OR TMIN */
/* --   1:  ITERATION COMPLETE; CONVERGED BY IMPLIED TEST */
/* --   2:  ITERATION COMPLETE; NEGLIGIBLE FUNCTION CHANGE */
/* --   3:  ITERATION COMPLETE; REACHED MAXIMUM # ITERATIONS */
/* --   4:  COULD NOT COMPLETE ITERATION */


/* --SET TIN (INITIAL TT), TMIN (SIMPLE INITIAL VALUES IF THIS IS FIRST */
/* --ITERATION WITH THIS NV) */

    vector<double> ths(my_data.param_NPV);
    vector<double> tht(my_data.param_NPV);

    /* Function Body */
    if (ifirst <= 0) {
   goto L10;
    }

/* --FIRST TIME */

    tmin = 0.;
    tin = .1;
    goto L30;

/* --NOT FIRST TIME */

L10:
    thpm = 1e30;
    l = -1;
    for (int i = 0; i < my_data.maxf1_.nt; ++i) {
   if (my_data.maxf2_.ist[i] > 2)
          continue;
   ++l;
   if (my_data.maxf2_.pdir[l] != 0.) {
            double d = dfn_(fabs(theta[i]), my_data.maxf1_.yota)/fabs(my_data.maxf2_.pdir[l]);
       thpm = min(thpm,d);
   }
    }
    tmin = .5 * thpm/my_data.maxf1_.epst;

    tin = 2. * my_data.maxf2_.fch/my_data.maxf2_.ptg;
    if (tin <= 0.) {
   tin = 1.;
    }
    tin = min(tin,1.);
    if (my_data.maxf2_.idif > 1) {
   tin = max(tin,1.);
    }

/* --ESTABLISH MAXIMUM TSTEP (TBND) TO SATISFY BOUNDARY RESTRICTIONS ON */
/* --INDEPENDENT PARAMETERS */

L30:
    tbnd = 1e5;
    l = -1;
    for (int i = 0; i < my_data.maxf1_.nt; ++i) {
   if (my_data.maxf2_.ist[i] > 2)
          continue;

   ++l;
   if (my_data.maxf2_.pdir[l] > 0.) {
            double d = (my_data.maxf1_.thu[i] - theta[i])/my_data.maxf2_.pdir[l];
       tbnd = min(tbnd,d);

            while(1)
            {
              if (theta[i] + tbnd * my_data.maxf2_.pdir[l] <= my_data.maxf1_.thu[i])
                break;

              tbnd = (1. - my_data.maxf1_.yota) * tbnd;
            }

   } else if (my_data.maxf2_.pdir[l] < 0.) {
       double d = (my_data.maxf1_.thl[i] - theta[i]) / my_data.maxf2_.pdir[l];
       tbnd = min(tbnd,d);

            while(1)
            {
              if (theta[i] + tbnd * my_data.maxf2_.pdir[l] >= my_data.maxf1_.thl[i])
                break;
              tbnd = (1. - my_data.maxf1_.yota) * tbnd;
            }
   }
    }

/* --MODIFY TIN IF NECESSARY */

    if (tin * (my_data.maxf1_.yota + 2.) >= tbnd) {
   tin = tbnd * (.5 - my_data.maxf1_.yota);
    }
    tt = tin;

/* --WRITE POSSIBLE TSTEP VALUES IN DETAIL FILE, IF ANY */

/* --BEGIN LINE SEARCH IN CHOSEN DIRECTION */

/* --SAVE INCOMING VALUES */

    ths = theta;
    fs = f;

/* --INITIALIZE IMPLIED BOUNDARY FLAG */

    my_data.maxf2_.impbnd = 0;

/* --COMPUTE THETA-1 */

L60:
    l = -1;
    for (int i = 0; i < my_data.maxf1_.nt; ++i) {
   if (my_data.maxf2_.ist[i] > 2)
          tht[i] = ths[i];
        else
        {
          ++l;
          tht[i] = ths[i] + tt * my_data.maxf2_.pdir[l];
        }
    }

/* --GET CORRESPONDING FUNCTIONAL VALUE F-1 */

    evaluate(tht, ft, nfe, lex);
    if (lex > 0) {
   goto L250;
    }
    if (ft <= f) {
   goto L180;
    }

/* --IMPROVING; TRY GOING FARTHER */

L90:
    tt2 = tt + tt;
    if (tt2 < tbnd) {
   goto L110;
    }

/* --CAN'T GO ANY FARTHER, SO STOP AT THETA-1 */

    my_data.maxf2_.tstep = tt;
    f = ft;
    theta = tht;
    goto L260;

/* --COMPUTE THETA-2 (FIXED PARAMETERS ALREADY ESTABLISHED; DEPENDENT ONES */
/* --WILL BE COMPUTED IN DEPAR VIA FUN) */

L110:
    l = -1;
    for (int i = 0; i < my_data.maxf1_.nt; ++i) {
   theta[i] = tht[i];
   if (my_data.maxf2_.ist[i] > 2)
          continue;
   ++l;
   tht[i] = theta[i] + tt * my_data.maxf2_.pdir[l];
    }

/* --GET CORRESPONDING FUNCTIONAL VALUE F-2 */

    evaluate(tht, ft2, nfe, lex);
    if (lex <= 0) {
   goto L130;
    }
    my_data.maxf2_.impbnd = 1;
    goto L140;
L130:
    if (ft2 > ft) {
   goto L150;
    }
    my_data.maxf2_.impbnd = 0;

/* --NO LONGER IMPROVING; STOP AT THETA-1 */

L140:
    my_data.maxf2_.tstep = tt;
    f = ft;
    goto L260;

/* --STILL IMPROVING; HOW FAST? */

L150:
    tt = tt2;
    diff1 = ft2 - ft;
    diff2 = ft - f;
    if (diff1 > diff2) {
   goto L170;
    }

/* --SLOWER; TOO SLOW? */

    if (diff1 / diff2 >= .75) {
   goto L170;
    }

/* --TOO SLOW; STOP AT THETA-2 */

    my_data.maxf2_.tstep = tt2;
    f = ft2;
    theta = tht;
    goto L260;

/* --IMPROVING FASTER OR AT LEAST FAST ENOUGH; */
/* --PROCEED WITH THETA-2 AS NEW THETA-1 */

L170:
    ft = ft2;
    goto L90;

/* --NOT IMPROVING; TRY NOT GOING SO FAR */

L180:
    my_data.maxf2_.impbnd = 0;

    // Added 2004.01.13:  When the steps are decreasing, sometimes the
    // change in thetas gets smaller than the numeric epsilon of the system.
    // At this point, the steps will keep going until the tt < tmin, but
    // there's no point in all that extra calculation, especially when tmin
    // is 0, as is the case where you've maximized in an earlier maxfun
    // operation and are finalizing using variable metric methods.

    numeric_difference = false;

    for (int i = 0; i < my_data.maxf1_.nt; ++i)
    {
        if(theta[i] != tht[i])
            numeric_difference = true;
    }

    my_data.maxf2_.impbnd = 0;
    if (tt > tmin && numeric_difference) {
        goto L190;
    }

/* --CAN'T CUT ANY FARTHER; QUIT WITHOUT COMPLETING ITERATION */

    ifl = 4;
    return 0;

/* --HALVE TT AND COMPUTE THETA-2 (FIXED PARAMETERS ALREADY ESTABLISHED; */
/* --DEPENDENT ONES WILL BE COMPUTED IN DEPAR VIA FUN) */

L190:
    tt *= .5;
    l = -1;
    for (int i = 0; i < my_data.maxf1_.nt; ++i) {
   if (my_data.maxf2_.ist[i] > 2)
          continue;
   ++l;
   tht[i] = ths[i] + tt * my_data.maxf2_.pdir[l];
    }

/* --GET CORRESPONDING FUNCTIONAL VALUE F-2 */

    evaluate(tht, ft2, nfe, lex);

    if (lex > 0) {
   goto L250;
    }
    if (ft2 <= f) {
   goto L220;
    }

/* --IMPROVED THIS TIME; STOP AT THETA-2 */

    my_data.maxf2_.tstep = tt;
    f = ft2;
    theta = tht;
    goto L260;

/* --STILL NOT IMPROVING; REDUCE TT FURTHER BY SCALING FACTOR */

L220:
    diff1 = ft + f - ft2 - ft2;
    scal = .1;
    if (diff1 < 0.) {
   scal = max(scal, (f - ft) * .5 / diff1 + 1.);
    }
    tt = scal * tt;

/* --COMPUTE NEW THETA-2 (FIXED PARAMETERS ALREADY ESTABLISHED; DEPENDENT */
/* --ONES WILL BE COMPUTED IN DEPAR VIA FUN) */

    l = -1;
    for (int i = 0; i < my_data.maxf1_.nt; ++i) {
   if (my_data.maxf2_.ist[i] > 2)
          continue;
   ++l;
   tht[i] = ths[i] + tt * my_data.maxf2_.pdir[l];
    }

/* --GET CORRESPONDING FUNCITONAL VALUE F-2 */

    evaluate(tht, ft, nfe, lex);
    if (lex > 0) {
   goto L250;
    }

/* --IF STILL NOT IMPROVED, GO BACK TO CUT TT MORE, WITH THETA-2 AS NEW */
/* --THETA-1 */

    if (ft <= f) {
   goto L180;
    }

/* --IMPROVED THIS TIME; STOP AT THETA-2 */

    my_data.maxf2_.tstep = tt;
    f = ft;
    theta = tht;
    goto L260;

/* --IF RUN INTO UNDEFINED FUNCTION BEFORE ACHIEVING ANY IMPROVEMENT, */
/* --HALVE TT AND START OVER */

L250:
    tt *= .5;
    my_data.maxf2_.impbnd = 1;
    goto L60;

/* --CHECK WHETHER TSTEP < TMIN */

L260:
    if (my_data.maxf2_.tstep >= tmin || my_data.maxf2_.tstep >= tin) {
   goto L270;
    }
    ifl = -1;
    goto L280;

/* --NORMAL TERMINATION */

L270:
    ifl = 0;

/* --FINISH UP THIS ITERATION */

L280:
    my_data.maxf2_.thpr = ths;
    my_data.maxf2_.fpr = fs;
    endit_();
    ifirst = 0;

/* --INCREMENT NUMBER OF ITERATIONS WITH THIS GRADIENT VECTOR */

    ++my_data.maxf2_.igage;

/* --WRITE OUT DETAILS OF THIS ITERATION IF DESIRED */
/* --CHECK FOR CONVERGENCE */

    vcnvch_(lex);

/* --SET FLAG (IF APPROPRIATE) AND RETURN */

    if (lex > 0) {
   ifl = lex;
    }
    return 0;

} /* lsrch_ */

int Maxfun::prepd_(int& lex)
{
    /* Local variables */
    int lexb;
    double sinit;

/* --PREPARE FOR DERIVATIVE COMPUTATION: */
/* --CHECK FOR CONVERGENCE TO BOUNDS AND WHETHER # PARAMETERS OK, */
/* --INITIALIZE INCREMENT STEPSIZE FACTORS FOR DERIVATIVE COMPUTATION */

/* --RETURNS FLAG LEX: */
/* --   0:  EVERYTHING OK */
/* --   1:  ALL INDEPENDENT PARAMETERS CONVERGED TO BOUNDS */
/* --   2:  TOO MANY PARAMETERS TO COMPUTE DERIVATIVES */


/* --INITIALIZE COUNTER OF PARAMETERS CONVERGED TO BOUNDS */

    /* Function Body */
    my_data.maxf2_.nb = 0;

/* --CHECK FOR CONVERGENCE TO BOUNDS */

    fixbnd_(lexb);
    if (lexb < 2) {
   goto L10;
    }
    lex = 1;
    return 0;

/* --CHECK FOR TOO MANY PARAMETERS */

L10:
    if (my_data.maxf2_.ni <= my_data.param_NPV) {
   goto L20;
    }
    lex = 2;
    return 0;

/* --INITIALIZE INCREMENT STEPSIZE FACTORS FOR 2ND DERIVATIVE COMPUTATION */

L20:
    if (my_data.maxf1_.ihit > 0) {
   sinit = my_data.maxf1_.epsd;
    } else {
   sinit = my_data.maxf1_.yota;
    }

    for (int i = 0; i < my_data.maxf1_.nt; ++i) {
   my_data.maxf2_.stp[i] = sinit;
    }

    lex = 0;
    return 0;

} /* prepd_ */

int Maxfun::fixbnd_(int& lex)
{
    /* Local variables */
    double thib;
    int lexb;
    double this__;
    int isti;
    double fs;
    int nbs;

/* --FIX PARAMETERS WHICH ARE CLOSE TO A BOUND */

/* --RETURNS FLAG LEX: */
/* --   0:  NOTHING NEW FIXED */
/* --   1:  FIXED AT LEAST ONE INDEPENDENT PARAMETER */
/* --   2:  ALL INDEPENDENT PARAMETERS NOW FIXED */


/* --SAVE INCOMING VALUE AND INITIALIZE EXIT FLAG */

    /* Function Body */
    nbs = my_data.maxf2_.nb;

/* --CHECK EACH VARYING INDEPENDENT PARAMETER */

    for (int i = 0; i < my_data.maxf1_.nt; ++i) {
   isti = my_data.maxf2_.ist[i];

   if (isti > 2)
          continue;

   fs = f;
   this__ = theta[i];
   thib = my_data.maxf1_.thu[i];
   bndchk_(this__, thib, lexb);
   if (lexb > 0)
       goto L10;

   thib = my_data.maxf1_.thl[i];
   bndchk_(this__, thib, lexb);
   if (lexb <= 0)
          continue;

/* --CLOSE TO BOUND; FIX */

L10:
   ++my_data.maxf2_.nb;
   theta[i] = thib;
   evaluate(theta, f, nfe, lexb);
   if (lexb > 0 || f < fs) {
           my_data.maxf2_.ist[i] = isti + 6;
           theta[i] = this__;
           f = fs;
   }
        else
          my_data.maxf2_.ist[i] = isti + 4;
    }

/* --FIX UP DEPENDENT PARAMETERS TO GO WITH FINAL THETA */
/* --(NEED NOT CHECK RETURN FLAG AS ALREADY TESTED THIS COMBINATION) */

    if (my_data.maxf2_.nd > 0) {
   depar(theta, lexb);
    }

/* --SET # PARAMETERS THAT MAY VARY */

    my_data.maxf2_.nv = my_data.maxf2_.ni - my_data.maxf2_.nb;

/* --SEE IF ANY PARAMETERS HAVE BEEN FIXED */

    if (my_data.maxf2_.nb > nbs) {
   goto L40;
    }
    lex = 0;
    return 0;

/* --SEE IF ANY LEFT UNFIXED */

L40:
    if (my_data.maxf2_.nv <= 0) {
   goto L50;
    }
    lex = 1;
    return 0;

L50:
    lex = 2;
    return 0;

} /* fixbnd_ */

int Maxfun::bndchk_(double& th, double& thb, int& lex)
{
    /* Local variables */
    double rth;


/* --TEST FOR CLOSENESS TO A BOUND */

/* --RETURNS FLAG LEX: */
/* --   0:  NOT CLOSE TO BOUND */
/* --   1:  CLOSE TO BOUND */


/* --PERFORM APPROPRIATE FORM OF TEST */

    if (fabs(th) >= my_data.maxf1_.epsd) {
   goto L10;
    }

    if (fabs(th - thb) <= my_data.maxf1_.epsd * my_data.maxf1_.epsd) {
   goto L30;
    }
    goto L20;

L10:
    rth = thb / th;
    if (rth >= 1. - my_data.maxf1_.epsd && rth <= my_data.maxf1_.epsd + 1.) {
   goto L30;
    }

/* --NOT CLOSE TO BOUND */

L20:
    lex = 0;
    return 0;

/* --CLOSE TO BOUND */

L30:
    lex = 1;
    return 0;

} /* bndchk_ */

int Maxfun::bcnvch_(int& lex)
{
    /* Local variables */
    double eloc;

/* --CHECK FOR CONVERGENCE -- IMPLIED AND REGULAR TESTS AND OTHER STUFF */
/* --APPROPRIATE FOR DIRECT SEARCH (BASIC) */

/* --RETURNS FLAG LEX: */
/* --   0:  FAILED CONVERGENCE TEST; CONTINUE ITERATING */
/* --   1:  PASSED CONVERGENCE TEST (CRITERION 1) */
/* --   2:  NEGLIGIBLE FUNCTION CHANGE (CONVERGENCE CRITERION 3) */
/* --   3:  FAILED CONVERGENCE TEST, BUT REACHED MAXIMUM # ITERATIONS */


/* --INITIALIZE INTERNAL EPSILON */

    /* Function Body */
    eloc = my_data.maxf1_.epsc1;

/* --LOOP THROUGH PARAMETERS, USING IMPLIED TEST FOR DEPENDENT ONES */

L10:
    for (int i = 0; i < my_data.maxf1_.nt; ++i) {
   switch ((int)(my_data.maxf2_.ist[i])) {
       case 1:  goto L20;
       case 2:  goto L20;
       case 3:  goto L30;
       case 4:  goto L40;
       case 5:  goto L40;
       case 6:  goto L40;
       case 7:  goto L40;
       case 8:  goto L40;
       case 9:  goto L40;
       case 10:  goto L40;
   }
L20:
   if (my_data.maxf2_.stp[i] > efn_(fabs(my_data.maxf2_.thpr[i]), eloc)) {
       goto L60;
   }
   goto L40;
L30:
   itest_(i, eloc, lex);
   if (lex <= 0) {
       goto L60;
   }
L40:
   ;
    }
/* --PASSED BASIC CONVERGENCE TEST */

    lex = 1;
    return 0;

/* --FAILED STANDARD CONVERGENCE TEST */

/* --TRY INCREASING EPSILON IF NOT ALREADY TRIED AND DIFMAX SMALL ENOUGH */

L60:
    if (my_data.maxf1_.epsc1 <= 0. || eloc > my_data.maxf1_.epsc1 || my_data.maxf2_.difmax > eloc*eloc) {
   goto L70;
    }
    eloc *= 10.;
    goto L10;

/* --COMPARE DIFMAX WITH USER-SPECIFIED EPSC3 */

L70:
    if (my_data.maxf2_.difmax >= my_data.maxf1_.epsc3) {
   goto L80;
    }
    lex = 2;
    return 0;

/* --GIVE UP ON CONVERGENCE; CHECK FOR MAXIMUM # ITERATIONS */

L80:
    if (my_data.maxf2_.it < my_data.maxf1_.maxit) {
   goto L90;
    }
    lex = 3;
    return 0;

/* --NOT THAT EITHER; CONTINUE ITERATING */

L90:
    lex = 0;
    return 0;

} /* bcnvch_ */

int Maxfun::vcnvch_(int& lex)
{
/* --CHECK FOR CONVERGENCE -- IMPLIED TEST APPROPRIATE FOR NEWTON-RAPHSON */
/* --OR VARIABLE METRIC METHODS */

/* --RETURNS FLAG LEX: */
/* --   0:  FAILED CONVERGENCE TEST; CONTINUE ITERATING */
/* --   1:  PASSED (IMPLIED) CONVERGENCE TEST (CRITERION 1) */
/* --   2:  NEGLIGIBLE FUNCTION CHANGE (CONVERGENCE CRITERION 3) */
/* --   3:  FAILED CONVERGENCE TEST, BUT REACHED MAXIMUM # ITERATIONS */


/* --DO IMPLIED TEST FOR EACH VARYING PARAMETER */

    /* Function Body */
    for (int i = 0; i < my_data.maxf1_.nt; ++i) {
   if (my_data.maxf2_.ist[i] > 3)
          continue;

   itest_(i, my_data.maxf1_.epsc1, lex);

   if (lex <= 0)
          break;
    }

/* --CONVERGED */
    if(lex > 0)
      return 0;

/* --DID NOT CONVERGE BY IMPLIED TEST */
/* --CHECK FOR NEGLIGIBLE CHANGE IN FUNCTION */

    if (my_data.maxf2_.fch >= my_data.maxf1_.epsc3)
    {
      if (my_data.maxf2_.it < my_data.maxf1_.maxit)
        return 0;
    }
    else
    {
      lex = 2;
      return 0;
    }

/* --DID NOT CONVERGE; CHECK FOR TOO MANY ITERATIONS */
/* --STOP AT MAXIMUM # ITERATIONS */

    lex = 3;
    return 0;

} /* vcnvch_ */

int Maxfun::itest_(int& i, double& eloc, int& lex)
{
    /* Local variables */
    double rth;

/* --CHECK FOR CONVERGENCE ON ONE PARAMETER -- IMPLIED TEST */

/* --RETURNS FLAG LEX: */
/* --   0:  FAILED IMPLIED CONVERGENCE TEST */
/* --   1:  PASSED IMPLIED CONVERGENCE TEST */

    /* Function Body */
    if (fabs(my_data.maxf2_.thpr[i]) > eloc) {
   goto L10;
    }
    if (fabs(my_data.maxf2_.cth[i]) > eloc * eloc) {
   goto L30;
    }
    goto L20;
L10:
    rth = theta[i] / my_data.maxf2_.thpr[i];
    if (rth < 1. - eloc || rth > eloc + 1.) {
   goto L30;
    }

/* --PASSED TEST */

L20:
    lex = 1;
    return 0;

/* --FAILED TEST */

L30:
    lex = 0;
    return 0;

} /* itest_ */

int Maxfun::endit_()
{
//    SAGE::MAXFUN::APIMaxFunction::output_iteration_end_results();

    Maxfun::iteration_ended = true;

    /* Local variables */
    double er;

/* --FINISH UP AN ITERATION */

    /* Function Body */
    my_data.maxf2_.erm = 0.;
    for (int i = 0; i < my_data.maxf1_.nt; ++i) {
   er = theta[i] - my_data.maxf2_.thpr[i];
   my_data.maxf2_.cth[i] = er;
   er = fabs(er);
   if (er > my_data.maxf2_.erm)
       my_data.maxf2_.erm = er;
    }

    my_data.maxf2_.fch = f - my_data.maxf2_.fpr;

    ++my_data.maxf2_.it;
    return 0;
} /* endit_ */


int Maxfun::augv_(int& ih)
{
    /* Local variables */
    double dthi;
    int isti;
    int j, k, l, m, n;
    double s;
    int istii, newit;
    int jj, kk, ll, mm, jm;
    double thi;
    int jmm;
    int lex;

#define u_ref(a_1,a_2) u[(a_2)*my_data.param_NPV+a_1]
#define v_ref(a_1,a_2) my_data.maxf2_.v[(a_2)*my_data.param_NPV+a_1]
#define df_ref(a_1,a_2) df[(a_2)*my_data.param_NPV+a_1]
#define av_ref(a_1,a_2) my_data.maxf2_.av[(a_2)*my_data.param_NP+a_1]
#define vn_ref(a_1,a_2) vn[(a_2)*my_data.param_NPV+a_1]


/* --AUGMENT THE VARIANCE-COVARIANCE MATRIX TO INCLUDE INDEPENDENT */
/* --PARAMETERS THAT HAVE CONVERGED TO BOUNDS AND DEPENDENT */
/* --PARAMETERS, PRINT IT, AND COMPUTE STANDARD DEVIATIONS */

/* --NOTE THAT VALUES COMPUTED HERE CORRESPOND TO THE FINAL ESTIMATES */
/* --THETA, WHILE THE REST OF THE MATRIX MAY NOT */


/* --FIRST CONSIDER INDEPENDENT PARAMETERS THAT HAVE CONVERGED TO BOUNDS
*/

    vector<int>    newv(my_data.param_NPV);
    vector<double> u(3*my_data.param_NPV);
    vector<double> prmuv(my_data.param_NPV);
    vector<double> df(my_data.param_NPV*my_data.param_NPV);
    vector<double> vn(my_data.param_NPV*my_data.param_NPV);
    vector<int>    jdp(my_data.param_NPV);
    vector<double> thm(my_data.param_NP);
    vector<double> thp(my_data.param_NP);

    int ireturn = 0; // due to JA for flagging round off error

    /* Function Body */
    if (my_data.maxf2_.nb > 0) {
   goto L20;
    }

/* --NO INDEPENDENT PARAMETERS CONVERGED TO BOUNDS (NI = NV); JUST */
/* --COPY V TO VN AND GO ON */

    for (kk = 0; kk < my_data.maxf2_.ni; ++kk)
   for (int k = 0; k < my_data.maxf2_.ni; ++k)
       vn_ref(k, kk) = v_ref(k, kk);
    goto L70;

/* --I INDEXES ALL NT PARAMETERS */
/* --K INDEXES NI INDEPENDENT PARAMETERS */
/* --L INDEXES NV ITERABLE (INDEPENDENT, NOT NEAR BOUND) PARAMETERS */

L20:
    k = -1;
    l = -1;
    for (int i = 0; i < my_data.maxf1_.nt; ++i) {
   isti = my_data.maxf2_.ist[i];

   if (isti >= 3 && isti <= 4)
            continue;

   ++k;

        if(isti > 2)
        {
          /* --PARAMETER I HAS CONVERGED TO A BOUND; INSERT ZEROS */

          for (kk = 0; kk < my_data.maxf2_.ni; ++kk) {
              vn_ref(k, kk) = 0.;
              vn_ref(kk, k) = 0.;
          }
          continue;
        }

/* --PARAMETER I IS INDEPENDENT AND VARYING */

   ++l;

/* --II, KK, LL PROVIDE SIMILAR INDEXING SCHEME FOR 2ND PARAMETER */

   kk = -1;
   ll = -1;
   for (int ii = 0; ii <= i; ++ii) {
       istii = my_data.maxf2_.ist[ii];
       if (istii >= 3 && istii <= 4)
              continue;

       ++kk;

       if (istii > 4)
              continue;

/* --PARAMETER II IS INDEPENDENT AND NOT FIXED */

       ++ll;

       vn_ref(k, kk) = v_ref(l, ll);
       vn_ref(kk, k) = v_ref(ll, l);
   }
    }

/* --NOW CONSIDER DEPENDENT PARAMETERS */

L70:
    if (my_data.maxf2_.nd > 0) {
   goto L90;
    }

/* --NO DEPENDENT PARAMETERS (NE = NI); JUST COPY VN TO AV AND GO ON TO */
/* --COMPUTE STANDARD DEVIATIONS */

    for (int jj = 0; jj < my_data.maxf2_.ne; ++jj)
   for (int j = 0; j < my_data.maxf2_.ne; ++j)
       av_ref(j, jj) = vn_ref(j, jj);
    goto L320;

/* --THERE ARE DEPENDENT PARAMETERS TO BE TAKEN CARE OF */

/* --FIRST MAJOR LOOP IS THROUGH INDEPENDENT PARAMETERS TO OBTAIN */
/* --DERIVATIVES OF ALL DEPENDENT PARAMETERS WITH RESPECT TO THEM */

/* --I INDEXES ALL NT PARAMETERS */
/* --J INDEXES ALL NE NON-FIXED (BY USER) PARAMETERS */
/* --K INDEXES NI INDEPENDENT PARAMETERS */
/* --M INDEXES ND DEPENDENT PARAMETERS */

L90:
    j = -1;
    k = -1;
    m = -1;

    for (int i = 0; i < my_data.maxf1_.nt; ++i) {
   isti = my_data.maxf2_.ist[i];
   if (isti == 4)
          continue;

   ++j;

        if(isti == 3)
        {
          /* --STORE INDEX (IN 1 TO NE SYSTEM) TO MTH DEPENDENT PARAMETER */

          ++m;
          jdp[m] = j;
          continue;
        }

/* --THIS IS AN INDEPENDENT PARAMETER */

   ++k;

   if (isti == 1) {
       goto L120;
   }

/* --FOR INDEPENDENT PARAMETERS WHICH HAVE CONVERGED TO A BOUND OR ARE NOT */
/* --INVOLVED IN FUNCTIONAL RELATIONSHIPS, SET TO ZERO ALL DERIVATIVES */
/* --WITH RESPECT TO THEM */

   for (mm = 0; mm < my_data.maxf2_.nd; ++mm)
       df_ref(mm, k) = 0.;
   goto L250;

/* --PREPARE TO COMPUTE DERIVATIVES WITH RESPECT TO PARAMETER I */

L120:
   thi = theta[i];

   for (mm = 0; mm < my_data.maxf2_.nd; ++mm) {
       prmuv[mm] = 999999.;
       newv[mm] = 1;
   }
   dthi = dfn_(fabs(thi), my_data.maxf2_.stp[i]);

   n = 0;

/* --OBTAIN NTH APPROXIMATION TO DERIVATIVE FOR EACH DEPENDENT PARAMETER */
/* --WHICH STILL NEEDS ANOTHER ITERATION */

L140:
        thp = theta;
        thm = theta;
   thp[i] = thi + dthi;
   thm[i] = thi - dthi;
   depar(thp, lex);
   if (lex > 0) {
       goto L170;
   }
   depar(thm, lex);
   if (lex > 0) {
       goto L170;
   }

/* --II INDEXES ALL NT PARAMETERS */
/* --MM INDEXES ND DEPENDENT PARAMETERS */

   mm = -1;
   for (int ii = 0; ii < my_data.maxf1_.nt; ++ii) {
       if (my_data.maxf2_.ist[ii] != 3)
              continue;

       ++mm;

       if (newv[mm] > 0)
      u_ref(mm, n) = (thp[ii] - thm[ii]) / (dthi + dthi);
   }

/* --PREPARE TO DO ANOTHER APPROXIMATION IF NECESSARY */

   if (ih <= 0) {
       goto L230;
   }
   ++n;
   if (n > 2) {
       goto L180;
   }

/* --PREPARE FOR ITERATION WITH SMALLER STEPSIZE */

L170:
   dthi *= .5;
   goto L140;

/* --FOR EACH DEPENDENT PARAMETER STILL BEING WORKED ON, FIT DERIVATIV
E */
/* --USING 3 (LATEST) APPROXIMATIONS */

/* --MM INDEXES ND DEPENDENT PARAMETERS */

L180:
   newit = 0;
   for (mm = 0; mm < my_data.maxf2_.nd; ++mm) {
       if (newv[mm] <= 0)
              continue;
             augv_fitder_(u_ref(mm, 0), u_ref(mm, 1), u_ref(mm, 2),
                    df_ref(mm, k), prmuv[mm], lex);

// Changed by JA in Aug. 2010, new subroutine reduces the
// discretization error in the first derivative
// Modified by JA, to flag round off error in the computation of first
// derivatives. Similar Feature was present in the original fortran code.

            ireturn = lex;

       if (lex - 1 < 0) {
      goto L210;
       } else if (lex - 1 == 0) {
      goto L200;
       } else {
      goto L190;
       }

/* --POSSIBLE ROUND-OFF ERROR DETECTED IN COMPUTATION OF DERIVATIVE */

L190:

/* --GO DIRECTLY HERE IF DERIVATIVE SUCCESSFULLY APPROXIMATED */

L200:
       newv[mm] = 0;
            continue;

/* --ANOTHER ITERATION NECESSARY FOR THIS DERIVATIVE */

L210:
       newit = 1;
   }

   if (newit == 0) {
       goto L250;
   }
   n = 2;
   goto L170;

/* --STORE 1ST APPROXIMATION (WHEN ONLY DOING ONE) FOR EACH DEPENDENT */
/* --PARAMETER */

L230:
   for (mm = 0; mm < my_data.maxf2_.nd; ++mm)
       df_ref(mm, k) = u_ref(mm, 0);

/* --COPY ROW, COLUMN FOR THIS INDEPENDENT PARAMETER TO OUTPUT */
/* --VARIANCE-COVARIANCE MATRIX */

/* --II INDEXES ALL NT PARAMETERS */
/* --JJ INDEXES ALL NE NON-FIXED PARAMETERS */
/* --KK INDEXES NI INDEPENDENT PARAMETERS */

L250:
   jj = -1;
   kk = -1;
   for (int ii = 0; ii <= i; ++ii) {

       istii = my_data.maxf2_.ist[ii];
       if (istii == 4)
              continue;

       ++jj;

       if (istii == 3)
              continue;

       ++kk;

       av_ref(j, jj) = vn_ref(k, kk);
       av_ref(jj, j) = vn_ref(kk, k);
   }
    }

/* --SECOND MAJOR LOOP IS THROUGH DEPENDENT PARAMETERS TO COMPUTE AND */
/* --STORE CORRESPONDING VARIANCES, COVARIANCES */

/* --M INDEXES ND DEPENDENT PARAMETERS */
/* --JM GIVES CORRESPONDING INDEX (IN 1 TO NE SYSTEM) FOR OUTPUT MATRIX */

    for (int m = 0; m < my_data.maxf2_.nd; ++m) {
   jm = jdp[m];

/* --FIRST DO COVARIANCES BETWEEN THIS DEPENDENT PARAMETER AND EACH */
/* --INDEPENDENT PARAMETER */

/* --I INDEXES ALL NT PARAMETERS */
/* --J INDEXES NE NON-FIXED PARAMETERS */
/* --K INDEXES NI INDEPENDENT PARAMETERS */

   j = -1;
   k = -1;

   for (int i = 0; i < my_data.maxf1_.nt; ++i) {
       isti = my_data.maxf2_.ist[i];
       if (isti == 4)
              continue;

       ++j;

       if (isti == 3)
              continue;

       ++k;
       s = 0.;
       for (kk = 0; kk < my_data.maxf2_.ni; ++kk) {
      s += df_ref(m, kk) * vn_ref(kk, k);
       }

       av_ref(jm, j) = s;
       av_ref(j, jm) = s;
   }

/* --THEN DO COVARIANCES BETWEEN THIS DEPENDENT PARAMETER AND OTHER */
/* --DEPENDENT PARAMETERS (INCLUDING ITS OWN VARIANCE) */

/* --MM INDEXES ND DEPENDENT PARAMETERS */
/* --JMM GIVES CORRESPONDING INDEX (IN 1 TO NE SYSTEM) FOR OUTPUT MATRIX */

   for (int mm = 0; mm <= m; ++mm) {
       jmm = jdp[mm];
       s = 0.;
       for (int kk = 0; kk < my_data.maxf2_.ni; ++kk)
      for (int k = 0; k < my_data.maxf2_.ni; ++k)
          s += df_ref(m, kk) * df_ref(mm, k) * vn_ref(kk, k);

       av_ref(jm, jmm) = s;
       av_ref(jmm, jm) = s;
   }
    }

/* --PRINT OUT VARIANCE-COVARIANCE MATRIX AND RELATED INFORMATION */

L320:

/* --COMPUTE STANDARD DEVIATIONS WHERE POSSIBLE */

    j = -1;
    for (int i = 0; i < my_data.maxf1_.nt; ++i)
    {
   if (my_data.maxf2_.ist[i] == 4)
        {
     my_data.maxf2_.stde[i] = 0.;
          continue;
        }

        ++j;
        s = av_ref(j, j);
        if (s >= 0.)
          my_data.maxf2_.stde[i] = sqrt(s);
        else
        {
          my_data.maxf2_.stde[i] = -sqrt(-s);
          my_data.maxf2_.ivfl = 2;
        }
    }
//    return 0;
    return ireturn;
} /* augv_ */

#undef vn_ref
#undef av_ref
#undef df_ref
#undef v_ref
#undef u_ref


int Maxfun::deriv1_()
{
    /* Local variables */
    double dthi, this__, thiy;
    int l;
    double fm = 0;
    double fp = 0;
    double si = 0;
    int lex;

/* --COMPUTE GRADIENT (VECTOR OF FIRST PARTIAL DERIVATIVES) AND ITS NORM
*/

/* --RETURNS FLAG IGFL: */
/* --   0:  DERIVATIVES SUCCESSFULLY COMPUTED */
/* --   1:  DERIVATIVES COULD NOT BE COMPUTED DUE TO UNDEFINED FUNCTION */


/* --INITIALIZE */

    vector<double> thy;

    /* Function Body */
    my_data.maxf2_.igage = -1;
    thy = theta;
    my_data.maxf2_.gtg = 0.;

/* --LOOP THROUGH NV INDEPENDENT VARYING PARAMETERS */

    l = -1;
    for (int i = 0; i < my_data.maxf1_.nt; ++i) {
   if (my_data.maxf2_.ist[i] > 2) continue;

   ++l;
   this__ = thy[i];
   si = my_data.maxf1_.yota;

L20:
   dthi = dfn_(fabs(this__), si);
   thiy = this__ + dthi;
   if (thiy > my_data.maxf1_.thu[i]) {
       goto L40;
   }
   thy[i] = thiy;
        evaluate(thy, fp, nfe, lex);
   if (lex <= 0) {
       goto L50;
   }

/* --RUNNING INTO IMPLIED BOUNDARY */

L30:
   my_data.maxf2_.impbnd = 1;

/* --RUNNING INTO TROUBLE; DECREASE INCREMENT */

L40:
   si *= .5;

   if (si < my_data.maxf1_.epsd * my_data.maxf1_.epsd / 8.) {
       goto L90;
   }
   goto L20;

/* --DIVIDE FORWARD, CENTRAL DIFFERENCE PATHS */

L50:
   if (my_data.maxf2_.idif > 1) {
       goto L60;
   }

/* --FORWARD DIFFERENCE */

   my_data.maxf2_.g[l] = (fp - f) / dthi;
   goto L70;

/* --CENTRAL DIFFERENCE */

L60:
   thiy = this__ - dthi;
   if (thiy < my_data.maxf1_.thl[i]) {
       goto L40;
   }
   thy[i] = thiy;
   evaluate(thy, fm, nfe, lex);
   if (lex > 0) {
       goto L30;
   }
   my_data.maxf2_.g[l] = (fp - fm) / (dthi + dthi);

L70:
/* Computing 2nd power */
   my_data.maxf2_.gtg += my_data.maxf2_.g[l] * my_data.maxf2_.g[l];
   thy[i] = this__;
    }

/* --FINISH NORM OF GRADIENT */

    my_data.maxf2_.gtg = sqrt(my_data.maxf2_.gtg);

/* --INDICATE HAVE G CORRESPONDING TO CURRENT THETA */

    my_data.maxf2_.igfl = 0;
    my_data.maxf2_.igage = 0;

    return 0;

/* --ERROR EXIT; PARAMETER I STUCK */

L90:
    my_data.maxf2_.igfl = 1;
    return 0;

} /* deriv1_ */

int Maxfun::deriv2_(int& ih, int& lex)
{
    /* Local variables */
    double athi, dthi, thim;
    int lexf;
    double thip, thii;
    int isti;
    int l, n;
    double dthii;
    int istii;
    double prmuh;
    double fm = 0;
    int ll, nl, nk;
    double si;
    int nm;
    double fp = 0;
    double e2d8;
    int nbd[4];
    double fmm = 0;
    double fmp = 0;
    double thi = 0;
    double sii = 0;
    double fpp = 0;
    double fpm = 0;

    double max_si, min_si, max_dthi, old_dthi;
    double hfunc;

#define h_ref(a_1,a_2)    my_data.maxf2_.h[(a_2)*my_data.param_NPV+a_1]
#define hn_ref(a_1,a_2,a_3) hn[((a_3)*my_data.param_NPV+(a_2))*my_data.param_NPV+a_1]
#define lrh_ref(a_1,a_2)    lrh[(a_2)*my_data.param_NPV+a_1]

/* --COMPUTE 2ND PARTIAL DERIVATIVES OF THE FUNCTION */

/* --RETURNS FLAG LEX: */
/* --   0:  NO PROBLEM WITH H */
/* --   1:  ROUND-OFF ERROR IN H */
/* --   2:  COULD NOT COMPUTE H */
/* --   3:  SOME ROWS COULD NOT BE COMPUTED */

/* --INDICATE TYPE OF APPROXIMATION USED FOR SECOND DERIVATIVES */

    vector<double> hn(3*my_data.param_NPV*my_data.param_NPV, 0.0);
    vector<int> lrh(my_data.param_NPV*my_data.param_NPV, 0);
    vector<double> thy(my_data.param_NP);
    vector<double> dthis(my_data.param_NP, numeric_limits<double>::quiet_NaN());

    /* Function Body */
/* --INITIALIZE */

    lex = 0;
    e2d8 = my_data.maxf1_.epsd * my_data.maxf1_.epsd / 8.;
    thy = theta;

    for (int ll = 0; ll < my_data.maxf2_.nv; ++ll)
   for (int l = 0; l < my_data.maxf2_.nv; ++l)
       lrh_ref(l, ll) = 0;

/* --LOOP THROUGH INDEPENDENT, VARYING PARAMETERS */

    l = -1;
    for (int i = 0; i < my_data.maxf1_.nt; ++i) {
   isti = my_data.maxf2_.ist[i];

   if (isti > 2)
            continue;

   ++l;
   thi = thy[i];
   athi = fabs(thi);
   si = my_data.maxf2_.stp[i];
   dthi = dfn_(athi, si);

        // Set initial minimum stepsize

        min_si = 0.0;

        // Set initial maximum stepsize.  This is defined such that it is
        // the minimum of 1.0 and the value for which thi +/- dthi * max_si / si
        // is always within the bounds.

        max_si = 1.0;

        thip = thi + dthi * max_si / si;

        if(thip > my_data.maxf1_.thu[i])
          max_si = (my_data.maxf1_.thu[i] - thi) * si / dthi;

        thim = thi - dthi * max_si / si;

        if(thim < my_data.maxf1_.thl[i])
          max_si = (thi - my_data.maxf1_.thl[i]) * si / dthi;

/* --CHECK WHETHER NEIGHBORING VALUES ARE OK WITH RESPECT TO BOUNDS,
*/
/* --DEPENDENT PARAMETERS */

L40:

   thip = thi + dthi;
   thim = thi - dthi;

   if (thip <= my_data.maxf1_.thu[i] && thim >= my_data.maxf1_.thl[i]) {
       goto L50;
   }

        // If we are greater than our bounds, we determine the maximum delta
        // theta possible and calculate a new delta theta that is 1 -
        // epsilon D from it.  si (step size) and max_si are changed
        // accordingly.

        max_dthi = min( my_data.maxf1_.thu[i] - thi,
                       -my_data.maxf1_.thl[i] + thi);

        old_dthi = dthi;

        dthi = (1.0 - my_data.maxf1_.epsd) * max_dthi;

        max_si = si * max_dthi / old_dthi;

        si = si * dthi / old_dthi;

   thip = thi + dthi;
   thim = thi - dthi;

   if (si < min_si + e2d8) {
       goto L84; // struck bound
   }

L50:
        // We no longer do this since we want to find any bounds that might
        // occur during evaluation.

   //if (isti == 2) {
   //    goto L80;
   //}
L60:
        // Calculate the function values for +/- delta theta.  If these
        // fail, there must be an implied bound.

   thy[i] = thip;
   evaluate(thy, fp, nfe, lexf);
   if (lexf > 0) {
       goto L70;
   }
   thy[i] = thim;
   evaluate(thy, fm, nfe, lexf);
   if (lexf <= 0) {
       goto L80;
   }

L70:
        // If we fail dependent parameter check, we're too far from the
        // current estimates.  In this case, lower the maximum stepsize to
        // the current stepsize, put the stepsize halfway to the lower bound
        // and calculate a new delta theta based on this change.

   my_data.maxf2_.impbnd = 1;

        max_si = si;

   si = (si + min_si) * 0.5;

        // As always, check for being too close to the minimum bound.

   if (si < min_si + e2d8) {
       goto L84; // struck bound
   }

   dthi = dthi * si / max_si;

   thim = thi - dthi;
   thip = thi + dthi;
   goto L60;

L80:
        // Test for a valid function return value.

        // Calculate the function value.

   hfunc = abs(fp - f - f + fm) / (dthi * dthi);

        if(hfunc >= 1.e-8) {
          goto L85;
        }

        // If the value is too small, we're not far enough away from the
        // current parameter values.  Set the minimum step size to the
        // current value, increase the step size halfway to the calculated
        // maximum and the delta theta proportionally.

        min_si = si;

        si = min(si*2.0, (si + max_si) * 0.5);

        // Check for being too close to the maximum bound.

        if (si > max_si - e2d8) {
            goto L84;  // struck bound
        }

        dthi = dthi * si / min_si;

   thim = thi - dthi;
   thip = thi + dthi;

   goto L40;

L84:
        // We have struck a bound.  This means that this value cannot have a
        // valid second derivative or standard error.  We make sure of this
        // by setting the delta theta to QNAN and then skipping over the
        // later calculations.

        h_ref(l,l) = numeric_limits<double>::quiet_NaN();
        dthis[i] = numeric_limits<double>::quiet_NaN();

        lex = 3;

        continue;

/* --RESTORE THY(I) AND SAVE CURRENT VALUE OF STEPSIZE FACTOR */

L85:
        // restore the thy(i) and store stepsize and delta theta

   my_data.maxf2_.stp[i] = si;
        dthis[i] = dthi;

   thy[i] = thi;

/* --COMPUTE 2ND PARTIAL DERIVATIVES FOR (I,II) AND (II,I) PAIRS */
/* --WHERE II < I */

   if (i == 0) {
       goto L360;
   }

   ll = -1;
   for (int ii = 0; ii < i; ++ii) {
       istii = my_data.maxf2_.ist[ii];
       if (istii > 2)
              continue;

       ++ll;

       thii = thy[ii];
       sii = my_data.maxf2_.stp[ii];
       dthii = dthis[ii];
       si = my_data.maxf2_.stp[i];
       dthi = dthis[i];

            // Skip this calculation if either quantity has a NaN as a
            // delta.  This indicates a parameter unable to calculate.

            if(SAGE::isnan(dthi) || SAGE::isnan(dthii)) continue;

/* --IF BOTH INDEPENDENT PARAMETERS ARE INVOLVED IN FUNCTIONAL */
/* --RELATIONSHIPS, CHECK ALL COMBINATIONS OF NEIGHBORING PARAMETE
R PAIRS */
/* --WITH RESPECT TO DEPENDENT PARAMETERS */

       if (isti == 2 || istii == 2) {
      goto L260;
       }
L90:
       nl = 1;
L100:
       for (nk = 0; nk < 4; ++nk)
      nbd[nk] = 0;

       thy[i] = thi + dthi;
       thy[ii] = thii + dthii;
       nm = 0;
       nk = 1;
L120:
       depar(thy, lexf);
       if (lexf <= 0) {
      goto L130;
       }
       my_data.maxf2_.impbnd = 1;
       nbd[nk - 1] = 1;
       ++nm;
L130:
       switch ((int)nk) {
      case 1:  goto L140;
      case 2:  goto L150;
      case 3:  goto L160;
      case 4:  goto L170;
       }
L140:
       thy[ii] = thii - dthii;
       nk = 2;
       goto L120;
L150:
       thy[i] = thi - dthi;
       thy[ii] = thii + dthii;
       nk = 3;
       goto L120;
L160:
       thy[ii] = thii - dthii;
       nk = 4;
       goto L120;
L170:
       if (nm <= 0) {
      goto L260;
       }
       switch ((int)nm) {
      case 1:  goto L180;
      case 2:  goto L220;
      case 3:  goto L230;
      case 4:  goto L230;
       }
L180:
       switch ((int)nl) {
      case 1:  goto L190;
      case 2:  goto L200;
      case 3:  goto L210;
      case 4:  goto L190;
       }
L190:
       si *= .5;
       if (si < e2d8) {
      goto L490;
       }
       dthi *= .5;
       nl = 2;
       goto L100;
L200:
       sii *= .5;
       if (sii < e2d8) {
      goto L480;
       }
       dthii *= .5;
       si += si;
       dthi += dthi;
       nl = 3;
       goto L100;
L210:
       si *= .5;
       if (si < e2d8) {
      goto L490;
       }
       dthi *= .5;
       nl = 4;
       goto L100;
L220:
       if (nbd[0] * nbd[1] > 0 || nbd[2] * nbd[3] > 0) {
      goto L250;
       }
       if (nbd[0] * nbd[2] > 0 || nbd[1] * nbd[3] > 0) {
      goto L240;
       }
L230:
       si *= .5;
       if (si < e2d8) {
      goto L490;
       }
       dthi *= .5;
L240:
       sii *= .5;
       if (sii < e2d8) {
      goto L480;
       }
       dthii *= .5;
       goto L90;
L250:
       si *= .5;
       if (si < e2d8) {
      goto L490;
       }
       dthi *= .5;
       goto L90;

/* --READY TO GO AHEAD WITH DERIVATIVE ESTIMATION */

L260:
       prmuh = 999999.;
       n = 0;

/* --COMPUTE NTH APPROXIMATION TO 2ND PARTIAL DERIVATIVE */

L270:
       thy[i] = thi + dthi;
       thy[ii] = thii + dthii;
       evaluate(thy, fpp, nfe, lexf);
       if (lexf > 0) {
      goto L280;
       }
       thy[i] = thi - dthi;

       evaluate(thy, fmp, nfe, lexf);
       if (lexf > 0) {
      goto L280;
       }
       thy[ii] = thii - dthii;
       evaluate(thy, fmm, nfe, lexf);
       if (lexf > 0) {
      goto L280;
       }
       thy[i] = thi + dthi;
       evaluate(thy, fpm, nfe, lexf);
       if (lexf <= 0) {
      goto L290;
       }
L280:
       my_data.maxf2_.impbnd = 1;
       if (n > 0) {
      goto L470;
       }
       si *= .5;
       if (si < e2d8) {
      goto L470;
       }
       sii *= .5;
       if (sii < e2d8) {
      goto L470;
       }
       dthi *= .5;
       dthii *= .5;
       goto L270;

L290:
       hn_ref(l, ll, n) = (fpp - fmp - fpm + fmm) / (dthi * 4. * dthii);

/* --STORE FIRST APPROXIMATION IF ONLY DOING ONE */

       if (ih > 0) {
      goto L300;
       }
       h_ref(l, ll) = hn_ref(l, ll, 0);
       goto L340;

/* --PREPARE TO DO ANOTHER APPROXIMATION IF APPROPRIATE */

L300:
       ++n;
       if (n <= 2) {
      goto L320;
       }

/* --FIT DERIVATIVES FROM 3 APPROXIMATIONS */

       fitder_(hn_ref(l, ll, 0), hn_ref(l, ll, 1), hn_ref(l, ll, 2),
                    h_ref(l, ll), prmuh, lexf);

/* --CHECK RESULTS OF 2ND DERIVATIVE ESTIMATION */

       if (lexf - 1 < 0) {
      goto L310;
       } else if (lexf - 1 == 0) {
      goto L340;
       } else {
      goto L330;
       }

/* --NEED TO DO ANOTHER ITERATION FOR 2ND DERIVATIVE */

L310:
       n = 2;

L320:
       si *= .5;
       dthi *= .5;
       sii *= .5;
       dthii *= .5;
       goto L270;

/* --ITERATION FINISHED BUT THERE IS ROUND-OFF ERROR IN COMPUTING */
/* --2ND DERIVATIVE */

L330:
       lrh_ref(l, ll) = 1;
       lex = 1;

/* --COPY 2ND DERIVATIVE TO ELEMENT (LL,L) OF SYMMETRIC MATRIX */

L340:
       h_ref(ll, l) = h_ref(l, ll);

/* --FINISH UP LOOP */

/* --RESTORE ORIGINAL VALUES OF ITH, IITH PARAMETERS */

       thy[i] = thi;
       thy[ii] = thii;
   }

/* --TAKE CARE OF (I,I) PAIR */

L360:
   si = my_data.maxf2_.stp[i];
   dthi = dthis[i];

        // Skip this calculation if delta theta is a NaN.  This indicates a
        // parameter unable to calculate.

        if(SAGE::isnan(dthi))
        {
          // Store a QNAN in the hessian matrix as well, so that we
          // know this is bad.

          h_ref(l,l) = numeric_limits<double>::quiet_NaN();

          continue;
        }

   prmuh = 999999.;
   n = 0;

/* --START HERE TO WORK ON NTH APPROXIMATION */

L370:
   thy[i] = thi + dthi;
   evaluate(thy, fp, nfe, lexf);
   if (lexf > 0) {
       goto L380;
   }
   thy[i] = thi - dthi;
   evaluate(thy, fm, nfe, lexf);
   if (lexf <= 0) {
       goto L390;
   }
L380:
   my_data.maxf2_.impbnd = 1;
   if (n > 0) {
       goto L490;
   }
   si *= .5;
   if (si < e2d8) {
       goto L490;
   }
   dthi *= .5;
   goto L370;

L390:
   hn_ref(l, l, n) = (fp - f - f + fm) / (dthi * dthi);

/* --STORE FIRST APPROXIMATION IF ONLY DOING ONE */

   if (ih > 0) {
       goto L400;
   }
   h_ref(l, l) = hn_ref(l, l, 0);
   goto L440;

/* --PREPARE TO DO ANOTHER APPROXIMATION IF NECESSARY */

L400:
   ++n;
   if (n <= 2) {
       goto L420;
   }

/* --FIT DERIVATIVES FROM 3 APPROXIMATIONS */

   fitder_(hn_ref(l, l, 0), hn_ref(l, l, 1), hn_ref(l, l, 2),
      h_ref(l, l), prmuh, lexf);

/* --CHECK RESULTS OF 2ND DERIVATIVE ESTIMATION */

   if (lexf - 1 < 0) {
       goto L410;
   } else if (lexf - 1 == 0) {
       goto L440;
   } else {
       goto L430;
   }

/* --NEED TO DO ANOTHER ITERATION FOR 2ND DERIVATIVE */

L410:
   n = 2;

L420:
   si *= .5;
   dthi *= .5;
   goto L370;

/* --ITERATION FINISHED BUT THERE IS ROUND-OFF ERROR IN COMPUTING */
/* --2ND DERIVATIVE */

L430:
   lrh_ref(l, l) = 1;
   lex = 1;

/* --RESTORE ORIGINAL VALUE OF I'TH PARAMETER */

L440:
   thy[i] = thi;
    }

/* --PRINT H AND WARNINGS ABOUT ROUND-OFF ERRORS IN DETAIL FILE, IF ANY */

    return 0;

/* --ERROR EXIT; PARAMETERS II AND I STUCK */

L470:
    goto L490;

/* --ERROR EXIT; PARAMETER II STUCK */

L480:

/* --ERROR EXIT; PARAMETER I STUCK */

L490:
    lex = 2;
    return 0;

} /* deriv2_ */

#undef lrh_ref
#undef hn_ref
#undef h_ref

int Maxfun::fitder_(double& d1, double& d2, double& d3,
                    double& dd, double& prmu, int& lex)
{
    /* Local variables */
    double rden, teps, r__, rmu;

/* --FITDER COMPUTES AN IMPROVED VALUE FOR A NUMERICAL DERIVATIVE */
/* --USING 3 ESTIMATES WHICH HAVE BEEN COMPUTED USING 3 DIFFERENT */
/* --STEP SIZES. */

/* --RETURNS FLAG LEX: */
/* --   0:  NEED ANOTHER ITERATION */
/* --   1:  SUCCESSFUL DERIVATIVE ESTIMATION */
/* --   2:  ROUND-OFF ERROR */


/* --SEE IF ALL ESTIMATES < 10**-8 */

    if (fabs(d1) < 1e-8 && fabs(d2) < 1e-8 && fabs(d3) < 1e-8) {
      dd = 0.;
      lex = 1;
      return 0;
    }

/* --COMPUTE RATIO R */

//    rden = d3 + d3 - d2;

    rden = (2*d2 - d1);
    if (rden != 0.) {
   goto L20;
    }
    r__ = 0.;
    goto L30;
L20:
//    r__ = (d1 - d2 * 6. + d3 * 8.) / 3.;

    r__ = ((d1 - d2*6.0 + d3*8.0)/3.0)/rden;

//    r__ = (d2 + d2 - d1) / rden;

/* --IF R "CLOSE" TO 1, FINISHED */

L30:
//    teps = my_data.maxf1_.epsd * 10.;

    teps = 0.5; // modification by JA, reduce tolerances


    if (r__ >= 1. - teps && r__ <= teps + 1.) {
   goto L50;
    }

/* --UNLESS R GETTING FARTHER FROM 1, DO ANOTHER ITERATION */

    rmu = fabs(r__ - 1.);
    if (rmu > prmu) {
   goto L40;
    }
    prmu = rmu;
    d1 = d2;
    d2 = d3;
    lex = 0;
    return 0;

L40:
    dd = d1;
    lex = 2;
    return 0;

L50:
    dd = (d1 - d2 * 6. + d3 * 8.) / 3.;
    lex = 1;
    return 0;

} /* fitder_ */

double Maxfun::dfn_(double ath, double sf)
{
    /* System generated locals */
    double ret_val;

/* --COMPUTES INCREMENT ("DELTA THETA") AS A FUNCTION OF ATH AND STEPSIZE  */
/* --FACTOR SF, WITH PARAMETER TAU */


    if (ath > .5) {
   ret_val = ath;
    } else {
   ret_val = .5;
    }
    ret_val = sf * ret_val;
    return ret_val;

} /* dfn_ */

double Maxfun::dfninv_(double ath, double delth)
{
    /* System generated locals */
    double ret_val;

/* --COMPUTES INCREMENT STEPSIZE FACTOR AS A FUNCTION OF ATH AND */
/* --INCREMENT DELTH, WITH PARAMETER TAU */

    if (ath > .5) {
   ret_val = delth / ath;
    } else {
   ret_val = delth / .5;
    }
    return ret_val;

} /* dfninv_ */

double Maxfun::efn_(double ath, double eloc)
{
    /* System generated locals */
    double ret_val;

/* --COMPUTES MINIMAL VALUE OF INCREMENT STEPSIZE FACTOR (DELTA) AS A */
/* --FUNCTION OF ATH, WITH PARAMETERS ELOC */

    ret_val = eloc;
    if (ath > .5) {
   return ret_val;
    }
    ret_val /= .5;
    if (ath > eloc) {
   ret_val *= ath;
    } else {
   ret_val *= eloc;
    }
    return ret_val;

} /* efn_ */

int Maxfun::mmult_(const vector<double>& a, int mra, int ma, int namb,
                   const vector<double>& b, int mrb, int nb,
                         vector<double>& c, int mrc)
{
    /* System generated locals */
    int a_dim1, b_dim1, c_dim1;

#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]
#define c_ref(a_1,a_2) c[(a_2)*c_dim1 + a_1]


/* --MULTIPLY MATRIX A (MA BY NAMB) BY MATRIX B (NAMB BY NB) AND STORE */
/* --THE RESULT IN MATRIX C (MA BY NB) */

/* --NOTE THAT INDICES I, J, K ARE LOCAL TO THIS GENERAL SUBROUTINE AND */
/* --MAY NOT CORRESPOND TO INDICES OF THE SAME NAMES IN OTHER PROGRAM */
/* --UNITS */

    /* Parameter adjustments */
    a_dim1 = mra;
    b_dim1 = mrb;
    c_dim1 = mrc;

    /* Function Body */
    for (int i = 0; i < ma; ++i) {
   for (int j = 0; j < nb; ++j) {
       double z = 0.;
       for (int k = 0; k < namb; ++k)
      z += a_ref(i, k) * b_ref(k, j);

       c_ref(i, j) = z;
   }
    }
    return 0;

} /* mmult_ */

#undef c_ref
#undef b_ref
#undef a_ref

// this method due to JA  for score testing

bool Maxfun::score_deriv1()
{
    /* Local variables */
    double dthi, this__, thiy;
    double fm = 0;
    double fp = 0;
    double si = 0;
    int lex;
    bool forward_deriv_done;
    bool backward_boundary_problem;

/* --INITIALIZE */
    vector<double> thy;

    /* Function Body */
    thy = theta; // replace by final estimate


    for (int i = 0; i < my_data.maxf1_.nt; ++i) {
        unsigned countfp = 0;
        unsigned countfm = 0;
   if (my_data.need_score[i] == 0) continue;

        forward_deriv_done = false;
        backward_boundary_problem = false;

        this__ = thy[i];
   si = my_data.maxf1_.yota;
        dthi = dfn_(fabs(this__), si);

L20:
        forward_deriv_done = false;
   dthi = dfn_(fabs(this__), si);
   thiy = this__ + dthi; // increment for the forward difference
   if (thiy > my_data.maxf1_.thu[i]) { // outside the boundary
       goto L40;
   }

   thy[i] = thiy;
   evaluate(thy, fp, nfe, lex);
        countfp++;
        if (lex <= 0) { // evaluated the forward difference with no trouble
            forward_deriv_done = true;
            goto L60;
   }

L40: // RUNNING INTO TROUBLE; DECREASE INCREMENT
   si *= .5; // shrink the step-size by half

   if (si < my_data.maxf1_.epsd * my_data.maxf1_.epsd / 8.) {
            if (!(forward_deriv_done)) goto L60; // never successfully computed a forward derivative
            if ( (forward_deriv_done) && (backward_boundary_problem) ) goto L80;
            if ( !(forward_deriv_done) && (backward_boundary_problem) ) goto L90; // BAD
   }
   goto L20; // try forward derivative again

L60: // setup for backward derivative
   thiy = this__ - dthi; // subtract for the backward difference
   if (thiy < my_data.maxf1_.thl[i]) { // outside lower boundary
        backward_boundary_problem = true;
       goto L40; // spilling outside lower boundary
   }

   thy[i] = thiy;
   evaluate(thy, fm, nfe, lex);
        countfm++;
        if (lex > 0) {
        backward_boundary_problem = true;
       goto L40;
   }
        backward_boundary_problem = false;
L80:
        if ( (forward_deriv_done) && !(backward_boundary_problem) )
         {
           my_data.score_vec[i] = (fp - fm)/(dthi + dthi);
           my_data.score_deriv_conv[i] = 0; // central derivative no problem
           my_data.score_stepsize[i] = 8.0*fabs(dthi);
         }

        if ( !(forward_deriv_done) && !(backward_boundary_problem) )
         {
           my_data.score_vec[i] = (f - fm)/dthi;
           my_data.score_deriv_conv[i] = 1; // very close to upper boundary, backward derivative
           my_data.score_stepsize[i] = fabs(4.0*dthi);
         }

        if ( (forward_deriv_done) && (backward_boundary_problem) )
         {
          my_data.score_vec[i] = (fp - f)/dthi;
          my_data.score_deriv_conv[i] = 2; // very close to lower boundary, forward derivative
          my_data.score_stepsize[i] = fabs(4.0*dthi);
         }
    } // end of loop over parameters

    return true;

/* --ERROR EXIT; PARAMETER I STUCK */

L90:
    return false;

} /* score_deriv1_ */


bool Maxfun::score_deriv2()
{
    /* Local variables */
    double dthi, thim;
    int lexf;
    double thip, thii;
    int n;
    double dthii;
    double prmuh;
    double fm = 0;
    int nl, nk;
    double si;
    int nm;
    double fp = 0;
    double e2d8;
    int nbd[4];
    int currenti, currentii;
    double fmm = 0;
    double fmp = 0;
    double thi = 0;
    double sii = 0;
    double fpp = 0;
    double fpm = 0;

    double max_si, min_si;

#define h_ref(a_1,a_2)    my_data.score_vec2[(a_2)*my_data.param_NPV+a_1]
#define hn_ref(a_1,a_2,a_3) hn[((a_3)*my_data.param_NPV+(a_2))*my_data.param_NPV+a_1]
#define lrh_ref(a_1,a_2)    lrh[(a_2)*my_data.param_NPV+a_1]

/* --COMPUTE 2ND PARTIAL DERIVATIVES OF THE FUNCTION */

/* --RETURNS FLAG LEX: */
/* --   0:  NO PROBLEM WITH H */
/* --   1:  ROUND-OFF ERROR IN H */
/* --   2:  COULD NOT COMPUTE H */
/* --   3:  SOME ROWS COULD NOT BE COMPUTED */

/* --INDICATE TYPE OF APPROXIMATION USED FOR SECOND DERIVATIVES */

    vector<double> hn(3*my_data.param_NPV*my_data.param_NPV, 0.0);
    vector<int> lrh(my_data.param_NPV*my_data.param_NPV, 0);
    vector<double> thy(my_data.param_NP);
    vector<double> dthis(my_data.param_NP, numeric_limits<double>::quiet_NaN());

    /* Function Body */
/* --INITIALIZE */

    int lex = 0;
    e2d8 = my_data.maxf1_.epsd *my_data.maxf1_.epsd /8.;
    thy = theta;
    int ih = 3; // to force interative fit

    for (int i = 0; i < my_data.maxf1_.nt; ++i) {
   if (my_data.need_score[i] == 0) continue; // not a parameter for which we need the score test
         currenti = my_data.score_deriv_conv[i];

        thi = thy[i];
        si = my_data.score_stepsize[i];
// si will be rescaled in order to be consistent with the first derivative
        dthi = dfn_(fabs(thi),si);

        // Set initial minimum stepsize

        min_si = 0.0;

        // Set initial maximum stepsize.

        max_si = 1.0;

/* --CHECK WHETHER NEIGHBORING VALUES ARE OK WITH RESPECT TO BOUNDS,
*/
/* --DEPENDENT PARAMETERS */

L40:
        if(currenti == 0)
         { // central first derivative
      thip = thi + dthi;
      thim = thi - dthi;
     if ((thip <= my_data.maxf1_.thu[i]) && (thim >= my_data.maxf1_.thl[i]) ) goto L60;
         }

        if(currenti == 1)
         { // backward first derivative close to upper boundary
      thip = thi - 2.0*dthi;
      thim = thi - dthi;
      if (thip >= my_data.maxf1_.thl[i]) goto L60;
         }

        if(currenti == 2)
         { // forward first derivative close to lower boundary
      thip = thi + 2.0*dthi;
      thim = thi + dthi;
      if (thip <= my_data.maxf1_.thu[i]) goto L60;
         }

// do a simple rescaling; we have step size information

        dthi = dthi*0.5;
        si = si*0.5;
        if (si < e2d8) { goto L84;} // struck implied bound
        goto L40;

L60:
        // Calculate the function values for thim thip.  If these
        // fail, there must be an implied bound.

   thy[i] = thip;
   evaluate(thy, fp, nfe, lexf);
   if (lexf > 0) {
       goto L70;
   }
   thy[i] = thim;
   evaluate(thy, fm, nfe, lexf);
   if (lexf <= 0) {
       goto L85;
   }

L70:
        si = si*0.5;

   if (si < e2d8) { goto L84;} // struck implied bound

   dthi = dthi*0.5;

        if(currenti == 0)
         { // central first derivative
      thip = thi + dthi;
      thim = thi - dthi;
         }

        if(currenti == 1)
         { // backward first derivative
      thip = thi - 2.0*dthi;
      thim = thi - dthi;
         }

        if(currenti == 2)
         { // forward first derivative
      thip = thi + 2.0*dthi;
      thim = thi + dthi;
         }

   goto L60;

L84:
        return false;

/* --RESTORE THY(I) AND SAVE CURRENT VALUE OF STEPSIZE FACTOR */

L85:

        dthis[i] = dthi;
   thy[i] = thi;

/* --COMPUTE 2ND PARTIAL DERIVATIVES FOR (I,II) AND (II,I) PAIRS */
/* --WHERE II < I */

   if (i == 0) { // diagonal element only
       goto L360;
   }

        for (int ii = 0; ii < i; ++ii) {
    if (my_data.need_score[ii] == 0) continue; // not a parameter for which we need the score test

          currentii = my_data.score_deriv_conv[ii];

// already known as ii < i
            thii = thy[ii];
            sii = my_data.score_stepsize[ii];
            dthii = dthis[ii];

            si = my_data.score_stepsize[i];
       dthi = dthis[i];

/*
this part of the code in deriv2 was done
only for at least one parameter in a functional relationship
it may not be strictly essential but is done here regardless, for safety's sake.
It is identical to what was done in deriv2, except for the extended cominations
of parameter arguments.
*/

L90:
       nl = 1;
L100:
       for (nk = 0; nk < 4; ++nk)
      nbd[nk] = 0;

             // initialize i and ii parameters

             // param i central derivative
             if( (currenti == 0) || (currenti == 2)) thy[i] = thi +dthi;

             // param i backward derivative
             if(currenti == 1) thy[i] = thi - dthi;

             // param ii central derivative
        if ( (currentii == 0) || (currentii == 2) ) thy[ii] = thii +dthii;

             // param ii backward derivative
             if (currentii == 1) thy[ii] = thii - dthii;

       nm = 0;
       nk = 1;
L120:
       depar(thy, lexf);
       if (lexf <= 0) {
      goto L130;
       }
       nbd[nk - 1] = 1;
       ++nm;
L130:
       switch ((int)nk) {
      case 1:  goto L140;
      case 2:  goto L150;
      case 3:  goto L160;
      case 4:  goto L170;
       }

L140:       // update ii parameters, i parameters held fixed
       if (currentii == 0) thy[ii] = thii -dthii;
       if (currentii == 1) thy[ii] = thii -2.0*dthii;
       if (currentii == 2) thy[ii] = thii +2.0*dthii;
       nk = 2;
       goto L120;

L150:       // completed first cycle , update i parameters , keep ii fixed
       if(currenti == 0) thy[i] = thi -dthi;
       if(currenti == 1) thy[i] = thi -2.0*dthi;
       if(currenti == 2) thy[i] = thi +2.0*dthi;

       nk = 3;
       goto L120;

L160:       // second cycle, i fixed, ii restored to starting set
       if( (currentii == 0) || (currentii == 2) ) thy[ii] = thii +dthii;
       if(currentii == 1) thy[ii] = thii -dthii;

       nk = 4;
       goto L120;
L170:
       if (nm <= 0) {
      goto L260;
       }
       switch ((int)nm) {
      case 1:  goto L180;
      case 2:  goto L220;
      case 3:  goto L230;
      case 4:  goto L230;
       }
L180:
       switch ((int)nl) {
      case 1:  goto L190;
      case 2:  goto L200;
      case 3:  goto L210;
      case 4:  goto L190;
       }
L190:
       si *= .5;
       if (si < e2d8) {
      goto L490;
       }
       dthi *= .5;
       nl = 2;
       goto L100;
L200:
       sii *= .5;
       if (sii < e2d8) {
      goto L480;
       }
       dthii *= .5;
       si += si;
       dthi += dthi;
       nl = 3;
       goto L100;
L210:
       si *= .5;
       if (si < e2d8) {
      goto L490;
       }
       dthi *= .5;
       nl = 4;
       goto L100;
L220:
       if (nbd[0] * nbd[1] > 0 || nbd[2] * nbd[3] > 0) {
      goto L250;
       }
       if (nbd[0] * nbd[2] > 0 || nbd[1] * nbd[3] > 0) {
      goto L240;
       }
L230:
       si *= .5;
       if (si < e2d8) {
      goto L490;
       }
       dthi *= .5;
L240:
       sii *= .5;
       if (sii < e2d8) {
      goto L480;
       }
       dthii *= .5;
       goto L90;
L250:
       si *= .5;
       if (si < e2d8) {
      goto L490;
       }
       dthi *= .5;
       goto L90;

/*
--READY TO GO AHEAD WITH DERIVATIVE ESTIMATION i
same code as in deriv_2 except for the extended set
of parameter combinations. The logic and the coding
for shrinking increments for dealing with boundaries
is unchanged.
*/

L260:
       prmuh = 999999.;
       n = 0;

/* --COMPUTE NTH APPROXIMATION TO 2ND PARTIAL DERIVATIVE */

L270:
            // initialise values for first parameter
       if ((currenti == 0) || (currenti == 2) ) thy[i] = thi +dthi;
       if (currenti == 1) thy[i] = thi -dthi;

            // initialise values for second parameter
            if ((currentii == 0) || (currentii == 2) ) thy[ii] = thii +dthii;
       if (currentii == 1) thy[ii] = thii -dthii;

       evaluate(thy, fpp, nfe, lexf);
       if (lexf > 0) {
      goto L280;
       }

            // first parameter second values, second parameter first values
       if (currenti == 0) thy[i] = thi -dthi;
       if (currenti == 1) thy[i] = thi -2.0*dthi;
       if (currenti == 2) thy[i] = thi +2.0*dthi;

       evaluate(thy, fmp, nfe, lexf);
       if (lexf > 0) {
      goto L280;
       }

            // second parameter second values, first parameter second values
            if(currentii == 0) thy[ii] = thii -dthii;
       if(currentii == 1) thy[ii] = thii -2.0*dthii;
       if(currentii == 2) thy[ii] = thii +2.0*dthii;

       evaluate(thy, fmm, nfe, lexf);
       if (lexf > 0) {
      goto L280;
       }

            // second parameter second values and first parameter first values
       if ((currenti == 0) || (currenti == 2) ) thy[i] = thi +dthi;
       if (currenti == 1) thy[i] = thi -dthi;

       evaluate(thy, fpm, nfe, lexf);
       if (lexf <= 0) {
      goto L290;
       }
L280:
            if (n > 0) {
      goto L470;
       }
       si *= .5;
       if (si < e2d8) {
      goto L470;
       }
       sii *= .5;
       if (sii < e2d8) {
      goto L470;
       }
       dthi *= .5;
       dthii *= .5;
       goto L270;

L290:
            if (currenti == 0){
             if (currentii == 0) hn_ref(i,ii,n) = (fpp - fmp - fpm + fmm) / (dthi * 4. * dthii); // from original code
             if (currentii == 1) hn_ref(i,ii,n) = (fpp - fmp +fmm - fpm)/(2.0*dthi*dthii);
             if (currentii == 2) hn_ref(i,ii,n) = (-fpp +fmp - fmm +fpm)/(2.0*dthi*dthii);
            }

            if (currenti == 1){
             if (currentii == 0) hn_ref(i,ii,n) = (fpp -fmp  +fmm - fpm)/(2.0*dthi*dthii);
             if (currentii == 1) hn_ref(i,ii,n) = (fpp -fmp  +fmm -fpm)/(dthi*dthii);
             if (currentii == 2) hn_ref(i,ii,n) = (-fpp +fmp -fmm +fpm)/(dthi*dthii);
            }

            if (currenti == 2){
             if (currentii == 0) hn_ref(i,ii,n) = (-fpp +fmp  -fmm +fpm)/(2.0*dthi*dthii);
             if (currentii == 1) hn_ref(i,ii,n) = (-fpp +fmp - fmm +fpm)/(dthi*dthii);
             if (currentii == 2) hn_ref(i,ii,n) = (fpp -fmp +fmm -fpm)/(dthi*dthii);
            }


/* --STORE FIRST APPROXIMATION IF ONLY DOING ONE */

       if (ih > 0) {
      goto L300;
       }
       h_ref(i, ii) = hn_ref(i, ii, 0);
       goto L340;

/* --PREPARE TO DO ANOTHER APPROXIMATION IF APPROPRIATE */

L300:
       ++n;
       if (n <= 2) {
      goto L320;
       }

/* --FIT DERIVATIVES FROM 3 APPROXIMATIONS */

       fitder_(hn_ref(i, ii, 0), hn_ref(i, ii, 1), hn_ref(i, ii, 2),
                    h_ref(i, ii), prmuh, lexf);


/* --CHECK RESULTS OF 2ND DERIVATIVE ESTIMATION */

       if (lexf - 1 < 0) {
      goto L310;
       } else if (lexf - 1 == 0) {
      goto L340;
       } else {
      goto L330;
       }

/* --NEED TO DO ANOTHER ITERATION FOR 2ND DERIVATIVE */

L310:
       n = 2;

L320:
       si *= .5;
       dthi *= .5;
       sii *= .5;
       dthii *= .5;
       goto L270;

/* --ITERATION FINISHED BUT THERE IS ROUND-OFF ERROR IN COMPUTING */
/* --2ND DERIVATIVE */

L330:
       lrh_ref(i, ii) = 1;
       lex = 1;

/* --COPY 2ND DERIVATIVE TO ELEMENT (LL,L) OF SYMMETRIC MATRIX */

L340:
       h_ref(ii, i) = h_ref(i, ii);

/* --FINISH UP LOOP */

/* --RESTORE ORIGINAL VALUES OF ITH, IITH PARAMETERS */

       thy[i] = thi;
       thy[ii] = thii;
   }

/* --TAKE CARE OF (I,I) PAIR */

L360:
   si = my_data.score_stepsize[i];
   dthi = dthis[i];

        prmuh = 999999.;
   n = 0;

/* --START HERE TO WORK ON NTH APPROXIMATION */

L370:
   if ( (currenti == 0) || (currenti == 2) ) thy[i] = thi +dthi;
   if (currenti == 1) thy[i] = thi -dthi;
   evaluate(thy, fp, nfe, lexf);
   if (lexf > 0) {
       goto L380;
   }

   if (currenti == 0) thy[i] = thi -dthi;
   if (currenti == 1) thy[i] = thi -2.0*dthi;
   if (currenti == 2) thy[i] = thi +2.0*dthi;

   evaluate(thy, fm, nfe, lexf);
   if (lexf <= 0) {
       goto L390;
   }
L380:
   if (n > 0) {
       goto L490;
   }
   si *= .5;
   if (si < e2d8) {
       goto L490;
   }
   dthi *= .5;
   goto L370;

L390:
   if (currenti == 0) hn_ref(i, i, n) = (fp - f - f + fm) / (dthi * dthi); // from original code
   if ( (currenti == 1) || (currenti == 2) ) hn_ref(i, i, n) = (f - 2.0*fp  + fm) / (dthi * dthi);

/* --STORE FIRST APPROXIMATION IF ONLY DOING ONE */

   if (ih > 0) {
       goto L400;
   }
   h_ref(i, i) = hn_ref(i, i, 0);
        goto L440;

/* --PREPARE TO DO ANOTHER APPROXIMATION IF NECESSARY */

L400:
   ++n;
   if (n <= 2) {
       goto L420;
   }

/* --FIT DERIVATIVES FROM 3 APPROXIMATIONS */

   fitder_(hn_ref(i, i, 0), hn_ref(i, i, 1), hn_ref(i, i, 2),
      h_ref(i, i), prmuh, lexf);


/* --CHECK RESULTS OF 2ND DERIVATIVE ESTIMATION */

   if (lexf - 1 < 0) {
       goto L410;
   } else if (lexf - 1 == 0) {
       goto L440;
   } else {
       goto L430;
   }

/* --NEED TO DO ANOTHER ITERATION FOR 2ND DERIVATIVE */

L410:
   n = 2;

L420:
   si *= .5;
   dthi *= .5;
   goto L370;

/* --ITERATION FINISHED BUT THERE IS ROUND-OFF ERROR IN COMPUTING */
/* --2ND DERIVATIVE */

L430:
   lrh_ref(i, i) = 1;
   lex = 1;

/* --RESTORE ORIGINAL VALUE OF I'TH PARAMETER */

L440:
   thy[i] = thi;
    }

/* --PRINT H AND WARNINGS ABOUT ROUND-OFF ERRORS IN DETAIL FILE, IF ANY */

    return true;

/* --ERROR EXIT; PARAMETERS II AND I STUCK */

L470:
    goto L490;

/* --ERROR EXIT; PARAMETER II STUCK */

L480:
    goto L490;

/* --ERROR EXIT; PARAMETER I STUCK */

L490:
    lex = 2;
    return false;

} /* score_deriv2_ */

#undef lrh_ref
#undef hn_ref
#undef h_ref

double Maxfun::calc_score()
{
    /* Local variables */
    int lex;

#define h_ref(a_1,a_2) my_data.maxf2_.score_vec2[(a_2)*my_data.param_NPV + a_1]
#define hn_ref(a_1,a_2) hn[(a_2)*my_data.param_NPV + a_1]

    vector<double> hn(my_data.param_NPV*my_data.param_NPV);
    unsigned count = 0;

    for (unsigned i = 0; i != my_data.score_vec.size() ; i++){
      if( !(isnan(my_data.score_vec[i])) ) count++;
    }

    vector<double> newscorevec;
    vector<double> newscorevec2;
    vector<double> inverse(count*count);
    newscorevec.clear();
    newscorevec2.clear();

    for (unsigned i = 0; i != my_data.score_vec.size() ; i++){
      if( !(isnan(my_data.score_vec[i])) ) newscorevec.push_back(my_data.score_vec[i]);
    }

   for(unsigned i = 0; i != my_data.score_vec2.size(); i++){
      if( !(isnan(my_data.score_vec2[i])) ) newscorevec2.push_back(-1.0*my_data.score_vec2[i]);
   }

//  should put check to ensure that newscorevec2 has dimensionality count^{2}

    bool getinverse = score_mxnvrt_(newscorevec2,count, count, inverse, count, lex);

    if (!(getinverse)) return numeric_limits<double>::quiet_NaN();

    double score = 0;

    for(unsigned i = 0; i != count; i++){
     for(unsigned j = 0; j != count; j++){
       score += newscorevec[i]*inverse[i*count + j]*newscorevec[j];
      }
    }


    return score;

} /* calc_score */
#undef hn_ref
#undef h_ref

bool Maxfun::score_mxnvrt_(vector<double>& a, int mra, int m, vector<double>& b,
            int mrb, int& lex)
{
    /* System generated locals */
    int a_dim1, b_dim1;

    /* Local variables */
    double p;


#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]


/* --INVERT POSITIVE DEFINITE MATRIX A (M BY M) AND PLACE INVERSE MATRIX */
/* --IN B (M BY M) */

/* --GAUSS-JORDAN ALGORITHM IS PERFORMED WITHOUT PERMUTATION OF ROWS: */
/* --A ZERO PIVOT ELEMENT AT ANY STEP CAUSES ERROR EXIT */

/* --RETURNS FLAG LEX: */
/* --   0:  NO ERROR */
/* --   1:  MATRIX CANNOT BE INVERTED BY THIS SUBROUTINE */
/* --   2:  INVALID ARGUMENT: ORDER OF MATRIX <= 0 */


/* --FIRST TAKE CARE OF SCALAR CASE */

    /* Parameter adjustments */
    a_dim1 = mra;
    b_dim1 = mrb;

    /* Function Body */
    if (m - 1 < 0) {
   goto L10;
    } else if (m - 1 == 0) {
   goto L20;
    } else {
   goto L30;
    }

/* --INVALID ORDER OF MATRIX */

L10:
     return false;

/* --SCALAR CASE:  M = 1 */

L20:
    if ( a_ref(0, 0) == 0.) {
   goto L100;
    }
    b_ref(0, 0) = 1. / a_ref(0, 0);
    goto L90;

/* --USUAL CASE:  M > 1 */

/* --INITIALIZE */

L30:
    for (int j = 0; j < m; ++j)
   for (int i = 0; i < m; ++i)
       b_ref(i, j) = a_ref(i, j);

/* --LOOP THROUGH PIVOT ELEMENTS */

    for (int k = 0; k < m; ++k) {

/* --GET PIVOT ELEMENT */

   p = b_ref(k, k);

/* --CHECK FOR ZERO PIVOT ELEMENT */

   if (p == 0.) {
       goto L100;
   }

/* --INVERT PIVOT ELEMENT */

   p = 1. / p;

/* --PROCESS K'TH COLUMN (EXCEPT K'TH ROW ELEMENT) */

   for (int i = 0; i < m; ++i)
       if (i != k)
      b_ref(i, k) = b_ref(i, k) * p;

/* --LOOP THROUGH OTHER COLUMNS */

   for (int j = 0; j < m; ++j) {
       if (j == k)
              continue;

/* --PROCESS J'TH COLUMN (EXCEPT K'TH ROW ELEMENT) */

       for (int i = 0; i < m; ++i)
      if (i != k)
          b_ref(i, j) = b_ref(i,j) - b_ref(i, k) * b_ref(k, j);
   }

/* --NOW TAKE CARE OF K'TH ROW */

   b_ref(k, k) = p;

   p = -p;
   for (int j = 0; j < m; ++j)
       if (j != k)
      b_ref(k, j) = p * b_ref(k, j);
    }

/* --NORMAL RETURN */

L90:
    return true;

/* --ERROR EXIT IF MATRIX CANNOT BE INVERTED */

L100:

    return false;

} /* score_mxnvrt_ */

// New subroutine for computing second derivatives in
// variance covariance matrices. Invoked by vcmx.
// Introduced by JA in bits and pieces in 2010.

int Maxfun::vcmx_deriv2_(int& ih, int& lex)
{
    /* Local variables */
    double athi, dthi, thim;
    int lexf;
    double thip, thii;
    int isti;
    int l, n;
    double dthii;
    int istii;
    double prmuh;
    double fm = 0;
    int ll, nl, nk;
    double si;
    int nm;
    double fp = 0;
    double e2d8;
    int nbd[4];
    double fmm = 0;
    double fmp = 0;
    double thi = 0;
    double sii = 0;
    double fpp = 0;
    double fpm = 0;

    double max_si, min_si, max_dthi, old_dthi;
    double hfunc;

#define h_ref(a_1,a_2)    my_data.maxf2_.h[(a_2)*my_data.param_NPV+a_1]
#define hn_ref(a_1,a_2,a_3) hn[((a_3)*my_data.param_NPV+(a_2))*my_data.param_NPV+a_1]
#define lrh_ref(a_1,a_2)    lrh[(a_2)*my_data.param_NPV+a_1]

/* --COMPUTE 2ND PARTIAL DERIVATIVES OF THE FUNCTION */

/* --RETURNS FLAG LEX: */
/* --   0:  NO PROBLEM WITH H */
/* --   1:  ROUND-OFF ERROR IN H */
/* --   2:  COULD NOT COMPUTE H */
/* --   3:  SOME ROWS COULD NOT BE COMPUTED */

/* --INDICATE TYPE OF APPROXIMATION USED FOR SECOND DERIVATIVES */

    vector<double> hn(4*my_data.param_NPV*my_data.param_NPV, 0.0);
    vector<int>    lrh(my_data.param_NPV*my_data.param_NPV, 0);
    vector<double> thy(my_data.param_NP);
    vector<double> dthis(my_data.param_NP, numeric_limits<double>::quiet_NaN());

    /* Function Body */
/* --INITIALIZE */

    lex = 0;
    e2d8 = my_data.maxf1_.epsd * my_data.maxf1_.epsd / 8.;
    thy = theta;

    for (int ll = 0; ll < my_data.maxf2_.nv; ++ll)
      for (int l = 0; l < my_data.maxf2_.nv; ++l)
        lrh_ref(l, ll) = 0;

/* --LOOP THROUGH INDEPENDENT, VARYING PARAMETERS */

    l = -1;
    for (int i = 0; i < my_data.maxf1_.nt; ++i)
    {
      isti = my_data.maxf2_.ist[i];

      if (isti > 2)
        continue;

      ++l;
      thi = thy[i];
      athi = fabs(thi);
      si = my_data.maxf2_.stp[i];
      dthi = dfn_(athi, si);

      // Set initial minimum stepsize

      min_si = 0.0;

      // Set initial maximum stepsize.  This is defined such that it is
      // the minimum of 1.0 and the value for which thi +/- dthi * max_si / si
      // is always within the bounds.

      max_si = 1.0;

      thip = thi + dthi * max_si / si;

      if(thip > my_data.maxf1_.thu[i])
          max_si = (my_data.maxf1_.thu[i] - thi) * si / dthi;

      thim = thi - dthi * max_si / si;

      if(thim < my_data.maxf1_.thl[i])
          max_si = (thi - my_data.maxf1_.thl[i]) * si / dthi;

/* --CHECK WHETHER NEIGHBORING VALUES ARE OK WITH RESPECT TO BOUNDS,
*/
/* --DEPENDENT PARAMETERS */

L40:

     thip = thi + dthi;
     thim = thi - dthi;

     if (thip <= my_data.maxf1_.thu[i] && thim >= my_data.maxf1_.thl[i])
       goto L50;

     // If we are greater than our bounds, we determine the maximum delta
     // theta possible and calculate a new delta theta that is 1 -
     // epsilon D from it.  si (step size) and max_si are changed
     // accordingly.

     max_dthi = min( my_data.maxf1_.thu[i] - thi,
                    -my_data.maxf1_.thl[i] + thi);

     old_dthi = dthi;

     dthi = (1.0 - my_data.maxf1_.epsd) * max_dthi;

     max_si = si * max_dthi / old_dthi;

     si = si * dthi / old_dthi;

     thip = thi + dthi;
     thim = thi - dthi;

     if (si < min_si + e2d8)
       goto L84; // struck bound

L50:
     // We no longer do this since we want to find any bounds that might
     // occur during evaluation.

     //if (isti == 2)
     //    goto L80;

L60:
     // Calculate the function values for +/- delta theta.  If these
     // fail, there must be an implied bound.

     thy[i] = thip;
     evaluate(thy, fp, nfe, lexf);
     if (lexf > 0)
       goto L70;

     thy[i] = thim;
     evaluate(thy, fm, nfe, lexf);
     if (lexf <= 0)
       goto L80;

L70:
     // If we fail dependent parameter check, we're too far from the
     // current estimates.  In this case, lower the maximum stepsize to
     // the current stepsize, put the stepsize halfway to the lower bound
     // and calculate a new delta theta based on this change.

     my_data.maxf2_.impbnd = 1;

     max_si = si;

     si = (si + min_si) * 0.5;

     // As always, check for being too close to the minimum bound.

     if (si < min_si + e2d8)
       goto L84; // struck bound

     dthi = dthi * si / max_si;

     thim = thi - dthi;
     thip = thi + dthi;
     goto L60;

L80:
     // Test for a valid function return value.

     // Calculate the function value.

     hfunc = abs(fp - f - f + fm) / (dthi * dthi);

     // New addition by JA as requested by RCE May 2010
     // We have a more refined check on the suitability of the choice of dthi

        if(hfunc >= 1.e-8) {

          if (!((fp < f) && (f > fm)) ){
            dthi = 2.0*dthi;
            goto L40;
           }

          double fm2 = 0;
          double fm4 = 0;
          double fp2 = 0;
          double fp4 = 0;
          double fm8 = 0;
          double fp8 = 0;



          double thim2 = thi -dthi/2;
          thy[i] = thim2;
     evaluate(thy,fm2,nfe,lexf);
          if (lexf < 0){
           thy[i] = thy[i] + dthi/2;
           dthi = 2.0*dthi;
           goto L40;
           }

          double thip2 = thi +dthi/2;
          thy[i] = thip2;
     evaluate(thy,fp2,nfe,lexf);
          if (lexf < 0) {
            thy[i] = thy[i] - dthi/2;
            dthi= 2.0*dthi;
            goto L40;
           }

          if (!((fp2 < f) && (f > fm2)) ){
            dthi = 2.0*dthi;
            goto L40;
           }

          double thim4 = thi -dthi/4;
     thy[i] = thim4;
     evaluate(thy,fm4,nfe,lexf);
          if (lexf < 0) {
            thy[i] = thy[i] + dthi/4;
            dthi= 2.0*dthi;
            goto L40;
            }

          double thip4 = thi +dthi/4;
          thy[i] = thip4;
     evaluate(thy,fp4,nfe,lexf);
          if (lexf < 0) {
            thy[i] = thy[i] -dthi/4;
            dthi= 2.0*dthi;
            goto L40;
           }

          if (!((fp4 < f) && (f > fm4)) ){
            dthi = 2.0*dthi;
            goto L40;
           }

          double thim8 = thi -dthi/8;
     thy[i] = thim8;
     evaluate(thy,fm8,nfe,lexf);
          if (lexf < 0) {
            thy[i] = thy[i] + dthi/8;
            dthi= 2.0*dthi;
            goto L40;
            }

          double thip8 = thi +dthi/8;
          thy[i] = thip8;
     evaluate(thy,fp8,nfe,lexf);
          if (lexf < 0) {
            thy[i] = thy[i] -dthi/8;
            dthi= 2.0*dthi;
            goto L40;
           }

          if (!((fp8 < f) && (f > fm8)) ){
            dthi = 2.0*dthi;
            goto L40;
           }
          goto L85;
        }

        // If the value is too small, we're not far enough away from the
        // current parameter values.  Set the minimum step size to the
        // current value, increase the step size halfway to the calculated
        // maximum and the delta theta proportionally.

        min_si = si;

        si = min(si*2.0, (si + max_si) * 0.5);

        // Check for being too close to the maximum bound.

        if (si > max_si - e2d8) {
            goto L84;  // struck bound
        }

        dthi = dthi * si / min_si;

     thim = thi - dthi;
     thip = thi + dthi;

     goto L40;

L84:
     // We have struck a bound.  This means that this value cannot have a
     // valid second derivative or standard error.  We make sure of this
     // by setting the delta theta to QNAN and then skipping over the
     // later calculations.

     h_ref(l,l) = numeric_limits<double>::quiet_NaN();
     dthis[i] = numeric_limits<double>::quiet_NaN();
     my_data.flat_dir[i] = my_data.label(i);

     lex = 3;

     continue;

/* --RESTORE THY(I) AND SAVE CURRENT VALUE OF STEPSIZE FACTOR */

L85:
     // restore the thy(i) and store stepsize and delta theta

     my_data.maxf2_.stp[i] = si;
     dthis[i] = dthi;

     thy[i] = thi;

/* --COMPUTE 2ND PARTIAL DERIVATIVES FOR (I,II) AND (II,I) PAIRS */
/* --WHERE II < I */

     if (i == 0)
       goto L360;

     ll = -1;
     for (int ii = 0; ii < i; ++ii)
     {
       istii = my_data.maxf2_.ist[ii];
       if (istii > 2)
         continue;

       ++ll;

       thii = thy[ii];
       sii = my_data.maxf2_.stp[ii];
       dthii = dthis[ii];
       si = my_data.maxf2_.stp[i];
       dthi = dthis[i];

       // Skip this calculation if either quantity has a NaN as a
       // delta.  This indicates a parameter unable to calculate.

       if(SAGE::isnan(dthi) || SAGE::isnan(dthii)) continue;

/* --IF BOTH INDEPENDENT PARAMETERS ARE INVOLVED IN FUNCTIONAL */
/* --RELATIONSHIPS, CHECK ALL COMBINATIONS OF NEIGHBORING PARAMETE
R PAIRS */
/* --WITH RESPECT TO DEPENDENT PARAMETERS */

       if (isti == 2 || istii == 2)
         goto L260;

L90:
       nl = 1;

L100:
       for (nk = 0; nk < 4; ++nk)
         nbd[nk] = 0;

       thy[i] = thi + dthi;
       thy[ii] = thii + dthii;
       nm = 0;
       nk = 1;

L120:
       depar(thy, lexf);
       if (lexf <= 0)
         goto L130;

       my_data.maxf2_.impbnd = 1;
       nbd[nk - 1] = 1;
       ++nm;
L130:
       switch ((int)nk)
       {
         case 1:  goto L140;
         case 2:  goto L150;
         case 3:  goto L160;
         case 4:  goto L170;
       }
L140:
       thy[ii] = thii - dthii;
       nk = 2;
       goto L120;
L150:
       thy[i] = thi - dthi;
       thy[ii] = thii + dthii;
       nk = 3;
       goto L120;
L160:
       thy[ii] = thii - dthii;
       nk = 4;
       goto L120;
L170:
       if (nm <= 0)
         goto L260;

       switch ((int)nm)
       {
         case 1:  goto L180;
         case 2:  goto L220;
         case 3:  goto L230;
         case 4:  goto L230;
       }
L180:
       switch ((int)nl)
       {
         case 1:  goto L190;
         case 2:  goto L200;
         case 3:  goto L210;
         case 4:  goto L190;
       }
L190:
       si *= .5;
       if (si < e2d8)
         goto L490;

       dthi *= .5;
       nl = 2;
       goto L100;
L200:
       sii *= .5;
       if (sii < e2d8)
         goto L480;

       dthii *= .5;
       si += si;
       dthi += dthi;
       nl = 3;
       goto L100;
L210:
       si *= .5;
       if (si < e2d8)
         goto L490;

       dthi *= .5;
       nl = 4;
       goto L100;
L220:
       if (nbd[0] * nbd[1] > 0 || nbd[2] * nbd[3] > 0)
         goto L250;

       if (nbd[0] * nbd[2] > 0 || nbd[1] * nbd[3] > 0)
         goto L240;

L230:
       si *= .5;
       if (si < e2d8)
         goto L490;

       dthi *= .5;
L240:
       sii *= .5;
       if (sii < e2d8)
         goto L480;

       dthii *= .5;
       goto L90;
L250:
       si *= .5;
       if (si < e2d8)
         goto L490;

       dthi *= .5;
       goto L90;

/* --READY TO GO AHEAD WITH DERIVATIVE ESTIMATION */

L260:

       prmuh = 999999.;
       n = 0;

/* --COMPUTE NTH APPROXIMATION TO 2ND PARTIAL DERIVATIVE */

L270:
       thy[i] = thi + dthi;
       thy[ii] = thii + dthii;
       evaluate(thy, fpp, nfe, lexf);
       if (lexf > 0) {
         goto L280;
       }
       thy[i] = thi -dthi;

       evaluate(thy, fmp, nfe, lexf);
       if (lexf > 0) {
         goto L280;
       }
       thy[ii] = thii - dthii;
       evaluate(thy, fmm, nfe, lexf);
       if (lexf > 0) {
         goto L280;
       }
       thy[i] = thi + dthi;
       evaluate(thy, fpm, nfe, lexf);
       if (lexf <= 0) {
         goto L290;
       }
L280:
       my_data.maxf2_.impbnd = 1;
       if (n > 0) {
         goto L470;
       }
       si *= .5;
       if (si < e2d8) {
         goto L470;
       }
       sii *= .5;
       if (sii < e2d8) {
         goto L470;
       }
       dthi *= .5;
       dthii *= .5;
       goto L270;

L290:
       hn_ref(l, ll, n) = (fpp - fmp - fpm + fmm) / (dthi * 4. * dthii);

/* --STORE FIRST APPROXIMATION IF ONLY DOING ONE */

       if (ih > 0) {
         goto L300;
       }
       h_ref(l, ll) = hn_ref(l, ll, 0);
       goto L340;

/* --PREPARE TO DO ANOTHER APPROXIMATION IF APPROPRIATE */

L300:
       ++n;
       if (n <= 3) {
         goto L320;
       }

/* --FIT DERIVATIVES FROM 3 APPROXIMATIONS */

            vcmx_fitder_(hn_ref(l, ll, 0), hn_ref(l, ll, 1), hn_ref(l, ll, 2), hn_ref(l,ll,3),
                    h_ref(l, ll), prmuh, lexf);

/* --CHECK RESULTS OF 2ND DERIVATIVE ESTIMATION */

       if (lexf - 1 < 0) {
         goto L310;
       } else if (lexf - 1 == 0) {
         goto L340;
       } else {
         goto L330;
       }

/* --NEED TO DO ANOTHER ITERATION FOR 2ND DERIVATIVE */

L310:
       n = 3;

L320:
       si *= .5;
       dthi *= .5;
       sii *= .5;
       dthii *= .5;
       if ( (si < e2d8) || (sii < e2d8) ) lexf = 1; // fix by JA
       if ( (si < e2d8) || (sii < e2d8) ) goto L340; // fix by JA
       goto L270;

/* --ITERATION FINISHED BUT THERE IS ROUND-OFF ERROR IN COMPUTING */
/* --2ND DERIVATIVE */

L330:
       lrh_ref(l, ll) = 1;
       lex = 1;

/* --COPY 2ND DERIVATIVE TO ELEMENT (LL,L) OF SYMMETRIC MATRIX */

L340:
       h_ref(ll, l) = h_ref(l, ll);

/* --FINISH UP LOOP */

/* --RESTORE ORIGINAL VALUES OF ITH, IITH PARAMETERS */

       thy[i] = thi;
       thy[ii] = thii;
     }

/* --TAKE CARE OF (I,I) PAIR */

L360:
   si = my_data.maxf2_.stp[i];
   dthi = dthis[i];

   // Skip this calculation if delta theta is a NaN.  This indicates a
   // parameter unable to calculate.

   if(SAGE::isnan(dthi))
   {
     // Store a QNAN in the hessian matrix as well, so that we
     // know this is bad.

     h_ref(l,l) = numeric_limits<double>::quiet_NaN();

     continue;
   }

   prmuh = 999999.;
   n = 0;

/* --START HERE TO WORK ON NTH APPROXIMATION */

L370:
   thy[i] = thi + dthi;
   evaluate(thy, fp, nfe, lexf);
   if (lexf > 0) {
       goto L380;
   }
   thy[i] = thi - dthi;
   evaluate(thy, fm, nfe, lexf);
   if (lexf <= 0) {
       goto L390;
   }
L380:
   my_data.maxf2_.impbnd = 1;
   if (n > 0) {
       goto L490;
   }
   si *= .5;
   if (si < e2d8) {
       goto L490;
   }
   dthi *= .5;
   goto L370;

L390:
   hn_ref(l, l, n) = (fp - f - f + fm) / (dthi * dthi);

/* --STORE FIRST APPROXIMATION IF ONLY DOING ONE */

   if (ih > 0) {
       goto L400;
   }
   h_ref(l, l) = hn_ref(l, l, 0);
   goto L440;

/* --PREPARE TO DO ANOTHER APPROXIMATION IF NECESSARY */

L400:
   ++n;
   if (n <= 3) {
       goto L420;
   }

/* --FIT DERIVATIVES FROM 3 APPROXIMATIONS */

    vcmx_fitder_(hn_ref(l, l, 0), hn_ref(l, l, 1), hn_ref(l, l, 2),hn_ref(l,l,3),
    h_ref(l, l), prmuh, lexf);

/* --CHECK RESULTS OF 2ND DERIVATIVE ESTIMATION */

   if (lexf - 1 < 0) {
       goto L410;
   } else if (lexf - 1 == 0) {
       goto L440;
   } else {
       goto L430;
   }

/* --NEED TO DO ANOTHER ITERATION FOR 2ND DERIVATIVE */

L410:
   n = 3;

L420:
   si *= .5;
        if (si < e2d8) lexf = 1; // temporary fix by JA
        if (si < e2d8) goto L440; // temporary fix by JA
        dthi *= .5;
   goto L370;

/* --ITERATION FINISHED BUT THERE IS ROUND-OFF ERROR IN COMPUTING */
/* --2ND DERIVATIVE */

L430:
   lrh_ref(l, l) = 1;
   lex = 1;

/* --RESTORE ORIGINAL VALUE OF I'TH PARAMETER */

L440:
   thy[i] = thi;
    }

/* --PRINT H AND WARNINGS ABOUT ROUND-OFF ERRORS IN DETAIL FILE, IF ANY */

    return 0;

/* --ERROR EXIT; PARAMETERS II AND I STUCK */

L470:
    goto L490;

/* --ERROR EXIT; PARAMETER II STUCK */

L480:

/* --ERROR EXIT; PARAMETER I STUCK */

L490:
    lex = 2;
    return 0;

} /* vcmx_deriv2_ */


#undef lrh_ref
#undef hn_ref
#undef h_ref

// Introduced by JA for use in augv
// Lowers the discretization error in computing the first derivative

int Maxfun::augv_fitder_(double& d1, double& d2, double& d3, double& dd, double& prmu, int& lex)
{
    /* Local variables */
    double rden, teps, r__, rmu;

/* --FITDER COMPUTES AN IMPROVED VALUE FOR A NUMERICAL DERIVATIVE */
/* --USING 3 ESTIMATES WHICH HAVE BEEN COMPUTED USING 3 DIFFERENT */
/* --STEP SIZES. */

/* --RETURNS FLAG LEX: */
/* --   0:  NEED ANOTHER ITERATION */
/* --   1:  SUCCESSFUL DERIVATIVE ESTIMATION */
/* --   2:  ROUND-OFF ERROR */


/* --SEE IF ALL ESTIMATES < 10**-8 */

    if (fabs(d1) < 1e-8 && fabs(d2) < 1e-8 && fabs(d3) < 1e-8) {
      dd = 0.;
      lex = 1;
      return 0;
    }

/* --COMPUTE RATIO R */

    rden = (4*d2/3.0 - d1/3.0);
    if (rden != 0.) {
   goto L20;
    }
    r__ = 0.;
    goto L30;
L20:

    r__ = (d1/45.0 - d2*4.0/9.0 + d3*64.0/45.0)/rden;


/* --IF R "CLOSE" TO 1, FINISHED */

L30:
//    teps = my_data.maxf1_.epsd * 10.;

    teps = 0.5; // modification by JA, reduce tolerances


    if (r__ >= 1. - teps && r__ <= teps + 1.) {
   goto L50;
    }

/* --UNLESS R GETTING FARTHER FROM 1, DO ANOTHER ITERATION */

    rmu = fabs(r__ - 1.);
    if (rmu > prmu) {
   goto L40;
    }
    prmu = rmu;
    d1 = d2;
    d2 = d3;
    lex = 0;
    return 0;

L40:
    dd = d1;
    lex = 2;
    return 0;

L50:
    dd = (d1/45.0 - d2*4.0/9.0 + d3*64.0/45.0);
    lex = 1;
    return 0;
}

// Introduced by JA on August 2010 on instructions of RCE
// To be used when computing variance covariance matrices
//
int Maxfun::vcmx_fitder_(double& d1, double& d2, double& d3, double& d4, double& dd, double& prmu, int& lex)
{
//  /* Local variables */
    double rden, teps, r__, rmu;

/* --FITDER COMPUTES AN IMPROVED VALUE FOR A NUMERICAL DERIVATIVE */
/* --USING 4 ESTIMATES WHICH HAVE BEEN COMPUTED USING 4 DIFFERENT */
/* --STEP SIZES. */

/* --RETURNS FLAG LEX: */
/* --   0:  NEED ANOTHER ITERATION */
/* --   1:  SUCCESSFUL DERIVATIVE ESTIMATION */
/* --   2:  ROUND-OFF ERROR */


/* --SEE IF ALL ESTIMATES < 10**-8 */

    if ((fabs(d1) < 1e-8) && (fabs(d2) < 1e-8) && (fabs(d3) < 1e-8) && (fabs(d4) < 1e-8)) {
      dd = 0.;
      lex = 1;
      return 0;
    }

/* --COMPUTE RATIO R */

    rden = (d1 -6.0*d2 +8.0*d3);
    if (rden != 0.) {
   goto L20;
    }
    r__ = 0.;
    goto L30;
L20:

    r__ = (d2 - d3*6.0 + d4*8.0)/rden;


/* --IF R "CLOSE" TO 1, FINISHED */

L30:

    teps = 0.5; // parameter can be tweaked


    if (r__ >= 1. - teps && r__ <= teps + 1.) {
   goto L50;
    }

/* --UNLESS R GETTING FARTHER FROM 1, DO ANOTHER ITERATION */

    rmu = fabs(r__ - 1.);
    if (rmu > prmu) {
   goto L40;
    }
    prmu = rmu;
    d1 = d2;
    d2 = d3;
    d3 = d4;
    lex = 0;
    return 0;

L40:
    dd = d1;
    lex = 2;
    return 0;

L50:
    dd = (d2 - d3*6. + d4*8.)/3.;
    lex = 1;
    return 0;

} /* vcmx_fitder_ */



// Local syntactic sugar
#define h_ref(a_1,a_2)  my_data.maxf2_.h[(a_2)*my_data.param_NPV + a_1]
#define hn_ref(a_1,a_2) hn[              (a_2)*my_data.param_NPV + a_1]
#define v_ref(a_1,a_2)  my_data.maxf2_.v[(a_2)*my_data.param_NPV + a_1]


/* ************************************************************************** **
** *************                Maxfun::vcmx_                    ************ **
** ************************************************************************** **
   COMPUTE VARIANCE-COVARIANCE MATRIX

   RETURNS FLAG IVFL:
      0:  NO PROBLEM WITH H
      1:  ROUND-OFF ERROR IN H
      3:  H COULD NOT BE INVERTED
      4:  H COULD NOT BE COMPUTED
----------------------------------------------------------------------------- */
int Maxfun::vcmx_(int& ih)
{
   int lex;
   vector<double> hn(my_data.param_NPV*my_data.param_NPV);

   my_data.maxf2_.ivfl = 0; // NO PROBLEM WITH H (default)

   // GET SECOND PARTIAL DERIVATIVES
   vcmx_deriv2_(ih, lex);

   if (lex - 1 < 0)
    goto SPD_OK;
   else if (lex - 1 == 0)
    goto SPD_ROUNDOFF_ERROR;
   else
    goto SPD_FAIL;

SPD_FAIL:
   my_data.maxf2_.ivfl = 4; // H COULD NOT BE COMPUTED
   return 0;

SPD_ROUNDOFF_ERROR:
   my_data.maxf2_.ivfl = 1; // ROUND-OFF ERROR IN H

SPD_OK:
   // PUT NEGATIVE OF MATRIX OF 2ND PARTIALS IN ANOTHER ARRAY
   for (int i = 0; i < my_data.maxf2_.nv; ++i)
      for (int j = 0; j < my_data.maxf2_.nv; ++j)
        hn_ref(i, j) = -h_ref(i, j);

#if 0
   // Print the two matrices
   int m = my_data.maxf2_.nv;
   cout << endl;
   cout << "Matrix hn (Before Inversion)" << endl;
   cout << "----------------------------" << endl;
   for (size_t i = 0; i < (size_t)m ; i++) {
      for (size_t j = 0; j < (size_t)m ; j++) {
         cout << setw(12) << hn_ref(i,j) << " ";
      }
      cout << endl;
   }
   cout << endl;
#endif

   mxnvrt_(hn, my_data.param_NPV, my_data.maxf2_.nv, my_data.maxf2_.v, my_data.param_NPV, lex);

   if (lex <= 0)
   {
      goto MATRIX_INV_OK;
   }

   my_data.maxf2_.ivfl  =  3; // H COULD NOT BE INVERTED
   my_data.maxf2_.ivage = -1;
   return 0;

MATRIX_INV_OK:

#if 0
   // Print the two matrices
   cout << endl;
   cout << "Matrix v (After Inversion)" << endl;
   cout << "----------------------------" << endl;
   for (size_t i = 0; i < (size_t)m ; i++) {
      for (size_t j = 0; j < (size_t)m ; j++) {
         cout << setw(12) << v_ref(i,j) << " ";
      }
      cout << endl;
   }
   cout << endl;
   cout << endl;
   cout << endl;
#endif

   // NOW HAVE NEW V; INDICATE THAT IT MATCHES CURRENT THETA
   my_data.maxf2_.ivage = 0;
   return 0;
} /* vcmx_ */

#undef hn_ref
#undef h_ref
#undef v_ref

/* ************************************************************************** **
** *************                Maxfun::mxnvrt_                  ************ **
** ************************************************************************** **
   INVERT POSITIVE DEFINITE MATRIX A (M BY M) AND PLACE INVERSE MATRIX
   IN B (M BY M).  GAUSS-JORDAN ALGORITHM IS PERFORMED WITHOUT PERMUTATION OF ROWS:
   A ZERO PIVOT ELEMENT AT ANY STEP CAUSES ERROR EXIT

   RETURNS FLAG LEX:
      0:  NO ERROR
      1:  MATRIX CANNOT BE INVERTED BY THIS SUBROUTINE
      2:  INVALID ARGUMENT: ORDER OF MATRIX <= 0

   NOTE THAT INDICES I, J, K ARE LOCAL TO THIS SUBROUTINE AND MAY NOT
   CORRESPOND TO INDICES OF THE SAME NAMES IN OTHER PROGRAM UNITS
----------------------------------------------------------------------------- */
int Maxfun::mxnvrt_(vector<double>& a, int a_dim, int m, vector<double>& b, int b_dim, int& lex)
{
  #undef  a_ref
  #define a_ref(x, y) a[y*a_dim + x]

  #undef  b_ref
  #define b_ref(x, y) b[y*b_dim + x]

#if 0
    cout << endl;
    cout << "lex = " << lex << endl;
    cout << "a before" << endl;
    cout << "--------------------------------" << endl;
    for (size_t i = 0; i < (size_t)m ; i++)
    {
       for (size_t j = 0; j < (size_t)m ; j++)
       {
          cout << setw(12) << a_ref(i,j) << " ";
       }
       cout << endl;
    }
    cout << endl;
    cout << "b before" << endl;
    cout << "--------------------------------" << endl;
    for (size_t i = 0; i < (size_t)m ; i++)
    {
       for (size_t j = 0; j < (size_t)m ; j++)
       {
          cout << setw(12) << b_ref(i,j) << " ";
       }
       cout << endl;
    }
#endif

   /* --------------------------------------------------------------------------
      Find all NaNs on the diagonal (indicating a variable whose second
      derivative can't be computed, and store a list.  Set the diagonal to 1.0
      and the row and column values to 0.0 for inversion.
   -------------------------------------------------------------------------- */
   std::list<size_t> invalid_rows;

   for(size_t i = 0; i < (size_t) m; ++i)
   {
      if( SAGE::isnan(a_ref(i,i)) )
      {
         invalid_rows.push_back(i);

         for(size_t j = 0; j < (size_t) m; ++j){
            if(i == j) a_ref(i,j) = 1.0;
            else       a_ref(i,j) = 0.0;
         }
      }
   }

   std::list<size_t>::const_iterator list_iter;

   if( m - 1 < 0 ) // Invalid matrix order
   {
      lex = 2;
      return 0;
   }
   else if( m - 1 == 0 ) // Scalar case
   {
      if( a_ref(0, 0) == 0.0 )
         goto ERROR_RETURN;

      b_ref(0, 0) = 1.0/a_ref(0, 0);
      goto NORMAL_RETURN;
   }
   else
   {
      //int return_code = do_LU_inverse(a, a_dim, m, b, b_dim, lex);
      int return_code = do_SVD_inverse(a, a_dim, m, b, b_dim, lex);

      if( return_code == 0 )
        goto ERROR_RETURN;
      else
        goto NORMAL_RETURN;
   }

NORMAL_RETURN:

#if 0
    cout << endl;
    cout << "a after" << endl;
    cout << "--------------------------------" << endl;
    for (size_t i = 0; i < (size_t)m ; i++)
    {
       for (size_t j = 0; j < (size_t)m ; j++)
       {
          cout << setw(12) << a_ref(i,j) << " ";
       }
       cout << endl;
    }
    cout << endl;
    cout << "b after" << endl;
    cout << "--------------------------------" << endl;
    for (size_t i = 0; i < (size_t)m ; i++)
    {
       for (size_t j = 0; j < (size_t)m ; j++)
       {
          cout << setw(12) << b_ref(i,j) << " ";
       }
       cout << endl;
    }
#endif

   // Set the diagonal elements of a to 0.0 and b to -1.0
   for( std::list<size_t>::const_iterator li  = invalid_rows.begin();
                                          li != invalid_rows.end(); ++li )
   {
      a_ref(*li, *li) =  0.0;
      b_ref(*li, *li) = -1.0;
   }

   lex = 0;
   return 0;

ERROR_RETURN:
   // Set the diagonal elements of a to 0.0 and b to -1.0
   for( std::list<size_t>::const_iterator li  = invalid_rows.begin();
                                          li != invalid_rows.end(); ++li )
   {
      a_ref(*li, *li) =  0.0;
      b_ref(*li, *li) = -1.0;
   }


   lex = 1;
   return 0;

   #undef b_ref
   #undef a_ref

} /* mxnvrt_ */

} // end of SAGE namespace
