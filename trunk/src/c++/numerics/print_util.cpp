//***************************************************************************
//* File:           print_util.cpp
//* Author:
//* History:        Version 0.0 Initial implementation
//* Modified by:    Alexandra Borosova, Jul. 2004 
//* Modification 1: "Fixing problem with asterisk."
//	Last change: 	Alexandra Borosova, 07/08/04
//* Copyright (c) 	2001 R.C. Elston
//* All Rights Reserved
//****************************************************************************

#include <math.h>
#include <string>
#include <algorithm>
#include "globals/config.h"
#include "LSF/parse_ops.h"
#include "numerics/print_util.h"

using namespace std;

namespace SAGE {

std::string fp(double d, size_t w, size_t p, char in)
{
  string invalid(in, 45);

  if( !finite(d) )
    return invalid.substr(0,w);

  string n = doub2str(d, w, p, ios::showpoint | ios::fixed);

//  if( n.size() > w )
//    n = doub2str(d, w, p, ios::scientific | ios::showpoint);

  if( n.size() > w)
    n = n.substr(0,w);

  return n;
}

std::string pval(double p, size_t w, int prec, int max_stars)
{
  // NOTE: Work in progress.  Should justify output to max_stars in all
  //       cases.  Right now it works only when max_stars == 2.

  if(prec < 0)
    prec = w - 3 - 2;

  string pv = fp(p,w-3,prec).substr(0,w);

  pv += " ";

  if( !finite(p) )
  {
    for( int i = 0; i < max_stars; ++i )
      pv += " ";
  }
  else if( p < 0.01 )
  {
    int stars = max_stars;

    if( p > 0 )
      stars = std::min( (int)-log10(p), max_stars );

    for( int i = 0; i < stars; ++i )
      pv += "*";

    if( stars < max_stars )
    {
      for( int i = stars; i < max_stars; ++i )
        pv += " ";
    }
  }
  else if( p < 0.05 )
  {
    pv += "*";

    for( int i = 1; i < max_stars; ++i )
      pv += " ";
  }
  else
  {
    for( int i = 0; i < max_stars; ++i )
      pv += " ";
  }

  return pv;
}

//* Modification 1: "Fixing problem with asterisk."
//	Last change: 	Alexandra Borosova, 07/08/04
std::string fp_scientific(double d, size_t w, size_t p, char in)
{
  string invalid(in, 45);

  if( !finite(d) )
    return invalid.substr(0,w);

//  string n = doub2str(d, w, p, ios::showpoint | ios::fixed);
//    cout << "   ALEX070604-1:inside fp n="  << n << endl;
//  if( n.size() > w )
    string n = doub2str(d, w, p, ios::scientific | ios::showpoint);
//    cout << "   ALEX070604-1:inside fp n="  << n << endl;

//  if( n.size() > w)
//    n = n.substr(0,w);
//    cout << "   ALEX070604-2:inside fp n="  << n << endl;

  return n;
}

std::string pval_scientific(double p, size_t w, int prec, int max_stars)
{
  // NOTE: Work in progress.  Should justify output to max_stars in all
  //       cases.  Right now it works only when max_stars == 2.

  if(prec < 0)
    prec = w - 3 - 2;

  string pv = fp_scientific(p,w-3,prec).substr(0,w);

  pv += " ";

  if( !finite(p) )
  {
    for( int i = 0; i < max_stars; ++i )
      pv += " ";
  }
  else if( p < 0.01 )
  {
    int stars = max_stars;

    if( p > 0 )
      stars = std::min( (int)-log10(p), max_stars );

    for( int i = 0; i < stars; ++i )
      pv += "*";

    if( stars < max_stars )
    {
      for( int i = stars; i < max_stars; ++i )
        pv += " ";
    }
  }
  else if( p < 0.05 )
  {
    pv += "*";

    for( int i = 1; i < max_stars; ++i )
      pv += " ";
  }
  else
  {
    for( int i = 0; i < max_stars; ++i )
      pv += " ";
  }

  return pv;
}

}
