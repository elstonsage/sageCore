#ifndef  RELTEST_PAIRS_H
#define  RELTEST_PAIRS_H

//==========================================================================
//  File:       putative_pair.h
//
//  Author:     Qing Sun & Yeunjoo Song
//
//  History:    Version 1.0   
//                      2.0  Took out from analysis
//                           & Updated to new libraries         yjs  Jul. 03
//
//  Notes:      This class defines the data structure for putative pairs.
//
//  Copyright (c) 2003 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "reltest/definitions.h"

namespace SAGE
{

namespace RELTEST
{

class putative_pair
{
  public:

    putative_pair()
     : my_s1(NULL), my_s2(NULL), my_Yj(0.), my_Yjp(0.), my_MIC(0.), my_type(UNRELATED) {};

    ~putative_pair() {};

    ind_id         first  (ind_id p)   { return my_s1=p; }
    ind_id         second (ind_id p)   { return my_s2=p; }

    ind_id         first  ()    const  { return my_s1; }
    ind_id         second ()    const  { return my_s2; }

    // get
    double         get_Yj ()    const  { return my_Yj;  }
    double         get_Yjp()    const  { return my_Yjp; }
    double         get_MIC()    const  { return my_MIC; }

    putative_type  get_pair_type()   const { return my_type; }
    
    //missing data rate : number of missng data markers/total markers
    double         get_s1_missing_rate() const { return my_missing_rate_s1; }
    double         get_s2_missing_rate() const { return my_missing_rate_s2; }
   
    // arithmetic operations
    double         add_Yj(double n)         { return (my_Yj += n);  }
    double         add_Yjp(double n)        { return (my_Yjp += n); }
    double         add_MIC(double n )       { return (my_MIC += n); }

    // set
    double         set_Yj(double n)         { return (my_Yj = n);  }
    double         set_Yjp(double n)        { return (my_Yjp = n); }

    double         set_s1_missing_rate(double d) { return (my_missing_rate_s1=d); }
    double         set_s2_missing_rate(double d) { return (my_missing_rate_s2=d); }

    putative_type  set_type(putative_type c)     { return (my_type=c); }
       
  private:

    ind_id         my_s1;     // first individual in pair
    ind_id         my_s2;     // second individual in pair

    double         my_Yj;     // mean allele_sharing statistic
    double         my_Yjp;    // parent-offspring sharing statistic
    double         my_MIC;    // marker information content

    putative_type  my_type;   // final classification

    //missing rate = # of missing data markers/total # markers
    double         my_missing_rate_s1;
    double         my_missing_rate_s2;
};

} // end of namespace RELTEST

} // end of namespace SAGE

#endif
