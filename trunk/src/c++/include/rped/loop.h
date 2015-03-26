#ifndef LOOPCHECKER_H
#define LOOPCHECKER_H

/////////////////////////////////////////////////////////////
// LoopChecker:  function object used to perform pedigree  //
//               loops (marriage loops, non-marriage loops)//
//               checking.                                 //
//                                                         //
// Author : Qing Sun                                       //
//                                                         //
// History: version 1.0, August 8, 1996                    //
//                                                         //
// Copyright (c) 1996 R.C. Elston                          //
//                                                         //
/////////////////////////////////////////////////////////////

#include <set>
#include "mped/sp.h"
#include "mped/mp.h"
#include "error/errorstream.h"
#include "rped/rped.h"

namespace SAGE {
namespace RPED {

class LoopChecker
{
 public:
   LoopChecker(cerrorstream &e = sage_cerr) : errors(e) 
   {
     reset();
   }
  ~LoopChecker(){}

   bool check_pedigree(const RefPedigree*, bool verbose = false);

   bool loops() const 
   { 
     return (my_marriage_loops > 0 || my_non_marriage_loops > 0); 
   }

   size_t marriage_loops()      const { return my_marriage_loops;     }
   size_t non_marriage_loops()  const { return my_non_marriage_loops; }
                                       
 private:
   bool check_mloop(const RefSubpedigree* p, bool verbose);
   bool check_nmloop(const RefSubpedigree* p);
   bool get_spouse_set(int, vector<uint>& ,const RefSubpedigree* p) const;
   void reset();
  
   size_t my_marriage_loops;
   size_t my_non_marriage_loops;

   cerrorstream errors;
};

} // End namespace RPED
} // End namespace SAGE

#endif
