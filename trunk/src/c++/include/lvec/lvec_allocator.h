#ifndef LVEC_ALLOCATOR_H
#define LVEC_ALLOCATOR_H

#include <memory>
#include "error/errorstream.h"
#include "error/errormanip.h"
#include <iomanip>

class lvec_allocator : public std::allocator<double>
{
public:

  pointer allocate (size_type n)
  {  
    double* tmp = (double*)(::operator new ((size_type)(n*sizeof (double))));

    // Note:  Change if we incorporate exceptions because we'd never
    // get here - GCW 4-15-99
    if (tmp == 0)
    {
      double ca = (double)n*sizeof(double)/(1<<20);
      double ta = (double)lvs/(1 << 20) + ca;
      SAGE::sage_cerr << SAGE::priority(SAGE::fatal)
                      << "Unable to allocate " << std::setprecision(2) << ca
                      << "MB of memory.  The total ammount of memory "
                      << "required for these analyses (>" 
                      << std::setprecision(2) << ta << "MB) is larger than "
                      << "your computer can currently provide.  Reducing "
                      << "the maximum pedigree size will dramatically "
                      << "reduce the amount of memory required." << std::endl;
      exit(1);
    }

    lvs += n * sizeof(double);

    return tmp; 
  }

  void deallocate (pointer p, size_type n)
  { 
       ::operator delete (p); 

       lvs -= n * sizeof(double);
  }

  static size_type lvec_size() { return lvs; }

private:

  static size_type lvs;
};

#endif
