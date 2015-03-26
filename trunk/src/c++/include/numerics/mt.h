#ifndef MT_H
#define MT_H

// NOTE:  This class implements a wonderful uniform pseudorandom number
//        generator and should serve as a basis for a general pseudorandom
//        number generation library that inludes more uniform generators that
//        will feed non-uniform generators as well as natural non-uniform
//        generators.   [kbj - 1999/11/13]

// Mersenne Twister (MT19937 or just MT) is a pseudorandom number generator
// developped by Makoto Matsumoto and Takuji Nishimura.  This implementation
// of MT generates 32bit uniform random numbers [0,2^32-1] internally to
// obtain both integer and real uniform random distributions.

#include <vector>
#include <limits>

namespace SAGE
{

class MersenneTwister
{
public:
  MersenneTwister(unsigned long = 0);

  // Interval R[  0,     1)
  double         uniform_real() const;
  // Interval R[  0, range)
  double         uniform_real(double range) const;
  // Interval R[min,   max)
  double         uniform_real(double min, double max) const;

  // Interval N[  0,2^32-1)
  unsigned long  uniform_integer() const;

  // Interval N[  0, range)
  unsigned long  uniform_integer(unsigned long range) const;

  // Interval N[min,   max)
  unsigned long  uniform_integer(unsigned long min,
                                 unsigned long max) const;

  // Seed may not be zero
  void reseed(unsigned long seed);

private:

  unsigned long genrand() const;

  mutable std::vector<unsigned long> mt;             // Working area
  mutable int k;
};

typedef MersenneTwister MT19937;
typedef MersenneTwister MT;

/////////////////////////////////////////////////////////////
///       Inline implementation of MersenneTwister        ///
/////////////////////////////////////////////////////////////

inline double MersenneTwister::uniform_real() const
{
  return ( (double)uniform_integer() / (unsigned long)0xffffffff );
}

inline double MersenneTwister::uniform_real(double range) const
{
  return uniform_real()*range;
}

inline double MersenneTwister::uniform_real(double min, double max) const
{
  if( min >= max )
    return std::numeric_limits<double>::quiet_NaN();
  return uniform_real()*(max-min) + min;
}

inline unsigned long
MersenneTwister::uniform_integer(unsigned long range) const
{
  return uniform_integer() % range;
}

inline unsigned long
MersenneTwister::uniform_integer(unsigned long min, unsigned long max) const
{
  return (uniform_integer() % (max-min)) + min;
}

}

#endif
