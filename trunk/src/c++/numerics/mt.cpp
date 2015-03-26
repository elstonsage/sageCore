#include "numerics/mt.h"
#include <ctime>

namespace SAGE
{

namespace {

unsigned long temper_u(unsigned long y) { return y >> 11; }
unsigned long temper_s(unsigned long y) { return y <<  7; }
unsigned long temper_t(unsigned long y) { return y << 15; }
unsigned long temper_l(unsigned long y) { return y >> 18; }
}

enum private_constants { N = 624, M = 397 };

/* Period parameters */  
#define MATRIX_A   0x9908b0df /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */

/* for tempering */   
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000

MersenneTwister::MersenneTwister(unsigned long seed)
{
  k = N+1;
  mt.resize(N);
  if(seed == 0)
    seed = time(NULL);
  reseed(seed);
}

void MersenneTwister::reseed(unsigned long seed) /* seed should not be 0 */
{
  /* setting initial seeds to mt[N] using     */
  /* the generator Line 25 of Table 1 in          */
  /* [KNUTH 1981, The Art of Computer Programming */
  /*    Vol. 2 (2nd Ed.), pp102]                  */

  mt[0]= seed & 0xffffffff;
  for (k=1; k<N; ++k)
	mt[k] = (69069 * mt[k-1]) & 0xffffffff;
}

unsigned long MersenneTwister::uniform_integer() const
{
  unsigned long y;
  static const unsigned long mag01[2]={0x0, MATRIX_A};
  /* mag01[x] = x * MATRIX_A  for x=0,1 */
  
  if(k == N)  /* generate N words at one time */
  { 
	int j;
	for (j = 0; j < N-M; ++j) 
	{
	  y = (mt[j]&UPPER_MASK)|(mt[j+1]&LOWER_MASK);
	  mt[j] = mt[j+M] ^ (y >> 1) ^ mag01[y & 0x1];
	}
	for (; j < N-1; ++j)
	{
	  y = (mt[j]&UPPER_MASK)|(mt[j+1]&LOWER_MASK);
	  mt[j] = mt[j+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1];
	}
	y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
	mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1];

	k = 0;
  }
  
  y = mt[k++];
  y ^= temper_u(y);
  y ^= temper_s(y) & TEMPERING_MASK_B;
  y ^= temper_t(y) & TEMPERING_MASK_C;
  y &= 0xffffffff; /* you may delete this line if word size = 32 */
  y ^= temper_l(y);

  return y;
}

}
