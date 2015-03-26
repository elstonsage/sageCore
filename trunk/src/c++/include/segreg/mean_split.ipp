#ifndef MEAN_SPLIT_H
#include "segreg/mean_split.h"
#endif

namespace SAGE
{
namespace SEGREG
{

inline bool mean_split::get_status() const { return my_status; }

inline double mean_split::get_y1bar() const { return y1bar;   }
inline double mean_split::get_y2bar() const { return y2bar;   }
inline double mean_split::get_S1()    const { return S1;      }
inline double mean_split::get_S2()    const { return S2;      }
inline double mean_split::get_p()     const { return p;       }
inline uint   mean_split::get_n()     const { return n;       }
inline uint   mean_split::get_n1()    const { return n1;      }
inline uint   mean_split::get_n2()    const { return n2;      }

inline double mean_split::get_sum()   const { return mean*n;  }
inline double mean_split::get_mean()  const { return mean;    }
inline double mean_split::get_var()   const { return var;     }
inline double mean_split::get_max()   const { return max_val; }
inline double mean_split::get_min()   const { return min_val; }

inline void mean_split::clear()
{
  // Initialize the storage

  y1bar = QNAN;
  y2bar = QNAN;
  S1    = QNAN;
  S2    = QNAN;
  p     = QNAN;
  n     = 0;
  n1    = 0;
  n2    = 0;

  mean    = QNAN;
  var     = QNAN;
  min_val = QNAN;
  max_val = QNAN;

  my_status = false;
}

inline void mean_split::calculate_values(uint tn1, vector<double>& values)
{
  n1 = tn1;
  n2 = n - n1;
  
  calculate_y1bar (values);
  calculate_y2bar (values);
  calculate_S1    (values);
  calculate_S2    (values);

}

//===================================================================================
inline void 
mean_split::calculate_y1bar(vector<double> & values)
{
  double summation = std::accumulate(values.begin(), values.begin() + n1, 0.0);

  y1bar = summation / n1;
}

//===================================================================================
inline void 
mean_split::calculate_y2bar(vector<double> & values)
{
  double summation = std::accumulate(values.begin() + n1, values.end(),0.0);

  y2bar = summation / n2;
}

//===================================================================================
inline void 
mean_split::calculate_S1(vector<double> & values)
{
  S1 = 0.0;

  for(size_t i = 0; i < n1; ++i)
    S1 += (values[i] - y1bar) * (values[i] - y1bar);
}

//===================================================================================
inline void 
mean_split::calculate_S2(vector<double> & values)
{
  S2 = 0.0;

  for(size_t i = n1; i < values.size(); ++i)
    S2 += (values[i] - y2bar) * (values[i] - y2bar);
}

//===================================================================================
// The vector isn't used in this case, but is passed for consistency.  The
// compiler should take care of this.
inline void 
mean_split::calculate_mean(vector<double>&)
{
  mean = (y1bar * n1 + y2bar * n2) / n;
}

//===================================================================================
inline void 
mean_split::calculate_var(vector<double> & values)
{
  var = 0.0;

  for(size_t i = 0; i < n; ++i)
    var += (values[i] - mean) * (values[i] - mean);

  var /= n;
}

}
}
