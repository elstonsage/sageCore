#ifndef MEAN_SPLIT_H
#define MEAN_SPLIT_H

#include "segreg/sub_model_base.h"
#include "fped/fped.h"
#include <string>
#include <numeric>

namespace SAGE
{
namespace SEGREG
{

/** This class partitions a set of ordered trait values into two
 *  sets such that
 *  \f$ {n_1 n_2 (\bar{y_1}-\bar{y_2})^2}/(S_1+S_2) \f$
 *  is maximized, where,
 *
 *  \f[\bar{y_1}=1/n_1 \sum _{i=1}^{n_1} y_i\f]
 *
 *  \f[\bar{y_2}=1/n_2 \sum _{i=n_1+1}^{n} y_i\f]
 *
 *  \f[n_2 = n-n_1\f]
 *
 *  \f[S_1= \sum _{i=1}^{n_1} (y_i-\bar{y_1})^2\f]
 *
 *  \f[S_1= \sum _{i=n_1+1}^{n} (y_i-\bar{y_2})^2\f]
 *
 *  \f[p=n_1/n\f]
 *
 *  This is used to determine initial estimates for means (
 *  \f$\bar{y_1}\f$, \f$\bar{y_2}\f$) and genotype frequencies
 *  (\f$1 - \sqrt {1-p}\f$).  See Appendix B of the Equ. Doc for more details.
 *  
 *  In some cases, these values can be uncalculatable (if, for example, all
 *  the traits are the same value, then the equation above is infinity).
 *  We use reasonable values in these cases. 
 */
class mean_split
{
public:
  mean_split() { clear(); }
  
  bool operator()(const FPED::Multipedigree& ped_data, const string & trait_name) 
  { 
    return calculate_trait_statistics(ped_data, trait_name); 
  }

  bool get_status() const; // Returns whether the previous mean_spit was valid;

  // The stuff from equations doc.

  double get_y1bar() const;
  double get_y2bar() const;
  double get_S1() const;
  double get_S2() const;
  double get_p() const;
  uint   get_n() const;
  uint   get_n1() const;
  uint   get_n2() const;

  // A few extra quantities to make things easy.

  double get_sum()  const;
  double get_mean() const;
  double get_var()  const;
  double get_max()  const;
  double get_min()  const;

private:
  void clear();

  bool calculate_trait_statistics(const FPED::Multipedigree& RMP, const string& trait_name);

  bool find_best_stats (vector<double>& values);

  void calculate_values(uint tn1, vector<double>& values);
  void calculate_y1bar (vector<double>& values);
  void calculate_y2bar (vector<double>& values);
  void calculate_S1    (vector<double>& values);
  void calculate_S2    (vector<double>& values);
  void calculate_mean  (vector<double>& values);
  void calculate_var   (vector<double>& values);

  double y1bar;
  double y2bar;
  double S1;
  double S2;
  double p;
  uint   n;
  uint   n1;
  uint   n2;

  double mean;
  double var;
  double max_val;
  double min_val;

  bool my_status;
};

}
}

#include "segreg/mean_split.ipp"

#endif
