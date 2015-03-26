#ifndef SIBPAL_MATRIX
#define SIBPAL_MATRIX

#include "sibpal/sib_cluster.h"

namespace SAGE   {
namespace SIBPAL {


struct sib_matrix_pattern
{
  sib_matrix_pattern()
  {
    common[0] = common[1] = common[2] = common[3] = std::numeric_limits<double>::quiet_NaN();
  }

  sib_matrix_pattern(double p2, double p1, double p0, double un = 0.0)
  {
    common[0] = p0;
    common[1] = p1;
    common[2] = p2;
    common[3] = un;
  }

  const double& operator[](size_t i) const  { return common[i]; }
        double& operator[](size_t i)        { return common[i]; }

private:  

  double common[4];
};


class sib_weights
{
  public:

    void weight_matrix(const sib_cluster& ship, matrix& W,
                       double p2, double p1, double p0, 
                       weight_status_type& status);

    void weight_matrix(const sib_cluster& ship, matrix& W,
                       const sib_matrix_pattern& p, 
                       weight_status_type& status);

    void weight_matrix_combined(const sib_cluster& ship, matrix& W, double c,
                                double p2_ff, double p1_ff, double p0_ff, 
                                double p2_hh, double p1_hh, double p0_hh, 
                                double p2_fh, double p1_fh, double p0_fh, 
                                weight_status_type& status);

    void weight_matrix_combined(const sib_cluster& ship, matrix& W, double c,
                                const sib_matrix_pattern& p_ff,
                                const sib_matrix_pattern& p_hh,
                                const sib_matrix_pattern& p_fh,
                                weight_status_type& status);
                            
    void normal_diff_weight_matrix(const sib_cluster& ship, matrix& W, double p, 
                                   weight_status_type& status);

    void normal_sum_weight_matrix(const sib_cluster& ship, matrix& W, double r, 
                                  weight_status_type& status);

    void normal_prod_weight_matrix(const sib_cluster& ship, matrix& W, double r,
                                   weight_status_type& status);
/*
    void weighted_sum_matrix(const sib_cluster& ship, matrix& W,
                             const sib_matrix_pattern& p1,
                             double ss1,
                             const sib_matrix_pattern& p2,
                             double ss2,
                             weight_status_type& status);

    void weighted_sum_matrix_full_rank(const sib_cluster& ship, matrix& W,
                             const sib_matrix_pattern& p1,
                             double ss1,
                             const sib_matrix_pattern& p2,
                             double ss2,
                             weight_status_type& status);
*/
  private:

    LUP lup;
    SVD svd;
    matrix W2;
};

void marker_matrix(const sib_cluster& ship, matrix& y, size_t m,
                   double w0 = 0.0, double w1 = 0.5, double w2 = 1.0);

sib_matrix_pattern magic_inverse(size_t s, const sib_matrix_pattern& p);

} //end of namespace SIBPAL
} //end of namespace SAGE

#endif
