#include <iomanip>
#include "segreg/RegPenetranceCommonAdjustments.h"

namespace SAGE
{
namespace SEGREG
{

void RegPenetranceCommonAdjustments::calculate_sib_values  ()
{
  if(my_max_sibs == 0) return;

  double Pbb = my_corr_sub_model->brother_brother_correlation();
  
  for(int informativity_counter = 0; informativity_counter < 4; ++informativity_counter)
  {
    bool   mother_informative = informativity_counter & 2, // 2 is 10 in binary
           father_informative = informativity_counter & 1; // 1 is 01 in binary

    double delta = my_corr_sub_model->delta(father_informative, mother_informative);

    double difference = Pbb - delta;

    my_prior_sib_coeffs[0].set_value(mother_informative, father_informative, 0.0);
    my_sib_var_adjs    [0].set_value(mother_informative, father_informative, 0.0);

    for(unsigned int j = 1; j <= my_max_sibs; ++j)
    {
      double Bj = (difference) / ((1.0 - Pbb) + j * difference);

      my_prior_sib_coeffs[j].set_value(mother_informative, father_informative, Bj);
                 
      my_sib_var_adjs[j].set_value(mother_informative, father_informative,
                                   - (Bj * Pbb * j));
    }
  }
}

void test_common_reg_corr_adjustments
   (ostream& o,
    model_class m,
    const residual_correlation_sub_model& r)
{
  if(m != model_A && m != model_D) return;
  
  unsigned sibsize = 0;
  
  if(m == model_D)
  {
    sibsize = 5;
  }

  RegPenetranceCommonAdjustments crca(m, r, sibsize);
  
  int err = crca.update();

  if(err)
  {
    o << "RegPenetranceCommonAdjustments failed to update properly.  Big ERROR!"
      << std::endl;
    exit(1);
  }
  
  o << "Testing RegPenetranceCommonAdjustments: " << ((m == model_A) ? "A" : "D") << std::endl;
  o << std::endl;
  
  o << "Important Correlations:" << std::endl
    << "=======================" << std::endl
    << "Pfm : " << r.father_mother_correlation()   << std::endl
    << "Pms : " << r.mother_son_correlation()      << std::endl
    << "Pfs : " << r.father_son_correlation()      << std::endl
    << "Pbb : " << r.brother_brother_correlation() << std::endl
    << std::endl;
  
  o << "Spousal Variance Adjustment: " << crca.get_spousal_var_adjustment()
    << std::endl << std::endl;
  
  o << "Moth Inf  Fath Inf  Par Var Adj" << std::endl
    << "========  ========  ===========" << std::endl;
    
  o << "false     false     "
    << setw(11) << crca.get_parental_var_adjustment(false,false) << std::endl;
    
  o << "false     true      "
    << setw(11) << crca.get_parental_var_adjustment(false,true) << std::endl;

  o << "true      false     "
    << setw(11) << crca.get_parental_var_adjustment(true,false) << std::endl;

  o << "true      true      "
    << setw(11) << crca.get_parental_var_adjustment(true,true) << std::endl;

  if(m == model_D)
  {
    o << std::endl;
    
    o << "Sib  FF Sib Coeff  FF Sib Var Adj  FT Sib Coeff  FT Sib Var Adj  TF Sib Coeff  TF Sib Var Adj  TT Sib Coeff  TT Sib Var Adj" << std::endl
      << "===  ============  ==============  ============  ==============  ============  ==============  ============  ==============" << std::endl;
    
    for(size_t j = 0; j <= sibsize; ++j)
    {
      o << setw(3)  << j << "  "
        << setw(12) << crca.get_prior_sib_coefficient (false, false, j) << "  "
        << setw(14) << crca.get_sibling_var_adjustment(false, false, j) << "  "
        << setw(12) << crca.get_prior_sib_coefficient (false, true,  j) << "  "
        << setw(14) << crca.get_sibling_var_adjustment(false, true,  j) << "  "
        << setw(12) << crca.get_prior_sib_coefficient (true,  false, j) << "  "
        << setw(14) << crca.get_sibling_var_adjustment(true,  false, j) << "  "
        << setw(12) << crca.get_prior_sib_coefficient (true,  true,  j) << "  "
        << setw(14) << crca.get_sibling_var_adjustment(true,  true,  j) << std::endl;
    }
  }

  o << std::endl;   
}

}
}
