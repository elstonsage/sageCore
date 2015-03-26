#ifndef RegPenetranceCommonAdjustments
#include "segreg/RegPenetranceCommonAdjustments.h"
#endif

namespace SAGE   {
namespace SEGREG {

inline
RegPenetranceCommonAdjustments::RegPenetranceCommonAdjustments()
  : my_model_class    (model_INVALID),
    my_corr_sub_model (NULL),
    my_max_sibs       (0),

    my_spousal_var_adj   (0),
    my_parental_var_adjs ()
{ }

inline
RegPenetranceCommonAdjustments::RegPenetranceCommonAdjustments
    (model_class                           m_class,
     const residual_correlation_sub_model& cor_sub_model,
     unsigned int                          max_sibship_size)
  : my_model_class    (m_class),
    my_corr_sub_model (&cor_sub_model),
    my_max_sibs       (max_sibship_size),
    
    my_spousal_var_adj   (0),
    my_parental_var_adjs ()
{
  my_prior_sib_coeffs.resize(my_max_sibs + 1);
  my_sib_var_adjs    .resize(my_max_sibs + 1);

  // Set fixed values for parental variance adjustments.  See documentation
  // for my_parental_var_adjs for details.
  
  my_parental_var_adjs.set_value(false, false, 0.0);
}

inline
RegPenetranceCommonAdjustments::RegPenetranceCommonAdjustments
    (const RegPenetranceCommonAdjustments& other)
{
  copy(other);
}

inline
RegPenetranceCommonAdjustments& RegPenetranceCommonAdjustments::operator=
    (const RegPenetranceCommonAdjustments & other)
{
  if(this != &other)
    copy(other);

  return *this;
}

inline
int
RegPenetranceCommonAdjustments::update()
{
  calculate_spousal_var_adjustment   ();
  calculate_parental_var_adjustments ();

  if(my_model_class == model_D)
  {
    calculate_sib_values();
  }

  return 0;
}

inline
double
RegPenetranceCommonAdjustments::get_spousal_var_adjustment() const
{
  return my_spousal_var_adj;
}

inline
double
RegPenetranceCommonAdjustments::get_parental_var_adjustment
    (bool mother_informative,
     bool father_informative) const
{
  return my_parental_var_adjs.get_value(mother_informative, father_informative);
}

inline
double
RegPenetranceCommonAdjustments::get_prior_sib_coefficient
    (bool          mother_informative,
     bool          father_informative,
     unsigned int  prior_sib_count    ) const
{
  return my_prior_sib_coeffs[prior_sib_count].get_value(mother_informative, father_informative);
}

inline
double
RegPenetranceCommonAdjustments::get_sibling_var_adjustment
    (bool          mother_informative,
     bool          father_informative,
     unsigned int  prior_sib_count    ) const
{
  return my_sib_var_adjs[prior_sib_count].get_value(mother_informative, father_informative);
}
                                       
inline
void
RegPenetranceCommonAdjustments::copy(const RegPenetranceCommonAdjustments & other)
{
  my_model_class    = other.my_model_class;
  my_corr_sub_model = other.my_corr_sub_model;
  my_max_sibs       = other.my_max_sibs;
  
  my_spousal_var_adj   = other.my_spousal_var_adj;
  my_parental_var_adjs = other.my_parental_var_adjs;
  my_prior_sib_coeffs  = other.my_prior_sib_coeffs;
  my_sib_var_adjs      = other.my_sib_var_adjs;
}

inline
void
RegPenetranceCommonAdjustments::calculate_spousal_var_adjustment   ()
{
  double Pfm = my_corr_sub_model->father_mother_correlation();

  my_spousal_var_adj = - Pfm * Pfm;
}

inline
void
RegPenetranceCommonAdjustments::calculate_parental_var_adjustments ()
{
  double alpha_m, alpha_f;

  double Pms = my_corr_sub_model->mother_son_correlation();
  double Pfs = my_corr_sub_model->father_son_correlation();

  // Set for mother = true, father = false
  
  alpha_m = my_corr_sub_model->alpha_mother(false, true);
  my_parental_var_adjs.set_value(true, false, - alpha_m * Pms);

  // Set for mother = false, father = true

  alpha_f = my_corr_sub_model->alpha_father(true, false);
  my_parental_var_adjs.set_value(false, true, - alpha_f * Pfs);

  // Set for mother = true, father = true

  alpha_m = my_corr_sub_model->alpha_mother(true, true);
  alpha_f = my_corr_sub_model->alpha_father(true, true);
  my_parental_var_adjs.set_value(true, true, - alpha_m * Pms - alpha_f * Pfs);
}

//lint -e{1401}  my_values not initialized, but that's not a problem here
inline
RegPenetranceCommonAdjustments::InformativeValueArray::InformativeValueArray()
{ }

inline
void
RegPenetranceCommonAdjustments::InformativeValueArray::set_value
    (bool   mother_informative,
     bool   father_informative,
     double value)
{
  my_values[calculate_index(mother_informative, father_informative)] = value;
}

inline
double
RegPenetranceCommonAdjustments::InformativeValueArray::get_value
    (bool mother_informative,
     bool father_informative) const
{
  return my_values[calculate_index(mother_informative, father_informative)];
}

inline
RegPenetranceCommonAdjustments::InformativeValueArray::ValueArray::size_type
  RegPenetranceCommonAdjustments::InformativeValueArray::calculate_index
    (bool mother_informative,
     bool father_informative) const
{
  //lint -e{732}
  ValueArray::size_type t = 2 * (mother_informative ? 1 : 0) 
                          + 1 * (father_informative ? 1 : 0);

  return t;
}

}
}
