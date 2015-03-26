//===================================================================
//
//  File:	SL_calculator.ipp
//
//  Author:	Stephen Gross
//
//===================================================================

#ifndef SL_CALCULATOR_H
#include "segreg/SL_calculator.h"
#endif

namespace SAGE {
namespace SEGREG {


//===================================================================
//  FPMM_SL::founder_SL2(...)
//===================================================================
inline log_double
FPMM_SL::founder_SL2(const penetrance_info&          ps_indiv) const
{
  log_double& v = my_founder_SL2_4[ps_indiv.member->subindex()]
                 (ps_indiv.genotype,ps_indiv.polygenotype);

  if(SAGE::isnan(v.get_double()))
    v = int_founder_SL2_4(ps_indiv);

  return v;

}

//===================================================================
//  FPMM_SL::founder_SL4(...)
//===================================================================
inline log_double
FPMM_SL::founder_SL4(const penetrance_info&          ps_indiv) const
{
  log_double& v = my_founder_SL2_4[ps_indiv.member->subindex()]
                 (ps_indiv.genotype,ps_indiv.polygenotype);

  if(SAGE::isnan(v.get_double()))
    v = int_founder_SL2_4(ps_indiv);

  return v;
}

//===================================================================
//  FPMM_SL::founder_SL6(...)
//===================================================================
inline log_double
FPMM_SL::founder_SL6(const penetrance_info&          ps_indiv,
                     const penetrance_info&          ps_spouse) const
{
//  return int_founder_SL6(ps_indiv,ps_spouse);

  for(vector<i_s_info>::iterator
      j  = my_founder_SL6[ps_indiv.member->subindex()].begin (); 
      j != my_founder_SL6[ps_indiv.member->subindex()].end   (); 
    ++j)
  {
    if(j->spouse_id == ps_spouse.member)
    {
      log_double& v = (*j)(ps_indiv .genotype,ps_indiv .polygenotype,
                           ps_spouse.genotype,ps_spouse.polygenotype);

      if(SAGE::isnan(v.get_double())) 
        v = int_founder_SL6(ps_indiv,ps_spouse);

      return v;
    }
  }
  my_founder_SL6[ps_indiv.member->subindex()].push_back(i_s_info(ps_spouse.member));
  return founder_SL6(ps_indiv,ps_spouse);

}

//===================================================================
//  FPMM_SL::nonfounder_SL2(...)
//===================================================================
inline log_double 
FPMM_SL::nonfounder_SL2(const penetrance_info&          ps_indiv,
                        const penetrance_info&          ps_mother,
                        const penetrance_info&          ps_father) const
{
//  return int_nonfounder_SL2(ps_indiv,ps_mother,ps_father);

  log_double& v = my_nonfounder_SL2_4[ps_indiv.member->subindex()]
                 (ps_indiv .genotype,ps_indiv .polygenotype,
                  ps_mother.genotype,ps_mother.polygenotype,
                  ps_father.genotype,ps_father.polygenotype);

  if(SAGE::isnan(v.get_double())) 
    v = int_nonfounder_SL2_4(ps_indiv,ps_mother,ps_father);

  return v;

}

//===================================================================
//  FPMM_SL::nonfounder_SL4(...)
//===================================================================
inline log_double 
FPMM_SL::nonfounder_SL4(const penetrance_info&          ps_indiv,
                        const penetrance_info&          ps_mother,
                        const penetrance_info&          ps_father) const
{
  // Because there are no familial correlations under an FPMM model,
  // SL2 and SL4 are equivalent.

  return nonfounder_SL2(ps_indiv, ps_mother, ps_father);
}

//===================================================================
//  FPMM_SL::nonfounder_SL6(...)
//===================================================================
inline log_double
FPMM_SL::nonfounder_SL6(const penetrance_info&          ps_indiv,
                        const penetrance_info&          ps_spouse) const
{
//  return int_nonfounder_SL6(ps_indiv,ps_spouse);

  for(vector<i_s_info>::iterator
      j  = my_nonfounder_SL6[ps_indiv.member->subindex()].begin ();
      j != my_nonfounder_SL6[ps_indiv.member->subindex()].end   (); 
    ++j)
  {
    if(j->spouse_id == ps_spouse.member)
    {
      log_double& v = (*j)(ps_indiv .genotype,ps_indiv .polygenotype,
                           ps_spouse.genotype,ps_spouse.polygenotype);

      if(SAGE::isnan(v.get_double())) 
        v = int_nonfounder_SL6(ps_indiv,ps_spouse);

      return v;        
    }
  }
  my_nonfounder_SL6[ps_indiv.member->subindex()].push_back(i_s_info(ps_spouse.member));
  return nonfounder_SL6(ps_indiv,ps_spouse);

}

//===================================================================
//  FPMM_SL::i_info::i_info()
//===================================================================
inline
FPMM_SL::i_info::i_info()
{
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < MAX_POLYGENOTYPE; ++j)
      data[i][j] = log_double(numeric_limits<double>::quiet_NaN());
}

//===================================================================
//  FPMM_SL::i_info::i_info() COPY Constructor
//===================================================================
inline
FPMM_SL::i_info::i_info(const i_info& other)
{
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < MAX_POLYGENOTYPE; ++j)
      data[i][j] = other.data[i][j];
}

//===================================================================
//  FPMM_SL::i_info::operator() (...)
//===================================================================
inline log_double &
FPMM_SL::i_info::operator() (int i, int j)
{
  return data[i][j];
}

//===================================================================
//  FPMM_SL::i_s_info::i_s_info()
//===================================================================
inline
FPMM_SL::i_s_info::i_s_info(member_const_pointer id)
{
  spouse_id = id;

  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < MAX_POLYGENOTYPE; ++j)
      for(int k = 0; k < 3; ++k)
        for(int l = 0; l < MAX_POLYGENOTYPE; ++l)
          data[i][j][k][l] = 
            log_double(numeric_limits<double>::quiet_NaN());
}

//===================================================================
//  FPMM_SL::i_s_info::i_s_info() COPY Constructor
//===================================================================
inline
FPMM_SL::i_s_info::i_s_info(const i_s_info& other)
{
  spouse_id = other.spouse_id;

  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < MAX_POLYGENOTYPE; ++j)
      for(int k = 0; k < 3; ++k)
        for(int l = 0; l < MAX_POLYGENOTYPE; ++l)
          data[i][j][k][l] = other.data[i][j][k][l];
}

//===================================================================
//  FPMM_SL::i_s_info::operator() (...)
//===================================================================
inline log_double &
FPMM_SL::i_s_info::operator() (int i, int j, int k, int l)
{
  return data[i][j][k][l];
}

//===================================================================
//  FPMM_SL::i_m_f_info::i_m_f_info()
//===================================================================
inline
FPMM_SL::i_m_f_info::i_m_f_info()
{
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < MAX_POLYGENOTYPE; ++j)
      for(int k = 0; k < 3; ++k)
        for(int l = 0; l < MAX_POLYGENOTYPE; ++l)
          for(int m = 0; m < 3; ++m)
            for(int n = 0; n < MAX_POLYGENOTYPE; ++n)
              data[i][j][k][l][m][n] = 
                log_double(numeric_limits<double>::quiet_NaN());
}

//===================================================================
//  FPMM_SL::i_m_f_info::i_m_f_info() COPY Constructor
//===================================================================
inline
FPMM_SL::i_m_f_info::i_m_f_info(const i_m_f_info& other)
{
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < MAX_POLYGENOTYPE; ++j)
      for(int k = 0; k < 3; ++k)
        for(int l = 0; l < MAX_POLYGENOTYPE; ++l)
          for(int m = 0; m < 3; ++m)
            for(int n = 0; n < MAX_POLYGENOTYPE; ++n)
              data[i][j][k][l][m][n] = other.data[i][j][k][l][m][n];
}

//===================================================================
//  FPMM_SL::i_m_f_info::operator() (...)
//===================================================================
inline log_double &
FPMM_SL::i_m_f_info::operator() (int i, int j, int k, int l, int m, int n)
{
  return data[i][j][k][l][m][n];
}

}}
