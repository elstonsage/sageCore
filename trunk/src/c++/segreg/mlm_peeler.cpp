//===========================================================================
//  File:       mlm_peeler.cpp
//
//  Author:     Kai He                          
//
//  History:
//                                                 
//  Copyright (c) 2001 R. C. Elston
//  All rights reserved
//===========================================================================

#include "segreg/mlm_peeler.h"
#include "mped/mp_utilities.h"

namespace SAGE
{

namespace SEGREG
{



//===========================================================================
//
//  internal_anterior_terminal(...) 
//
//===========================================================================
/// Calculates \f$ant(i,u_i)\f$ if \f$i\f$ is a founder
/**
 *  Calculates:
 *
 *  \f[ant(i,u_i) = \Psi(u_i) \frac{e^{\theta_u(i)y_i}}{1 + e^{\theta_u(i)}}\f]
 */
const mlm_peeler::result_type&
mlm_peeler::internal_anterior_terminal
(const member_type& i, const data_type& g, result_type& r)
{
  double pop_freq = mod.freq_sub_model.prob(g.get_index());
  double penetrance = mcc->get_penetrance(i, g.get_index());

  return r = pop_freq * penetrance;
}

//===========================================================================
//
//  internal_anterior(...)
//
//===========================================================================
/// Calculates \f$ant(i,u_i)\f$ if \f$i\f$ is a nonfounder
const mlm_peeler::result_type&
mlm_peeler::internal_anterior(const member_type& i, const data_type& i_data, result_type& r)
{
  result_type my_anterior   (0.0);

  member_const_pointer m = i.get_mother();
  member_const_pointer f = i.get_father();

  const TypeDescription& tdesc = my_lelements.get_type_description();
  
  for(TypeDescription::StateIterator m_data = tdesc.begin();
      m_data != tdesc.end(); ++m_data)
  {
    result_type anterior1    = anterior              (*m,     *m_data);
    result_type posterior1   = posterior_except_mate (*m, *f, *m_data);

    result_type mother_prod = anterior1 * posterior1;

    for(TypeDescription::StateIterator f_data = tdesc.begin();
        f_data != tdesc.end(); ++f_data)
    {
      result_type anterior2    = anterior              (*f,     *f_data);
      result_type posterior2   = posterior_except_mate (*f, *m, *f_data);

      result_type father_prod = anterior2 * posterior2;

      double rho = family_rho(*m, *f, *m_data, *f_data);
      
      // Get i's transmission and penetrance

      double transm = my_lelements.get_transmission(i_data,*m_data,*f_data);

      double i_pen = mcc->get_penetrance(i, i_data.get_index());

      // Calculate sibling likelihoods.

      result_type sib_prob(1.0);

      member_type::sibling_const_iterator sib = i.sibling_begin();

      for( ; sib != i.sibling_end(); ++sib)
      {
        result_type sib_sum(0.0);

        for(TypeDescription::StateIterator sib_data = tdesc.begin();
            sib_data != tdesc.end(); ++sib_data)
        {
          double sib_transm = mod.transm_sub_model.prob(sib_data->get_index(),
                                                        m_data->get_index(),
                                                        f_data->get_index());
          double sib_pen = mcc->get_penetrance(*sib, sib_data->get_index());
          
          result_type sib_post = posterior(*sib, *sib_data);

          sib_sum += (sib_transm * sib_pen * sib_post);
        }

        sib_prob *= sib_sum;
      }

      my_anterior += mother_prod * father_prod * rho * transm * i_pen * sib_prob;

    }
  }
  
  if(SAGE::isnan(my_anterior.get_double())) r = 0.0;
  else                                r = my_anterior;

  return r;

}

//==========================================================================
//
//  internal_posterior_with_mate(...)
//
//==========================================================================
/// Calculates \f$post(i,s,u_i)\f$
const mlm_peeler::result_type&
mlm_peeler::internal_posterior_with_mate
(const member_type& i, const member_type& sp, const data_type& i_data, result_type& r)
{
  result_type my_posterior   (0.0);

  const TypeDescription& tdesc = my_lelements.get_type_description();
  
  for(TypeDescription::StateIterator sp_data = tdesc.begin();
      sp_data != tdesc.end(); ++sp_data)
  {
    result_type anterior1     = anterior               (sp,  *sp_data);
    result_type posterior1    = posterior_except_mate  (sp,i,*sp_data);
    
    result_type spouse_prod = anterior1 * posterior1;

    double rho = family_rho(i, sp, i_data, *sp_data);

    
    // Calculate child likelihoods.

    result_type child_prob(1.0);

    member_type::offspring_const_iterator child = i.offspring_begin(sp);

    for( ; child != i.offspring_end(); ++child)
    {
      result_type child_sum(0.0);

      for(TypeDescription::StateIterator ch_data = tdesc.begin();
          ch_data != tdesc.end(); ++ch_data)
      {
        double child_transm = my_lelements.get_transmission(*ch_data,i_data,*sp_data);
        double child_pen = mcc->get_penetrance(*child, ch_data->get_index());
        
        result_type child_post = posterior(*child, *ch_data);
        
        child_sum += (child_transm * child_pen * child_post);
      }

      child_prob *= child_sum;
    }

    my_posterior += (spouse_prod * rho * child_prob);

  }

  //cout<<"num_post="<<num<<endl;

  if(SAGE::isnan(my_posterior.get_double())) r = 0.0;
  else                                 r = my_posterior;

  return r;

}

//===========================================================================
//
//  internal_posterior(...) 
//
//===========================================================================
const mlm_peeler::result_type&
mlm_peeler::internal_posterior(const member_type& ind, const data_type& g, result_type& r)
{
  log_double posterior_loop_total(1.0);
  
  member_type::mate_const_iterator mate_loop = ind.mate_begin();

  for( ; mate_loop != ind.mate_end(); ++mate_loop)
  {
      posterior_loop_total *= posterior_with_mate(ind, mate_loop->mate(), g);
  }
   
  if(SAGE::isnan(posterior_loop_total.get_double())) r = 0.0;
  else                                         r = posterior_loop_total;
       
  return r;

}

//==========================================================================
//
//  internal_posterior_except_mate(...)
//  
//==========================================================================
const mlm_peeler::result_type&
mlm_peeler::internal_posterior_except_mate
(const member_type& ind, const member_type& mate, const data_type& g, result_type& r)
{
  log_double posterior_loop_total(1.0);
  
  member_type::mate_const_iterator mate_loop = ind.mate_begin();

  for( ; mate_loop != ind.mate_end(); ++mate_loop)
  {                         
      if(mate_loop->mate().index() != mate.index())
         posterior_loop_total *= posterior_with_mate(ind,mate_loop->mate(),g);
  }
  
  if(SAGE::isnan(posterior_loop_total.get_double())) r = 0.0;
  else                                         r = posterior_loop_total;
  
  return r;
}

//==========================================================================
//
//  internal_posterior_terminal(...)
//
//==========================================================================
const mlm_peeler::result_type&
mlm_peeler::internal_posterior_terminal(const member_type& ind, const data_type& g, result_type& r)
{

  return r = 1.0;

} 

}
}
