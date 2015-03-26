//=============================================================================
// File:    covariance_calculator.cpp
//
// Author:  Yeunjoo Song
//
// History: Initial implementation                                   yjs Jan 10
//
// Notes:   This source file implements a calculator of 
//          allele-sharing covariance matrix.
//
// Copyright (c) 2010   R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "relpal/covariance_calculator.h"

namespace SAGE   {
namespace RELPAL {

covariance_calculator::covariance_calculator()
                     : my_valid(false)
{}

covariance_calculator::~covariance_calculator()
{}

void
covariance_calculator::compute(FPED::SubpedigreeConstPointer sp)
{
  if( sp == NULL || !sp->member_count() )
  {
    my_valid = false;
    return;
  }

  my_sped = sp;

  compute_kinship_2p();
  compute_kinship_3p();
  //compute_kinship_4p();
  compute_covariance();

  return;
}

void
covariance_calculator::get_covariance_matrix(const vector<mem_pointer>& mem, trimatrix& cov_matrix) const
{
  size_t m = mem.size();

  cov_matrix.resize(m*m, QNAN);

  for( size_t i = 0; i < m; ++i )
  {
    size_t a = mem[i]->subindex();

    for( size_t j = 0; j <= i; ++j )
    {
      size_t b = mem[j]->subindex();

      for( size_t k = 0; k < m; ++k )
      {
        size_t c = mem[k]->subindex();

        for( size_t l = 0; l <= k; ++l )
        {
          size_t d = mem[l]->subindex();
          size_t row = a*mem[i]->subpedigree()->member_count() + b;
          size_t col = c*mem[i]->subpedigree()->member_count() + d;
#if 0
          if( my_cov_matrix(row, col) )
          {
            cout << i << "(" << mem[i]->name() << ") a = " << a << " ";
            cout << j << "(" << mem[j]->name() << ") b = " << b << " ";
            cout << k << "(" << mem[k]->name() << ") c = " << c << " ";
            cout << l << "(" << mem[l]->name() << ") d = " << d << "  = ";
            cout << my_cov_matrix(row, col) << endl;
          }
#endif
            cov_matrix(i*m + j, k*m + l) = cov_matrix(i*m + j, l*m + k)
          = cov_matrix(j*m + i, k*m + l) = cov_matrix(j*m + i, l*m + k) = my_cov_matrix(row, col);
        }
      }
    }
  }
#if 0
  print_trimatrix(cov_matrix, cout, "new cov_matrix");
#endif

  return;
}

void
covariance_calculator::get_covariance_matrix(const vector<mem_pointer>& mem,
                                             const ibd_state_info&      is_info,
                                             size_t                     mi,
                                             trimatrix&                 cov_matrix) const
{
  size_t m = mem.size();

  cov_matrix.resize(m*m, 0.0);

  const vector<id_pair>&    p_ids    = is_info.pair_ids;
  const a_marker_ibd_state& i_states = is_info.ibd_states[mi];

  CorrelationInfo  mean_share_info;
  mean_share_info.resize(p_ids.size());

  for( size_t i = 0; i < i_states.size(); ++i )
  {
    double l_prob = i_states[i].first;
    const vector<size_t>& a_state = i_states[i].second;

    vector<double>   a_state_d(a_state.size());
    for( size_t s = 0; s < a_state.size(); ++s )
      a_state_d[s] = (double)a_state[s]; 

    mean_share_info.add(a_state_d, l_prob/1.0);
  }

  for( size_t i = 0; i < m; ++i )
  {
    for( size_t j = 0; j <= i; ++j )
    {
      for( size_t k = 0; k < m; ++k )
      {
        for( size_t l = 0; l <= k; ++l )
        {
          size_t row_m = is_info.find_pair_ids(mem[i], mem[j]);
          size_t col_m = is_info.find_pair_ids(mem[k], mem[l]);

          double cov_m = 0.0;

          if( row_m < is_info.pair_ids.size() && col_m < is_info.pair_ids.size() )
            cov_m = mean_share_info.covariance(row_m, col_m);
#if 0
          cout << i << "(" << mem[i]->name() << ") a = " << a << " ";
          cout << j << "(" << mem[j]->name() << ") b = " << b << " ";
          cout << k << "(" << mem[k]->name() << ") c = " << c << " ";
          cout << l << "(" << mem[l]->name() << ") d = " << d << "  = ";
          cout << cov_m << endl;
#endif
            cov_matrix(i*m + j, k*m + l) = cov_matrix(i*m + j, l*m + k)
          = cov_matrix(j*m + i, k*m + l) = cov_matrix(j*m + i, l*m + k) = cov_m;
        }
      }
    }
  }
#if 0
  print_trimatrix(cov_matrix, cout, "cov_m_matrix");
#endif

  return;
}

void
covariance_calculator::compute_kinship_2p()
{
  int m = my_sped->member_count();

  // 2p kinship coefficient
  my_kinships_2p.resize(m, QNAN);

  ostringstream names;

  for( int a = 0; a < m; ++a )
  {
    names << my_sped->member_index(a).name() << " ";

    for( int b = 0; b <= a; ++b )
    {
      const member& mid1 = my_sped->member_index(a);
      const member& mid2 = my_sped->member_index(b);

      if( mid1.parent1() || mid2.parent1() )
      {
        int ap1 = -1, ap2 = -1, bp1 = -1, bp2 = -1;

        if( mid1.parent1() )
        {
          ap1 = mid1.parent1()->subindex();
          ap2 = mid1.parent2()->subindex();
        }
        if( mid2.parent1() )
        {
          bp1 = mid2.parent1()->subindex();
          bp2 = mid2.parent2()->subindex();
        }

        if( a == b )
        {
          my_kinships_2p(a,b) = (1.0 + my_kinships_2p(ap1, ap2)) / 2.0;
        }
        else if( max(max(ap1, ap2), max(bp1, bp2)) == max(ap1, ap2) )
        {
          my_kinships_2p(a,b) = (my_kinships_2p(ap1, b) + my_kinships_2p(ap2, b)) / 2.0;
        }
        else
        {
          my_kinships_2p(a,b) = (my_kinships_2p(bp1, a) + my_kinships_2p(bp2, a)) / 2.0;
        }
      }
      else
      {
        if( a == b )
        {
          my_kinships_2p(a,b) = 0.5;
        }
        else
        {
          my_kinships_2p(a,b) = 0.0;
        }
      }
    }
  }

#if 0
  cout << names.str() << endl;
  print_trimatrix(my_kinships_2p, cout, "2p");
#endif

  return;
}

void
covariance_calculator::compute_kinship_3p()
{
  int m = my_sped->member_count();

  // 3p kinship coefficient
  my_kinships_3p.resize(0);

  for( int a = 0; a < m; ++a )
  {
    trimatrix kinships_bc(a+1, QNAN);
    my_kinships_3p.push_back(kinships_bc);
  }

  //cout << "my_kinships_3p size = " << my_kinships_3p.size() << endl;

  for( int a = 0; a < m; ++a )
  {
    trimatrix my_kinships_bc(a+1, QNAN);

    ostringstream names_bc;

    for( int b = 0; b <= a; ++b )
    {
      names_bc << my_sped->member_index(b).name() << " ";

      for( int c = 0; c <= b; ++c )
      {
        const member& mid1 = my_sped->member_index(a);
        const member& mid2 = my_sped->member_index(b);
        const member& mid3 = my_sped->member_index(c);
#if 0
        cout << a << "(" << mid1.name() << ") ";
        cout << b << "(" << mid2.name() << ") ";
        cout << c << "(" << mid3.name() << "), ";
#endif
        double k3_val = QNAN;

        if( mid1.parent1() || mid2.parent1() || mid3.parent1() )
        {
          int ap1 = -1, ap2 = -1, bp1 = -1, bp2 = -1, cp1 = -1, cp2 = -1;

          if( mid1.parent1() )
          {
            ap1 = mid1.parent1()->subindex();
            ap2 = mid1.parent2()->subindex();

            //cout << "ap1 = " << ap1 << " ap2 = " << ap2 << " ";
          }
          if( mid2.parent1() )
          {
            bp1 = mid2.parent1()->subindex();
            bp2 = mid2.parent2()->subindex();

            //cout << "bp1 = " << bp1 << " bp2 = " << bp2 << " ";
          }
          if( mid3.parent1() )
          {
            cp1 = mid3.parent1()->subindex();
            cp2 = mid3.parent2()->subindex();

            //cout << "cp1 = " << cp1 << " cp2 = " << cp2 << " ";
          }

          if( a == b )
          {
            if( b == c ) // a=b=c
            {
              //cout << "case a=b=c " << flush;

              k3_val = (1.0 + 3.0*my_kinships_2p(ap1,ap2)) / 4.0;
            }
            else // a=b, c
            {
              //cout << "case a=b, c " << flush;

              if( max(max(ap1,ap2), max(cp1,cp2)) == max(ap1,ap2) )
                k3_val = (my_kinships_2p(a,c) + get_3p_value(ap1, ap2, c)) / 2.0;
              else
                k3_val = (my_kinships_2p(c,a) + get_3p_value(cp1, cp2, a)) / 2.0;
            }
          }
          else if( b == c ) // a, b=c
          {
            //cout << "case a, b=c " << flush;

            if( max(max(ap1,ap2), max(bp1,bp2)) == max(ap1,ap2) )
              k3_val = (my_kinships_2p(a,b) + get_3p_value(ap1, ap2, b)) / 2.0;
            else
              k3_val = (my_kinships_2p(b,a) + get_3p_value(bp1, bp2, a)) / 2.0;
          }
          else // a, b, c
          {
            //cout << "case a, b, c " << flush;

            if(    max(max(ap1,ap2), max(bp1,bp2)) == max(ap1,ap2)
                && max(max(ap1,ap2), max(cp1,cp2)) == max(ap1,ap2) )
            {
              k3_val =  (get_3p_value(ap1, b, c)
                       + get_3p_value(ap2, b, c)) / 2.0;
            }
            else if(    max(max(bp1,bp2), max(ap1,ap2)) == max(bp1,bp2)
                     && max(max(bp1,bp2), max(cp1,cp2)) == max(bp1,bp2) )
            {
              k3_val =  (get_3p_value(bp1, a, c)
                       + get_3p_value(bp2, a, c)) / 2.0;
            }
            else
            {
              k3_val =  (get_3p_value(cp1, a, b)
                       + get_3p_value(cp2, a, b)) / 2.0;
            }
          }
        }
        else
        {
          //cout << "case all founders " << flush;

          if( a == b && b == c )
          {
            k3_val = 0.25;
          }
          else
          {
            k3_val = 0.0;
          }
        }

        if( k3_val < 1.0e-15 )
          k3_val = 0.0;
#if 0
        cout << " = " << k3_val << endl;
#endif
        my_kinships_3p[a](b,c) = k3_val;
      }
    }
#if 0
    cout << "a = " << my_sped->member_index(a).name() << endl
         << names_bc.str() << endl;
    print_trimatrix(my_kinships_3p[a], cout, "kinship_bc");
#endif
  }

  return;
}

void
covariance_calculator::compute_kinship_4p()
{
  int m = my_sped->member_count();

  // 4p kinship coefficient
  my_kinships_4p.resize(0);

  for( int a = 0; a < m; ++a )
  {
    vector< trimatrix > kinships_bcd;

    for( int b = 0; b <= a; ++b )
    {
      trimatrix kinships_cd(b+1, QNAN);
      kinships_bcd.push_back(kinships_cd);
    }

    my_kinships_4p.push_back(kinships_bcd);
  }

  //cout << "my_kinships_4p size = " << my_kinships_4p.size() << endl;

  for( int a = 0; a < m; ++a )
  {
    for( int b = 0; b <= a; ++b )
    {
      ostringstream names_cd;

      for( int c = 0; c <= b; ++c )
      {
        names_cd << my_sped->member_index(c).name() << " ";

        for( int d = 0; d <= c; ++d )
        {
          const member& mid1 = my_sped->member_index(a);
          const member& mid2 = my_sped->member_index(b);
          const member& mid3 = my_sped->member_index(c);
          const member& mid4 = my_sped->member_index(d);
#if 0
          cout << a << "(" << mid1.name() << ") ";
          cout << b << "(" << mid2.name() << ") ";
          cout << c << "(" << mid3.name() << ") ";
          cout << d << "(" << mid4.name() << "), ";
#endif
          double k4_val = QNAN;

          if( mid1.parent1() || mid2.parent1() || mid3.parent1() || mid4.parent1() )
          {
            int ap1 = -1, ap2 = -1, bp1 = -1, bp2 = -1, cp1 = -1, cp2 = -1, dp1 = -1, dp2 = -1;

            if( mid1.parent1() )
            {
              ap1 = mid1.parent1()->subindex();
              ap2 = mid1.parent2()->subindex();

              //cout << "ap1 = " << ap1 << " ap2 = " << ap2 << " ";
            }
            if( mid2.parent1() )
            {
              bp1 = mid2.parent1()->subindex();
              bp2 = mid2.parent2()->subindex();

              //cout << "bp1 = " << bp1 << " bp2 = " << bp2 << " ";
            }
            if( mid3.parent1() )
            {
              cp1 = mid3.parent1()->subindex();
              cp2 = mid3.parent2()->subindex();

              //cout << "cp1 = " << cp1 << " cp2 = " << cp2 << " ";
            }
            if( mid4.parent1() )
            {
              dp1 = mid4.parent1()->subindex();
              dp2 = mid4.parent2()->subindex();

              //cout << "dp1 = " << dp1 << " dp2 = " << dp2 << " ";
            }

            if( a == b )
            {
              if( b == c )
              {
                if( c == d ) // a=b=c=d
                {
                  //cout << "case a=b=c=d " << flush;

                  k4_val = (1.0 + 7.0*my_kinships_2p(ap1, ap2)) / 8.0;
                }
                else // a=b=c, d
                {
                  //cout << "case a=b=c, d " << flush;

                  if( max(max(ap1,ap2), max(dp1,dp2)) == max(ap1,ap2) )
                  {
                    k4_val = (my_kinships_2p(a,d) + 3.0*get_3p_value(ap1, ap2, d)) / 4.0;
                  }
                  else
                  {
                    k4_val = (my_kinships_2p(d,a) + 3.0*get_3p_value(dp1, dp2, a)) / 4.0;
                  }
                }
              }
              else if( c == d ) // a=b, c=d
              {
                //cout << "case a=b, c=d " << flush;

                if( max(max(ap1,ap2), max(cp1,cp2)) == max(ap1,ap2) )
                {
                  k4_val = (my_kinships_2p(a,c) + 3.0*get_3p_value(ap1, ap2, c)) / 4.0;
                }
                else
                {
                  k4_val = (my_kinships_2p(c,a) + 3.0*get_3p_value(cp1, cp2, a)) / 4.0;
                }
              }
              else // a=b, c, d
              {
                //cout << "case a=b, c, d " << flush;

                if(    max(max(ap1,ap2), max(cp1,cp2)) == max(ap1,ap2)
                    && max(max(ap1,ap2), max(dp1,dp2)) == max(ap1,ap2) )
                {
                  k4_val = (my_kinships_3p[a](c,d) + get_4p_value(ap1, ap2, c, d)) / 2.0;
                }
                else if(    max(max(cp1,cp2), max(ap1,ap2)) == max(cp1,cp2)
                         && max(max(cp1,cp2), max(dp1,dp2)) == max(cp1,cp2) )
                {
                  k4_val = (my_kinships_3p[a](c,d) + get_4p_value(cp1, cp2, a, d)) / 2.0;
                }
                else
                {
                  k4_val = (my_kinships_3p[a](c,d) + get_4p_value(dp1, dp2, a, c)) / 2.0;
                }
              }
            }
            else if( b == c )
            {
              if( c == d ) // a, b=c=d
              {
                //cout << "case a, b=c=d " << flush;

                if( max(max(ap1,ap2), max(bp1,bp2)) == max(ap1,ap2) )
                {
                  k4_val = (my_kinships_2p(a,b) + 3.0*get_3p_value(ap1, ap2, b)) / 4.0;
                }
                else
                {
                  k4_val = (my_kinships_2p(b,a) + 3.0*get_3p_value(bp1, bp2, a)) / 4.0;
                }
              }
              else // a, b=c, d
              {
                //cout << "case a, b=c, d " << flush;

                if(    max(max(ap1,ap2), max(bp1,bp2)) == max(ap1,ap2)
                    && max(max(ap1,ap2), max(dp1,dp2)) == max(ap1,ap2) )
                {
                  k4_val = (my_kinships_3p[a](b,d) + get_4p_value(ap1, ap2, b, d)) / 2.0;
                }
                else if(    max(max(bp1,bp2), max(ap1,ap2)) == max(bp1,bp2)
                         && max(max(bp1,bp2), max(dp1,dp2)) == max(bp1,bp2) )
                {
                  k4_val = (my_kinships_3p[a](b,d) + get_4p_value(bp1, bp2, a, d)) / 2.0;
                }
                else
                {
                  k4_val = (my_kinships_3p[a](b,d) + get_4p_value(dp1, dp2, a, b)) / 2.0;
                }
              }
            }
            else if( c == d ) // a, b, c=d
            {
              //cout << "case a, b, c=d " << flush;

              if(    max(max(ap1,ap2), max(bp1,bp2)) == max(ap1,ap2)
                  && max(max(ap1,ap2), max(cp1,cp2)) == max(ap1,ap2) )
              {
                k4_val = (my_kinships_3p[a](b,c) + get_4p_value(ap1, ap2, b, c)) / 2.0;
              }
              else if(    max(max(bp1,bp2), max(ap1,ap2)) == max(bp1,bp2)
                       && max(max(bp1,bp2), max(cp1,cp2)) == max(bp1,bp2) )
              {
                k4_val = (my_kinships_3p[a](b,c) + get_4p_value(bp1, bp2, a, c)) / 2.0;
              }
              else
              {
                k4_val = (my_kinships_3p[a](b,c) + get_4p_value(cp1, cp2, a, b)) / 2.0;
              }
            } 
            else // a, b, c, d
            {
              //cout << "case a, b, c, d " << flush;

              if(    max(max(ap1,ap2), max(bp1,bp2)) == max(ap1,ap2)
                  && max(max(ap1,ap2), max(cp1,cp2)) == max(ap1,ap2)
                  && max(max(ap1,ap2), max(dp1,dp2)) == max(ap1,ap2) )
              {
                k4_val =  (get_4p_value(ap1, b, c, d)
                         + get_4p_value(ap2, b, c, d)) / 2.0;
              }
              else if(    max(max(bp1,bp2), max(ap1,ap2)) == max(bp1,bp2)
                       && max(max(bp1,bp2), max(cp1,cp2)) == max(bp1,bp2)
                       && max(max(bp1,bp2), max(dp1,dp2)) == max(bp1,bp2) )
              {
                k4_val =  (get_4p_value(bp1, b, c, d)
                         + get_4p_value(bp2, b, c, d)) / 2.0;

              }
              else if(    max(max(cp1,cp2), max(ap1,ap2)) == max(cp1,cp2)
                       && max(max(cp1,cp2), max(bp1,bp2)) == max(cp1,cp2)
                       && max(max(cp1,cp2), max(dp1,dp2)) == max(cp1,cp2) )
              {
                k4_val =  (get_4p_value(cp1, a, b, d)
                         + get_4p_value(cp2, a, b, d)) / 2.0;
              }
              else
              {
                k4_val =  (get_4p_value(dp1, a, b, c)
                         + get_4p_value(dp2, a, b, c)) / 2.0;
              }
            }
          }
          else
          {
            //cout << "case all founders " << flush;

            if( a == b && b == c && c == d )
            {
              k4_val = 0.125;
            }
            else
            {
              k4_val = 0.0;
            }
          }

          if( k4_val < 1.0e-15 )
            k4_val = 0.0;
#if 0
          cout << " = " << k4_val << endl;
#endif
          my_kinships_4p[a][b](c,d) = k4_val;
        }
      }
#if 0
      cout << "a = " << my_sped->member_index(a).name() << " "
           << "b = " << my_sped->member_index(b).name() << " c,d = "
           << names_cd.str() << endl;
      print_trimatrix(my_kinships_4p[a][b], cout, "kinship_cd");
#endif
    }
  }

  return;
}

void
covariance_calculator::compute_covariance()
{
  int m = my_sped->member_count();
#if 0
  cout << endl << "*** Start covariance matrix for " << my_sped->name() << " ***" << endl;
#endif
  my_cov_matrix.resize(m*m, QNAN);

  ostringstream names_ab;
  for( int a = 0; a < m; ++a )
  {
    for( int b = 0; b <= a; ++b )
    {
      names_ab << my_sped->member_index(a).name() << ","
               << my_sped->member_index(b).name() << " ";

      for( int c = 0; c < m; ++c )
      {
        for( int d = 0; d <= c; ++d )
        {
          const member& mid1 = my_sped->member_index(a);
          const member& mid2 = my_sped->member_index(b);
          const member& mid3 = my_sped->member_index(c);
          const member& mid4 = my_sped->member_index(d);
#if 0
          cout << a << "(" << mid1.name() << ") ";
          cout << b << "(" << mid2.name() << ") ";
          cout << c << "(" << mid3.name() << ") ";
          cout << d << "(" << mid4.name() << "), ";
#endif
          double cov_val = QNAN;

          if( is_constant_pair(&mid1, &mid2) || is_constant_pair(&mid3, &mid4) )
          {
#if 0
            cout << "case constant " << flush;
#endif
            cov_val = 0.0;
          }
          else
          {
#if 0
            cout << "case not constant " << flush;
#endif
            int ap1 = -1, ap2 = -1, bp1 = -1, bp2 = -1, cp1 = -1, cp2 = -1;

            if( mid1.parent1() )
            {
              ap1 = mid1.parent1()->subindex();
              ap2 = mid1.parent2()->subindex();
#if 0
              cout << "ap1 = " << ap1 << " ap2 = " << ap2 << " ";
#endif
            }
            if( mid2.parent1() )
            {
              bp1 = mid2.parent1()->subindex();
              bp2 = mid2.parent2()->subindex();
#if 0
              cout << "bp1 = " << bp1 << " bp2 = " << bp2 << " ";
#endif
            }
            if( mid3.parent1() )
            {
              cp1 = mid3.parent1()->subindex();
              cp2 = mid3.parent2()->subindex();
#if 0
              cout << "cp1 = " << cp1 << " cp2 = " << cp2 << " ";
#endif
            }

            if( a == c )
            {
              double pbqd = my_cov_matrix(ap1*m+b, ap2*m+d);
              double qbpd = my_cov_matrix(ap2*m+b, ap1*m+d);
              double pbd  = get_3p_value(ap1, b, d);
              double qbd  = get_3p_value(ap2, b, d);
              double pb   = my_kinships_2p(ap1, b);
              double pd   = my_kinships_2p(ap1, d);
              double qb   = my_kinships_2p(ap2, b);
              double qd   = my_kinships_2p(ap2, d);
#if 0
              cout << "case ab,ad" << flush;
              cout << " pbqd = " << pbqd
                   << " qbpd = " << qbpd
                   << " pbd = " << pbd
                   << " qbd = " << qbd
                   << " pb = " << pb
                   << " pd = " << pd
                   << " qb = " << qb
                   << " qd = " << qd
                   << " " << flush;
#endif
              cov_val = ((pbqd+qbpd) / 4.0) + pbd + qbd - (pb*pd) - (qb*qd);
            }
            else if( a == d )
            {
              double pbqd = my_cov_matrix(ap1*m+b, ap2*m+c);
              double qbpd = my_cov_matrix(ap2*m+b, ap1*m+c);
              double pbd  = get_3p_value(ap1, b, c);
              double qbd  = get_3p_value(ap2, b, c);
              double pb   = my_kinships_2p(ap1, b);
              double pd   = my_kinships_2p(ap1, c);
              double qb   = my_kinships_2p(ap2, b);
              double qd   = my_kinships_2p(ap2, c);
#if 0
              cout << "case ab,ca" << flush;
              cout << " pbqd = " << pbqd
                   << " qbpd = " << qbpd
                   << " pbd = " << pbd
                   << " qbd = " << qbd
                   << " pb = " << pb
                   << " pd = " << pd
                   << " qb = " << qb
                   << " qd = " << qd
                   << " " << flush;
#endif
              cov_val = ((pbqd+qbpd) / 4.0) + pbd + qbd - (pb*pd) - (qb*qd);
            }
            else if( b == c )
            {
              double pbqd = my_cov_matrix(bp1*m+a, bp2*m+d);
              double qbpd = my_cov_matrix(bp2*m+a, bp1*m+d);
              double pbd  = get_3p_value(bp1, a, d);
              double qbd  = get_3p_value(bp2, a, d);
              double pb   = my_kinships_2p(bp1, a);
              double pd   = my_kinships_2p(bp1, d);
              double qb   = my_kinships_2p(bp2, a);
              double qd   = my_kinships_2p(bp2, d);
#if 0
              cout << "case ab,bd" << flush;
              cout << " pbqd = " << pbqd
                   << " qbpd = " << qbpd
                   << " pbd = " << pbd
                   << " qbd = " << qbd
                   << " pb = " << pb
                   << " pd = " << pd
                   << " qb = " << qb
                   << " qd = " << qd
                   << " " << flush;
#endif
              cov_val = (pbqd + qbpd) / 4.0 + pbd + qbd - (pb*pd) - (qb*qd);
            }
            else if( b == d )
            {
              double pbqd = my_cov_matrix(bp1*m+a, bp2*m+c);
              double qbpd = my_cov_matrix(bp2*m+a, bp1*m+c);
              double pbd  = get_3p_value(bp1, a, c);
              double qbd  = get_3p_value(bp2, a, c);
              double pb   = my_kinships_2p(bp1, a);
              double pd   = my_kinships_2p(bp1, c);
              double qb   = my_kinships_2p(bp2, a);
              double qd   = my_kinships_2p(bp2, c);
#if 0
              cout << "case ab,cb" << flush;
              cout << "bp1*m+a = " << bp1*m+a
                   << "bp2*m+c = " << bp2*m+c
                   << "bp2*m+a = " << bp2*m+a
                   << "bp1*m+c = " << bp1*m+c << flush;
              cout << " pbqd = " << pbqd
                   << " qbpd = " << qbpd
                   << " pbd = " << pbd
                   << " qbd = " << qbd
                   << " pb = " << pb
                   << " pd = " << pd
                   << " qb = " << qb
                   << " qd = " << qd
                   << " " << flush;
#endif
              cov_val = ((pbqd+qbpd) / 4.0) + pbd + qbd - (pb*pd) - (qb*qd);
            }
            else // a, b, c, d
            {
#if 0
              cout << "case ab, cd " << flush;
#endif
              if( a > c )
              {
                double pbcd = my_cov_matrix(ap1*m+b, c*m+d);
                double qbcd = my_cov_matrix(ap2*m+b, c*m+d);

                cov_val = (pbcd + qbcd) / 2.0;
              }
              else
              {
                double pbcd = my_cov_matrix(cp1*m+d, a*m+b);
                double qbcd = my_cov_matrix(cp2*m+d, a*m+b);

                cov_val = (pbcd + qbcd) / 2.0;
              }
            }
          }

          if( fabs(cov_val) < 1.0e-15 )
            cov_val = 0.0;
#if 0
//          if( cov_val )
//          {
//            cout << a << "(" << mid1.name() << ") ";
//            cout << b << "(" << mid2.name() << ") ";
//            cout << c << "(" << mid3.name() << ") ";
//            cout << d << "(" << mid4.name() << ") ";
            cout << " = " << cov_val << endl;
//          }
#endif
            my_cov_matrix(a*m + b, c*m + d) = my_cov_matrix(a*m + b, d*m + c)
          = my_cov_matrix(b*m + a, c*m + d) = my_cov_matrix(b*m + a, d*m + c) = cov_val;
        }
      }
    }
  }
#if 0
  // print
  cout << names_ab.str() << endl;
  print_trimatrix(my_cov_matrix, cout, "my_cov_matrix");
#endif

  return;
}

double
covariance_calculator::get_3p_value(int a, int b, int c)
{
  int is[] = {a, b, c};
  vector<int> ivs(is, is+3);
  sort(ivs.begin(), ivs.end());

  double val = my_kinships_3p[ivs[2]](ivs[1],ivs[0]);

  return val;
}

double
covariance_calculator::get_4p_value(int a, int b, int c, int d)
{
  int is[] = {a, b, c, d};
  vector<int> ivs(is, is+4);
  sort(ivs.begin(), ivs.end());

  double val = my_kinships_4p[ivs[3]][ivs[2]](ivs[1],ivs[0]);

  return val;
}

bool
covariance_calculator::is_constant_pair(mem_pointer mid1, mem_pointer mid2) const
{
  // self pair
  if( mid1 == mid2 )
    return true;

  // both founders, unrelated
  if( !(mid1->parent1() || mid2->parent1()) )
    return true;

  // parent-offspring pair
  if(    mid1->parent1() == mid2 || mid1->parent2() == mid2
      || mid2->parent1() == mid1 || mid2->parent2() == mid1 )
    return true;

  // spouse pair, assumed to be unrelated
  FPED::MateConstIterator s1 = mid1->mate_begin();
  for( ; s1 != mid1->mate_end(); ++s1 )
  {
    if( &(s1->mate()) == mid2 )
      return true;
  }

  // mid1 is founder, if mid1 not an ancestor of mid2, unreated
  if( !mid1->parent1() )
  {
    if( !(is_ancestor(mid1,  mid2->parent1()) || is_ancestor(mid1,  mid2->parent2())) )
      return true;
  }

  // mid2 is founder, if mid2 not an ancestor of mid1, unrelated
  if( !mid2->parent1() )
  {
    if( !(is_ancestor(mid1,  mid2->parent1()) || is_ancestor(mid1,  mid2->parent2())) )
      return true;
  }

  return false;
}

bool
covariance_calculator::is_ancestor(mem_pointer mid1, mem_pointer p) const
{
  if( mid1 == p )
    return true;
  else if( !p )
    return false;

  if( is_ancestor(mid1, p->parent1()) )
    return true;
  
  return is_ancestor(mid1, p->parent2());
}

} // end of namespace RELPAL
} // end of namespace SAGE

