//==========================================================================
//  File:      ibd_analysis.cpp
//
//  Author:
//
//  History:   Initial implementation.
//             Updated to new libraries.                         yjs Jul. 03
//
//  Notes:
//
//  Copyright (c) 2003 R.C. Elston
//    All Rights Reserved
//==========================================================================

#include "lvec/dgraph.h"
#include "ibd/ibd_analysis.h"
#include "numerics/corinfo.h"

namespace SAGE
{

ibd_analysis::ibd_analysis()
            : my_valid(false)
{}

bool ibd_analysis::compute(const lvector& mkr, const meiosis_map& _mm, bool ibd_state_out)
{
  if( !_mm.built() || !mkr.is_valid() )
    return my_valid = false;

  my_meiosis_map = _mm;

  size_t member_count = my_meiosis_map.get_subpedigree()->member_count();
  size_t pair_count   = member_count * member_count;
  size_t r_pair_count = (member_count*(member_count-1))/2;

  // Allocate our data structures
  my_ibd_sharing.resize(0);
  my_ibd_sharing.resize(pair_count, 0.0);

  my_sib_ibd_sharing.resize(0);
  my_sib_ibd_sharing.resize(pair_count, 0.0);

  my_ibd_state.resize(0);

  descent_graph dg(0, my_meiosis_map);

  double total = mkr.total();

#if 0
  cout << "start ibd_analysis::compute()..." << endl
       << "total = " << total << endl;
  dg.print_allele_vector();
#endif

  if( total == 0.0 )
    return my_valid = false;

  for( lvector::const_iterator i = mkr.begin(); i != mkr.end(); ++i )
  {
    if( *i == 0.0 )
      continue;

    lvector::equivalence_class e = i - mkr.begin();

    dg.move_to(e, my_meiosis_map);

#if 0
  cout << "equivalence_class e = " << e << " : " << mkr[e] << endl;
  dg.print_allele_vector();
#endif

    vector<size_t>   a_state(r_pair_count);

    for( size_t j = 0, s = 0; j < member_count-1; ++j )
    {
      for( size_t k = j+1; k < member_count; ++k, ++s )
      {
        size_t sharing = dg.sharing(j,k);
#if 0
  cout << "j = " << j << ":" << my_meiosis_map.member(j)->name()
       << ", k = " << k << ":" << my_meiosis_map.member(k)->name()
       << ", dg.sharing(j,k) = " << sharing << endl;
#endif
        if(    my_meiosis_map.is_x_linked()
            && is_brother_brother(my_meiosis_map, j, k) )
          --sharing;

        switch( sharing )
        {
          case 0 : my_prob_share(j, k, 0) += *i; a_state[s] = sharing; break;
          case 2 : my_prob_share(j, k, 2) += *i; a_state[s] = sharing; break;

          // Added for m or paternal bit split. - yjs Jun. 2002
          case 1 : if( is_sib(my_meiosis_map, j, k) )
                   {
                     if( dg.mother_sharing(j,k) )
                     {
                       my_sib_prob_share(j, k, 0) += *i;
                       a_state[s] = 3;
                     }

                     if( dg.father_sharing(j,k) )
                       if(    !my_meiosis_map.is_x_linked()
                           || !is_brother_brother(my_meiosis_map, j, k) )
                       {
                         my_sib_prob_share(j, k, 2) += *i;
                         a_state[s] = 4;
                       }
#if 0
  cout << "  dg.mother_sharing(j,k) = " << dg.mother_sharing(j,k)
       << ", dg.father_sharing(j,k) = " << dg.father_sharing(j,k) << endl;
#endif
                   }
                   else if( is_hsib(my_meiosis_map, j, k) )
                   {
                     if( is_maternal_hsib(my_meiosis_map, j, k) && dg.mother_sharing(j,k) )
                     {
                       my_sib_prob_share(j, k, 0) += *i;
                       a_state[s] = 3;
                     }

                     if( is_paternal_hsib(my_meiosis_map, j, k) && dg.father_sharing(j,k) )
                       if(    !my_meiosis_map.is_x_linked()
                           || !is_brother_brother(my_meiosis_map, j, k) )
                       {
                         my_sib_prob_share(j, k, 2) += *i;
                         a_state[s] = 4;
                       }
#if 0
  cout << "  dg.mother_sharing(j,k) = " << dg.mother_sharing(j,k)
       << ", dg.father_sharing(j,k) = " << dg.father_sharing(j,k) << endl;
#endif
                   }
                   else
                     a_state[s] = sharing;

                   break;

          default: break;
        }
      }
    }

    if( ibd_state_out )
    {
      my_ibd_state.push_back(make_pair((*i)/total, a_state));

#if 0
      cout << setw(5) << doub2str((*i), 5) << "   ";
      for( size_t s = 0; s < a_state.size(); ++s )
        cout << setw(5) << a_state[s] << " ";
      cout << endl;
#endif
    }
  }

  for( size_t j = 0; j < member_count-1; ++j )
  {
    for( size_t k = j+1; k < member_count; ++k )
    {
      my_prob_share(j, k, 0) /= total;
      my_prob_share(j, k, 2) /= total;

      my_sib_prob_share(j, k, 0) /= total;
      my_sib_prob_share(j, k, 2) /= total;
    } 
  }

#if 0
  CorrelationInfo  mean_share_info;

  mean_share_info.resize(r_pair_count);

  cout << "      ";
  for( size_t j = 0; j < member_count-1; ++j )
    for( size_t k = j+1; k < member_count; ++k )
      cout << "(" << setw(1) << j << "," << setw(1) << k << ") ";
  cout << endl;

  for( lvector::const_iterator i = mkr.begin(); i != mkr.end(); ++i )
  {
    if( *i == 0.0 ) continue;

    lvector::equivalence_class e = i - mkr.begin();
    dg.move_to(e, my_meiosis_map);

    vector<double>   a_state(r_pair_count);

    for( size_t j = 0, s = 0; j < member_count-1; ++j )
    {
      for( size_t k = j+1; k < member_count; ++k, ++s )
      {
        int sharing = dg.sharing(j,k);

        assert( j != k || sharing == 2 );
        assert( sharing == dg.sharing(k,j) );   

        a_state[s] = sharing;
      }
    }

    double weight = *i;

    cout << setw(5) << doub2str(weight, 5) << "   ";
    for( size_t s = 0; s < a_state.size(); ++s )
      cout << setw(5) << a_state[s] << " ";
    cout << endl;

    mean_share_info.add(a_state, weight/total);
  }

  cout << "***Reduced Cov(phi)" << endl;
  cout << "      ";
  for( size_t i = 0; i < member_count-1; ++i )
    for( size_t j = i+1; j < member_count; ++j )
      cout << "(" << setw(1) << i << "," << setw(1) << j << ") ";
  cout << endl;

  for( size_t i = 0, r = 0; i < member_count-1; ++i )
    for( size_t j = i+1; j < member_count ; ++j, ++r )
    {
      cout << "(" << setw(1) << i << "," << setw(1) << j << ") ";

      for(size_t k = 0, c = 0; k < member_count-1; ++k)
        for(size_t l = k+1; l < member_count; ++l, ++c)
          cout << setw(5) << doub2str(mean_share_info.covariance(r,c), 5) << " ";
      cout << endl;
    }
  cout << endl;
#endif

  return my_valid = true;
}

} // end of namespace SAGE
