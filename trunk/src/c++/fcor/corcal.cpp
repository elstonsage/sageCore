//****************************************************************************
//* File:      corcal.cpp                                                    *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   0.0 Initial implementation                         yjs Oct 99 *
//*            1.0 Revised to implement MainCorCal                yjs May 01 *
//*                                                                          *
//* Notes:     This source file computes  correlations for each type of      *
//*            pairset in PairSetData.                                       *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "fcor/corcal.h"

namespace SAGE {
namespace FCOR {

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     CorrelationCal                                               ~
// ~                                                                         ~
// ~ Purpose:   Calculate correlations of pairsetdata.                       ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CorrelationCal::CorrelationCal()
{
  my_valid        = false;
  my_parser       = NULL;
  my_pairset      = NULL;
  my_pairset_info = NULL;
}

CorrelationCal::CorrelationCal(const FcorParser& p)
{
  my_valid        = false;
  my_parser       = &p;
  my_pairset      = NULL;
  my_pairset_info = NULL;
}

void
CorrelationCal::compute_correlations(const pairset_vector&         pairsets,
                                     const pairset_info_vector&    pinfos,
                                           pairset_result_vector&  results)
{
  my_pairset      = &pairsets;
  my_pairset_info = &pinfos;

  compute_correlations(results);
}

void
CorrelationCal::compute_correlations(const pairset_vector&         pairsets,
                                     const pairset_info_vector&    pinfos,
                                     const weight_matrix_vector&   weights,
                                           pairset_result_vector&  results)
{
  my_pairset      = &pairsets;
  my_pairset_info = &pinfos;

  compute_correlations(weights, results);
}

void
CorrelationCal::compute_correlations(pairset_result_vector& results)
{
  if( my_parser == NULL || my_pairset == NULL)
    return;

  bool se = my_parser->get_analysis_options().standard_error;
  size_t trait_count = my_parser->get_trait_count();

  my_corinfo.resize(0);
  results.resize(my_pairset->size(), pairset_result(trait_count, WEIGHT_COUNT + 1, se));

  for( size_t r = 0; r < my_pairset->size(); ++r )
  {
    // For multiple weights
    //
    corinfo_by_weight_type cor;
    cor.resize(WEIGHT_COUNT);
    for( size_t i = 0; i < cor.size(); ++i )
      cor[i].resize(trait_count * 2);

    if( (*my_pairset_info)[r].total_pair_count > 2 )
    {
      compute_pairset_correlation((*my_pairset)[r], cor);
      set_correlation(results[r], cor);
    }

    my_corinfo.push_back(cor);
  }

  build_pair_to_corinfo_map();

  my_valid  = true;
}

void
CorrelationCal::compute_correlations(const weight_matrix_vector&   weights,
                                           pairset_result_vector&  results)
{
  if( my_parser == NULL || my_pairset == NULL)
    return;

  const vector<name_index_type>& trait_info = my_parser->get_trait_list();
  size_t t_count = trait_info.size();

  my_corinfo.resize(0);

  for( size_t r = 0; r < my_pairset->size(); ++r )
  {
    corinfo_by_weight_type cor;
    cor.resize(1);  // only one weight.
    cor[0].resize(t_count * 2);

    if( (*my_pairset_info)[r].total_pair_count > 2 )
    {
      const pairset_by_pedigree_type& pairset = (*my_pairset)[r];

      //compute_pairset_correlation((*my_pairset)[r], cor);
      for( size_t p = 0; p < pairset.size(); ++p )
      {
        const weight_matrix& ped_weight = weights[r][p];

        if( !ped_weight.first.rows() )
          continue;

        // Build triangle weight matrix.
        TriangleMatrix<double> tri_weight;
        tri_weight.resize(t_count * 2);

        //cout << "tri weight" << endl;
        for( size_t t1 = 0; t1 < tri_weight.size(); ++t1 )
        {
          for( size_t t2 = t1; t2 < tri_weight.size(); ++t2 )
          {
            tri_weight(t1, t2) = ped_weight.first(t1%t_count, t2%t_count);
            //cout << " " << tri_weight(t1, t2);
          }
          //cout << endl;
        }
        //cout << endl;

        for( size_t i = 0; i < pairset[p].size(); ++i )
        {
          vector<double> traits(t_count * 2);

          const_pedigree_member_pair m_pair = pairset[p][i].member_pair;

          for( size_t t = 0; t < t_count; ++t )
          {
            traits[t]
            = m_pair.first ->pedigree()->info().trait(m_pair.first ->index(), trait_info[t].second);
          }

          for( size_t t = 0; t < t_count; ++t )
          {
            traits[t + t_count]
            = m_pair.second->pedigree()->info().trait(m_pair.second->index(), trait_info[t].second);
          }

          cor[0].add(traits, tri_weight);
        }
      }

      //set_correlation(results[r], cor);
      for( size_t t1 = 0; t1 < t_count; ++t1 )
        for( size_t t2 = 0; t2 < t_count; ++t2 )
          results[r].corr(t1, t2).correlation[WEIGHT_COUNT + 1] = cor[0].correlation(1, 0);
    }

    my_corinfo.push_back(cor);
  }

  build_pair_to_corinfo_map();

  my_valid  = true;
}


void
CorrelationCal::set_correlation(pairset_result& pr, const vector<CorrelationInfo>& cor)
{
  size_t trait_count = my_parser->get_trait_count();

  for( size_t t1 = 0; t1 < trait_count; ++t1 )
    for( size_t t2 = 0; t2 < trait_count; ++t2 )
      for( size_t w = 0; w < WEIGHT_COUNT; ++w )
      {
        pr.corr(t1, t2).count          = cor[w].count(t1 + trait_count, t2);
        pr.corr(t1, t2).correlation[w] = cor[w].correlation(t1 + trait_count, t2);
      }
}

void
CorrelationCal::build_pair_to_corinfo_map()
{
  size_t ped_count = my_parser->get_multi_pedigree()->pedigree_count();
  my_corinfo_map.resize(ped_count);

  for( size_t p = 0; p < ped_count; ++p )
    my_corinfo_map[p].resize(my_parser->get_multi_pedigree()->pedigree_index(p).member_count());

  for( size_t r = 0; r < my_pairset->size(); ++r )
  {
    pairset_by_pedigree_const_iterator ped_begin = (*my_pairset)[r].begin();
    pairset_by_pedigree_const_iterator ped_end   = (*my_pairset)[r].end();

    for( ; ped_begin != ped_end; ++ped_begin )
    {
      pairset_const_iterator pair_begin = ped_begin->begin();
      pairset_const_iterator pair_end   = ped_begin->end();

      for( ; pair_begin != pair_end; ++pair_begin )
      {
        size_t p = pair_begin->member_pair.first->pedigree()->index();

        size_t i = pair_begin->member_pair.first->index();
        size_t j = pair_begin->member_pair.second->index();

        my_corinfo_map[p](i, j) = r;
      }
    }
  }
}

void
CorrelationCal::compute_pairset_correlation(const pairset_by_pedigree_type& pairset,
                                                  corinfo_by_weight_type&   corinfo) const
{
  const vector<name_index_type>& trait_info = my_parser->get_trait_list();

  size_t t_count = trait_info.size();

  for( size_t p = 0; p < pairset.size(); ++p )
  {
    for( size_t i = 0; i < pairset[p].size(); ++i )
    {
      vector<double> traits(t_count * 2);

      const_pedigree_member_pair m_pair = pairset[p][i].member_pair;

      for( size_t t = 0; t < t_count; ++t )
      {
        traits[t]
        = m_pair.first ->pedigree()->info().trait(m_pair.first ->index(), trait_info[t].second);
      }

      for( size_t t = 0; t < t_count; ++t )
      {
        traits[t + t_count]
        = m_pair.second->pedigree()->info().trait(m_pair.second->index(), trait_info[t].second);
      }

      for( size_t w = 0; w < WEIGHT_COUNT; ++w )
        corinfo[w].add(traits, pairset[p][i].pair_weight[w]);

#if 0
  for( size_t t = 0; t < traits.size(); ++t )
    cout << setw(10) << traits[t] << ", ";
  cout << pairset[p][i].pair_weight[PAIR_WISE] << ", "
       << pairset[p][i].pair_weight[UNIFORM] << ", "
       << pairset[p][i].pair_weight[MEAN] << endl;
#endif
    }
  }

  return;
}

void
CorrelationCal::view_corinfo(ostream& out, bool see_class) const
{
  size_t t_count = my_parser->get_trait_count();

  for( size_t r = 0; r < my_corinfo.size(); ++r )
  {
    out << (*my_pairset_info)[r].gname << endl;

    const vector<CorrelationInfo>& cor = my_corinfo[r];

    for( size_t w = 0; w < cor.size(); ++w )
    {
      out << "w = " << w << endl;

      size_t t1 = 0;
      if( !see_class )
        t1 = t_count;

      for( ; t1 < cor[w].size(); ++t1 )
      {
        for( size_t t2 = 0; t2 <= t1; ++t2 )
        {
          if( !see_class && t2 >= t_count )
            break;

          out << setw(5) << cor[w].count(t1, t2) << " ";
          out << setprecision(7);
          out.setf(ios_base::fixed,ios_base::floatfield);
          out << setw(10) << cor[w].correlation(t1, t2) << " : "
              << setw(10) << cor[w].covariance(t1, t2) << "  ";
        }
        out << endl;
      }
      out << endl;
    }
    out << endl;
  }

  return;
}

// end of CorrelationCal Implementation

} // end of namespace FCOR
} // end of namespace SAGE
