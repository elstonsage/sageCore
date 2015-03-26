// ---------------------------------------------------------------------------
// Inline Implementation of SubCalBaseM
// ---------------------------------------------------------------------------

inline
size_t
SubCalBaseM::pedigree_count() const
{
  return my_pedigree_count;
}

inline
size_t
SubCalBaseM::trait_count() const
{
  return my_trait_count;
}

inline
size_t
SubCalBaseM::pair_count(const pairset_by_pedigree_type* pairset, size_t ped) const
{
  return (*pairset)[ped].size();
}

inline                                                                                  
void
SubCalBaseM::pedigree_weight(const vector< Matrix2D< vector<double> > >& weight,
                             Matrix2D< vector<internal_real_type> >&     sum_ped) const
{                                                                                       
  for( size_t t1 = 0; t1 < sum_ped.rows(); ++t1 )                                  
    for( size_t t2 = 0; t2 < sum_ped.cols(); ++t2 )                                  
      for( size_t i = 0; i < weight.size(); ++i )
        for( size_t w = 0; w < weight[i](t1, t2).size(); ++w )
          sum_ped(t1, t2)[w] += weight[i](t1, t2)[w];                                                
}

inline
vector< pair<double, size_t> >
SubCalBaseM::get_corinfo(const pairset_by_pedigree_type* pairset1,
                         const pairset_by_pedigree_type* pairset2,
                         size_t i, size_t j, size_t p, pair_type p_t,
                         size_t t1, size_t t2) const

{
  size_t p1   = size_t(-1);
  size_t p2   = size_t(-1);
  size_t mem1 = size_t(-1);
  size_t mem2 = size_t(-1);

  if( p_t == XX || p_t == XY )
  {
    p1   = (*pairset1)[p][i].member_pair.second->pedigree()->index();
    mem1 = (*pairset1)[p][i].member_pair.second->index();
  }

  if( p_t == YX || p_t == YY )
  {
    p1   = (*pairset1)[p][i].member_pair.first->pedigree()->index();
    mem1 = (*pairset1)[p][i].member_pair.first->index();
  }

  if( p_t == UX || p_t == UY )
  {
    p1   = (*pairset1)[p][i].member_pair.second->pedigree()->index();
    mem1 = (*pairset1)[p][i].member_pair.second->index();
  }

  if( p_t == VX || p_t == VY )
  {
    p1   = (*pairset1)[p][i].member_pair.first->pedigree()->index();
    mem1 = (*pairset1)[p][i].member_pair.first->index();
  }

  if( p_t == XX || p_t == YX || p_t == UX || p_t == VX )
  {
    p2   = (*pairset2)[p][j].member_pair.second->pedigree()->index();
    mem2 = (*pairset2)[p][j].member_pair.second->index();
  }

  if( p_t == XY || p_t == YY || p_t == UY || p_t == VY )
  {
    p2   = (*pairset2)[p][j].member_pair.first->pedigree()->index();
    mem2 = (*pairset2)[p][j].member_pair.first->index();
  }

  vector< pair<double, size_t> > cor_cnt_pair;

  if( p1 != p2 )
  {
    for( size_t w = 0; w < my_weight_count; ++w )
      cor_cnt_pair.push_back(make_pair(std::numeric_limits<double>::quiet_NaN(), size_t(-1)));

    return cor_cnt_pair;
  }

  size_t weight_start = 0;
  size_t weight_end   = my_weight_count;

//  if( my_weight_count < WEIGHT_COUNT )
//  {
//    weight_start = 0 + WEIGHT_COUNT;
//    weight_end   = my_weight_count + WEIGHT_COUNT;
//  }

  size_t cor_index = my_correlation->get_corinfo_map()[p1](mem1, mem2); 

  for( size_t w = weight_start; w < weight_end; ++w )
  {
    double pair_cor   = my_correlation->get_corinfo()[cor_index][w].correlation(t1, t2);
    size_t pair_count = my_correlation->get_corinfo()[cor_index][w].count(t1, t2);

      if( (*my_correlation->get_pairset_info())[cor_index].type == INTRA )
        pair_count /= 2;

    cor_cnt_pair.push_back(make_pair(pair_cor, pair_count));
  }

  return cor_cnt_pair;
}

// end of SubCalBaseM Implementation

