#include "palbase/pair_filter.h"

namespace SAGE    {
namespace PALBASE {

bool
pair_filter::valid(const rel_pair& sp, bool filter_sex, double cp) const
{
#if 0
  cout << "Valid ? "
       << sp.rels().pair.first->pedigree()->name() << " ("
       << sp.rels().pair.first->name()             << ","
       << sp.rels().pair.second->name()            << ")"
       << endl;
#endif

  if( !filter_any() ) return true;

  if( filter_sex && !valid_sex_info(sp) )
    return false;

  if(    !valid_marker(sp) || !valid_trait(sp, cp)
      || !valid_subset(sp) || !valid_covariate(sp) )
    return false;

  return true;
}

bool
pair_filter::valid_marker(const rel_pair& sp) const
{
  if( !my_markers.size() )
    return true;

  double prior_f0 = sp.prior_prob_share(0);
  double prior_f2 = sp.prior_prob_share(2);
  bool   informative = false;
  bool   informative_filter = false;

  marker_set::const_iterator i;
  for( i = my_markers.begin(); i != my_markers.end(); ++i )
  {
    const double f0 = sp.prob_share(i->marker, 0);
    const double f2 = sp.prob_share(i->marker, 2);

    if( !finite(f0) || !finite(f2) )
      return false;

    if( i->tolerance < INF )
    {
      informative_filter = true;

      double df0 = fabs( prior_f0 - f0 );
      double df2 = fabs( prior_f2 - f2 );
      double deltaf  = std::max(df0, df2);

      if( deltaf > i->tolerance)
        informative = true;
    }
  }
  
  if( informative_filter && !informative )
    return false;

  return true;
}

bool
pair_filter::valid_trait(const rel_pair& sp, double cutpoint) const
{
  if( !my_traits.size() )
    return true;

  const size_t ped = sp.pedigree_number();
  const size_t i1  = sp.rels().pair.first->index();
  const size_t i2  = sp.rels().pair.second->index();

  trait_set::const_iterator i;
  for( i = my_traits.begin(); i != my_traits.end(); ++i )
  {
    const size_t t = i->trait;
    double trait1 = sp.pair_data()->trait(ped, i1, t);
    double trait2 = sp.pair_data()->trait(ped, i2, t);

#if 0
    cout << "Valid trait ? "
         << sp.rels().first->pedigree()->name() << " ("
         << sp.rels().first->name()             << ","
         << sp.rels().second->name()            << " = "
         << trait1 << ", " << trait2            << "), min="
         << i->min << ", max=" << i->max << endl;
#endif

    if( !finite(trait1) || !finite(trait2) )
      return false;

    if(    trait1 < i->min || trait1 > i->max 
        || trait2 < i->min || trait2 > i->max )
      return false;

    const RPED::RefTraitInfo& info = sp.pair_data()->fped_info().trait_info(t);

    if( info.type() == RPED::RefTraitInfo::binary_trait )
    {
      size_t a = (size_t) (trait1 + trait2);
      if(    trait1 < 0 || trait1 > 1 
          || trait2 < 0 || trait2 > 1 
          || a > 2      || !i->affection[a] )
        return false;
    }  
    else
    {
      size_t aff = 0;
      if( trait1 > cutpoint )
        ++aff;
      if( trait2 > cutpoint )
        ++aff;
      if( !i->affection[aff] )
        return false;      
    }
  }

  return true;
}

bool
pair_filter::valid_subset(const rel_pair& sp) const
{
  if( !my_subsets.size() )
    return true;

  const size_t ped = sp.pedigree_number();
  const size_t i1  = sp.rels().pair.first->index();
  const size_t i2  = sp.rels().pair.second->index();

  trait_set::const_iterator i;
  for( i = my_subsets.begin(); i != my_subsets.end(); ++i )
  {
    const size_t t = i->trait;
    double trait1 = sp.pair_data()->trait(ped, i1, t);
    double trait2 = sp.pair_data()->trait(ped, i2, t);

#if 0
    cout << "Valid subset ? "
         << sp.rels().first->pedigree()->name() << " ("
         << sp.rels().first->name()             << ","
         << sp.rels().second->name()            << " = "
         << trait1 << ", " << trait2            << ") ";
#endif

    if( !finite(trait1) || !finite(trait2) )
    {
#if 0
      cout << "No1" << endl;
#endif
      return false;
    }

    if(    trait1 < i->min || trait1 > i->max 
        || trait2 < i->min || trait2 > i->max )
    {
#if 0
      cout << "No2" << endl;
#endif
      return false;
    }

    const RPED::RefTraitInfo& info = sp.pair_data()->fped_info().trait_info(t);

    if( info.type() == RPED::RefTraitInfo::binary_trait )
    {
      size_t a = (size_t) (trait1 + trait2);
      if(    trait1 < 0 || trait1 > 1 
          || trait2 < 0 || trait2 > 1 
          || a > 2      || !i->affection[a] )
      {
#if 0
        cout << "No3" << endl;
#endif
        return false;
      }
    }  
    else
    {
      size_t aff = 0;
      if(trait1 > 0)
        ++aff;
      if(trait2 > 0)
        ++aff;
      if( !i->affection[aff] )
      {
#if 0
        cout << "No4" << endl;
#endif
        return false;
      }
    }
  }

#if 0
  cout << "Yes!" << endl;
#endif

  return true;
}

bool
pair_filter::valid_covariate(const rel_pair& sp) const
{
  if( !my_covariates.size() && !my_pair_covariates.size() )
    return true;

  const size_t ped = sp.pedigree_number();
  const size_t i1  = sp.rels().pair.first->index();
  const size_t i2  = sp.rels().pair.second->index();

  trait_set::const_iterator i;
  for(i = my_covariates.begin(); i != my_covariates.end(); ++i)
  {
    const size_t t = i->trait;
    double trait1 = sp.pair_data()->trait(ped, i1, t);
    double trait2 = sp.pair_data()->trait(ped, i2, t);

    if( !finite(trait1) || !finite(trait2) )
      return false;

    if(    trait1 < i->min || trait1 > i->max 
        || trait2 < i->min || trait2 > i->max )
      return false;
  }

  for( i = my_pair_covariates.begin(); i != my_pair_covariates.end(); ++i )
  {
    const size_t t = i->trait;
    double trait = sp.pair_data()->get_pair_covariate(sp.pair_number(), t);

    if( !finite(trait) )
      return false;

    if( trait < i->min || trait > i->max )
      return false;
  }

  return true;
}

bool
pair_filter::valid_sex_info(const rel_pair& sp) const
{
  if( sp.rels().pair.first->is_sex_unknown() )
    return false;

  if( sp.rels().pair.second->is_sex_unknown() )
    return false;

  return true;
}

bool
pair_filter::is_concordant_aff_pair(const rel_pair& sp, double cutpoint) const
{
  const size_t ped = sp.pedigree_number();
  const size_t i1  = sp.rels().pair.first->index();
  const size_t i2  = sp.rels().pair.second->index();

  trait_set::const_iterator i;
  for( i = my_traits.begin(); i != my_traits.end(); ++i )
  {
    const size_t t = i->trait;
    double trait1 = sp.pair_data()->trait(ped, i1, t);
    double trait2 = sp.pair_data()->trait(ped, i2, t);

    if( !finite(trait1) || !finite(trait2) )
      return false;

    if(    trait1 < i->min || trait1 > i->max 
        || trait2 < i->min || trait2 > i->max )
      return false;

    bool is_con = false;

    bool cv1 = (trait1 > cutpoint) ? 1 : 0;
    bool cv2 = (trait2 > cutpoint) ? 1 : 0;

    if( !cv1 || !cv2 )
      return false;

    is_con = (cv1 == cv2) ? true : false;

    if( !is_con )
      return false;
  }

  return true;
}

bool
pair_filter::is_concordant_unaff_pair(const rel_pair& sp, double cutpoint) const
{
  const size_t ped = sp.pedigree_number();
  const size_t i1  = sp.rels().pair.first->index();
  const size_t i2  = sp.rels().pair.second->index();

  trait_set::const_iterator i;
  for( i = my_traits.begin(); i != my_traits.end(); ++i )
  {
    const size_t t = i->trait;

    double trait1 = sp.pair_data()->trait(ped, i1, t);
    double trait2 = sp.pair_data()->trait(ped, i2, t);

    if( !finite(trait1) || !finite(trait2) )
      return false;

    if(    trait1 < i->min || trait1 > i->max 
        || trait2 < i->min || trait2 > i->max )
      return false;

    bool cv1 = (trait1 > cutpoint) ? 1 : 0;
    bool cv2 = (trait2 > cutpoint) ? 1 : 0;

    if( cv1 || cv2 )
      return false;
#if 0
    cout << sp.rels().first->pedigree()->name() << " ("
         << sp.rels().first->name()             << ","
         << sp.rels().second->name()            << " = "
         << trait1 << ", " << trait2            << ", min="
         << i->min << ", max=" << i->max << ", name = " << info.name() << endl;
#endif
  }

  return true;
}

bool
pair_filter::is_discordant_pair(const rel_pair& sp, double cutpoint) const
{
  const size_t ped = sp.pedigree_number();
  const size_t i1  = sp.rels().pair.first->index();
  const size_t i2  = sp.rels().pair.second->index();

  trait_set::const_iterator i;
  for( i = my_traits.begin(); i != my_traits.end(); ++i )
  {
    const size_t t = i->trait;

    double trait1 = sp.pair_data()->trait(ped, i1, t);
    double trait2 = sp.pair_data()->trait(ped, i2, t);

    if( !finite(trait1) || !finite(trait2) )
      return false;

    if(    trait1 < i->min || trait1 > i->max 
        || trait2 < i->min || trait2 > i->max )
      return false;

    bool cv1 = (trait1 > cutpoint) ? 1 : 0;
    bool cv2 = (trait2 > cutpoint) ? 1 : 0;

    if( cv1 == cv2 )
      return false;
#if 0
    cout << sp.rels().first->pedigree()->name() << " ("
         << sp.rels().first->name()             << ","
         << sp.rels().second->name()            << " = "
         << trait1 << ", " << trait2            << ", min="
         << i->min << ", max=" << i->max << ", name = " << info.name() << endl;
#endif
  }

  return true;
}

void
pair_filter::add_marker(size_t m, double tolerance)
{
  marker_iterator i = my_markers.find(m);
  if( i != my_markers.end() )
  {
    if( i->tolerance > tolerance )
    {
      i->tolerance = tolerance;
    }
  }
  else
    my_markers.insert( filter_marker(m, tolerance) );
}

void
pair_filter::add_trait(size_t t, double min, double max)
{
  trait_iterator i = my_traits.find(t);
  if( i != my_traits.end() )
  {
    if( i->min > min )
    {
      i->min = min;
    } 
    if( i->max < max )
    {
      i->max = max;
    } 
  }
  else
    my_traits.insert( filter_trait(t, min, max) );
}

void
pair_filter::add_trait(size_t t, size_t a, double min, double max)
{
  bool aa[3];
  aa[0] = aa[1] = aa[2] = false;
  aa[a] = true;

  add_trait(t, aa, min, max);
}

void
pair_filter::add_trait(size_t t, const bool* a, double min, double max)
{
  trait_iterator i = my_traits.find(t);
  if( i != my_traits.end() )
  {
    if( i->min > min )
    {
      i->min = min;
    } 
    if( i->max < max )
    {
      i->max = max;
    } 
    for(size_t j = 0; a && j < 3; ++j)
      if( i->affection[j] != a[j] )
      {
        i->affection[j] = a[j];
      }
  }
  else
    my_traits.insert( filter_trait(t, a, min, max) );
}

void
pair_filter::add_subset(size_t t, double min, double max)
{
  trait_iterator i = my_subsets.find(t);
  if( i != my_subsets.end() )
  {
    if( i->min > min )
    {
      i->min = min;
    } 
    if( i->max < max )
    {
      i->max = max;
    } 
  }
  else
    my_subsets.insert( filter_trait(t, min, max) );
}

void
pair_filter::add_subset(size_t t, size_t a, double min, double max)
{
  bool aa[3];
  aa[0] = aa[1] = aa[2] = false;
  aa[a] = true;

  add_subset(t, aa, min, max);
}

void
pair_filter::add_subset(size_t t, const bool* a, double min, double max)
{
  trait_iterator i = my_subsets.find(t);
  if( i != my_subsets.end() )
  {
    if( i->min > min )
    {
      i->min = min;
    } 
    if( i->max < max )
    {
      i->max = max;
    } 
    for(size_t j = 0; a && j < 3; ++j)
      if( i->affection[j] != a[j] )
      {
        i->affection[j] = a[j];
      }
  }
  else
    my_subsets.insert( filter_trait(t, a, min, max) );
}

void
pair_filter::add_covariate(size_t t, double min, double max)
{
  trait_iterator i = my_covariates.find(t);
  if( i != my_covariates.end() )
  {
    if( i->min > min )
    {
      i->min = min;
    } 
    if( i->max < max )
    {
      i->max = max;
    } 
  }
  else
    my_covariates.insert( filter_trait(t, min, max) );
}

void
pair_filter::add_pair_covariate(size_t t, double min, double max)
{
  trait_iterator i = my_pair_covariates.find(t);
  if( i != my_pair_covariates.end() )
  {
    if( i->min > min )
    {
      i->min = min;
    } 
    if( i->max < max )
    {
      i->max = max;
    } 
  }
  else
    my_pair_covariates.insert( filter_trait(t, min, max) );
}

} // end namespace PALBASE
} // end namespace SAGE
