#include "palbase/relative_pairs.h"

namespace SAGE    {
namespace PALBASE {

relative_pairs::relative_pairs()
{
  my_valid_distance_exist = false;
  my_fsib_pair_count = my_hsib_pair_count = 0;
  my_mm_pair_count = my_mf_pair_count = my_ff_pair_count = 0;

  invalidate();
  invalidate_pair_info();
}

relative_pairs::~relative_pairs()
{}

void
relative_pairs::prebuild()
{
  my_prior_ibd.resize( get_fped().pedigree_count() );
  for( size_t p=0; p < get_fped().pedigree_count(); ++p )
  {
    my_prior_ibd[p].compute( get_fped().pedigree_index(p) );
  }
}

void
relative_pairs::build(bool build_sib_cluster, bool build_x_type)
{
  if( !my_pairs.size() )
    return;

  // pair_num in my_pairs, my_relmap, & my_ibd_pairs is all wrong now after sorting.
  //
  sort(my_pairs.begin(), my_pairs.end(), ibd_pair_info_less());

  vector< ibd_probability_info > temp_ibd_probs;
  temp_ibd_probs.resize(my_pairs.size());

  for( size_t i = 0; i < my_pairs.size(); ++i )
  {
    temp_ibd_probs[i] = my_ibd_probs[my_pairs[i].index];

    my_pairs[i].index = i;

    if( find_pair(my_pairs[i].pair.first, my_pairs[i].pair.second) < my_pairs.size() )
      my_relmap[member_pair(my_pairs[i].pair.first, my_pairs[i].pair.second)] = i;
    else if( find_pair(my_pairs[i].pair.second, my_pairs[i].pair.first) < my_pairs.size() )
      my_relmap[member_pair(my_pairs[i].pair.second, my_pairs[i].pair.first)] = i;
    else
      cout << "Impossible";

    if( my_pairs[i].type == pair_generator::SIBSIB )
      ++my_fsib_pair_count;
    else if( my_pairs[i].type == pair_generator::HALFSIB )
      ++my_hsib_pair_count;
    else
      continue;

    if( my_pairs[i].pair_sex == MM )
      ++my_mm_pair_count;
    else if( my_pairs[i].pair_sex == MF )
      ++my_mf_pair_count;
    else if( my_pairs[i].pair_sex == FF )
      ++my_ff_pair_count;
  }

  my_ibd_probs = temp_ibd_probs;

  if( build_sib_cluster )
  {
    my_sibship_cluster.resize(0);

    for( size_t i = 0; i < my_pairs.size(); ++i )
    {
      if(    my_pairs[i].type != pair_generator::SIBSIB
          && my_pairs[i].type != pair_generator::HALFSIB )
        continue;

      mem_pointer dad1 = my_pairs[i].pair.first->parent1();
      mem_pointer mom1 = my_pairs[i].pair.first->parent2();

      mem_pointer dad2 = my_pairs[i].pair.second->parent1();
      mem_pointer mom2 = my_pairs[i].pair.second->parent2();

      bool into_cluster = false;

      for( size_t ci = 0; ci < my_sibship_cluster.size(); ++ci )
      {
        const set<mem_pointer>& p_set = my_sibship_cluster[ci].parent_set;

        if(    is_cluster_parent(dad1, p_set)
            || is_cluster_parent(mom1, p_set)
            || is_cluster_parent(dad2, p_set)
            || is_cluster_parent(mom2, p_set) )
        {
          my_sibship_cluster[ci].parent_set.insert(dad1);
          my_sibship_cluster[ci].parent_set.insert(mom1);
          my_sibship_cluster[ci].parent_set.insert(dad2);
          my_sibship_cluster[ci].parent_set.insert(mom2);

          if( my_pairs[i].type == pair_generator::SIBSIB )
            my_sibship_cluster[ci].fsib_pairs.push_back(i);
          else
            my_sibship_cluster[ci].hsib_pairs.push_back(i);

          into_cluster = true;
        }
      }

      if( !into_cluster )
      {
        set<mem_pointer> p_set;
    
        p_set.insert(dad1);
        p_set.insert(mom1);
        p_set.insert(dad2);
        p_set.insert(dad2);

        sibcluster_info new_cluster;
        new_cluster.parent_set = p_set;

        if( is_fsib(my_pairs[i].pair.first, my_pairs[i].pair.second) )
          new_cluster.fsib_pairs.push_back(i);
        else
          new_cluster.hsib_pairs.push_back(i);

        my_sibship_cluster.push_back(new_cluster);
      }
    }

    for( size_t ci = 0; ci < my_sibship_cluster.size(); ++ci )
    {
      const set<mem_pointer>& p_set   = my_sibship_cluster[ci].parent_set;
      const vector<size_t>&   f_pairs = my_sibship_cluster[ci].fsib_pairs;
      
      if( p_set.size() > 2 && f_pairs.size() > 1 )
      {
        map< id_pair, vector<size_t> >& fmap = my_sibship_cluster[ci].full_sibship_map;
        
        for( size_t j = 0; j < f_pairs.size(); ++j )
        {
          mem_pointer dad = my_pairs[f_pairs[j]].pair.first->parent1();
          mem_pointer mom = my_pairs[f_pairs[j]].pair.first->parent2();
                    
          if( dad->is_female() )
          {
            dad = my_pairs[f_pairs[j]].pair.first->parent2();
            mom = my_pairs[f_pairs[j]].pair.first->parent1();
          }
                                            
          fmap[make_pair(mom, dad)].push_back(f_pairs[j]);
        }
      }
    }
  }

  if( build_x_type )
  {
    for( size_t i = 0; i < my_pairs.size(); ++i )
    {
      my_pairs[i].x_type = get_pair_x_type(my_pairs[i].pair.first, my_pairs[i].pair.second);
    }
  }

  validate();
}

void
relative_pairs::build_from_pedfile(const RPED::RefMultiPedigree& mped, bool build_sib_cluster, bool build_x_type)
{
  for( size_t p = 0; p < mped.pedigree_count(); ++p )
  {
    unsigned int type_wanted = pair_generator::EVERY_MASK;

    const RPED::RefPedigree& rp = mped.pedigree_index(p);

    const FPED::Pedigree*    fp = get_fped().pedigree_find(rp.name());

    pair_generator pair_gen(const_cast<RPED::RefPedigree*>(&rp), type_wanted);

    pair_generator::iterator pi = pair_gen.begin();
    for( ; pi != pair_gen.end(); ++pi )
    {
      pair_generator::relative_pair& rel_pair = *pi;

      mem_pointer member_one = fp->member_find(rel_pair.member_one()->name());
      mem_pointer member_two = fp->member_find(rel_pair.member_two()->name());

      if( !member_one || !member_two )
        continue;

      add_pair(member_one, member_two, rel_pair.type());
    }
  }

  build(build_sib_cluster, build_x_type);
}

size_t
relative_pairs::add_marker(const string& name, double dist, gmodel_type mt)
{
  ibd_marker_info m_info(name, dist, mt);

  my_markers.push_back(m_info);

  if( !my_valid_distance_exist && !SAGE::isnan(dist) )
    my_valid_distance_exist = true;

  return my_markers.size() - 1;
}

size_t
relative_pairs::add_pair(mem_pointer i1, mem_pointer i2, pair_type ppt)
{
  pair_type pt = get_pair_type(i1, i2);
  pair_sex_type s_type = get_pair_sex_type(i1, i2);

  invalidate();

#if 0
  cout << pair_count() << " (" << i1->pedigree()->name() << ") "
       << i1->index() << "(" << i1->name() << "):"
       << i2->index() << "(" << i2->name() << ") " << pt << " " << flush
       << endl;
#endif
#if 0
  if( pt == pair_generator::GRANDP || pt == pair_generator::AVUNC )
  {
    // make sure they are in generation order
    size_t i1_p1 = (i1->parent1())? i1->parent1()->index() : 0;
    size_t i1_p2 = (i1->parent2())? i1->parent2()->index() : 0;
    size_t i2_p1 = (i2->parent1())? i2->parent1()->index() : 0;
    size_t i2_p2 = (i2->parent2())? i2->parent2()->index() : 0;

    size_t mi1 = std::max(i1_p1, i1_p2);
    size_t mi2 = std::max(i2_p1, i2_p2);

    if( mi2 == std::max(mi1, mi2) )
      swap(i1, i2);
  }
  else if( i1->name() > i2->name() )
    swap(i1, i2);

  cout << i1->index() << "(" << i1->name() << "):"
       << i2->index() << "(" << i2->name() << ")" << endl;
#endif

  if( i1->name() > i2->name() )
    swap(i1, i2);

  size_t pid = find_pair(i1,i2);

  if( pid < pair_count() )
  {
    if( my_ibd_probs[pid].f0s.size() != marker_count() )
    {
      for( size_t ii = my_ibd_probs[pid].f0s.size(); ii < marker_count(); ++ii )
      {
        my_ibd_probs[pid].f0s.push_back(QNAN);
        my_ibd_probs[pid].f1mp.push_back(QNAN);
        my_ibd_probs[pid].f2s.push_back(QNAN);
      }
    }

    return pid;
  }

  pid = my_pairs.size();

  my_pairs.push_back( rel_pair_data(i1, i2, pt, s_type, pid) );
  my_relmap.insert( rel_map::value_type(member_pair(i1,i2), pid) );

  my_ibd_probs.push_back(ibd_probability_info(marker_count()));

  return pid;
}

size_t
relative_pairs::find_pair(const mem_pointer& i1, const mem_pointer& i2) const
{
#if 0
  cout << i1->name() << " "
       << i2->name() << " " << flush;
#endif
  if( !i1 || !i2 || i1 == i2 )
    return (size_t)-1;

  member_pair mp1(i1, i2);
  member_pair mp2(i2, i1);

  rel_map::const_iterator j = my_relmap.find( mp1 );

  if( j == my_relmap.end() )
  {
    rel_map::const_iterator jj = my_relmap.find( mp2 );

    if( jj == my_relmap.end() )
      return (size_t)-1;
    else
      return jj->second;
  }

  return j->second;
}

double
relative_pairs::get_average_marker_distance() const
{
  if( !valid_distance_exist() )
    return numeric_limits<double>::quiet_NaN();

  if( !my_markers.size() )
    return numeric_limits<double>::quiet_NaN();

  size_t skip_cnt = 0;
  double dis_sum  = 0.0;
  double prev_dis = 0.0;

  for( size_t i = 0; i < my_markers.size(); ++i )
  {
    if( SAGE::isnan(my_markers[i].distance) )
    {
      skip_cnt += 1;
      continue;
    }

    dis_sum += (my_markers[i].distance - prev_dis) ;
    prev_dis = my_markers[i].distance;
  }

  return dis_sum / double(my_markers.size() - skip_cnt - 1);
}

pair_type
relative_pairs::get_pair_type(const mem_pointer i1, const mem_pointer i2) const
{
  if( is_fsib(i1, i2) )
    return pair_generator::SIBSIB;
  else if( is_hsib(i1, i2) )
    return pair_generator::HALFSIB;
  else if( is_grandp(i1, i2) )
    return pair_generator::GRANDP;
  else if( is_avunc(i1, i2) )
    return pair_generator::AVUNC;
  else if( is_cousin(i1, i2) )
    return pair_generator::COUSIN;

  return pair_generator::EVERY;
}

pair_sex_type
relative_pairs::get_pair_sex_type(const mem_pointer i1, const mem_pointer i2) const
{
  if( is_mm_pair(i1, i2) )
    return MM;
  else if( is_mf_pair(i1, i2) )
    return MF;
  else if( is_ff_pair(i1, i2) )
    return FF;

  return UNKNOWN;
}

pair_x_type
relative_pairs::get_pair_x_type(const mem_pointer i1, const mem_pointer i2) const
{
  if( is_fsib(i1, i2) )
  {
    if( is_mm_pair(i1, i2) )
      return M_M;
    else if( is_mf_pair(i1, i2) )
      return M_F;
    else if( is_ff_pair(i1, i2) )
      return F_F;
    else
      return INVALID;
  }
  else if( is_paternal_hsib(i1, i2) )
  {
    return INVALID;
  }
  else if( is_hsib(i1, i2) )
  {
    if( is_mm_pair(i1, i2) )
      return M_F_M;
    else if( is_mf_pair(i1, i2) )
      return M_F_F;
    else if( is_ff_pair(i1, i2) )
      return F_F_F;
    else
      return INVALID;
  }
  else if( is_p_grandp(i1, i2) )
  {
    return INVALID;
  }
  else if( is_grandp(i1, i2) )
  {
    if( is_mm_pair(i1, i2) )
      return MFM;
    else if( is_mf_pair(i1, i2) )
    {
      if( is_first_ind_grandp(i1, i2) )
      {
        if( i1->is_male() )
          return MFF;
        else if( i1->is_female() )
          return FFM;
      }
      else
      {
        if( i2->is_male() )
          return MFF;
        else if( i2->is_female() )
          return FFM;
      }
    }
    else if( is_ff_pair(i1, i2) )
      return FFF;
    else
      return INVALID;
  }
  else if( is_m_avunc(i1, i2) )
  {
    if( is_mm_pair(i1, i2) )
      return M_FM;
    else if( is_mf_pair(i1, i2) )
    {
      if( is_first_ind_avunc(i1, i2) )
      {
        if( i1->is_male() )
          return M_FF;
        else if( i1->is_female() )
          return F_FM;
      }
      else
      {
        if( i2->is_male() )
          return M_FF;
        else if( i2->is_female() )
          return F_FM;
      }
    }
    else if( is_ff_pair(i1, i2) )
      return F_FF;
    else
      return INVALID;
  }
  else if( is_avunc(i1, i2) )
  {
    if( is_mm_pair(i1, i2) )
      return INVALID;
    else if( is_mf_pair(i1, i2) )
    {
      if( is_first_ind_avunc(i1, i2) )
      {
        if( i1->is_male() )
          return M_MF;
      }
      else
      {
        if( i2->is_male() )
          return M_MF;
      }
    }
    else if( is_ff_pair(i1, i2) )
      return F_MF;
    else
      return INVALID;
  }
  else if( is_p_p_cousin(i1, i2) )
  {
    if( is_ff_pair(i1, i2) )
      return FM_MF;
    else
      return INVALID;
  }
  else if( is_m_m_cousin(i1, i2) )
  {
    if( is_mm_pair(i1, i2) )
      return MF_FM;
    else if( is_mf_pair(i1, i2) )
      return MF_FF;
    else if( is_ff_pair(i1, i2) )
      return FF_FF;
    else
      return INVALID;
  }
  else if( is_cousin(i1, i2) )
  {
    if( is_ff_pair(i1, i2) )
      return FM_FF;
    else if( is_p_m_cousin(i1, i2) && i1->is_male() && i2->is_female() )
      return FM_FM;
    else
      return INVALID;
  }

  return INVALID;
}

//
//------------------------------------------------------------------
//

int
relative_pairs::set_pair_covariate(size_t pn, size_t c,
                                  const string&          value,
                                  const pair_pheno_info& info)
{
  if( c >= pair_covariate_count() || pn >= pair_count() )
    return 3;

  int code = 0;

  double d = str2doub(value);

  if( !finite(d) )
    code = 2;

  string smiss = info.get_string_missing_code();
  double nmiss = info.get_numeric_missing_code();

  if( value == smiss || (finite(nmiss) && d == nmiss) )
  {
    code = 1;
    d = numeric_limits<double>::quiet_NaN();
  }

  my_pair_covariates[c][pn] = d;

  return code;
}

int
relative_pairs::set_pair_weight(size_t pn, size_t w,
                               const string&          value,
                               const pair_pheno_info& info)
{
  if( w >= pair_weight_count() || pn >= pair_count() )
    return 3;

  int code = 0;

  double d = str2doub(value);

  if( !finite(d) )
    code = 2;

  string smiss = info.get_string_missing_code();
  double nmiss = info.get_numeric_missing_code();

  if( value == smiss || (finite(nmiss) && d == nmiss) )
  {
    code = 1;
    d = numeric_limits<double>::quiet_NaN();
  }
  else if( d < 0.0 || d > 1.0 )
  {
    code = 2;
    d = numeric_limits<double>::quiet_NaN();
  }

  my_pair_weights[w][pn] = d;

  return code;
}

void
relative_pairs::dump_pairs(ostream &out) const
{
  out << "========================" << endl
      << "  Relative Pairs Dump  " << endl
      << "========================" << endl << endl;

  out << " Pairs:" << endl;

  out << "fsib_pair count = " << my_fsib_pair_count << endl;
  out << "hsib_pair count = " << my_hsib_pair_count << endl;

  for( size_t i = 0; i < my_pairs.size(); ++i )
  {
    if( i && i % 5 == 0 )
      out << endl;

    out << " " << i << " " << my_pairs[i].index
        << " " << my_pairs[i].pair.first->pedigree()->name()
        << "(" << my_pairs[i].pair.first->name()
        << "," << my_pairs[i].pair.second->name()
        << ")";

    if( my_pairs[i].type == pair_generator::SIBSIB )
      out << "F";
    else if( my_pairs[i].type == pair_generator::HALFSIB )
      out << "H";
    else if( my_pairs[i].type == pair_generator::GRANDP )
      out << "G";
    else if( my_pairs[i].type == pair_generator::AVUNC )
      out << "A";
    else if( my_pairs[i].type == pair_generator::COUSIN )
      out << "C";
    else
      out << "?";
    out << " ";

    if( my_pairs[i].pair_sex == MM )
      out << "MM";
    else if( my_pairs[i].pair_sex == MF )
      out << "MF";
    else if( my_pairs[i].pair_sex == FF )
      out << "FF";
    else
      out << "??";
  }

  out << endl << " Markers:" << endl;

  out << "marker count = " << my_markers.size() << endl;

  for( size_t i = 0; i < my_markers.size(); ++i )
  {
    if( i && i % 5 == 0 )
      out << endl;

    out << " " << i << " " << my_markers[i].name
        << " " << my_markers[i].type
        << "(" << my_markers[i].distance
        << ")";
  }

  out << endl;

  return;
}

void
relative_pairs::dump_sibship_cluster(ostream &out) const
{
  out << "========================" << endl
      << "  Sibship Cluster Dump  " << endl
      << "========================" << endl << endl;

  out << "fsib_pair count = " << my_fsib_pair_count << endl;
  out << "hsib_pair count = " << my_hsib_pair_count << endl;

  out << "sibship_cluster :" << endl;
  for( size_t ci = 0; ci < my_sibship_cluster.size(); ++ci )
  {
    const set<mem_pointer>&               p_set = my_sibship_cluster[ci].parent_set;
    const vector<size_t>&                 fsibs = my_sibship_cluster[ci].fsib_pairs;
    const vector<size_t>&                 hsibs = my_sibship_cluster[ci].hsib_pairs;
    const map< id_pair, vector<size_t> >& f_map = my_sibship_cluster[ci].full_sibship_map;

    out << "Cluster No. " << ci+1
        << ", parents = " << p_set.size()
        << ", sib pair = " << fsibs.size() + hsibs.size()
        << endl;

    set<mem_pointer>::const_iterator pi = p_set.begin();
    for( ; pi != p_set.end(); ++pi )
    {
      out << (*pi)->name() << " ";
    }
    out << endl;

    for( size_t i = 0; i < fsibs.size(); ++i )
    {
      size_t si = fsibs[i];
      out << " "
          << si
          << "(" << my_pairs[si].pair.first->name()
          << "," << my_pairs[si].pair.second->name()
          << ")";
    }
    out << endl;

    for( size_t i = 0; i < hsibs.size(); ++i )
    {
      size_t si = hsibs[i];
      out << " "
          << si
          << "(" << my_pairs[si].pair.first->name()
          << "," << my_pairs[si].pair.second->name()
          << ")";
    }
    out << endl;

    map< id_pair, vector<size_t> >::const_iterator mi = f_map.begin();
    for( ; mi != f_map.end(); ++mi )
    {
      out << "parents: " << mi->first.first->name() << "," << mi->first.second->name() << endl;
      const vector<size_t>& full_sibs = mi->second;

      for( size_t fi = 0; fi < full_sibs.size(); ++fi )
      {
        size_t si = full_sibs[fi];
        out << " " << si
            << "(" << my_pairs[si].pair.first->name()
            << "," << my_pairs[si].pair.second->name()
            << ")";
      }
      out << endl;
    }
  }
}

void
relative_pairs::dump_pairs_by_pedigree(ostream &out) const
{
  return;
}

} // end namespace PALBASE
} // end namespace SAGE
