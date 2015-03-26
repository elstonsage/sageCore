// ===========================================
// Inline Implementation of basic_storage_ibd
// ===========================================

inline
basic_storage_ibd::basic_storage_ibd(const meiosis_map& mmap)
                 : my_built(false)
{
  my_pedigree = mmap;
}
inline
basic_storage_ibd::basic_storage_ibd(const meiosis_map&  mmap,
                                     const region_type&  region,
                                     bool                use_intervals)
                 : my_built(false)
{
  my_pedigree = mmap;

  gmodel_type m_type = MLOCUS::AUTOSOMAL;
  if( region.is_x_linked() )
    m_type = MLOCUS::X_LINKED;

  if( use_intervals )
  {
    string m_name = "";
    double m_dist = 0.0;

    for( size_t i = 0; i < region.locus_count() - 1; ++i )
    {
      m_name = region.locus(i).name();
      m_dist = region.locus(i).point_location(0);

      add_marker(m_name, m_dist, m_type);

      for( size_t j = 1; j < region.locus(i).interval_point_count(1); ++j )
      {
        string loc_name("");
        if( region.name().size() > 0 )
          loc_name = region.name() + '_';

	loc_name += doub2str(region.locus(i).point_location(j), 0, 1, ios::showpoint | ios::fixed);

	double loc_dist = region.locus(i).point_location(j);

        add_marker(loc_name, loc_dist, m_type);
      }
    }

    m_name = region.locus(region.locus_count() - 1).name();
    m_dist = region.locus(region.locus_count() - 1).point_location(0);

    add_marker(m_name, m_dist, m_type);
  }
  else
  {
    for( size_t i = 0; i < region.locus_count(); ++i )
    {
      string m_name = region.locus(i).name();
      double m_dist = region.locus(i).point_location(0);

      add_marker(m_name, m_dist, m_type);
    }
  }
}

inline
basic_storage_ibd::~basic_storage_ibd()
{}

inline void
basic_storage_ibd::build()
{
  my_ibd_probs.resize(pair_count());

  for( size_t i = 0; i < pair_count(); ++i )
  {
    my_ibd_probs[i].f0s.resize(marker_count(), QNAN);
    my_ibd_probs[i].f1mp.resize(marker_count(), QNAN);
    my_ibd_probs[i].f2s.resize(marker_count(), QNAN);
  }

  my_ibd_states.resize(marker_count());

  my_built = true;
}

inline bool
basic_storage_ibd::built() const
{
  return my_built;
}

inline bool
basic_storage_ibd::has_pedigree()
{
  return my_pedigree.built();
}

inline size_t
basic_storage_ibd::pair_count() const
{
  return my_pairs.size();
}

inline const id_pair
basic_storage_ibd::get_pair(size_t i) const
{
  return my_pairs[i].pair;
}

inline id_pair
basic_storage_ibd::get_pair(size_t i)
{
  return my_pairs[i].pair;
}

inline const id_pair 
basic_storage_ibd::get_pair(const std::string &ped,
                            const std::string &i1,
                            const std::string &i2,
                            error_t &e) const
{
  if(!test_pedigree(ped))
  {
    e = bad_pedigree;
    
    return std::make_pair((mem_pointer) NULL, (mem_pointer) NULL);
  }

  mem_pointer id1 = NULL, id2 = NULL;
  
  for(size_t i = 0; i < my_pedigree.get_subpedigree()->member_count(); ++i)
  {
    if(my_pedigree.member(i)->name() == i1) id1 = my_pedigree.member(i);
    if(my_pedigree.member(i)->name() == i2) id2 = my_pedigree.member(i);
  }
  
  if(!id1)
  {
    if(!id2)
    {
      e = bad_ind_both;

      return std::make_pair((mem_pointer) NULL, (mem_pointer) NULL);
    }
    else
    {
      e = bad_ind1;

      return std::make_pair((mem_pointer) NULL, id2);
    }
  }
  else
  {
    if(!id2)
    {
      e = bad_ind2;

      return std::make_pair(id1, (mem_pointer) NULL);
    }
    else
    {
      e = no_error;

      return std::make_pair(id1, id2);
    }
  }  
}

inline id_pair 
basic_storage_ibd::get_pair(const std::string &ped,
                            const std::string &i1,
                            const std::string &i2,
                            error_t &e)
{
  if(!my_pedigree.built() || ped != my_pedigree.get_pedigree()->name())
  {
    e = bad_pedigree;
    
    return std::make_pair((mem_pointer) NULL, (mem_pointer) NULL);
  }

  mem_pointer id1 = NULL, id2 = NULL;
  
  for(size_t i = 0; i < my_pedigree.get_subpedigree()->member_count(); ++i)
  {
    if(my_pedigree.member(i)->name() == i1) id1 = my_pedigree.member(i);
    if(my_pedigree.member(i)->name() == i2) id2 = my_pedigree.member(i);
  }
  
  if(!id1)
  {
    if(!id2)
    {
      e = bad_ind_both;

      return std::make_pair((mem_pointer) NULL, (mem_pointer) NULL);
    }
    else
    {
      e = bad_ind1;

      return std::make_pair((mem_pointer) NULL, id2);
    }
  }
  else
  {
    if(!id2)
    {
      e = bad_ind2;

      return std::make_pair(id1, (mem_pointer) NULL);
    }
    else
    {
      e = no_error;

      return std::make_pair(id1, id2);
    }
  }  
}

inline bool
basic_storage_ibd::get_pair(size_t i, std::string &ped, std::string &i1, std::string &i2) const
{
  if(!built() || !test_pedigree() || i > my_pairs.size()) return false;
  
  ped = my_pedigree.get_pedigree()->name();
  i1 = my_pairs[i].pair.first->name();
  i2 = my_pairs[i].pair.second->name();
  
  return true;
}

inline bool
basic_storage_ibd::use_pair(size_t i) const
{
  return i < pair_count();
}

inline bool
basic_storage_ibd::use_pair(const mem_pointer, const mem_pointer) const
{
  return true;
}

inline bool
basic_storage_ibd::valid_pair(size_t i) const
{
  return (built() && i < pair_count());
}

inline bool
basic_storage_ibd::invalidate_pair(size_t i) const
{
  return !valid_pair(i);
}

inline size_t
basic_storage_ibd::add_pair(mem_pointer i1, mem_pointer i2, pair_type pt)
{
  if(!test_pedigree() || i1->subindex() > my_pedigree.get_subpedigree()->member_count()
                      || i2->subindex() > my_pedigree.get_subpedigree()->member_count())
    return (size_t) -1;

  my_pairs.push_back(ibd_pair_info(i1, i2, pt));

  return pair_count() - 1;
}

inline size_t
basic_storage_ibd::pair_index(const mem_pointer i1, const mem_pointer i2) const
{
  if(!built()) return (size_t) -1;
  
  // Pairs are sorted, but we'll not worry about that just now.
  for(size_t i = 0; i < pair_count(); ++i)
    if((my_pairs[i].pair.first == i1 &&
        my_pairs[i].pair.second == i2) ||
       (my_pairs[i].pair.first == i2 &&
        my_pairs[i].pair.second == i1)) return i;

  return (size_t) -1;
}

inline bool
basic_storage_ibd::set_ibd(size_t i, size_t m, double f0, double f2)
{
  if(i > pair_count() || m > marker_count()) return false;

  if(!finite(f0) || !finite(f2) 
                 || f0 < -20 * EPS
                 || f2 < -20 * EPS
                 || f0 + f2 > 1.0 + 20 * EPS)
  {
    my_ibd_probs[i].f0s[m] = QNAN;
    my_ibd_probs[i].f2s[m] = QNAN;
  }
  else
  {
    if(f0 > 1.0) f0 = 1.0;
    if(f2 > 1.0) f2 = 1.0;
    if(f0 < 0.0) f0 = 0.0;
    if(f2 < 0.0) f2 = 0.0;

    my_ibd_probs[i].f0s[m] = f0;
    my_ibd_probs[i].f2s[m] = f2;
  }
  
  return true;
}

inline bool
basic_storage_ibd::set_ibd(size_t i, const std::vector<double> &f0,
                                     const std::vector<double> &f2)
{
  if(i > pair_count() || marker_count() != f0.size()
                      || marker_count() != f2.size()) return false;

  for(size_t m = 0; m < marker_count(); ++m)
    set_ibd(i, m, f0[m], f2[m]);

  return true;
}

inline bool
basic_storage_ibd::set_ibd(size_t i, size_t m, double f0, double f1mp, double f2)
{
  if(i > pair_count() || m > marker_count()) return false;

  if(!finite(f0) || !finite(f2) 
                 || f0 < -20 * EPS
                 || f2 < -20 * EPS
                 || f0 + f2 > 1.0 + 20 * EPS)
  {
    my_ibd_probs[i].f0s[m]  = QNAN;
    my_ibd_probs[i].f1mp[m] = QNAN;
    my_ibd_probs[i].f2s[m]  = QNAN;
  }
  else
  {
    if(f0 > 1.0) f0 = 1.0;
    if(f2 > 1.0) f2 = 1.0;
    if(f0 < 0.0) f0 = 0.0;
    if(f2 < 0.0) f2 = 0.0;

    my_ibd_probs[i].f0s[m]  = f0;
    my_ibd_probs[i].f1mp[m] = f1mp;
    my_ibd_probs[i].f2s[m]  = f2;
  }
  
  return true;
}

inline bool
basic_storage_ibd::set_ibd(size_t i, const std::vector<double> &f0,
                                     const std::vector<double> &f1mp,
                                     const std::vector<double> &f2)
{
  if(i > pair_count() || marker_count() != f0.size()
                      || marker_count() != f1mp.size()
                      || marker_count() != f2.size()) return false;

  for(size_t m = 0; m < marker_count(); ++m)
    set_ibd(i, m, f0[m], f1mp[m], f2[m]);

  return true;
}

inline bool
basic_storage_ibd::get_ibd(size_t i, size_t m, double &f0, double &f2) const
{
  if(i > pair_count() || m > marker_count()) return false;

  f0 = my_ibd_probs[i].f0s[m];
  f2 = my_ibd_probs[i].f2s[m];

  return true;
}

inline bool
basic_storage_ibd::get_ibd(size_t i, std::vector<double> &f0,
                                     std::vector<double> &f2) const
{
  if(i > pair_count()) return false;

  f0 = my_ibd_probs[i].f0s;
  f2 = my_ibd_probs[i].f2s;

  return true;
}

inline bool
basic_storage_ibd::get_ibd(size_t i, size_t m, double &f0, double &f1mp, double &f2) const
{
  if(i > pair_count() || m > marker_count()) return false;

  f0   = my_ibd_probs[i].f0s[m];
  f1mp = my_ibd_probs[i].f1mp[m];
  f2   = my_ibd_probs[i].f2s[m];

  return true;
}

inline bool
basic_storage_ibd::get_ibd(size_t i, std::vector<double> &f0,
                                     std::vector<double> &f1mp,
                                     std::vector<double> &f2) const
{
  if(i > pair_count()) return false;

  f0   = my_ibd_probs[i].f0s;
  f1mp = my_ibd_probs[i].f1mp;
  f2   = my_ibd_probs[i].f2s;

  return true;
}

inline const sped_pointer
basic_storage_ibd::get_subped(size_t i) const
{
  return my_pairs[i].pair.first->subpedigree();
}

inline sped_pointer
basic_storage_ibd::get_subped(size_t i)
{
  return my_pairs[i].pair.first->subpedigree();
}

inline bool
basic_storage_ibd::set_ibd_state(size_t m, const a_marker_ibd_state& i_state)
{
  if( m > marker_count() )
    return false;

  my_ibd_states[m] = i_state;
  
  return true;
}

inline bool
basic_storage_ibd::get_ibd_state(size_t m, a_marker_ibd_state& i_state) const
{
  if( m > marker_count() )
    return false;

  i_state = my_ibd_states[m];
  
  return true;
}

inline bool
basic_storage_ibd::set_ibd_state(const sped_pointer sp, const ibd_state_info& i_state)
{
  my_ibd_states = i_state.ibd_states;

  return true;
}

inline bool
basic_storage_ibd::get_ibd_state(const sped_pointer sp, ibd_state_info& i_state) const
{
  i_state.ibd_states = my_ibd_states;

  return true;
}

inline bool
basic_storage_ibd::test_pedigree(const std::string& s) const
{
  if(!my_pedigree.built() || !my_pedigree.get_subpedigree()
                          || (s.size() && my_pedigree.get_pedigree()->name() != s))
    return false;

  return true;
}
