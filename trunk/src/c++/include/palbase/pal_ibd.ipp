// ===========================================
// Inline Implementation of pair_analysis_ibd
// ===========================================

inline
pair_analysis_ibd::pair_analysis_ibd(relative_pairs &p)
                 : my_pairs(&p)
{}

inline
pair_analysis_ibd::~pair_analysis_ibd() {}

inline void
pair_analysis_ibd::build()
{
  my_pairs->build();
}

inline bool
pair_analysis_ibd::built() const
{
  return my_pairs->built();
}

inline bool
pair_analysis_ibd::has_pedigree()
{
  return true;
}

inline size_t
pair_analysis_ibd::add_marker(const string &name, double dist, gmodel_type mt)
{
  IBD::add_marker(name, dist, mt);

  return my_pairs->add_marker(name, dist, mt);
}

inline size_t
pair_analysis_ibd::pair_count() const
{
  return my_pairs->pair_count();
}

inline const id_pair
pair_analysis_ibd::get_pair(size_t i) const
{
  if( i >= my_pairs->pair_count() )
    return id_pair(NULL,NULL);
  
  return make_pair(my_pairs->rels(i).pair.first, my_pairs->rels(i).pair.second);
}

inline id_pair
pair_analysis_ibd::get_pair(size_t i)
{
  if( i >= my_pairs->pair_count() )
    return id_pair(NULL,NULL);
  
  return make_pair(my_pairs->rels(i).pair.first, my_pairs->rels(i).pair.second);
}

inline const id_pair
pair_analysis_ibd::get_pair(const std::string &ped,
                            const std::string &i1,
                            const std::string &i2,
                            error_t &e) const
{
  FPED::PedigreePointer p = my_pairs->get_fped().pedigree_find(ped);

  if( !p )
  {
    e = bad_pedigree;
    return id_pair(NULL,NULL);
  }
  
  FPED::MemberPointer m1 = p->member_find(i1);
  FPED::MemberPointer m2 = p->member_find(i2);

  if( m1 && m2 )
  {
    e = no_error;
    return id_pair( m1, m2 );
  }
  else if( !m1 && !m2 )
    e = bad_ind_both;
  else if( !m1 )
    e = bad_ind1;
  else if( !m2 )
    e = bad_ind2;
    
  return id_pair(NULL,NULL);
}
        
inline id_pair
pair_analysis_ibd::get_pair(const std::string &ped,
                            const std::string &i1,
                            const std::string &i2,
                            error_t &e)
{
  FPED::PedigreePointer p = my_pairs->get_fped().pedigree_find(ped);

  if( !p )
  {
    e = bad_pedigree;
    return id_pair(NULL,NULL);
  }
  
  FPED::MemberPointer m1 = p->member_find(i1);
  FPED::MemberPointer m2 = p->member_find(i2);

  if( m1 && m2 )
  {
    e = no_error;
    return id_pair( &(*m1), &(*m2) );
  }
  else if( !m1 && !m2 )
    e = bad_ind_both;
  else if(!m1)
    e = bad_ind1;
  else if(!m2)
    e = bad_ind2;
    
  return id_pair(NULL,NULL);
}
    
inline bool
pair_analysis_ibd::get_pair(size_t i, std::string &ped, std::string &i1, std::string &i2) const
{
  if( i >= my_pairs->pair_count() )
  {
    ped = i1 = i2 = "";
    return false;
  }
        
  ped = my_pairs->rels(i).pair.first->pedigree()->name();
  i1  = my_pairs->rels(i).pair.first->name();
  i2  = my_pairs->rels(i).pair.second->name();

  return true;
}

inline bool
pair_analysis_ibd::use_pair(size_t i) const
{
  if( i >= my_pairs->pair_count() )
    return false;
    
  // move this into the pairs
  return true;
}

inline bool
pair_analysis_ibd::use_pair(const mem_pointer i1, const mem_pointer i2) const
{
  if( i1 == NULL || i2 == NULL )
    return false;

  // move this into the pairs
  return true;
}

inline bool
pair_analysis_ibd::valid_pair(size_t i) const
{
  if( i >= my_pairs->pair_count() )
    return false;

  return true;
}

inline bool
pair_analysis_ibd::invalidate_pair(size_t i) const
{
  // provide real functionality in pairs
  return false;
}

inline size_t
pair_analysis_ibd::add_pair(mem_pointer i1, mem_pointer i2, pair_type pt)
{
  if( i1 == NULL || i2 == NULL )
    return (size_t)-1;

  return my_pairs->add_pair(i1, i2, pt);
}
    
inline size_t
pair_analysis_ibd::pair_index(const mem_pointer i1, const mem_pointer i2) const
{
  return my_pairs->find_pair(i1, i2);
}

inline bool
pair_analysis_ibd::set_ibd(size_t i, size_t m, double f0, double f2)
{
  if( i >= pair_count() || m >= marker_count() )
  {
    cout << "Invalid pair_count (" << i << "," << pair_count() << "), "
         <<       "marker_count (" << m << "," << marker_count()
         << ")" << endl;
    return false;
  }
  
  if( !finite(f0) || !finite(f2) || f0 < 0 || f2 < 0 || f0+f2 > 1.01 )
    f0 = f2 = QNAN;
    
  my_pairs->set_f0(i,m,f0);
  my_pairs->set_f2(i,m,f2);

  return true;
}
    
inline bool
pair_analysis_ibd::set_ibd(size_t i, const std::vector<double> &f0,
                                     const std::vector<double> &f2)
{
  if( i >= pair_count() || f0.size() != f2.size()
                        || f0.size() != marker_count() )
    return false;
    
  size_t mc = marker_count();

  for( size_t m = 0; m < mc; ++m )
  {
    double f0m = f0[m];
    double f2m = f2[m];

    if( !finite(f0m) || !finite(f2m) || f0m < 0 || f2m < 0 || f0m+f2m > 1.01 )
      f0m = f2m = QNAN;

    my_pairs->set_f0(i,m,f0m);
    my_pairs->set_f2(i,m,f2m);
  }

  return true;
}

inline bool
pair_analysis_ibd::set_ibd(size_t i, size_t m, double f0, double f1mp, double f2)
{
  if( i >= pair_count() || m >= marker_count() )
  {
    cout << "Invalid pair_count (" << i << "," << pair_count() << "), "
         <<       "marker_count (" << m << "," << marker_count() 
         << ")" << endl;

    return false;
  }
  
  if( !finite(f0) || !finite(f2) || f0 < 0 || f2 < 0 || f0+f2 > 1.01 )
    f0 = f2 = QNAN;
    
  my_pairs->set_f0(i,m,f0);
  my_pairs->set_f1mp(i,m,f1mp);
  my_pairs->set_f2(i,m,f2);

  return true;
}
    
inline bool
pair_analysis_ibd::set_ibd(size_t i, const std::vector<double> &f0,
                                     const std::vector<double> &f1mp,
                                     const std::vector<double> &f2)
{
  if( i >= pair_count() || f0.size()   != f2.size()
                        || f1mp.size() != f2.size()
                        || f0.size()   != marker_count() )
    return false;
    
  size_t mc = marker_count();

  for( size_t m = 0; m < mc; ++m )
  {
    double f0m = f0[m];
    double f1  = f1mp[m];
    double f2m = f2[m];

    if( !finite(f0m) || !finite(f2m) || f0m < 0 || f2m < 0 || f0m+f2m > 1.01 )
      f0m = f2m = QNAN;

    my_pairs->set_f0(i,m,f0m);
    my_pairs->set_f1mp(i,m,f1);
    my_pairs->set_f2(i,m,f2m);
  }

  return true;
}

inline bool
pair_analysis_ibd::get_ibd(size_t i, size_t m, double &f0, double &f2) const
{
  f0 = f2 = QNAN;

  if( i >= pair_count() || m >= marker_count() )
    return false;
  
  f0 = my_pairs->prob_share(i,m,0);
  f2 = my_pairs->prob_share(i,m,2);

  return true;
}
    
inline bool
pair_analysis_ibd::get_ibd(size_t i, std::vector<double> &f0,
                                     std::vector<double> &f2) const
{
  f0.resize(0);
  f0.resize( marker_count(), QNAN );
  f2.resize(0);
  f2.resize( marker_count(), QNAN );

  if( i >= pair_count() )
    return false;

  size_t mc = marker_count();
  
  for( size_t m = 0; m < mc; ++m )
  {
    f0[m] = my_pairs->prob_share(i,m,0);
    f2[m] = my_pairs->prob_share(i,m,2);
  }

  return true;
}

inline bool
pair_analysis_ibd::get_ibd(size_t i, size_t m, double &f0, double &f1mp, double &f2) const
{
  f0 = f1mp = f2 = QNAN;

  if( i >= pair_count() || m >= marker_count() )
    return false;
  
  f0   = my_pairs->prob_share(i,m,0);
  f1mp = my_pairs->prob_share(i,m,3);
  f2   = my_pairs->prob_share(i,m,2);

  return true;
}
    
inline bool
pair_analysis_ibd::get_ibd(size_t i, std::vector<double> &f0,
                                     std::vector<double> &f1mp,
                                     std::vector<double> &f2) const
{
  f0.resize(0);
  f0.resize( marker_count(), QNAN );
  f1mp.resize(0);
  f1mp.resize( marker_count(), QNAN );
  f2.resize(0);
  f2.resize( marker_count(), QNAN );

  if( i >= pair_count() )
    return false;

  size_t mc = marker_count();
  
  for( size_t m = 0; m < mc; ++m )
  {
    f0[m]   = my_pairs->prob_share(i,m,0);
    f1mp[m] = my_pairs->prob_share(i,m,3);
    f2[m]   = my_pairs->prob_share(i,m,2);
  }

  return true;
}

inline const sped_pointer
pair_analysis_ibd::get_subped(size_t i) const
{
  return my_pairs->get_subped(i);
}

inline sped_pointer
pair_analysis_ibd::get_subped(size_t i)
{
  return my_pairs->get_subped(i);
}

inline bool
pair_analysis_ibd::set_ibd_state(size_t m, const a_marker_ibd_state& i_state)
{
  return false;
}

inline bool
pair_analysis_ibd::get_ibd_state(size_t m, a_marker_ibd_state& i_state) const
{
  return false;
}

inline bool
pair_analysis_ibd::set_ibd_state(const sped_pointer sp, const ibd_state_info& i_info)
{
  return my_pairs->set_ibd_state(sp, i_info);
}

inline bool
pair_analysis_ibd::get_ibd_state(const sped_pointer sp, ibd_state_info& i_info) const
{
  return false;
}
