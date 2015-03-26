// ===========================================
// Inline Implementation of sim_storage_ibd
// ===========================================

inline
sim_storage_ibd::sim_storage_ibd(const meiosis_map&         ped,
                                 const mcmc_parameters&     par,
                                 vector<sim_relative_pair>& pairs)
               : my_pedigree(ped), my_params(par),
                 my_pair_ibds(pairs), my_built(false)
{}

inline
sim_storage_ibd::~sim_storage_ibd()
{}

inline void
sim_storage_ibd::build()
{
  my_built = true;
}

inline bool
sim_storage_ibd::built() const
{
  return my_built;
}

inline bool
sim_storage_ibd::has_pedigree()
{
  return my_pedigree.built();
}

inline size_t
sim_storage_ibd::pair_count() const
{
  return my_pair_ibds.size();
}

inline const id_pair
sim_storage_ibd::get_pair(size_t i) const
{
  return make_pair(my_pair_ibds[i].get_first_ind(),
                   my_pair_ibds[i].get_second_ind());
}

inline id_pair
sim_storage_ibd::get_pair(size_t i)
{
  return make_pair(my_pair_ibds[i].get_first_ind(),
                   my_pair_ibds[i].get_second_ind());
}

inline const id_pair 
sim_storage_ibd::get_pair(const std::string &ped, const std::string &i1,
                          const std::string &i2, error_t &e) const
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
sim_storage_ibd::get_pair(const std::string &ped, const std::string &i1,
                          const std::string &i2, error_t &e)
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
sim_storage_ibd::get_pair(size_t i, std::string &ped, std::string &i1, std::string &i2) const
{
//cout << "get_pair(i = " << i << ")" << endl;
 
  if(!built() || !test_pedigree() || i > my_pair_ibds.size()) return false;

//cout << my_pair_ibds[i].get_first_ind() << " "
//     << my_pair_ibds[i].get_second_ind() << endl;

  ped = my_pair_ibds[i].get_first_ind()->pedigree()->name();

//cout << ped << endl;

  i1  = my_pair_ibds[i].get_first_ind()->name();

//cout << i1 << endl;

  i2  = my_pair_ibds[i].get_second_ind()->name();

//cout << i2 << endl;
  
  return true;
}

inline bool
sim_storage_ibd::use_pair(size_t i) const
{
  return i < pair_count();
}

inline bool
sim_storage_ibd::use_pair(const mem_pointer, const mem_pointer) const
{
  return true;
}

inline bool
sim_storage_ibd::valid_pair(size_t i) const
{
  return (built() && i < pair_count());
}

inline bool
sim_storage_ibd::invalidate_pair(size_t i) const
{
  return !valid_pair(i);
}

inline size_t
sim_storage_ibd::add_pair(mem_pointer i1, mem_pointer i2, pair_type pt)
{
//cout << "sim_storage_ibd::add_pair(" << i1->name()
//     << ", " << i2->name() << ", " << p_type << ")" << endl;

  if(!test_pedigree() || i1->subindex() > my_pedigree.get_subpedigree()->member_count()
                      || i2->subindex() > my_pedigree.get_subpedigree()->member_count())
    return (size_t) -1;

  sim_relative_pair a_pair(marker_count(), i1, i2, pt);

  my_pair_ibds.push_back(a_pair);

//cout << "added, " << my_pair_ibds[pair_count() - 1].get_first_ind()->name()
//     << my_pair_ibds[pair_count() - 1].get_second_ind()->name()
//     << ", my_pair_ibds.size() = " << my_pair_ibds.size() << endl;

  return pair_count() - 1;
}

inline size_t
sim_storage_ibd::pair_index(const mem_pointer i1, const mem_pointer i2) const
{
  if(!built()) return (size_t) -1;
  
  // Pairs are sorted, but we'll not worry about that just now.
  for(size_t i = 0; i < pair_count(); ++i)
    if((my_pair_ibds[i].get_first_ind() == i1 &&
        my_pair_ibds[i].get_second_ind() == i2) ||
       (my_pair_ibds[i].get_first_ind() == i2 &&
        my_pair_ibds[i].get_second_ind() == i1)) return i;

  return (size_t) -1;
}

inline bool
sim_storage_ibd::set_ibd(size_t i, size_t m, double f0, double f2)
{
  if(i > pair_count() || m > marker_count()) return false;

  size_t steps = 0;

  if(!finite(f0) || !finite(f2) 
                 || f0 < -20 * EPS
                 || f2 < -20 * EPS
                 || f0 + f2 > 1.0 + 20 * EPS)
  {
    my_pair_ibds[i].set_values(m, QNAN, QNAN, steps);
  }
  else
  {
    if(f0 > 1.0) f0 = 1.0;
    if(f2 > 1.0) f2 = 1.0;
    if(f0 < 0.0) f0 = 0.0;
    if(f2 < 0.0) f2 = 0.0;

    my_pair_ibds[i].set_values(m, f0, f2, steps);
  }
  
  return true;
}

inline bool
sim_storage_ibd::set_ibd(size_t i, const std::vector<double> &f0, const std::vector<double> &f2)
{
  if(i > pair_count() || marker_count() != f0.size()
                      || marker_count() != f2.size()) return false;

  for(size_t m = 0; m < marker_count(); ++m)
    set_ibd(i, m, f0[m], f2[m]);

  return true;
}

inline bool
sim_storage_ibd::set_ibd(size_t i, size_t m, double f0, double f1mp, double f2)
{
  if(i > pair_count() || m > marker_count()) return false;

  size_t steps = 0;

  if(!finite(f0) || !finite(f2) 
                 || f0 < -20 * EPS
                 || f2 < -20 * EPS
                 || f0 + f2 > 1.0 + 20 * EPS)
  {
    my_pair_ibds[i].set_values(m, QNAN, QNAN, QNAN, steps);
  }
  else
  {
    if(f0 > 1.0) f0 = 1.0;
    if(f2 > 1.0) f2 = 1.0;
    if(f0 < 0.0) f0 = 0.0;
    if(f2 < 0.0) f2 = 0.0;

    my_pair_ibds[i].set_values(m, f0, f2, f1mp, steps);
  }
  
  return true;
}

inline bool
sim_storage_ibd::set_ibd(size_t i, const std::vector<double> &f0,
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
sim_storage_ibd::get_ibd(size_t i, size_t m, double &f0, double &f2) const
{
  if(i > pair_count() || m > marker_count()) return false;

  f0 = my_pair_ibds[i].get_value(m, 0);
  f2 = my_pair_ibds[i].get_value(m, 2);

  return true;
}

inline bool
sim_storage_ibd::get_ibd(size_t i, std::vector<double> &f0,
                                   std::vector<double> &f2) const
{
  if(i > pair_count() || marker_count() != f0.size()
                      || marker_count() != f2.size()) return false;

  for(size_t m = 0; m < marker_count(); ++m)
  {
    f0[m] = my_pair_ibds[i].get_value(m, 0);
    f2[m] = my_pair_ibds[i].get_value(m, 2);
  }

  return true;
}

inline bool
sim_storage_ibd::get_ibd(size_t i, size_t m, double &f0, double &f1mp, double &f2) const
{
  if(i > pair_count() || m > marker_count()) return false;

  f0   = my_pair_ibds[i].get_value(m, 0);
  f1mp = my_pair_ibds[i].get_value(m, 1);
  f2   = my_pair_ibds[i].get_value(m, 2);

  return true;
}

inline bool
sim_storage_ibd::get_ibd(size_t i, std::vector<double> &f0,
                                   std::vector<double> &f1mp,
                                   std::vector<double> &f2) const
{
  if(i > pair_count() || marker_count() != f0.size()
                      || marker_count() != f1mp.size()
                      || marker_count() != f2.size()) return false;

  for(size_t m = 0; m < marker_count(); ++m)
  {
    f0[m]   = my_pair_ibds[i].get_value(m, 0);
    f1mp[m] = my_pair_ibds[i].get_value(m, 1);
    f2[m]   = my_pair_ibds[i].get_value(m, 2);
  }

  return true;
}

inline const sped_pointer
sim_storage_ibd::get_subped(size_t i) const
{
  return my_pair_ibds[i].get_first_ind()->subpedigree();
}

inline sped_pointer
sim_storage_ibd::get_subped(size_t i)
{
  return my_pair_ibds[i].get_first_ind()->subpedigree();
}

inline bool
sim_storage_ibd::set_ibd_state(size_t m, const a_marker_ibd_state& i_state)
{
  return false;
}

inline bool
sim_storage_ibd::get_ibd_state(size_t m, a_marker_ibd_state& i_state) const
{
  return false;
}

inline bool
sim_storage_ibd::set_ibd_state(const sped_pointer sp, const ibd_state_info& i_info)
{
  return false;
}

inline bool
sim_storage_ibd::get_ibd_state(const sped_pointer sp, ibd_state_info& i_info) const
{
  return false;
}

inline bool
sim_storage_ibd::test_pedigree(const std::string& s) const
{
  if(!my_pedigree.built() || !my_pedigree.get_subpedigree()
                          || (s.size() && my_pedigree.get_pedigree()->name() != s))
    return false;

  return true;
}

