//////////////////////////////////////////////////////////////////////////
//                 Implementation of rel_pair (Inline)                  //
//////////////////////////////////////////////////////////////////////////

inline void
rel_pair::set_pair(size_t n)
{
  pair_num = n;
}

inline bool
rel_pair::is_fsib_pair() const
{
  if(!data)
    return false;

  return data->rels(pair_num).type == pair_generator::SIBSIB;
}

inline bool
rel_pair::is_hsib_pair() const
{
  if(!data)
    return false;

  return data->rels(pair_num).type == pair_generator::HALFSIB;
}

inline bool
rel_pair::is_mm_pair() const
{
  if(!data)
    return false;

  return data->rels(pair_num).pair_sex == MM;
}

inline bool
rel_pair::is_mf_pair() const
{
  if(!data)
    return false;

  return data->rels(pair_num).pair_sex == MF;
}

inline bool
rel_pair::is_ff_pair() const
{
  if(!data)
    return false;

  return data->rels(pair_num).pair_sex == FF;
}

inline size_t
rel_pair::pair_number() const
{
  if(!data)
    return (size_t)-1;

  return pair_num;
}

inline size_t
rel_pair::pedigree_number() const
{
  if(!data)
    return (size_t)-1;

  return data->pedigree_number(pair_num);
}

inline const rel_pair_data&
rel_pair::rels() const
{
  return data->rels(pair_num);
}

inline FPED::PedigreeConstPointer
rel_pair::pedigree() const
{
  if(!data)
    return NULL;

  return rels().pair.first->pedigree();
}

inline relative_pairs*
rel_pair::pair_data() const
{
  return data;
}

inline size_t
rel_pair::marker_count() const
{
  if(!data)
    return 0;

  return data->marker_count();
}

inline size_t
rel_pair::marker_find(const string &name) const
{
  if(!data)
    return (size_t)-1;

  return data->marker_find(name);
}

inline string
rel_pair::marker_name(size_t i) const
{
  if(!data)
    return string();

  return data->marker_name(i);
}

inline double
rel_pair::avg_share(size_t m) const
{
  if(!data)
    return std::numeric_limits<double>::quiet_NaN();

  return data->avg_share(pair_num,m);
}

inline double
rel_pair::prob_share(size_t m, size_t n) const
{
  if(!data)
    return std::numeric_limits<double>::quiet_NaN();

  return data->prob_share(pair_num,m,n);
}

inline double
rel_pair::weighted_share(size_t m, double w0,double w1,double w2) const
{
  if(!data)
    return std::numeric_limits<double>::quiet_NaN();

  return data->weighted_share(pair_num, m, w0, w1, w2);
}

inline double
rel_pair::prior_avg_share() const
{
  if(!data)
    return std::numeric_limits<double>::quiet_NaN();

  return data->prior_avg_share(pair_num);
}

inline double
rel_pair::prior_prob_share(size_t n) const
{
  if(!data)
    return std::numeric_limits<double>::quiet_NaN();

  return data->prior_prob_share(pair_num,n);
}

inline double
rel_pair::prior_weighted_share(double w0,double w1,double w2) const
{
  if(!data)
    return std::numeric_limits<double>::quiet_NaN();

  return data->prior_weighted_share(pair_num, w0, w1, w2);
}

inline bool
rel_pair::operator==(const rel_pair &rhs) const
{
  if(!data || !rhs.data)
    return false;

  return (data == rhs.data && pair_num == rhs.pair_num);
}

inline bool
rel_pair_less::operator()(rel_pair p1, rel_pair p2) const
{
  if( p1.rels().pair.first->pedigree()->name() == p2.rels().pair.first->pedigree()->name() )
    if( p1.rels().pair.first->name() == p2.rels().pair.first->name() )
      return p1.rels().pair.second->name() < p2.rels().pair.second->name();
    else
      return p1.rels().pair.first->name() < p2.rels().pair.first->name();
  else
    return p1.rels().pair.second->pedigree()->name() < p2.rels().pair.second->pedigree()->name();
}
