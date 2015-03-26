//////////////////////////////////////////////////////////////////////////
//                 Implementation of relative_pairs (Inline)            //
//////////////////////////////////////////////////////////////////////////

inline bool
relative_pairs::valid()  const
{
  return my_valid;
}

inline void
relative_pairs::invalidate()
{
  my_valid = false;
}

inline void
relative_pairs::validate()
{
  my_valid = true;
}

inline void
relative_pairs::set_fped(FPED::Multipedigree& fp)
{
  my_fped = &fp;
}

inline const FPED::Multipedigree&
relative_pairs::get_fped() const
{
  return *my_fped;
}

inline FPED::Multipedigree&
relative_pairs::get_fped()
{
  return *my_fped;
}

inline const FPED::FilteredMultipedigreeInfo&
relative_pairs::fped_info() const
{
  return get_fped().info();
}

inline FPED::FilteredMultipedigreeInfo&
relative_pairs::fped_info()
{
  return get_fped().info();
}

inline const FPED::FilteredPedigreeInfo&
relative_pairs::ped_info(size_t p) const
{
  return get_fped().pedigree_index(p).info();
}

inline FPED::FilteredPedigreeInfo&
relative_pairs::ped_info(size_t p)
{
  return get_fped().pedigree_index(p).info();
}

inline size_t
relative_pairs::pedigree_count()  const
{
  return get_fped().pedigree_count();
}

inline string
relative_pairs::pedigree_name(size_t p)  const
{
  return get_fped().pedigree_index(p).name();
}

inline
size_t
relative_pairs::pedigree_number(size_t n) const
{
  return my_pairs[n].pair.first->pedigree()->index();
}

inline size_t
relative_pairs::member_count(size_t p)  const
{
  return ped_info(p).member_count();
}

inline string
relative_pairs::member_name(size_t p, size_t i)  const
{
  return get_fped().pedigree_index(p).member_index(i).name();
}

inline size_t
relative_pairs::marker_count()  const
{
  return my_markers.size();
}

inline size_t
relative_pairs::marker_find(const std::string &name)  const
{
  for( size_t i = 0; i < marker_count(); ++i )
    if( marker_name(i) == name )
      return i;

  return my_markers.size();
}

inline string
relative_pairs::marker_name(size_t m)  const
{
  return my_markers[m].name;
}

inline double
relative_pairs::marker_distance(size_t m)  const
{
  return my_markers[m].distance;
}

inline gmodel_type
relative_pairs::marker_genotype_model(size_t m)  const
{
  return my_markers[m].type;
}

inline bool
relative_pairs::valid_distance_exist() const
{
  return my_valid_distance_exist;
}

inline size_t
relative_pairs::trait_count()  const
{
  return fped_info().trait_count();
}

inline string
relative_pairs::trait_name(size_t t)  const
{
  return fped_info().trait_info(t).name();
}

inline size_t
relative_pairs::trait_find(const std::string &name)  const
{
  return fped_info().trait_find(name);
}

inline size_t
relative_pairs::pair_count() const
{
  return my_pairs.size();
}

inline size_t
relative_pairs::fsib_pair_count() const
{
  return my_fsib_pair_count;
}

inline size_t
relative_pairs::hsib_pair_count() const
{
  return my_hsib_pair_count;
}

inline size_t
relative_pairs::mm_pair_count() const
{
  return my_mm_pair_count;
}

inline size_t
relative_pairs::mf_pair_count() const
{
  return my_mf_pair_count;
}

inline size_t
relative_pairs::ff_pair_count() const
{
  return my_ff_pair_count;
}

inline double
relative_pairs::trait(size_t p, size_t i, size_t t)  const
{
  return ped_info(p).trait(i,t);
}

inline double
relative_pairs::trait_missing_code(size_t t)  const
{
  return fped_info().trait_info(t).numeric_missing_code();
}

inline bool
relative_pairs::trait_missing(size_t p, size_t i, size_t t)  const
{
  return ped_info(p).trait_missing(i,t);
}

inline double
relative_pairs::trait(const MPED::member_base& m, size_t t)  const
{
  return trait( m.pedigree()->index(), m.index(), t);
}

inline bool
relative_pairs::trait_missing(const MPED::member_base& m, size_t t)  const
{
  return trait_missing( m.pedigree()->index(), m.index(), t);
}

inline const rel_pair_data&
relative_pairs::rels(size_t i) const
{
  return my_pairs[i];
}

inline void
relative_pairs::set_f0(size_t i, size_t m, double f0)
{
  my_ibd_probs[i].f0s[m] = f0;
}

inline void
relative_pairs::set_f2(size_t i, size_t m, double f2)
{
  my_ibd_probs[i].f2s[m] = f2;
}

inline void
relative_pairs::set_f1mp(size_t i, size_t m, double f1)
{
  my_ibd_probs[i].f1mp[m] = f1;
}

inline bool
relative_pairs::built() const
{
  return valid();
}

inline bool
relative_pairs::is_x_linked(size_t m) const
{
  if( m >= my_markers.size() )
    return false;

  return my_markers[m].type == MLOCUS::X_LINKED;
}

inline bool
relative_pairs::x_linked_marker_exist() const
{
  for( size_t i = 0; i < my_markers.size(); ++i )
    if( my_markers[i].type == MLOCUS::X_LINKED )
      return true;

  return false;
}

inline bool
relative_pairs::autosomal_marker_exist() const
{
  for( size_t i = 0; i < my_markers.size(); ++i )
    if( my_markers[i].type == MLOCUS::AUTOSOMAL )
      return true;

  return false;
}

inline double
relative_pairs::prob_share(size_t i, size_t m, size_t n) const
{
  double f0   = my_ibd_probs[i].f0s[m];
  double f2   = my_ibd_probs[i].f2s[m];
  double f1mp = my_ibd_probs[i].f1mp[m];

  if( finite(f0) != finite(f2) )
    cout << "problem with " << rels(i).pair.first->pedigree()->name() << ": "
                            << rels(i).pair.first->name()  << ","
                            << rels(i).pair.second->name() << " at marker "
                            << marker_name(m) << endl;

  if( !finite (f0) || !finite(f2) )
    return QNAN;

  double f1 = 1.0-f0-f2;

  switch( n )
  {
    case 0: return f0;
    case 1: return f1;
    case 2: return f2;
    case 3: return f1mp;
  }

  return QNAN;
}

inline double
relative_pairs::avg_share(size_t i, size_t m)  const
{
  double f1 = prob_share(i,m,1);
  double f2 = prob_share(i,m,2);

  if( !finite(f1) || !finite(f2) )
    return QNAN;

  return f2 + 0.5*f1;
}

inline double
relative_pairs::weighted_share(size_t i, size_t m,
                              double w0, double w1, double w2) const
{
  double f0 = prob_share(i,m,0);
  double f2 = prob_share(i,m,2);

  if( !finite(f0) || !finite(f2) )
    return QNAN;

  double f1 = 1.0 - f0 - f2;

  return w0*f0 + w1*f1 + w2*f2;
}

inline double
relative_pairs::prior_prob_share(size_t p, size_t i1, size_t i2, size_t n) const
{
  if( p >= my_prior_ibd.size() || n > 2 )
    return QNAN;

  return my_prior_ibd[p].f(i1,i2,n);
}

inline double
relative_pairs::prior_prob_share(size_t i, size_t n) const
{
  if( i >= pair_count() || n > 2 )
    return QNAN;

  size_t p  = my_pairs[i].pair.first->pedigree()->index();
  size_t i1 = my_pairs[i].pair.first->index();
  size_t i2 = my_pairs[i].pair.second->index();

  return my_prior_ibd[p].f(i1,i2,n);
}

inline double
relative_pairs::prior_avg_share(size_t p, size_t i1, size_t i2)  const
{
  double f1 = prior_prob_share(p, i1, i2, 1);
  double f2 = prior_prob_share(p, i1, i2, 2);

  if( !finite(f1) || !finite(f2) )
    return QNAN;

  return f2 + 0.5*f1;
}

inline double
relative_pairs::prior_avg_share(size_t i)  const
{
  double f1 = prior_prob_share(i,1);
  double f2 = prior_prob_share(i,2);

  if( !finite(f1) || !finite(f2) )
    return QNAN;

  return f2 + 0.5*f1;
}

inline double
relative_pairs::prior_weighted_share(size_t i, double w0, double w1, double w2) const
{
  double f0 = prior_prob_share(i,0);
  double f2 = prior_prob_share(i,2);

  if( !finite(f0) || !finite(f2) )
    return QNAN;

  double f1 = 1.0 - f0 - f2;

  return w0*f0 + w1*f1 + w2*f2;
}

inline sibship_cluster_iterator
relative_pairs::sibship_cluster_begin()
{
  return my_sibship_cluster.begin();
}

inline sibship_cluster_iterator
relative_pairs::sibship_cluster_end()
{
  return my_sibship_cluster.end();
}

inline sibship_cluster_const_iterator
relative_pairs::sibship_cluster_begin() const
{
  return my_sibship_cluster.begin();
}

inline sibship_cluster_const_iterator
relative_pairs::sibship_cluster_end() const
{
  return my_sibship_cluster.end();
}

//
//------------------------------------------------------------------
//

inline bool
relative_pairs::valid_pair_info() const
{
  return my_valid_pair_info;
}

inline void
relative_pairs::invalidate_pair_info()
{
  my_valid_pair_info = false;
  resize_pair_covariates(0);
  resize_pair_weights(0);
  my_pair_covariate_info.resize(0);
  my_pair_weight_info.resize(0);
}

inline void
relative_pairs::validate_pair_info()
{
  my_valid_pair_info = true;
}

inline size_t
relative_pairs::pair_covariate_count() const
{
  return my_pair_covariate_count;
}

inline size_t
relative_pairs::pair_weight_count() const
{
  return my_pair_weight_count;
}

inline size_t
relative_pairs::pair_covariate_info_count() const
{
  return my_pair_covariate_info.size();
}

inline size_t
relative_pairs::pair_weight_info_count() const
{
  return my_pair_weight_info.size();
}

inline double
relative_pairs::get_pair_covariate(size_t pn, size_t c) const
{
  return my_pair_covariates[c][pn];
}

inline double
relative_pairs::get_pair_weight(size_t pn, size_t w) const
{
  return my_pair_weights[w][pn];
}

inline
pair_pheno_info&
relative_pairs::pair_covariate_info(size_t c)
{
  return my_pair_covariate_info[c];
}

inline
const pair_pheno_info&
relative_pairs::pair_covariate_info(size_t c) const
{
  return my_pair_covariate_info[c];
}

inline
pair_pheno_info&
relative_pairs::pair_weight_info(size_t w)
{
  return my_pair_weight_info[w];
}

inline
const pair_pheno_info&
relative_pairs::pair_weight_info(size_t w) const
{
  return my_pair_weight_info[w];
}
    
inline string
relative_pairs::get_pair_covariate_name(size_t c) const
{
  return pair_covariate_info(c).get_name();
}

inline string
relative_pairs::get_pair_weight_name(size_t w) const
{
  return pair_weight_info(w).get_name();
}

inline size_t
relative_pairs::pair_covariate_find(const string& name) const
{
  string cn = toUpper(name);
  for( size_t c = 0; c < pair_covariate_info_count(); ++c )
    if( toUpper(pair_covariate_info(c).get_name()) == cn )
      return c;
      
  return (size_t)-1;
}
    
inline size_t
relative_pairs::pair_weight_find(const string& name) const
{
  string wn = toUpper(name);
  for( size_t w = 0; w < pair_weight_info_count(); ++w )
    if( toUpper(pair_weight_info(w).get_name()) == wn )
      return w;
      
  return (size_t)-1;
}

inline size_t
relative_pairs::add_pair_covariate(const string& name, pair_pheno_info::info_use usage,
                                   double value)
{
  size_t c = pair_covariate_find(name);
  if( c < pair_covariate_info_count() )
      return c;

  my_pair_covariate_info.push_back( pair_pheno_info(name, usage, value) );
  return pair_covariate_info_count() - 1;
}

inline size_t
relative_pairs::add_pair_weight(const string& name, pair_pheno_info::info_use usage,
                                double value)
{
  size_t c = pair_weight_find(name);
  if( c < pair_weight_info_count() )
      return c;

  my_pair_weight_info.push_back( pair_pheno_info(name, usage, value) );
  return pair_weight_info_count() - 1;
}

inline void
relative_pairs::resize_pair_covariates(size_t c)
{
  my_pair_covariate_count = c;
  my_pair_covariates.resize(c);
  
  double qnan = numeric_limits<double>::quiet_NaN();
  for( size_t i = 0; i < c; ++i )
    my_pair_covariates[i].resize(pair_count(), qnan);

  if( c )
    validate_pair_info();
}

inline void
relative_pairs::resize_pair_weights(size_t w)
{
  my_pair_weight_count = w;
  my_pair_weights.resize(w);
  
  double qnan = numeric_limits<double>::quiet_NaN();
  for( size_t i = 0; i < w; ++i )
    my_pair_weights[i].resize(pair_count(), qnan);

  if( w )
    validate_pair_info();
}

inline const sped_pointer
relative_pairs::get_subped(size_t i) const
{
  return my_pairs[i].pair.first->subpedigree();
}

inline sped_pointer
relative_pairs::get_subped(size_t i)
{
  return my_pairs[i].pair.first->subpedigree();
}

inline bool
relative_pairs::set_ibd_state(const sped_pointer p, const ibd_state_info& i_info)
{
  my_ibd_states[p] = i_info;

  return true;
}

inline bool
relative_pairs::get_ibd_state(const sped_pointer p, ibd_state_info& i_info) const
{
  if( my_ibd_states.find(p) != my_ibd_states.end() )
  {
    i_info = my_ibd_states.find(p)->second;

    return true;
  }

  return false;
}

