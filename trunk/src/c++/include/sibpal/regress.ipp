////////////////////////////////////////////////////////////////////////////
//             Implementation of Trait Regression (Inline)                //
////////////////////////////////////////////////////////////////////////////

//
//------------------------------------------------------------------------
//

inline string
RegressionVariant::name() const
{
  return my_name;
}

inline bool
RegressionVariant::apply_weighting() const
{
  return true;
}

inline const map<size_t, pair<double, double> >&
RegressionVariant::get_w() const
{
  return w;
}

inline void
RegressionVariant::clear_w()
{
  w.clear();
}

inline void
RegressionVariant::set_name(const string& n)
{
  my_name = n;
}

//
//------------------------------------------------------------------------
//

inline void
TraitRegression::validate()
{
  my_valid_regression = true;
}

inline void
TraitRegression::set_regression_method(RegressionVariant* r)
{
  my_regression_method.reset(r);
}

inline void
TraitRegression::set_regression_method(regression_variant_pointer r)
{
  my_regression_method = r;
}

inline void
TraitRegression::disown_regression_method()
{
  my_regression_method.reset();
}

inline void
TraitRegression::set_model(const regression_model& p)
{
  invalidate();
  reset();
  my_model = p;
}

inline void
TraitRegression::set_sib_clusters(const vector<sib_cluster>& sc)
{
  my_sib_clusters = sc;
}

inline void
TraitRegression::set_filter(const pair_filter& pf)
{
  my_filter = pf;
}

inline void
TraitRegression::set_use_pairs(const pair<bool, bool>& up)
{
  my_use_fsib = up.first;
  my_use_hsib = up.second;
}

inline void
TraitRegression::set_pair_counts(const pair<size_t, size_t>& pc)
{
  my_fsib_pair_count = pc.first;
  my_hsib_pair_count = pc.second;
}

inline void
TraitRegression::set_use_x_pairs(const vector<bool>& x)
{
  my_use_x_pair = x;
}

inline void
TraitRegression::set_fix_x_pairs(const vector<bool>& x)
{
  my_fix_x_pair = x;
}

inline void
TraitRegression::set_x_pair_counts(const vector<size_t>& pc)
{
  my_x_pair_count = pc;
}

inline void
TraitRegression::set_x_types(const vector<pair_type>& pt)
{
  my_x_types = pt;
}

inline void
TraitRegression::set_reg_results(const regression_results& r)
{
  my_reg_results = r;
}

//--------------------------------------------------------------

inline relative_pairs&
TraitRegression::get_pairs() const
{
  return pairs;
}

inline RegressionVariant*
TraitRegression::regression_method()
{
  return my_regression_method.get();
}

inline const RegressionVariant*
TraitRegression::regression_method() const
{
  return my_regression_method.get();
}

inline regression_model&
TraitRegression::get_model()
{
  return my_model;
}

inline const regression_model&
TraitRegression::get_model() const
{
  return my_model;
}

inline const vector<sib_cluster>&
TraitRegression::get_sib_clusters() const
{
  return my_sib_clusters;
}

inline const pair_filter&
TraitRegression::get_filter() const
{
  return my_filter;
}

inline pair<bool, bool>
TraitRegression::get_use_pairs() const
{
  return make_pair(my_use_fsib, my_use_hsib);
}

inline pair<size_t, size_t>
TraitRegression::get_pair_counts() const
{
  return make_pair(my_fsib_pair_count, my_hsib_pair_count);
}

inline const vector<bool>&
TraitRegression::get_use_x_pairs() const
{
  return my_use_x_pair;
}

inline const vector<bool>&
TraitRegression::get_fix_x_pairs() const
{
  return my_fix_x_pair;
}

inline const vector<size_t>&
TraitRegression::get_x_pair_counts() const
{
  return my_x_pair_count;
}

inline const vector<pair_type>&
TraitRegression::get_x_types() const
{
  return my_x_types;
}

inline regression_results&
TraitRegression::get_reg_results()
{
  return my_reg_results;
}

inline const regression_results&
TraitRegression::get_reg_results() const
{
  return my_reg_results;
}

inline regression_results&
TraitRegression::get_reg_intercept_results()
{
  return my_reg_results_intercept;
}

inline const regression_results&
TraitRegression::get_reg_intercept_results() const
{
  return my_reg_results_intercept;
}

inline bool
TraitRegression::built() const
{
  if( !get_model().valid() )
    return false;

  return my_built;
}

inline bool
TraitRegression::valid() const
{
  return my_valid_regression;
}

inline size_t
TraitRegression::pair_count() const
{
  return my_fsib_pair_count + my_hsib_pair_count;
}

inline size_t
TraitRegression::fsib_pair_count() const
{
  return my_fsib_pair_count;
}

inline size_t
TraitRegression::hsib_pair_count() const
{
  return my_hsib_pair_count;
}

inline size_t
TraitRegression::valid_pair_type_count() const
{
  if( my_use_fsib && my_use_hsib )
    return 2;
  else if( my_use_fsib || my_use_hsib )
    return 1;
  else
    return 0;
}

inline size_t
TraitRegression::valid_intercept_count() const
{
  if( my_use_fsib )
  {
    if( get_model().is_x_linked() )
      return my_x_intercept[MM].size();
    else
      return my_fsib_intercept.size();
  }

  return my_hsib_intercept.size();
}

inline size_t
TraitRegression::mm_pair_count() const
{
  if( !my_x_pair_count.size() )
    return 0;

  return my_x_pair_count[MM];
}

inline size_t
TraitRegression::mf_pair_count() const
{
  if( !my_x_pair_count.size() )
    return 0;

  return my_x_pair_count[MF];
}

inline size_t
TraitRegression::ff_pair_count() const
{
  if( !my_x_pair_count.size() )
    return 0;

  return my_x_pair_count[FF];
}

inline size_t
TraitRegression::valid_trait_count() const
{
   return my_valid_trait_count;
}

inline size_t
TraitRegression::valid_parameter_count() const
{
   return my_valid_parameter_count;
}

inline size_t
TraitRegression::parameter_count() const
{
  return get_model().get_parameter_count() + 1;
}
/*
inline double
TraitRegression::rmse() const
{
  double ss = get_residual_variance();

  if(finite(ss) && ss > 0.0)
    return sqrt(ss);

  return numeric_limits<double>::quiet_NaN();
}
*/
inline size_t
TraitRegression::get_svd_return_code() const
{
  return my_gls.get_svd_return_code();
}

inline cerrorstream&
TraitRegression::error_stream()
{
  return errors;
}

inline void
TraitRegression::set_error_stream(cerrorstream& e)
{
  errors = e;
}

//----------------------------------------------------------------

inline double
TraitRegression::member_trait(const MPED::member_base* m, size_t t) const
{
  return pairs.trait(*m, t);
}

inline pair<double,double>
TraitRegression::get_trait_means(const sib_pair& spair) const
{
  if( get_model().is_x_linked() )
  {
    return get_trait_means_x(spair);
  }
   
  double m = get_model().get_trait().get_mean();
  
  return pair<double,double>(m, m);
}

inline pair<double,double>
TraitRegression::get_sibship_means(const sib_pair& spair) const
{
  if( get_model().is_x_linked() )
  {
    return get_sibship_means_x(spair);
  }
   
  double m1 = get_model().get_trait().trait_all_sibs_info.mean();
  double m2 = get_model().get_trait().trait_all_sibs_info.mean();

  sib_mean_map::const_iterator s1 = my_sib_mean_map.find(spair.rels().pair.first);

  if( s1 != my_sib_mean_map.end() )
    m1 = s1->second.first;

  sib_mean_map::const_iterator s2 = my_sib_mean_map.find(spair.rels().pair.second);

  if( s2 != my_sib_mean_map.end() )
    m2 = s2->second.first;

  return pair<double,double>(m1, m2);
}

inline pair<double,double>
TraitRegression::get_blup_means(const sib_pair& spair) const
{
  if( get_model().is_x_linked() )
  {
    return get_blup_means_x(spair);
  }
   
  double m1 = get_model().get_trait().trait_all_sibs_info.mean();
  double m2 = get_model().get_trait().trait_all_sibs_info.mean();

  sib_mean_map::const_iterator s1 = my_sib_mean_map.find(spair.rels().pair.first);

  if( s1 != my_sib_mean_map.end() )
    m1 = s1->second.second;

  sib_mean_map::const_iterator s2 = my_sib_mean_map.find(spair.rels().pair.second);

  if( s2 != my_sib_mean_map.end() )
    m2 = s2->second.second;

  return pair<double,double>(m1, m2);
}

inline pair<double,double>
TraitRegression::get_trait_means_x(const sib_pair& spair) const
{
  double m1;
  
  if( spair.rels().pair.first->is_male() )
  {
    m1 = get_model().get_trait().fixed_male_mean;

    if( !finite(m1) )
      m1 = get_model().get_trait().trait_male_sibs_info.mean();
  }
  else
  {   
    m1 = get_model().get_trait().fixed_female_mean;

    if( !finite(m1) )
      m1 = get_model().get_trait().trait_female_sibs_info.mean();
  }
   
  double m2;

  if( spair.rels().pair.second->is_male() )
  {
    m2 = get_model().get_trait().fixed_male_mean;

    if( !finite(m2) )
      m2 = get_model().get_trait().trait_male_sibs_info.mean();
  }
  else
  {   
    m2 = get_model().get_trait().fixed_female_mean;

    if( !finite(m2) )
      m2 = get_model().get_trait().trait_female_sibs_info.mean();
  }
   
  return pair<double,double>(m1, m2);
}

inline pair<double,double>
TraitRegression::get_sibship_means_x(const sib_pair& spair) const
{
  double m1;

  if( spair.rels().pair.first->is_male() )
    m1 = get_model().get_trait().trait_male_sibs_info.mean();
  else
    m1 = get_model().get_trait().trait_female_sibs_info.mean();

  double m2;

  if( spair.rels().pair.second->is_male() )
    m2 = get_model().get_trait().trait_male_sibs_info.mean();
  else
    m2 = get_model().get_trait().trait_female_sibs_info.mean();

  sib_mean_map::const_iterator s1;

  if( spair.rels().pair.first->is_male() )
    s1 = my_male_sib_mean_map.find(spair.rels().pair.first);
  else
    s1 = my_female_sib_mean_map.find(spair.rels().pair.first);

  if( s1 != my_male_sib_mean_map.end() || s1 != my_female_sib_mean_map.end() )
    m1 = s1->second.first;

  sib_mean_map::const_iterator s2;

  if( spair.rels().pair.second->is_male() )
    s2 = my_male_sib_mean_map.find(spair.rels().pair.second);
  else
    s2 = my_female_sib_mean_map.find(spair.rels().pair.second);

  if( s2 != my_male_sib_mean_map.end() || s2 != my_female_sib_mean_map.end() )
    m2 = s2->second.first;

  return pair<double,double>(m1, m2);
}

inline pair<double,double>
TraitRegression::get_blup_means_x(const sib_pair& spair) const
{
  double m1;

  if( spair.rels().pair.first->is_male() )
    m1 = get_model().get_trait().trait_male_sibs_info.mean();
  else
    m1 = get_model().get_trait().trait_female_sibs_info.mean();

  double m2;

  if( spair.rels().pair.second->is_male() )
    m2 = get_model().get_trait().trait_male_sibs_info.mean();
  else
    m2 = get_model().get_trait().trait_female_sibs_info.mean();

  sib_mean_map::const_iterator s1;

  if( spair.rels().pair.first->is_male() )
    s1 = my_male_sib_mean_map.find(spair.rels().pair.first);
  else
    s1 = my_female_sib_mean_map.find(spair.rels().pair.first);

  if( s1 != my_male_sib_mean_map.end() || s1 != my_female_sib_mean_map.end() )
    m1 = s1->second.second;

  sib_mean_map::const_iterator s2;

  if( spair.rels().pair.second->is_male() )
    s2 = my_male_sib_mean_map.find(spair.rels().pair.second);
  else
    s2 = my_female_sib_mean_map.find(spair.rels().pair.second);

  if( s2 != my_male_sib_mean_map.end() || s2 != my_female_sib_mean_map.end() )
    m2 = s2->second.second;

  return pair<double,double>(m1, m2);
}

//----------------------------------------------------------------

inline bool 
TraitRegression::params_dom_equal(const independent_variable& total_param,
                                  const independent_variable& param) const
{
  if( total_param.covariates != param.covariates )
    return false;
  if( total_param.markers.size() != param.markers.size() )
    return false;
  for( size_t i = 0; i < total_param.markers.size(); ++i )
    if( total_param.markers[i].marker_index != param.markers[i].marker_index )
      return false;
  return true;
}

inline pair<double,double>
TraitRegression::pair_traits(const sib_pair&        spair,
                             const trait_parameter& param) const
{
  const double t1 = member_trait( spair.rels().pair.first,  param.trait_index );
  const double t2 = member_trait( spair.rels().pair.second, param.trait_index );

  return pair<double,double>(t1,t2);
}

inline bool
TraitRegression::simulating() const
{
  return my_simulating;
}

inline const sib_pair
TraitRegression::permutation_pair(const sib_pair& pr) const
{
 if( my_group_permutation_vector.size() && my_simulation_map.size() )
  {
    const size_t pair_number       = pr.pair_number();
    assert(pair_number < my_simulation_map.size());
    const pair_index& index        = my_simulation_map[pair_number];

    assert(index.sibship_size   < my_group_permutation_vector.size() );
    assert(index.sibship_number < my_group_permutation_vector[index.sibship_size].size());
    assert(index.sib_number     < my_group_permutation_vector[index.sibship_size][index.sibship_number].size());

    size_t new_pair_number = my_group_permutation_vector[index.sibship_size][index.sibship_number][index.sib_number];

    assert(new_pair_number < pairs.pair_count());

    return sib_pair(new_pair_number, &pairs);
  }
  return pr;
}
