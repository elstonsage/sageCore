////////////////////////////////////////////////////////////////////////////
//             Implementation of lodpal_params (Inline)                   //
////////////////////////////////////////////////////////////////////////////

//
//--------------------------------------------------------------------------
//

inline
void
parameter_estimate::set_value(double v)
{
  my_value = v;
}

inline
void
parameter_estimate::set_initial_value(double v)
{
  my_initial_value = v;
}

inline
void
parameter_estimate::set_first_derivative(double d)
{
  my_first_derivative = d;
}

inline
void
parameter_estimate::set_stderr(double d)
{
  my_stderr = d;
}

inline
void
parameter_estimate::fix_value()
{
  my_fixed = true;
}

inline
void
parameter_estimate::fix_value(double v)
{
  my_value = v;
  my_fixed = true;
}

inline
void
parameter_estimate::release_value()
{
  my_fixed = false;
}

inline
double
parameter_estimate::value() const
{
  return my_value;
}

inline
double
parameter_estimate::initial_value() const
{
  return my_initial_value;
}

inline
double
parameter_estimate::first_derivative() const
{
  return my_first_derivative;
}

inline
double
parameter_estimate::get_stderr() const
{
  return my_stderr;
}

inline
bool
parameter_estimate::fixed() const
{
  return my_fixed;
}

//
//--------------------------------------------------------------------------
//

inline
bool
trait_parameter::operator==(const trait_parameter& t) const
{
  if( trait == (size_t)-1 || t.trait == (size_t)-1 )
    return false;
  return trait == t.trait;
}

inline
bool
trait_parameter::operator!=(const trait_parameter& t) const
{
  return !( (*this) == t );
}

//
//--------------------------------------------------------------------------
//

inline
bool
autosomal_model_type::operator==(const autosomal_model_type& m) const
{
  return    model == m.model
         && constraint == m.constraint
         && fixed == m.fixed
         && parent_of_origin == m.parent_of_origin
         && alpha == m.alpha;
}

inline
bool
autosomal_model_type::operator!=(const autosomal_model_type& m) const
{
  return !( (*this) == m);
}

//
//--------------------------------------------------------------------------
//

inline
bool
x_linkage_model_type::operator==(const x_linkage_model_type& m) const
{
  return   (lambda1_equal == m.lambda1_equal && lambda2_fixed == m.lambda2_fixed)
         && alpha == m.alpha;
}

inline
bool
x_linkage_model_type::operator!=(const x_linkage_model_type& m) const
{
  return !( (*this) == m);
}

//
//--------------------------------------------------------------------------
//

inline
bool
marker_type::operator==(const marker_type& m) const
{
  return marker == m.marker;
}

inline
bool
marker_type::operator!=(const marker_type& m) const
{
  return !( (*this) == m);
}

//
//--------------------------------------------------------------------------
//

inline
bool
covariate_type::operator==(const covariate_type& c) const
{
  return covariate == c.covariate && power == c.power &&
         operation == c.operation && adjust == c.adjust &&
         adjust_value == c.adjust_value;
}

inline
bool
covariate_type::operator!=(const covariate_type& c) const
{
  return !( (*this) == c );
}

//
//--------------------------------------------------------------------------
//

inline
bool
lodpal_weight_type::operator==(const lodpal_weight_type& w) const
{
  return weight == w.weight &&
         operation == w.operation;
}

inline
bool
lodpal_weight_type::operator!=(const lodpal_weight_type& w) const
{
  return !( (*this) == w );
}

//
//--------------------------------------------------------------------------
//

inline
bool
independent_variable::operator==(const independent_variable& c) const
{
    return covariates == c.covariates && markers == c.markers;
}

inline
bool
independent_variable::operator!=(const independent_variable& c) const
{
  return !( (*this) == c );
}

//
//--------------------------------------------------------------------------
//

inline
bool
lodpal_parameters::valid() const
{
  return my_valid;
}

inline
bool
lodpal_parameters::use_pair_cache() const
{
  return my_use_pair_cache;
}

inline
void
lodpal_parameters::set_use_pair_cache(bool c)
{
  my_use_pair_cache = c;
}

inline
bool
lodpal_parameters::skip_uninformative_pairs() const
{
  return my_skip_uninformative_pairs;
}

inline
void
lodpal_parameters::set_skip_uninformative_pairs(bool s)
{
  if(s != my_skip_uninformative_pairs)
  {
    invalidate();
    my_skip_uninformative_pairs = s;
  }
}

inline
bool
lodpal_parameters::turn_off_default() const
{
  return my_turn_off_default;
}

inline
void
lodpal_parameters::set_turn_off_default(bool d)
{
  my_turn_off_default = d;
}

inline
bool
lodpal_parameters::multipoint() const
{
  return my_multipoint;
}

inline
void
lodpal_parameters::set_multipoint(bool m)
{
  my_multipoint = m;
}

inline
bool
lodpal_parameters::print_lambda() const
{
  return my_print_lambda;
}

inline
void
lodpal_parameters::set_print_lambda(bool m)
{
  my_print_lambda = m;
}

inline
bool
lodpal_parameters::sib_pairs_only() const
{
  return my_sib_pairs_only;
}

inline
void
lodpal_parameters::set_sib_pairs_only(bool s)
{
  my_sib_pairs_only = s;
}

inline
bool
lodpal_parameters::autosomal_marker_exist() const
{
  return my_autosomal_marker_exist;
}

inline
void
lodpal_parameters::set_autosomal_marker_exist(bool s)
{
  my_autosomal_marker_exist = s;
}

inline
bool
lodpal_parameters::x_linked_marker_exist() const
{
  return my_x_linked_marker_exist;
}

inline
void
lodpal_parameters::set_x_linked_marker_exist(bool s)
{
  my_x_linked_marker_exist = s;
}

inline
bool
lodpal_parameters::use_mm_pair() const
{
  return my_use_mm_pair;
}

inline
bool
lodpal_parameters::use_mf_pair() const
{
  return my_use_mf_pair;
}

inline
bool
lodpal_parameters::use_ff_pair() const
{
  return my_use_ff_pair;
}

inline
size_t
lodpal_parameters::used_pair_type() const
{
  return my_total_used_pair_type;
}

inline
void
lodpal_parameters::set_x_linkage_pair_type(bool mm, bool mf, bool ff)
{
  my_use_mm_pair = mm;
  my_use_mf_pair = mf;
  my_use_ff_pair = ff;

  my_total_used_pair_type = 0;
  if( my_use_mm_pair )
    my_total_used_pair_type += 1;
  if( my_use_mf_pair )
    my_total_used_pair_type += 1;
  if( my_use_ff_pair )
    my_total_used_pair_type += 1;
}

inline
size_t
lodpal_parameters::diagnostic_marker() const
{
  return my_diagnostic_marker;
}

inline
void
lodpal_parameters::set_diagnostic_marker(size_t m)
{
  my_diagnostic_marker = m;
}

inline
size_t
lodpal_parameters::max_marker_name_size() const
{
  return my_max_marker_name_size;
}

inline
void
lodpal_parameters::set_max_marker_name_size(size_t m)
{
  my_max_marker_name_size = m;
}

inline
size_t
lodpal_parameters::trait_count() const
{
  return my_traits.size();
}

inline
size_t
lodpal_parameters::subset_count() const
{
  return my_subsets.size();
}

inline
size_t
lodpal_parameters::marker_count() const
{
  return my_parameters.markers.size();
}

inline
size_t
lodpal_parameters::covariate_count() const
{
  return my_parameters.covariates.size();
}

inline
size_t
lodpal_parameters::parameter_count() const
{
  return marker_count() + covariate_count();
}

inline
lodpal_parameters::marker_iterator
lodpal_parameters::marker_begin()
{
   return my_parameters.markers.begin();
}

inline
lodpal_parameters::marker_iterator
lodpal_parameters::marker_end()
{
   return my_parameters.markers.end();
}

inline
lodpal_parameters::covariate_iterator
lodpal_parameters::covariate_begin()
{
   return my_parameters.covariates.begin();
}

inline
lodpal_parameters::covariate_iterator
lodpal_parameters::covariate_end()
{
   return my_parameters.covariates.end();
}

inline
lodpal_parameters::marker_const_iterator
lodpal_parameters::marker_begin() const
{
   return my_parameters.markers.begin();
}

inline
lodpal_parameters::marker_const_iterator
lodpal_parameters::marker_end() const
{
   return my_parameters.markers.end();
}

inline
lodpal_parameters::covariate_const_iterator
lodpal_parameters::covariate_begin() const
{
   return my_parameters.covariates.begin();
}

inline
lodpal_parameters::covariate_const_iterator
lodpal_parameters::covariate_end() const
{
   return my_parameters.covariates.end();
}

inline
const trait_parameter&
lodpal_parameters::trait_parameters(size_t t) const
{
  return my_traits[t];
}

inline
const trait_parameter&
lodpal_parameters::subset_parameters(size_t t) const
{
  return my_subsets[t];
}

inline
const marker_type&
lodpal_parameters::marker_parameters(size_t m) const
{
  return my_parameters.markers[m];
}

inline
const covariate_type&
lodpal_parameters::covariate_parameters(size_t c) const
{
  return my_parameters.covariates[c];
}

inline
const independent_variable&
lodpal_parameters::parameters() const
{
  return my_parameters;
}

inline
const lodpal_weight_type&
lodpal_parameters::weight_parameter() const
{
  return my_weight;
}

inline
const autosomal_model_type&
lodpal_parameters::autosomal_model() const
{
  return my_autosomal_model;
}

inline
const x_linkage_model_type&
lodpal_parameters::x_linkage_model() const
{
  return my_x_linkage_model;
}

inline
trait_parameter&
lodpal_parameters::trait_parameters(size_t t)
{
  return my_traits[t];
}

inline
trait_parameter&
lodpal_parameters::subset_parameters(size_t t)
{
  return my_subsets[t];
}

inline
marker_type&
lodpal_parameters::marker_parameters(size_t m)
{
  return my_parameters.markers[m];
}

inline
covariate_type&
lodpal_parameters::covariate_parameters(size_t c)
{
  return my_parameters.covariates[c];
}

inline
independent_variable&
lodpal_parameters::parameters()
{
  return my_parameters;
}

inline
lodpal_weight_type&
lodpal_parameters::weight_parameter()
{
  return my_weight;
}

inline
autosomal_model_type&
lodpal_parameters::autosomal_model()
{
  return my_autosomal_model;
}

inline
x_linkage_model_type&
lodpal_parameters::x_linkage_model()
{
  return my_x_linkage_model;
}

inline
void
lodpal_parameters::clear_traits()
{
  if(my_traits.size())
  {
    my_traits.clear();
    invalidate();
  }
}

inline
lodpal_parameters::tib_value
lodpal_parameters::set_trait(const trait_parameter& t)
{
  clear_traits();
  return add_trait(t);
}

inline
lodpal_parameters::tib_value
lodpal_parameters::add_trait(size_t t, trait_parameter::trait_type tt, double cp)
{
  return add_trait( trait_parameter(t, tt, cp) );
}

inline
lodpal_parameters::tib_value
lodpal_parameters::set_trait(size_t t, trait_parameter::trait_type tt, double cp)
{
  invalidate();
  clear_traits();
  return add_trait(t, tt, cp);
}

inline
void
lodpal_parameters::clear_subsets()
{
  if(my_subsets.size())
  {
    my_subsets.clear();
    invalidate();
  }
}

inline
lodpal_parameters::tib_value
lodpal_parameters::set_subset(const trait_parameter& t)
{
  clear_subsets();
  return add_subset(t);
}

inline
lodpal_parameters::tib_value
lodpal_parameters::add_subset(size_t t)
{
  return add_subset( trait_parameter(t) );
}

inline
lodpal_parameters::tib_value
lodpal_parameters::set_subset(size_t t)
{
  invalidate();
  clear_subsets();
  return add_subset(t);
}

inline
void
lodpal_parameters::clear_parameters()
{
  clear_markers();
  clear_covariates();
}

inline
void
lodpal_parameters::clear_markers()
{
  my_parameters.markers.clear();
  invalidate();
}

inline
void
lodpal_parameters::clear_covariates()
{
  my_parameters.covariates.clear();
  invalidate();
}

inline
lodpal_parameters::mib_value
lodpal_parameters::add_marker(size_t m, marker_type::inheritance_type in, double b1, double b2)
{
  return add_parameter( marker_type(m, in, b1, b2) );
}

inline
lodpal_parameters::cib_value
lodpal_parameters::add_covariate(size_t c, covariate_type::cov_ad ad, covariate_type::cov_op op,
                                 double ad_val, double power, double d1, double d2)
{
  return add_parameter( covariate_type(c, ad, op, ad_val, power, d1, d2) );
}

inline
void
lodpal_parameters::clear_weight()
{
  my_weight = lodpal_weight_type();
  invalidate();
}

inline
bool
lodpal_parameters::set_weight(const lodpal_weight_type& w)
{
  clear_weight();
  return add_weight(w);
}

inline
bool
lodpal_parameters::add_weight(size_t w, lodpal_weight_type::weight_op ww)
{
  return add_weight( lodpal_weight_type(w, ww) );
}

inline
bool
lodpal_parameters::set_weight(size_t w, lodpal_weight_type::weight_op ww)
{
  invalidate();
  clear_weight();
  return add_weight(w, ww);
}

inline
void
lodpal_parameters::clear_autosomal_model()
{
  my_autosomal_model = autosomal_model_type();
  invalidate();
}

inline
bool
lodpal_parameters::set_autosomal_model(const autosomal_model_type& m)
{
  clear_autosomal_model();
  return add_autosomal_model(m);
}

inline
bool
lodpal_parameters::add_autosomal_model(autosomal_model_type::model_type      m,
                                       autosomal_model_type::constraint_type cont,
                                       autosomal_model_type::poo_fixed_type  fixed,
                                       bool poo_test, double alpha)
{
  return add_autosomal_model( autosomal_model_type(m, cont, fixed, poo_test, alpha) );
}

inline
bool
lodpal_parameters::set_autosomal_model(autosomal_model_type::model_type      m,
                                       autosomal_model_type::constraint_type cont,
                                       autosomal_model_type::poo_fixed_type  fixed,
                                       bool poo_test, double alpha)
{
  invalidate();
  clear_autosomal_model();
  return add_autosomal_model(m, cont, fixed, poo_test, alpha);
}

inline
void
lodpal_parameters::clear_x_linkage_model()
{
  my_x_linkage_model = x_linkage_model_type();
  invalidate();
}

inline
bool
lodpal_parameters::set_x_linkage_model(const x_linkage_model_type& m)
{
  clear_x_linkage_model();
  return add_x_linkage_model(m);
}

inline
bool
lodpal_parameters::add_x_linkage_model(bool l1, bool l2, double al)
{
  return add_x_linkage_model( x_linkage_model_type(l1, l2, al) );
}

inline
bool
lodpal_parameters::set_x_linkage_model(bool l1, bool l2, double al)
{
  invalidate();
  clear_x_linkage_model();
  return add_x_linkage_model(l1, l2, al);
}
