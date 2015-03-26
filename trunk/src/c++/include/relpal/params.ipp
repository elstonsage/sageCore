////////////////////////////////////////////////////////////////////////////
//  Implementation of params.h (Inline)
////////////////////////////////////////////////////////////////////////////

inline bool
dependent_variable::operator==(const dependent_variable& t) const
{
  if( trait_index == (size_t)-1 || t.trait_index == (size_t)-1 )
    return false;

  return trait_index == t.trait_index;
}

inline bool
dependent_variable::operator!=(const dependent_variable& t) const
{
  return !( (*this) == t );
}

//
//------------------------------------------------------------------------
//

inline bool
covariate_type::operator==(const covariate_type& c) const
{
  return    covariate_index == c.covariate_index
         && adj_trait_index == c.adj_trait_index
         && test_variable   == c.test_variable;
}

inline bool
covariate_type::operator!=(const covariate_type& c) const
{
  return !( (*this) == c );
}

inline bool
covariate_type::operator<(const covariate_type& c) const
{
  if( covariate_index < c.covariate_index )
    return true;
  if( covariate_index > c.covariate_index )
    return false;
  return test_variable < c.test_variable;
}

//
//------------------------------------------------------------------------
//

inline bool
marker_type::operator==(const marker_type& m) const
{
  return marker_index == m.marker_index && test_variable == m.test_variable;
}

inline bool
marker_type::operator!=(const marker_type& m) const
{
  return !( (*this) == m);
}

inline bool
marker_type::operator<(const marker_type& m) const
{
  if( marker_index < m.marker_index )
    return true;
  if( marker_index > m.marker_index )
    return false;
  return (test_variable < m.test_variable);
}

//
//------------------------------------------------------------------------
//

inline bool
interaction_type::operator==(const interaction_type& c) const
{
    return    type == c.type
           && covariates == c.covariates
           && markers == c.markers
           && test_variable == c.test_variable
           && batch == c.batch;
}

inline bool
interaction_type::operator!=(const interaction_type& c) const
{
  return !( (*this) == c );
}

inline bool
interaction_type::operator<(const interaction_type& c) const
{
  if(type < c.type)
    return true;
  if(type > c.type)
    return false;

  if(covariates < c.covariates)
    return true;
  if(covariates > c.covariates)
    return false;

  return (markers < c.markers);
}

//
//------------------------------------------------------------------------
//

inline bool
independent_variable::operator==(const independent_variable& c) const
{
    return    type == c.type
           && covariates == c.covariates
           && markers == c.markers
           && t1 == c.t1
           && t2 == c.t2
           && test_variable == c.test_variable;
}

inline bool
independent_variable::operator!=(const independent_variable& c) const
{
  return !( (*this) == c );
}

inline bool
independent_variable::operator<(const independent_variable& c) const
{
  if(type < c.type)
    return true;
  if(type > c.type)
    return false;

  if(covariates < c.covariates)
    return true;
  if(covariates > c.covariates)
    return false;

  return (markers < c.markers);
}

//
//------------------------------------------------------------------------
//

inline bool
filtering_type::operator==(const filtering_type& t) const
{
  if( subset_index == (size_t)-1 || t.subset_index == (size_t)-1 )
    return false;

  return subset_index == t.subset_index;
}

inline bool
filtering_type::operator!=(const filtering_type& t) const
{
  return !( (*this) == t );
}
