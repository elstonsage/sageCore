/// Static functions for use by the iterative models.
///
/// These functions are to make the choice between discrete and continuous
/// models a little easier.
///
/// NOTE: They currently only distinguish between discrete (binary) and
///       continuous models.  They do not, at present, know anything
///       about age of onset variables which are both.

//@{

inline bool is_model_continuous(const model& m)
{
  return m.get_primary_trait_type() == pt_CONTINUOUS ||
         (m.get_primary_trait_type() == pt_ONSET &&
          m.ons_sub_model.t_option() == onset_sub_model::t_A);
}

inline bool is_model_binary(const model& m)
{
  return m.get_primary_trait_type() == pt_BINARY ||
         (m.get_primary_trait_type() == pt_ONSET &&
          m.ons_sub_model.t_option() == onset_sub_model::t_S);
}

inline bool is_model_onset(const model& m)
{
  return m.get_primary_trait_type() == pt_ONSET;
}

/// get_type_model() returns the mean or susceptibility sub_model as
/// necessary based on model type

inline
const genotype_specific_mean_susc_sub_model& get_type_model(const model& m)
{
  return m.type_dependent_sub_model();
}

/// get_type_model() returns the mean or susceptibility sub_model as
/// necessary based on model type

inline
genotype_specific_mean_susc_sub_model& get_type_model(model& m)
{
  return m.type_dependent_sub_model();
}

inline
string singular(const model& m)
{
  const string mean_singular("mean");
  const string susc_singular("susc.");

  if(is_model_continuous(m)) return mean_singular;
  else                       return susc_singular;
}

inline
string plural(const model& m)
{
  const string mean_plural("means");
  const string susc_plural("suscs.");

  if(is_model_continuous(m)) return mean_plural;
  else                       return susc_plural;
}

//@}


