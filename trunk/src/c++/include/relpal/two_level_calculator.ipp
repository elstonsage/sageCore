////////////////////////////////////////////////////////////////////////////
//             Implementation of two_level_calculator.h (Inline)          //
////////////////////////////////////////////////////////////////////////////

inline void
relpal_least_square::reset()
{
  reset(parameters);
}

inline void
relpal_score::reset()
{
  reset(parameters);
}

inline bool
relpal_score::is_naive_var() const
{
  return my_naive_var;
}

inline bool
relpal_score::is_sandwich_var() const
{
  return my_sandwich_var;
}

inline bool
relpal_score::is_alternative_var() const
{
  return my_alternative_var;
}

inline bool
relpal_score::is_IBD_var() const
{
  return my_IBD_var;
}

inline void
relpal_score::set_naive_var(bool v)
{
  my_naive_var = v;
}

inline void
relpal_score::set_sandwich_var(bool v)
{
  my_sandwich_var = v;
}

inline void
relpal_score::set_alternative_var(bool v)
{
  my_alternative_var = v;
}

inline void
relpal_score::set_IBD_var(bool v)
{
  my_IBD_var = v;
}
