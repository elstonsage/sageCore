//=================================================
// Inline Implemantation: starting_state_generator
//=================================================

inline
int starting_state_generator::set_indicator_count(int i)
{
  return my_indicator_count = i;
}

inline
int starting_state_generator::get_indicator_count() const
{
  return my_indicator_count;
}

inline
int starting_state_generator::set_random_trial_count(int i)
{
  return my_trial_count = i;
}

inline
int starting_state_generator::get_random_trial_count() const
{
  return my_trial_count;
}

inline
bool starting_state_generator::is_valid() const
{
  return my_valid;
}
