// ------------------------------------------------------
// Inline Implementation of mcmc_parameters
// ------------------------------------------------------

inline bool mcmc_parameters::is_multipoint() const  
{
  return my_multipoint;
}

inline size_t mcmc_parameters::get_max_tunnel() const
{
  return  my_max_tunnel;  
}

inline double mcmc_parameters::get_single_marker_tunneling() const
{
  return my_single_marker;  
}

inline double mcmc_parameters::get_local_marker() const
{
  return my_local_marker;
}

inline double mcmc_parameters::get_local_individual() const
{
  return my_local_individual;
}

inline double mcmc_parameters::get_transition_weight(size_t i) const
{
  return my_T[i];
}

inline long mcmc_parameters::get_dememorization_step() const  
{
  return my_dememorization_step;
}

inline long mcmc_parameters::get_simulation_step()const  
{
  return my_simulation_step;
}

inline long mcmc_parameters::get_batch_count()const  
{
  return  my_batch_count;
}

inline bool mcmc_parameters::get_use_factor() const
{
  return my_use_factor;
}

inline double mcmc_parameters::get_dememorization_factor() const
{
  return my_dememorization_factor;
}

inline double mcmc_parameters::get_simulation_factor() const
{
  return my_simulation_factor;
}

inline double mcmc_parameters::get_batch_factor() const
{
  return my_batch_factor;
}

inline unsigned long mcmc_parameters::get_random_seed() const
{
  return my_random_seed;
}

inline bool mcmc_parameters::set_multipoint(bool m) 
{
  my_multipoint = m;

  return true;
}

inline bool mcmc_parameters::set_max_tunnel(size_t p)
{
  my_max_tunnel = p;

  return true;
}

inline bool mcmc_parameters::set_single_marker_tunneling(double p)
{
  my_single_marker = p;

  return true;
}

inline bool mcmc_parameters::set_local_marker(double p)
{
  my_local_marker = p;

  return true;
}

inline bool mcmc_parameters::set_local_individual(double p)
{
  my_local_individual = p;

  return true;
}

inline bool mcmc_parameters::set_transition_weight(size_t i, double p)
{
  if(p < 0.0 || p > 1.0 || i > 2)
    return false;

  my_T[i] = p;

  return true;
} 

inline
bool mcmc_parameters::set_dememorization_step(long p ) 
{
  if(p<0)
    return false;
  
  my_dememorization_step = p;

  return true;
}

inline
bool mcmc_parameters::set_simulation_step(long p)
{
  if(p<=0) return false;

  my_simulation_step=p;

  return true;

}    

inline
bool mcmc_parameters::set_batch_count(long p)
{
  if(p<=0) return false;

  my_batch_count=p;

  return true;

}    

inline
bool mcmc_parameters::set_use_factor(bool p) 
{
  my_use_factor = p;

  return true;
}

inline
bool mcmc_parameters::set_dememorization_factor(double p ) 
{
  if(p<0) return false;

  my_dememorization_factor = p;

  return true;
}

inline
bool mcmc_parameters::set_simulation_factor(double p)
{
  if(p<=0) return false;

  my_simulation_factor = p;

  return true;
}    

inline
bool mcmc_parameters::set_batch_factor(double p)
{
  if(p<=0) return false;

  my_batch_factor = p;

  return true;
}    

inline
bool mcmc_parameters::set_random_seed(unsigned long p)
{
  my_random_seed = p;

  return true;
}

inline void mcmc_parameters::normalize_weights()
{
  double sum = my_T[0] + my_T[1] + my_T[2];
  
  if(sum < numeric_limits<double>::epsilon())
  {
    // Possibly produce error message?
  
    my_T[0] = 1.0;
    my_T[1] = 0.0;
    my_T[2] = 0.0;
    
    return;
  }

  // Possibly produce warning if sum differs significantly from 1.0?
  
  my_T[0] /= sum;
  my_T[1] /= sum;
  my_T[2] /= sum;
}

