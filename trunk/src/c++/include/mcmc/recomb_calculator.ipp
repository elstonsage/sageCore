//=================================================
// Inline Implemantation: recombination_calculator 
//=================================================

inline
double
recombination_calculator::operator() (size_t m) const
{
  return log_recombination(m);
}

inline
double
recombination_calculator::log_recombination(size_t m) const
{
  const bit_field& b1 = my_data.get_indicator(m);
  const bit_field& b2 = my_data.get_indicator(m+1);

  return log_recombination(m, b1, b2);
}

inline
double
recombination_calculator::operator() () const
{
  return log_recombination();
}

inline
double
recombination_calculator::log_recombination() const
{
  double ratio = 0.0;

  for(size_t i = 0; i < my_ped_region.get_region().locus_count()-1; ++i)
    ratio += log_recombination(i);

  return ratio;
}

inline
double
recombination_calculator::operator() (size_t m, const bit_field& b1, const bit_field& b2) const
{
  return log_recombination(m, b1, b2);
}

inline
double
recombination_calculator::log_recombination(size_t m, const bit_field& b1, const bit_field& b2) const
{
  double ratio = 0.0;

  for(size_t i = 0; i < b1.size(); ++i)
  {
    if(b1[i] == b2[i]) ratio += my_log_one_minus_thetas[m];
    else               ratio += my_log_thetas[m];
  }

  return ratio;
}

inline
double
recombination_calculator::log_recombination_ratio(size_t m, const bit_field& bl, const bit_field& br) const
{
  double ratio = 0.0;

  const bit_field& left  = my_data.get_indicator(m);
  const bit_field& right = my_data.get_indicator(m+1);

  double theta_ratio = my_theta_ratios[m];

  for(size_t i = 0; i < left.size(); ++i)
  {
    if(bl[i] ^ br[i])
    {
      bool l = left[i];
      bool r = right[i];

      if(l == r) ratio -= theta_ratio;
      else       ratio += theta_ratio;
    }
  }

  return ratio;  
}


