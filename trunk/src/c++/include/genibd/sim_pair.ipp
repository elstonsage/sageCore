// ======================================
// Inline functions of sim_relative_pair
// ======================================

inline
fmember_const_pointer
sim_relative_pair::get_first_ind() const
{
  return my_member_one;
}

inline
fmember_const_pointer
sim_relative_pair::get_second_ind() const
{
  return my_member_two;
}

inline
bool
sim_relative_pair::is_sib() const
{
  return (my_pair_type == pair_generator::SIBSIB);
}

inline
bool
sim_relative_pair::is_hsib() const
{
  return (my_pair_type == pair_generator::HALFSIB);
}

inline
bool
sim_relative_pair::is_brother_brother() const
{
/*
  if( !is_sib() && !is_hsib() )
    return false;

  if(    my_member_one->is_male()
      && my_member_two->is_male() )
    return true;

  return false;
*/
  return SAGE::is_brother_brother(my_member_one, my_member_two);
}

inline
bool
sim_relative_pair::is_sister_sister() const
{
  if( !is_sib() && !is_hsib() )
    return false;

  if(    my_member_one->is_female()
      && my_member_two->is_female() )
    return true;

  return false;
}

inline
void
sim_relative_pair::set_last_sharing(size_t m, size_t i)
{
  my_data[m].last_sharing = i;
}

inline
void
sim_relative_pair::set_last_sib_sharing(size_t m, size_t i)
{
  my_data[m].last_sib_sharing = i;
}

inline
size_t
sim_relative_pair::get_last_sharing(size_t m) const
{
  return my_data[m].last_sharing;
}

inline
size_t
sim_relative_pair::get_last_sib_sharing(size_t m) const
{
  return my_data[m].last_sib_sharing;
}
  
inline
void
sim_relative_pair::set_values(size_t m, double f0, double f2, size_t steps)
{
  marker_data& d = my_data[m];

  d.f0 = size_t(steps * f0);
  d.f2 = size_t(steps * f2);

  d.total = steps;

  d.last_sharing = (size_t) - 1;
}

inline
void
sim_relative_pair::set_values(size_t m, double f0, double f1mp, double f2, size_t steps)
{
  marker_data& d = my_data[m];

  d.f0  = size_t(steps * f0);
  d.f2  = size_t(steps * f2);

  double f1  = 1. - (f0 + f2);
  double f1m = (f1 + f1mp) / 2.0;
  double f1p = (f1 - f1mp) / 2.0;

  d.f1m = size_t(steps * f1m);
  d.f1p = size_t(steps * f1p);

  d.total = steps;

  d.last_sharing     = (size_t) - 1;
  d.last_sib_sharing = (size_t) - 1;
}

inline
double
sim_relative_pair::get_value(size_t m, size_t i) const
{
  if( my_data[m].total == 0 )
    return std::numeric_limits<double>::quiet_NaN();

  switch(i)
  {
    case 0 : return (double) my_data[m].f0 / (double) my_data[m].total;
    case 2 : return (double) my_data[m].f2 / (double) my_data[m].total;
    case 1 : if( !(is_sib() || is_hsib()) )
               return QNAN;

             return (double) my_data[m].f1m / (double) my_data[m].total
                  - (double) my_data[m].f1p / (double) my_data[m].total;

    default: size_t t = my_data[m].total - my_data[m].f0 - my_data[m].f2;
             return (double) t             / (double) my_data[m].total;
  }
}

inline
void
sim_relative_pair::increment_value(size_type m, size_t increment, bool x_linked)
{
  size_t last_sharing = my_data[m].last_sharing;

#if 0
  cout << my_member_one->name() << ":" << my_member_two->name()
       << " sharing " << last_sharing
       << ", last_sib_sharing = " << my_data[m].last_sib_sharing << endl;
#endif

  switch(last_sharing)
  {
    case 0 : my_data[m].f0 += increment; break;
    case 2 : my_data[m].f2 += increment; break;

    case 1 : if( is_sib() )
               if( my_data[m].last_sib_sharing == 0 )
                 my_data[m].f1m += increment;
               else
                 if( !x_linked || is_sister_sister() )
                   my_data[m].f1p += increment;

             else if( is_hsib() )
               if( my_data[m].last_sib_sharing == 0 )
                 my_data[m].f1m += increment;
               else
                 if( !x_linked || is_sister_sister() )
                   my_data[m].f1p += increment;

             break;

    default: return;
  }

  my_data[m].total += increment;

#if 0
 cout << my_data[m].f0 << ", "  << my_data[m].f1m << ", "
      << my_data[m].f1p << ", " << my_data[m].f2 << ", "
      << my_data[m].total << endl;
#endif
}
