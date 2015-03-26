//=======================================
// Inline Implemantation: mcmc_simulator
//=======================================

inline
MersenneTwister&
mcmc_simulator::get_random() const
{
  return my_random_generator;
}

inline
double
mcmc_simulator::get_random_real() const
{
  return my_random_generator.uniform_real();
}

inline
int
mcmc_simulator::get_random_uniform_int(size_t t) const
{
  return my_random_generator.uniform_integer(t);
}

//selecting a locus to meet P(locus1->locus2)=P(locus2->locus1)
//according to the weight of selecting neighboring loci
inline
void
mcmc_simulator::select_a_locus(transition& t, bool& neighbor)  const
{   
  bool b = t.rule == ILLEGAL_SWITCH_RULE;

  // Should we pull a neighboring marker?
  neighbor = !b && get_random_real() < my_parameters->get_local_marker();

  // Less than 3 markers

  if( my_total_loci == 1 )
  {
    t.marker = 0;
    return;
  }

  if( my_total_loci <= 3 )
  {
    t.marker = get_random_uniform_int(my_total_loci);
    return;
  }

  // more than three markers.

  int last_marker = my_last_switch_marker;

  int new_marker = -1;

  if( neighbor )  //select neighboring loci
  { 
    double d = get_random_real();

         if(d < 1.0/3.0) new_marker = last_marker - 1;
    else if(d < 2.0/3.0) new_marker = last_marker;
    else                 new_marker = last_marker + 1;  
  }

  if( new_marker < 0 || new_marker >= (int) my_total_loci )
  {
    new_marker = get_random_uniform_int(my_total_loci);
  }

  t.marker = new_marker;
}

inline
void
mcmc_simulator::select_local_pivot_person(transition& t) const
{
  //cout << "Select Local!" << endl;

  const local_individual_type& c = my_local_individuals[my_last_pivot_person];

  t.person = c.local[get_random_uniform_int(c.local.size())];
}

inline
void
mcmc_simulator::select_pivot_person(transition& t) const
{
  int rand_num;

  if( get_random_real() > my_marker_individuals[t.marker].typed_probability ) //choose an untyped
  {
    size_t untyped = my_marker_individuals[t.marker].untyped_individual.size();
    rand_num       = get_random_uniform_int(untyped);
    t.person       = my_marker_individuals[t.marker].untyped_individual[rand_num];
  }
  else
  {
    size_t typed  = my_marker_individuals[t.marker].typed_individual.size();
    rand_num      = get_random_uniform_int(typed);
    t.person      = my_marker_individuals[t.marker].typed_individual[rand_num];
  }
}

inline
void
mcmc_simulator::select_trans_rule(transition& t) const
{
  bool composite = false;
  
  const FPED::Member& mem = my_pedigree.get_subpedigree().member_index(t.person);

  if(    mem.is_founder()
      || (mem.offspring_count() && get_random_real() < 0.2) )
    composite = true;

  if( !composite )
  {
    double d = get_random_real();

    if     ( d < 1.0/3.0 )    t.rule = T01;
    else if( d < 2.0/3.0 )    t.rule = T02;
    else                      t.rule = T03;
  }
  else
  {
    double d = get_random_real();

    if( d > 0.4 )             t.rule = T1;
    else
    {
      if( d > 0.2 )           t.rule = T2a;
      else                    t.rule = T2b;

      // If we're doing a T2, we must change to a nuclear family member for
      // the pivot person

      size_t n = my_nuclear_families[t.person].size();

      t.person = my_nuclear_families[t.person][get_random_uniform_int(n)];
    }
  }
}

inline
void
mcmc_simulator::select_extent(transition& t) const
{
  double d = get_random_real();

  // Choose our extent (either all left, all right or just current marker)
  // We may want to make the extent probability a user modifiable parameter,
  // but for now, it's fixed at 0.2

  if     (d < 0.2) t.extent = -1;
  else if(d < 0.8) t.extent =  0;
  else             t.extent =  1;
}

inline
void
mcmc_simulator::change_bit(size_t i, bool mother)
{
  const FPED::Member& mem = my_pedigree.get_subpedigree().member_index(i);
  
  if( mother )
  {
    my_transition_pattern[my_pedigree.get_mother_meiosis(mem)].flip();
  }
  else
  {
    my_transition_pattern[my_pedigree.get_father_meiosis(mem)].flip();
  }
}

inline
void
mcmc_simulator::apply_bit_pattern_at_marker(int marker)
{
  my_data->get_indicator(marker) ^= my_bits_changed[marker];
}

inline
double
mcmc_simulator::log_recomb_prob_ratio() const
{
  // Calculate the inheritance_likelihood_ratio
  double recomb_prob_ratio = 0.0;

  for( size_t i = 0; i < my_markers_used.size() - 1; ++i )
  {
    if( my_markers_used[i] || my_markers_used[i+1] )
    {
      const bit_field& b1 = my_bits_changed[i];
      const bit_field& b2 = my_bits_changed[i + 1];

      recomb_prob_ratio += my_recomb_like->log_recombination_ratio(i, b1, b2);
    }
  }

  return recomb_prob_ratio;
}
/*
inline double mcmc_simulator::calculate_inheritance_around_marker(int marker)
{
  double probability = 1.0;

  if(marker != 0)
    probability *= calculate_inheritance_at_interval(marker - 1);

  if(marker != my_total_loci - 1)
    probability *= calculate_inheritance_at_interval(marker);

  return probability;
}

inline double mcmc_simulator::calculate_inheritance_at_interval(int interval)
{
  double probability = 1.0;

  bool marker1_bit;
  bool marker2_bit;

  double theta = my_ped_region.get_region().locus(interval).locus_theta(1);

  for(int j = 0; j < my_pedigree.individual_count(); ++j)
  {
    if(my_pedigree.founder(j)) continue;

    marker1_bit = my_data->father_bit(j,interval);
    marker2_bit = my_data->father_bit(j,interval+1);

    if(marker1_bit != marker2_bit) probability *= theta;
    else                           probability *= (1.0 - theta);

    marker1_bit = my_data->mother_bit(j,interval);
    marker2_bit = my_data->mother_bit(j,interval+1);

    if(marker1_bit != marker2_bit) probability *= theta;
    else                           probability *= (1.0 - theta);
  }

  return probability;
}
*/
