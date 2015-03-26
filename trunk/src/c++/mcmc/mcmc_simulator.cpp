//==========================================================================
//  File:       mcmc_simulator.cpp
//
//  Author:     Qing Sun
//
//  History:    Version 0.90
//              Updated to new libraries                       - yjs Sep. 04
//
//  Notes:      This file implements MCMC simulator which walks on the space
//              based on inheritance vetors and data structures defined in
//              MCMC_DataStruct.h file.(See Soble and Lang, Am.,J. Hum
//              Genet. 58:1323-1337,1996)
//
//  Copyright (c) 1998 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "numerics/geometric.h"
#include "lvec/fixed_bit_calculator.h"
#include "mcmc/mcmc_simulator.h"

namespace SAGE
{

namespace MCMC
{

mcmc_simulator::mcmc_simulator
    (const McmcMeiosisMap&   ped,
     const pedigree_region&  ped_region,
     mcmc_parameters*        param,
     cerrorstream&           err)
  : my_pedigree(ped),
    my_ped_region(ped_region),
  
    my_starting_state(NULL), 
    my_random_generator(),
    my_loop(false),
    my_current_marker(0),
    my_transition_steps(0),
    errors(err)
{
  if( param == NULL )
  {
    errors << priority(fatal)
           << " Can't initiate mcmc_simulator." << endl;
    exit(1);
  }

  if(param->get_random_seed())
  {
    my_random_generator.reseed(param->get_random_seed() +
       ped.get_subpedigree().member_index(0).mpindex());
  }
  
  my_parameters = param;

  my_multipoint          = param->is_multipoint();
  my_max_tunnel          = param->get_max_tunnel();
  my_single_marker_ratio = param->get_single_marker_tunneling();

  my_total_loci          = my_ped_region.inheritance_model_count();
  my_data                = new mcmc_data_accessor(my_pedigree, my_total_loci);
  my_total_useful_loci   = get_valid_markers();

  my_last_switch_marker  = my_random_generator.uniform_integer(my_total_loci);
  my_last_pivot_person   = my_random_generator.uniform_integer(my_pedigree.get_individual_count());

  my_log_theta_ratio.resize(my_total_loci-1);
  my_marker_individuals.resize(my_total_loci);

  for(size_t i=0;i < my_total_loci; ++i )
  {
    my_marker_individuals[i].typed_individual   =
    my_marker_individuals[i].untyped_individual = ind_vector();

    if( i != my_total_loci - 1 )
    {
      double theta = ped_region.get_region().locus(i).locus_theta(1);

      // This is BAD, but needed right now.
      if(theta == 0.0) theta = numeric_limits<double>::epsilon();
      if(theta == 1.0) theta = (double) 1.0 - numeric_limits<double>::epsilon();

      my_log_theta_ratio[i] = log(theta) - log1p(-theta);
    }
  }

  my_nuclear_families.resize(my_pedigree.get_individual_count());

  my_markers_used.resize(my_total_loci);
  my_bits_changed.resize(my_total_loci, bit_field(ped.get_meiosis_count(), false));

  my_transition_pattern.resize(ped.get_meiosis_count());

  // build inheritance calculator
  my_recomb_like = new recombination_calculator(my_ped_region, *my_data);

  // build Founder_graph
  my_marker_like = new marker_likelihood_calculator(my_ped_region, my_pedigree, *my_data);

#if 0
  my_data->dump_accessor(cout);
  my_recomb_like->dump_recombination_calculator(cout);
  my_marker_like->dump_graphs(cout);
#endif
}  

mcmc_simulator::~mcmc_simulator()
{ }

//This function is called from main(). 
bool
mcmc_simulator::start(ostream& out)
{
 //add application specific code here.

 return false;
}

bool
mcmc_simulator::initialize_analysis(ostream& info)
{
  if( my_total_useful_loci == 0 )
  {
    errors << SAGE::priority(SAGE::information)
           << " There is no valid marker data available for pedigree "
           << my_pedigree.get_subpedigree().pedigree()->name() << " . Skip the pedigree"
           << endl;

    info << " There is no valid marker data avaible for pedigree "
         << my_pedigree.get_subpedigree().pedigree()->name() << " . Skip the pedigree"
         << endl;

    return false;
  }

  //allocate memories for typed/untyped vectors.
  init_storage();
  update_storage();

  //generate inheritance vector
  my_starting_state = new starting_state_generator(my_pedigree, my_ped_region,
                                                   *my_marker_like, *my_recomb_like,
                                                   get_random());

  if( !my_starting_state->is_valid() )
    return false;

  // Generate the remappings, if used

  bool use_remapping = true;

#if 0
  cout << "my_starting_state is valid!!" << endl;
  ped.dump_map();
  my_pedigree.dump_map(cout);
#endif

//  size_t ind_count = ped.individual_count();

  my_bit_remappers.resize(0, bit_remap_type(0));
  my_bit_remappers.resize(my_total_loci, bit_remap_type(0));

  for( size_t t = 0; t < my_total_loci; ++t )
  {
    if( !my_data->is_valid_locus(t) )
      continue;

    my_bit_remappers[t] = bit_remap_type(my_data->get_indicator(0).size());

    my_bit_remappers[t].clear();

    // Create the bit patterns.
    fixed_bit_calculator fbc(my_pedigree.get_subpedigree(), my_ped_region[t]);

    const fixed_bit_container& fb = fbc.get_fixed_bit_container();

#if 0
    cout << "marker " << t << endl;

    fb.dump(cout);
#endif

    const FPED::Subpedigree& sped = my_pedigree.get_subpedigree();

    for(FPED::MemberConstIterator i = sped.member_begin();
        i != sped.member_end(); ++i )
    {
      if( i->is_founder() )
        continue;

      size_t mpos = my_pedigree.get_mother_meiosis(*i);
      size_t fpos = my_pedigree.get_father_meiosis(*i);

      bool mf = fb.mother_fixed(*i);
      bool ff = fb.father_fixed(*i);

      if( !mf || !use_remapping ) my_bit_remappers[t].remap_bit(mpos, mpos);
      if( !ff || !use_remapping ) my_bit_remappers[t].remap_bit(fpos, fpos);

      if( (mf && ff) || !use_remapping )
        continue;

      for(FPED::SiblingConstIterator j = i->sibling_begin();
          j != i->sibling_end(); ++j)
      {
        if( !mf && fb.mother_synchronized(*i, *j) )
        {
          size_t sib_mpos = my_pedigree.get_mother_meiosis(*j);

          my_bit_remappers[t].remap_bit(mpos, sib_mpos);
        }
        if( !ff && fb.father_synchronized(*i, *j) )
        {
          size_t sib_fpos = my_pedigree.get_father_meiosis(*j);

          my_bit_remappers[t].remap_bit(fpos, sib_fpos);
        }
      }
    }
  }

  return true;
}

int
mcmc_simulator::get_valid_markers()
{
  int useful_loci = my_total_loci;

  for(size_t i = 0 ; i < my_total_loci; ++i )
  {
//    const inheritance_model& imodel = my_ped_region[i];

    if( !my_ped_region.model_consistent(i) )
    {
      --useful_loci;
    }
    else
      my_data->set_valid_locus(i);
  }

  return useful_loci;
}
 
void
mcmc_simulator::init_storage()
{
  // T2 transitions may only be reasonable on families with grandchildren. 
  // Currently, all nuclear families are valid choices.  We may want to test
  // the other option eventually.

  for( size_t i = 0; i < my_pedigree.get_subpedigree().family_count(); ++i )
  {
    const FPED::Family& fam = my_pedigree.get_subpedigree().family_index(i);

    size_t m = fam.parent1()->subindex();
    size_t f = fam.parent2()->subindex();

    if( fam.parent1()->is_male() )
      std::swap(m, f);

    FPED::OffspringConstIterator oi = fam.offspring_begin();

    my_nuclear_families[m].push_back(oi->subindex());
    my_nuclear_families[f].push_back(oi->subindex());
  }

  int untyped, typed;

  for(size_t m = 0; m < my_total_loci; m++ )
  {
    untyped = typed = 0;

    for( size_t j = 0; j < my_pedigree.get_individual_count(); ++j )
    {
#if 0
  cout << "j = " << j << ", p_id = " << my_pedigree[j]->index()
       << ", sp_id = " << my_pedigree[j]->subindex()
       << ", name = " << my_pedigree[j]->name() << " : ";
#endif

      if( my_data->is_geno_miss(j, m) )
      {
        ++untyped;
#if 0
  cout << "u ";
#endif
      }
      else 
      {
        ++typed;
#if 0
  cout << "y ";
#endif
      }                               
#if 0
  cout << endl;
#endif
    }

    //allocate memory
    typed_individual_type& m_ind = my_marker_individuals[m];

    if( typed )
    { 
      m_ind.typed_individual.resize(typed);
    } 

    if( untyped )
    {
      m_ind.untyped_individual.resize(untyped);
    }

    double t = typed;
    double u = untyped * SM_UTP_weight;

    m_ind.typed_probability = t / (t + u);
  } 

  init_local_individuals();

#if 0
  dump_nuc_fams(cout);
  dump_marker_inds(cout);
#endif
}

void
mcmc_simulator::update_storage()
{
  int ti,uti;
  
  for(size_t m = 0; m < my_total_loci; m++ )
  {
    ti = uti = 0;
    for( size_t j = 0; j != my_pedigree.get_individual_count(); ++j )
    {
      if( my_data->is_geno_miss(j, m) )
      {
        my_marker_individuals[m].untyped_individual[uti++] = j;
      }
      else
      {
        my_marker_individuals[m].typed_individual[ti++] = j;
      }
    }
  }

#if 0
  cout << "mcmc_simulator::update_storage()" << endl;
  dump_marker_inds(cout);
#endif
}

void
mcmc_simulator::init_local_individuals()
{
  my_local_individuals.resize(my_pedigree.get_individual_count());

  const FPED::Subpedigree& sped = my_pedigree.get_subpedigree();

  for(FPED::MemberConstIterator i = sped.member_begin();
      i != sped.member_end(); ++i )
  {
    local_individual_type& c = my_local_individuals[i->subindex()];
 
    // self

    c.local.push_back(i->subindex());

    // parents

    if( i->is_nonfounder() )
    {
      c.local.push_back(my_pedigree.get_mother(*i)->subindex());
      c.local.push_back(my_pedigree.get_father(*i)->subindex());
    }

    // mates

    FPED::MateConstIterator mi = i->mate_begin();
    for( ; mi != i->mate_end(); ++mi )
    {
      c.local.push_back(mi->mate().subindex());
    }

    // siblings

    FPED::SiblingConstIterator si = i->sibling_begin();
    for( ; si != i->sibling_end(); ++si )
    {
      c.local.push_back(si->subindex());
    }

    // offspring

    mi = i->mate_begin();
    for( ; mi != i->mate_end(); ++mi )
    {
      FPED::OffspringConstIterator oi = i->offspring_begin(mi);
      for( ; oi != i->offspring_end(); ++oi )
      {
        c.local.push_back(oi->subindex());
      }
    }

#if 0
    // testing output
    
    cout << "Individual: " << my_pedigree[i]->name() << endl;;
    
    if(c.local.size())
    {
      cout << "T0: ";
      for(size_type t = 0; t < c.local.size(); ++t)
        cout << my_pedigree[c.local[t]]->name() << " ";
      cout << endl;
    }
#endif
  }
}

bool
mcmc_simulator::create_starting_state()
{
  //find a legal state to start the simulation.
  my_starting_state->generate_starting_state(my_data->get_indicator());

  for( size_t i = 0; i != my_total_loci; ++i )
  {
#if 0
    double logl = my_marker_like->log_likelihood(i);
  cout << "i = " << i << " : " << logl << endl;
#endif
  }

  // Test the starting state.

  if( !my_marker_like->valid_likelihood() )
  {
    my_marker_like->dump_graphs(cout);
    my_data->dump_accessor(cout);

    exit(0);
  }

  return true;
}

bool
mcmc_simulator::dememorize(dot_formatter& out, int demem_steps)
{
  for( int i = 0; i < demem_steps; ++i )
  {
    start_step();

    generate_transition();

    double rho = apply_transition();

    if( get_random_real() > rho )
    { 
      reject_transition(); 
    }

    out.trigger();

    end_step();
  }

  return true;
}

// Should be moved so can be inlined in derived classes
void
mcmc_simulator::start_step()
{
  // This does maintenance to make sure the mcmc_simulator is ready for
  // a transition

  for( size_t i = 0; i < my_total_loci; ++i )
    if(my_markers_used[i])
      my_bits_changed[i].clear();

  std::fill(my_markers_used.begin(), my_markers_used.end(), 0);

  my_old_last_switch = my_last_switch_marker; // Keep in case we reject.
  my_old_last_pivot  = my_last_pivot_person;

#if 0
  cout << "my_old_last_switch = " << my_old_last_switch << endl;
  cout << "my_old_last_pivot  = " << my_old_last_pivot  << endl;
  cout << "end mcmc_simulator::start_step()..." << endl; 
#endif
}


void
mcmc_simulator::generate_transition(bool do_tunnel)
{
//  cout << "My Multipoint" << my_multipoint << endl;
  
  static const float q = 2.0/3.0;

  bool single_marker = !my_multipoint || my_random_generator.uniform_real() < my_single_marker_ratio;

//  cout << single_marker << endl;
  
  if( do_tunnel && my_max_tunnel > 1 )
  {
    my_transition_steps = geometric(q, get_random()) + 1;

    if( my_transition_steps > my_max_tunnel ) 
      my_transition_steps = my_transition_steps % my_max_tunnel;
  }
  else
    my_transition_steps = 1;

  bool local_pivot = false;

  transition t;

  if( single_marker )
  {
    if( !my_multipoint )  //singlepoint
    {
      t.marker = my_current_marker;
    }
    else //multipoint
    {
      select_a_locus(t, local_pivot);

      my_last_switch_marker = t.marker;   // Store for next transition.
    }
  }

  for( int tunnel = 0; tunnel < my_transition_steps; ++tunnel )
  {
    my_transition_pattern.clear();

    if( !single_marker )
    {
      select_a_locus(t, local_pivot);

      my_last_switch_marker = t.marker;   // Store for next transition.
    }
    else
      local_pivot = get_random_real() < my_parameters->get_local_individual();

    if( local_pivot ) select_local_pivot_person(t);
    else              select_pivot_person(t);

    // local transitioning! - Note certain this is necessary anymore.  Never
    // have been really certain.
    //if(get_random_real() < 0.2)
    //{
    //  if(t.rule < T1 && my_pedigree.founder(t.person)) select_trans_rule(t);
    //  else if(!my_pedigree.parent(t.person))           select_trans_rule(t);
    //  else if(t.rule == ILLEGAL_SWITCH_RULE)           select_trans_rule(t);
    //}
    //else
      //choosing a transition rule(T0, T1 or T2) and pivot person
      select_trans_rule(t);

    my_last_pivot_person = t.person;   // Store for next transition.

    switch(t.rule)
    {
      case T01 :
      case T02 :
      case T03 : T0_switch(t);
                 break;

      case T1  : T1_switch(t);
                 break;

      case T2a :
      case T2b : T2_switch(t);
                 break;

      default  : cout<<"illegal switch rule.\n";
                 exit(1);
    }

    select_extent(t);
    
    // Apply transition to local_bits

    size_t begin = size_t(-1);
    size_t end = size_t(-1);
  
    switch(t.extent)
    {
      case -1 : begin = 0;        end = t.marker + 1;  break;
      case  0 : begin = t.marker; end = t.marker + 1;  break;
      case  1 : begin = t.marker; end = my_total_loci; break;
    }

    for( size_t i = begin; i < end; ++i )
    {
      if( my_data->is_valid_locus(i) )
        my_bits_changed[i] ^= my_bit_remappers[i](my_transition_pattern);

      my_markers_used[i] = 1;
    }
  } //tunneling through
}

double
mcmc_simulator::apply_transition()
{
  double marker_likelihood_ratio = 0.0;

#ifdef MCMC_DEBUG
  double calc_inheritance_ratio = 1.0;

  if(!my_marker_like->valid_likelihood())
  {
    my_marker_like->dump_graphs(cout);

    exit(0);
  }
#endif

  bool valid_transition = false;

  // Rebuild the markers
  for( size_t i = 0; i < my_total_loci; ++i )
  {
    double new_ratio = 0.0;

    if( !my_markers_used[i] )
      continue;

    new_ratio -= my_marker_like->log_likelihood(i);

#ifdef MCMC_DEBUG
    calc_inheritance_ratio /= calculate_inheritance_around_marker(i);
#endif

    apply_bit_pattern_at_marker(i);

    new_ratio += my_marker_like->log_likelihood(i);

    if( SAGE::isnan(new_ratio) )
    {
      apply_bit_pattern_at_marker(i);

      my_bits_changed[i].clear();
      my_markers_used[i] = 0;
    }    
    else
    {
      marker_likelihood_ratio += new_ratio;
      valid_transition = true;
    }

#ifdef MCMC_DEBUG
    // Finish the ratio
    calc_inheritance_ratio *= calculate_inheritance_around_marker(i);
#endif
  }

  if( !valid_transition )
    return 0.0;

  if( !my_multipoint )
  {
    return exp(marker_likelihood_ratio);
  }
  else
  {
    double recombination_probability_ratio = log_recomb_prob_ratio();

#ifdef MCMC_DEBUG
    if(exp(inheritance_likelihood_ratio) < calc_inheritance_ratio * 0.99 ||
       exp(inheritance_likelihood_ratio) > calc_inheritance_ratio * 1.01)
    {
      cout << "Inheritance Ratios: "
           << exp(inheritance_likelihood_ratio) << ' ' 
           << calc_inheritance_ratio << endl;

      print_descent_pattern(switch_marker);
    }
#endif

    return exp(marker_likelihood_ratio + recombination_probability_ratio);
  }
}

void
mcmc_simulator::reject_transition()
{
  my_last_switch_marker = my_old_last_switch;
  my_last_pivot_person  = my_old_last_pivot;

  for( size_t i = 0; i != my_total_loci; ++i )
  {
    if( !my_markers_used[i] )
      continue;

    apply_bit_pattern_at_marker(i);
  }
}

void
mcmc_simulator::end_step()
{}

void
mcmc_simulator::T0_switch(const transition& t)
{
  switch(t.rule)
  {
    case T01 : change_bit(t.person, true);  break;
    case T02 : change_bit(t.person, false); break;
    case T03 : change_bit(t.person, true);
               change_bit(t.person, false);
    default  : ; // Do nothing
  }
}

void
mcmc_simulator::T1_switch(const transition& t)
{
  const FPED::Member& mem = my_pedigree.get_subpedigree().member_index(t.person);

  FPED::MateConstIterator mi = mem.mate_begin();

  if( mi != mem.mate_end() )
  {
    FPED::OffspringConstIterator oi = mem.offspring_begin(mi);

    bool b = mem.is_female();

    for( ; oi != mem.offspring_end(); ++oi )
    {
      change_bit(oi->subindex(), b);
    }
  }
}

void
mcmc_simulator::T2_switch(const transition& t)
{
//  size_t m = my_pedigree.mother_index(t.person);
//  size_t f = my_pedigree.father_index(t.person);
  const FPED::Member& mem = my_pedigree.get_subpedigree().member_index(t.person);

  FPED::SiblingConstIterator si = mem.sibling_begin();

  for( ; si != mem.sibling_end(); ++ si )
  {
    transition trans;

    trans.person = si->subindex();
    trans.marker = t.marker;

    bool l = my_data->mother_bit(si->subindex(), t.marker);
    bool r = my_data->father_bit(si->subindex(), t.marker);

    if( (t.rule == T2a && l != r) || (t.rule == T2b && l == r) )
    {
      trans.rule = T03;

      T0_switch(trans);
    }

    trans.rule = T1;

    T1_switch(trans);
  }
}

void
mcmc_simulator::print_descent_pattern(size_t switch_marker) const
{
  cout << "PATTERN:" << endl;

  for(size_t q = 0; q != my_pedigree.get_individual_count(); ++q )
  {
    const FPED::Member& mem = my_pedigree.get_subpedigree().member_index(q);

    if(mem.is_founder()) continue;

    cout << mem.name() << " : ";

    if( switch_marker != 0 )
    {
      cout << my_data->mother_bit(q,switch_marker-1) << '/'
           << my_data->father_bit(q,switch_marker-1) << ' ';
    }
    else
      cout << "    ";

    cout << my_data->mother_bit(q,switch_marker) << '/'
         << my_data->father_bit(q,switch_marker) << ' ';


    if( switch_marker != my_total_loci-1 )
    {
      cout << my_data->mother_bit(q,switch_marker+1) << '/'
           << my_data->father_bit(q,switch_marker+1) << endl;
    }
    else
      cout << endl;
  }
}

void
mcmc_simulator::dump_nuc_fams(ostream& o) const
{
  o << "Nuclear Families:" << endl;
  for( size_t f = 0; f < my_nuclear_families.size(); ++f )
  {
//    o << f << " - " << my_pedigree[f].name() << " : ";

//    for( size_t fi = 0; fi < my_nuclear_families[f].size(); ++fi )
//      o << my_pedigree[my_nuclear_families[f][fi]].name() << " ";

//    o << endl;
  }
  o << endl;
}

void
mcmc_simulator::dump_marker_inds(ostream& o) const
{
/*
  o << "Marker individuals:" << endl;
  for( size_t m = 0; m < my_marker_individuals.size(); ++m )
  {
    o << "m = " << m << endl
      << "  typed_individual.size   = " << my_marker_individuals[m].typed_individual.size();

    if( my_marker_individuals[m].typed_individual.size() )
      for( size_t t = 0; t < my_marker_individuals[m].typed_individual.size(); ++t )
        o << " " << my_marker_individuals[m].typed_individual[t]
          << ":" << my_pedigree[my_marker_individuals[m].typed_individual[t]].name();

    o << endl;

    o << "  untyped_individual.size = " << my_marker_individuals[m].untyped_individual.size();

    if( my_marker_individuals[m].untyped_individual.size() )
      for( size_t t = 0; t < my_marker_individuals[m].untyped_individual.size(); ++t )
        o << " " << my_marker_individuals[m].untyped_individual[t]
          << ":" << my_pedigree[my_marker_individuals[m].untyped_individual[t]].name();

    o << endl
      << "  typed_probability       = " << my_marker_individuals[m].typed_probability
      << endl;
  }
  o << endl;
*/
}

void mcmc_simulator::dump_storage() const
{
  dump_nuc_fams(cout);
  dump_marker_inds(cout);
}

} // end of namespace MCMC

} // end of namespace SAGE
