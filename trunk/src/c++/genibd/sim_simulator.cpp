//==========================================================================
//  File:       sim_simulator.cpp
//
//  Author:     Qing Sun
//
//  History:    Version 0.90
//              Updated to new libraries.                       - yjs Oct 04
//
//  Notes:      This header defines interfaces to IBD simulation object.
//
//  Copyright (c) 1998 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "genibd/sim_simulator.h"

namespace SAGE
{

namespace GENIBD
{

ibd_mcmc_simulator::ibd_mcmc_simulator(const mcmc_meiosis_map& ped,
                                       const pedigree_region&  ped_region,
                                       mcmc_parameters*        param,
                                       cerrorstream&           err)
                  : mcmc_simulator(ped, ped_region, param, err)
{
  my_relative_pairs = NULL;

  my_trees.resize(my_total_loci);
  for( size_t i = 0; i < my_total_loci; ++i )
  {
    // Do not size the pattern here!
    my_trees[i].tree.resize(my_pedigree.get_individual_count() * 2);
  }
}

ibd_mcmc_simulator::~ibd_mcmc_simulator()
{}

bool
ibd_mcmc_simulator::start(ostream& info)
{
  if( !initialize_analysis(info) )
    return false;

  //check pair set.
  if( !my_relative_pairs || my_relative_pairs->empty() ) 
  {
    errors << SAGE::priority(SAGE::information)
           << " No pairs available in IBD_SIM analysis." << endl; 

    info << " No pairs available in IBD_SIM analysis." << endl;

    return false;
  }

  // Initialize output

  cerrorstream out1(cout);
  cerrorstream out2(info);

  out1.prefix("      ");
  out2.prefix("      ");

  cerrormultistream outstream;

  outstream.insert(out1);
  outstream.insert(out2);

  if(my_multipoint)     do_multipoint_analysis(outstream);
  else                  do_singlepoint_analysis (outstream);

  return true;
}

void
ibd_mcmc_simulator::do_singlepoint_analysis(cerrorstream& out_file)
{
  int rej_counter;

  long batch_count = my_parameters->get_batch_count();
  long demem_count = my_parameters->get_dememorization_step();
  long sim_count   = my_parameters->get_simulation_step();

  if( my_parameters->get_use_factor() )
  {
    size_t base_amount = my_pedigree.get_individual_count();

    double log_base = log((double) base_amount);

    batch_count = std::max((size_t) 100, (size_t) (my_parameters->get_batch_factor() * log_base));
    demem_count = (size_t)my_parameters->get_dememorization_factor() * base_amount;
    sim_count   = (size_t)my_parameters->get_simulation_factor() * base_amount;
  }

  out_file << "Simulation will do "
           << batch_count << " batches using "
           << demem_count << " dememorization steps and "
           << sim_count   << " simulation steps for each marker." 
           << endl << endl;

  out_file << "Generating starting state..." << endl;

  out_file.set_raw_mode();

  size_t total_steps = batch_count * (demem_count + sim_count)
                                   * my_ped_region.get_region().locus_count();

  time_dot_formatter timer(out_file);

  timer.set_prefix_width(40);

  timer.set_trigger_count(total_steps + 1);

  timer.set_time_check_count(min( (int)1000, (int)total_steps/1000));

  for( unsigned long i = 0; i < (unsigned long)batch_count; ++i )
  {
    create_starting_state();

    if( !i )
      timer.trigger();  // Initalize everything

    //simulate on one marker at a time
    for( size_t m = 0; m < my_total_useful_loci; ++m ) 
    {
      if( !my_data->is_valid_locus(m) )
        continue;

      my_current_marker = m;

      dememorize(timer, demem_count);

      // simulation phase

      set_init_ibd_sc();

      timer.trigger();

      rej_counter = 0;

      for( size_t s = 1; s != (size_t)sim_count; ++s )
      {
        start_step();

        generate_transition();

        double rho = apply_transition();

        if( get_random_real() > rho )
        {  
          ++rej_counter;
          reject_transition();   
        }
        else
        {
          singlepoint_ibd_sharing(s);
          rej_counter = 0;
        }

        timer.trigger();

        end_step();
      }

      //last update
      if( rej_counter != 0 )
      {
        singlepoint_ibd_sharing(sim_count - 1);
      }
    }
  }
}

void
ibd_mcmc_simulator::do_multipoint_analysis(SAGE::cerrorstream& out_file)
{
  size_t num_of_rejected, rej_counter = 0;

  long batch_count = my_parameters->get_batch_count();
  long demem_count = my_parameters->get_dememorization_step();
  long sim_count   = my_parameters->get_simulation_step();

  if( my_parameters->get_use_factor() )
  {
    size_t base_amount = my_pedigree.get_individual_count() * my_ped_region.get_region().locus_count();

    double log_base = log((double) base_amount);

    batch_count = std::max((size_t) 100, (size_t) (my_parameters->get_batch_factor() * log_base));
    demem_count = (size_t)my_parameters->get_dememorization_factor() * base_amount;
    sim_count   = (size_t)my_parameters->get_simulation_factor() * base_amount;
  }

  out_file << "Simulation will do " << batch_count << " batches using "
           << demem_count << " dememorization steps and "
           << sim_count   << " simulation steps." << endl << endl;

  out_file << "Generating starting state..." << endl;

  out_file.set_raw_mode();

  size_t total_steps = batch_count * (demem_count + sim_count);

  time_dot_formatter timer(out_file);

  timer.set_prefix_width(40);

  timer.set_trigger_count(total_steps + 1);

  timer.set_time_check_count(std::min( (int)1000, (int)total_steps/1000));

  for( unsigned long i = 0; i < (unsigned long)batch_count; ++i )
  {
    create_starting_state();

    if( !i )
      timer.trigger();  // Initalize everything

    //dememorization phase
 
    dememorize(timer, demem_count);

    // simulation phase

    num_of_rejected = set_init_ibd_sc();

    timer.trigger();

    rej_counter = 0;

    for( size_t s = 1; s < (size_t)sim_count; ++s )
    {
      start_step();

      generate_transition();

      double rho = apply_transition();

      if( get_random_real() > rho )
      {
        ++num_of_rejected; 
        ++rej_counter;
        reject_transition();
      }
      else  
      {
         multipoint_ibd_sharing(s);
         rej_counter = 0;
      }

      timer.trigger();

      end_step();
    }

    if( rej_counter!=0 )
    {
      // Make sure we update all the markers
      for( size_t m = 0; m != my_total_loci; ++m )
        ibd_sharing(sim_count - 1, m);
    }
  }
}

//initiate pairs ibd(f0,f2,fc) counts. This function is called only once
// when first time updating pairs allele ibd counts.
int
ibd_mcmc_simulator::set_init_ibd_sc()
{
  int reject_flag = 0;
 
  start_step();

  generate_transition();

  double rho = apply_transition();

  if( get_random_real() > rho )
  {  
    reject_transition(); 
    reject_flag = 1;
  }

  for( size_t i = 0; i < my_total_loci; ++i )
  {
    my_trees[i].pattern.resize(0);
    my_trees[i].last_used = 0;  
    ibd_sharing(1, i);
    my_trees[i].last_used = 0;  
  }

  end_step();

  return reject_flag;
}

void
ibd_mcmc_simulator::ibd_sharing(int current, int mid)
{
  const bit_field& b    = my_data->get_indicator(mid);
  vector<size_t>&  tree = my_trees[mid].tree;

  size_t offset = current - my_trees[mid].last_used;

#if 0
  cout << "ibd_sharing(current = " << current
       << ", mid = " << mid << ", bit_field = ";
  MCMC::print_bit_field(cout, b);
  cout << ", offset = " << offset << endl;
#endif

  if( !offset )
    return;

  my_trees[mid].last_used = current;

  if( !my_trees[mid].pattern.size() || my_trees[mid].pattern != b )
  {
    my_trees[mid].pattern = b;

    size_t count = 0;

    for( size_t i = 0; i != my_pedigree.get_individual_count(); ++i )
    {
      if( my_pedigree.get_subpedigree().member_index(i).is_founder() )
      {
        tree[i*2]   = count++;
        tree[i*2+1] = count++;
      }
      else
      {
        const FPED::Member& mem = my_pedigree.get_subpedigree().member_index(i);
        
        size_t mloc = my_pedigree.get_mother(mem)->subindex();

        if( b[my_pedigree.get_mother_meiosis(mem)] ) tree[i*2] = tree[mloc*2+1];
        else                                         tree[i*2] = tree[mloc*2];

        size_t floc = my_pedigree.get_father(mem)->subindex();

        if( b[my_pedigree.get_father_meiosis(mem)] ) tree[i*2+1] = tree[floc*2+1];
        else                                         tree[i*2+1] = tree[floc*2];
      }
    }
  }

  for( size_t i = 0; i != my_relative_pairs->size(); ++i )
  {
    sim_relative_pair& r = (*my_relative_pairs)[i];

    if( offset != 0 )
      r.increment_value(mid, offset, my_pedigree.is_x_linked());

    size_t loc1 = r.get_first_ind()->subindex();
    size_t loc2 = r.get_second_ind()->subindex();

    int id11 = tree[loc1*2];
    int id12 = tree[loc1*2+1];
    int id21 = tree[loc2*2];
    int id22 = tree[loc2*2+1];

    int p = (id11 == id21) + (id11 == id22) +
            (id12 == id21) + (id12 == id22);

    if( p > 2 ) p = 2;

    if(    my_pedigree.is_x_linked()
        && r.is_brother_brother() )
      --p;

    r.set_last_sharing(mid, p);

    // Added for maternal & paternal bit split. - yjs Jun. 2002
    //
    if( p == 1 )
    {
      if( r.is_sib() )
      {
        if( id11 == id21 )              // mother sharing
        {
          r.set_last_sib_sharing(mid, 0);
        }
        else if( id12 == id22 )         // father sharing
        {
          if(    !my_pedigree.is_x_linked()
              || r.is_sister_sister() )
          {
            r.set_last_sib_sharing(mid, 1);
          }
        }
      }
      else if( r.is_hsib() )
      {
        if( id11 == id21 )              // mother sharing
        {
          r.set_last_sib_sharing(mid, 0);
        }
        else if( id12 == id22 )         // father sharing
        {
          if(    !my_pedigree.is_x_linked()
              || r.is_sister_sister() )
          {
            r.set_last_sib_sharing(mid, 1);
          }
        }
      }
    } // End of addition
  }

  return;
}

} // end of namespace GENIBD

} // end of namespace SAGE
