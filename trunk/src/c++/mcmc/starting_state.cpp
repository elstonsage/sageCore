//==========================================================================
//  File:     starting_state.cpp
//
//  Author:   Geoff Wedig
//
//  History:  0.1 Initial Implementation
//            1.0 Updated to new libraries                       yjs Jun. 04
//
//  Notes:
//
//  Copyright (c) 2004 R. C. Elston
//  All Rights Reserved
//==========================================================================

#include "mcmc/starting_state.h"

namespace SAGE
{

namespace MCMC
{

starting_state_generator::starting_state_generator
    (const McmcMeiosisMap&               mmap,
     const pedigree_region&              ped_r,
     const marker_likelihood_calculator& m_like,
     const recombination_calculator&     re_cal,
     MersenneTwister&                    ra_gen)
  : my_ped_region(ped_r),
    my_pedigree(mmap),
    my_marker_like(m_like),
    my_recomb_calculator(re_cal),
    my_random_generator(ra_gen)
{
  // XXX - Add changes to check compatibility? 0008.31
  if( mmap.get_subpedigree().name() != ped_r.get_subpedigree().name() )
  {
    my_valid = false;
    return;
  }

  my_valid = true;

  set_indicator_count();
  set_random_trial_count();
}

starting_state_generator::~starting_state_generator()
{ }

double
starting_state_generator::generate_starting_state(mcmc_data_accessor::indicator_type& state)
{
#if 0
  cout << "starting_state_generator::generate_starting_state()..." << endl;
#endif

  if( !my_valid ) 
    return numeric_limits<double>::quiet_NaN();

  int marker_count = my_ped_region.get_region().locus_count();

  genotype_eliminator gelim;
  gelim.set_subpedigree(my_ped_region.get_subpedigree());

  my_indicators.resize(0); 
  my_indicators.resize(marker_count, marker_indicators(my_indicator_count));

  // First marker

  for( int test = 0; test < my_indicator_count; ++test )
  {
    create_marker_state(my_indicators[0][test], 0, gelim);
  }

  for( int marker = 1; marker < marker_count; ++marker )
  {
    for( int test = 0; test < my_indicator_count; ++test )
    {
      create_marker_state(my_indicators[marker][test], marker, gelim);
      choose_previous_marker(my_indicators[marker][test], marker);
    }
  }

  // pick which last marker to take set values return.
  //
  vector<double> values(my_indicator_count, numeric_limits<double>::quiet_NaN());

  for( int i = 0; i < my_indicator_count; ++i )
  {
    indicator_data prev = my_indicators[marker_count-1][i];

    double prob = prev.log_likelihood;

    if( SAGE::isnan(prob) ) continue;

    values[i] = prob;
  }

  int best = choose_one(values);

  vector<int>& best_vector = my_indicators[marker_count-1][best].previous_markers;

  best_vector.push_back(best);

  for( int marker = 0; marker < marker_count; ++marker )
  {
    state[marker] = my_indicators[marker][best_vector[marker]].my_indicator;

#if 0
  cout << "marker = " << marker << " : ";
  print_bit_field(cout, state[marker]);
  cout << endl;
#endif
  }

#if 0
  cout << "end of starting_state_generator::generate_starting_state()..." << endl;
#endif

  return values[best];
}

void
starting_state_generator::create_marker_state(indicator_data&      id,
                                              int                  marker,
                                              genotype_eliminator& gelim)
{
#if 0
  cout << "starting_state_generator::create_marker_state(marker = "
       << marker << ")..." << endl;
#endif

  id.my_indicator.resize(my_pedigree.get_meiosis_count());

  allele_vector sorted_ped_allele(my_pedigree.get_individual_count());

  const MLOCUS::inheritance_model* org_imodel = my_ped_region.get_region().locus(marker).locus();

  assert(org_imodel != NULL);  // Should never happen

  if(    !my_ped_region.model_consistent(marker) )
  {
    assign_random_bit(id);
    return;
  }

  if( !random_genotype_reduction(sorted_ped_allele, my_ped_region[marker], marker, gelim) )
  {
    id.log_likelihood = numeric_limits<double>::quiet_NaN();
    return;
  }

  assign_bit(sorted_ped_allele, id); 

  id.log_likelihood = my_marker_like.log_likelihood(marker, id.my_indicator);

#if 0
  print_bit_field(cout, id.my_indicator);
  cout << " id.log_likelihood = " << id.log_likelihood << endl;
#endif

  assert(!SAGE::isnan(id.log_likelihood));
  if( SAGE::isnan(id.log_likelihood) )
  {
    cout << marker << endl;
    my_marker_like.dump_graphs(cout);
    exit(1);
  }

#if 0
  cout << "end of starting_state_generator::create_marker_state()..." << endl;
#endif
}

void
starting_state_generator::choose_previous_marker(indicator_data& id, int marker)
{
#if 0
  cout << "starting_state_generator::choose_previous_marker(marker = "
       << marker << ")" << endl;
#endif

  vector<double> values(my_indicator_count, numeric_limits<double>::quiet_NaN());

  for( int i = 0; i < my_indicator_count; ++i )
  {
    indicator_data prev = my_indicators[marker-1][i];

    double prob = prev.log_likelihood;

    if(SAGE::isnan(prob)) continue;

    double recomb = my_recomb_calculator(marker-1, prev.my_indicator, id.my_indicator);

    values[i] = prob + recomb;
  }

  int best = choose_one(values);

  id.log_likelihood += values[best];

  id.previous_markers.reserve(marker-1);
  id.previous_markers = my_indicators[marker-1][best].previous_markers;
  id.previous_markers.push_back(best);

#if 0
  cout << "end of starting_state_generator::choose_previous_marker()..." << endl;
#endif
}

int
starting_state_generator::choose_one(vector<double>& values) const
{
  // Find smallest - 'Center' on that at 0.0 (all now positive).  Sum is sum of all values

  double smallest = numeric_limits<double>::infinity();

  for(size_t i = 0; i < values.size(); ++i)
    if(!SAGE::isnan(values[i]) && values[i] < smallest)
      smallest = values[i];

  // Calculate sum of exp(centered values).  All values are at least 1.

  vector<double> sums(values.size());

  double sum = 0.0;

  for(size_t i = 0; i < values.size(); ++i)
    if(!SAGE::isnan(values[i]))
    {
      sum += exp(values[i] - smallest);
     
      sums[i] = sum;
    }

  double rand_num = my_random_generator.uniform_real(sum);

  return std::lower_bound(sums.begin(), sums.end(), rand_num) - sums.begin();
}

void
starting_state_generator::assign_bit(const allele_vector& alvec, indicator_data& id)
{
#if 0
  cout << "starting_state_generator::assign_bit()..." << endl;
#endif

  size_t     mother_loc, father_loc;
  CommonAlleleSet::AlleleID  ffa, fsa, mfa, msa, cfa, csa;

  vector<bool> flipped(my_pedigree.get_individual_count(), false);

  // for each child, assign inh_vec bits
  for( size_t i = 0; i < my_pedigree.get_individual_count(); ++i)
  {
    // We don't care about founders
    if( my_pedigree.get_subpedigree().member_index(i).is_founder() )
      continue;

    // Get the individual
    const FPED::Member& mem = my_pedigree.get_subpedigree().member_index(i);

    // Get mother and father positions.
    mother_loc = my_pedigree.get_mother(mem)->subindex();
    father_loc = my_pedigree.get_father(mem)->subindex();

#if 0
  cout << "child = " << my_pedigree[i]->name()
       << ", mother = " << my_pedigree[mother_loc]->name()
       << ", father = " << my_pedigree[father_loc]->name() << endl;
#endif

    // Get parental alleles
    mfa = alvec[mother_loc].first;
    msa = alvec[mother_loc].second;
    ffa = alvec[father_loc].first;
    fsa = alvec[father_loc].second;

    // If parents flipped, flip them.
    if( flipped[mother_loc] )
      std::swap(mfa, msa);

    if( flipped[father_loc] )
      std::swap(ffa, fsa);

#if 0
    cout << mfa << '/' << msa << " X " << ffa << '/' << fsa << endl;
#endif

    cfa = alvec[i].first;
    csa = alvec[i].second;

#if 0
  cout << my_pedigree[i]->name() << " ";
#endif
    // Sort phase of child - after this block, cfa comes from mother, csa
    // from father.
    {
      bool phase1 = (mfa == cfa || msa == cfa) && (ffa == csa || fsa == csa);
      bool phase2 = (mfa == csa || msa == csa) && (ffa == cfa || fsa == cfa);

      if((phase2 && !phase1) || (phase2 && my_random_generator.uniform_integer(2) == 1))
      {
        std::swap(cfa, csa);

        flipped[i] = true;
      }
    }

#if 0
  cout << cfa << '/' << csa << " --> ";
#endif

    // Generate Mother's bit
    size_t k = my_pedigree.get_mother_meiosis(mem);

    if(mfa == msa)          // If mother homozygous
    {
      id.my_indicator[k] = my_random_generator.uniform_integer(2);
    }
    else                    // If mother heterozygous, set based upon which allele matches
    {
      id.my_indicator[k] = (mfa != cfa);
    }

    // Generate Father's bit
    if(ffa == fsa)          // If father homozygous
    {
      id.my_indicator[k+1] = my_random_generator.uniform_integer(2);
    }
    else                    // If father heterozygous, set based upon which allele matches
    {
      id.my_indicator[k+1] = (ffa != csa);
    }

#if 0
  cout << id.my_indicator[k] << id.my_indicator[k+1] << endl;
#endif
  }

#if 0
  cout << "end starting_state_generator::assign_bit()..." << endl;
#endif
}

void
starting_state_generator::assign_random_bit(indicator_data& id)
{
#if 0
  cout << "starting_state_generator::assign_random_bit().." << endl;
#endif

  for(size_t i = 0; i < id.my_indicator.size(); ++i )
    id.my_indicator[i] = my_random_generator.uniform_integer(2);

  id.log_likelihood = 0.0;

#if 0
  cout << "end starting_state_generator::assign_random_bit().." << endl;
#endif
}

bool
starting_state_generator::random_genotype_reduction(allele_vector&           sorted_ped_allele,
                                                    const MLOCUS::inheritance_model& model,
                                                    int                      mid,
                                                    genotype_eliminator&     gelim )
{
#if 0
  cout << "Begin genotype elimination for marker " << model.gmodel().name() << endl;
  cout << "  Eliminating..." << endl;
#endif

  bool valid_random_reduction = false;

  for( size_t attempt = 0; attempt < my_trial_count * my_pedigree.get_individual_count(); ++attempt ) 
  {
    MLOCUS::inheritance_model imodel(model);

    // Generate our set of individuals needing elimination
    //
    vector<pair<size_t, vector<size_t> > > overtyped;

    for( size_t i = 0; i < my_pedigree.get_individual_count(); ++i )
    {
#if 0
  cout << my_pedigree[i]->name() << " : ";
#endif
      if( imodel.unphased_penetrance_count(i + 1) > 1 )
      {
        vector<size_t> genos;
        MLOCUS::inheritance_model::unphased_penetrance_iterator upi =
            imodel.unphased_penetrance_begin(i + 1);

        for( ; upi != imodel.unphased_penetrance_end(i + 1); ++upi )
          genos.push_back(upi.geno_id());

        overtyped.push_back(make_pair(i, genos));
#if 0
  cout << "overtyped ";
#endif
      }
#if 0
  cout << endl;
#endif
    }

#if 0
  cout << "initial overtyped.size() = " << overtyped.size() << endl;
  for(size_t s = 0; s < overtyped.size(); ++s)
  {
    cout << overtyped[s].first << ":";
    for(size_t g = 0; g < overtyped[s].second.size(); ++g)
      cout << overtyped[s].second[g] << " ";
    cout << endl;
  }
  cout << endl;
#endif

    while( overtyped.size() )
    {
      size_t pick_ind = my_random_generator.uniform_integer(overtyped.size());

      size_t ind_index                     = overtyped[pick_ind].first;
      const vector<size_t>& overtype_genos = overtyped[pick_ind].second;

      size_t pick_geno = my_random_generator.uniform_integer(overtype_genos.size());

#if 0
  cout << "ind_index = " << ind_index
       << ", count = " << imodel.unphased_penetrance_count(ind_index + 1)
       << ", pick_index = " << pick_geno << endl;
#endif

      for(size_t g_id = 0; g_id < overtype_genos.size(); ++g_id )
      {
        if( g_id == pick_geno )
          continue;

        imodel.remove_unphased_penetrance(ind_index + 1, overtype_genos[g_id], true);
      }

#if 0
  cout << "after remove_unphased_penetrance" << endl;
  imodel.print_info_sparse_matrix();
#endif

      // XXX! We can speed this up by only looking at the families of
      // individual indiv.
      //
      size_t gelim_return = gelim.process(imodel, mid, genotype_eliminator::genotype, false);

#if 0
  cout << "after gelim  return code = " << gelim_return << endl;
  imodel.print_info_sparse_matrix();
#endif

      if( gelim_return != 0 )
        break;

      // Compress the genotypes
      overtyped.resize(0);

      for( size_t j = 0; j < my_pedigree.get_individual_count(); ++j )
      {
        if( imodel.unphased_penetrance_count(j + 1) > 1 )
        {
          vector<size_t> new_genos;
          MLOCUS::inheritance_model::unphased_penetrance_iterator new_upi =
              imodel.unphased_penetrance_begin(j + 1);

          for( ; new_upi != imodel.unphased_penetrance_end(j + 1); ++new_upi )
            new_genos.push_back(new_upi.geno_id());

          overtyped.push_back(make_pair(j, new_genos));
        }
      }

#if 0
  cout << "after remove: overtyped.size() = " << overtyped.size() << endl;
#endif

    }

    bool valid_pattern = true;

    for( size_t i = 0; i < my_pedigree.get_individual_count(); ++i )
    {
      if( imodel.unphased_penetrance_count(i + 1) != 1 )
      {
        valid_pattern = false;
        break;
      }
    }

    if( valid_pattern )
    {
#if 0
  cout << "  Assigning genotypes:" << endl;
#endif

      for( size_t i = 0; i < my_pedigree.get_individual_count(); ++i )
      {
        MLOCUS::inheritance_model::unphased_penetrance_iterator upi =
            imodel.unphased_penetrance_begin(i + 1);

        pair<CommonAlleleSet::AlleleID, CommonAlleleSet::AlleleID> pa;

        pa.first  = upi.unphased_geno().allele1();
        pa.second = upi.unphased_geno().allele2();

        sorted_ped_allele[i] = pa;

#if 0
  cout <<"    id = " << my_pedigree[i]->name()
       << " a1 = " << upi.unphased_geno().allele1().name()
       << " a2 = " << upi.unphased_geno().allele2().name() << endl;
#endif
      }

      valid_random_reduction = true;

      break;
    }
  }

#if 0
  cout << "End genotype elimination, return " << valid_random_reduction << endl;
#endif

  return valid_random_reduction;
}

} // end of namespace MCMC

} // end of namespace SAGE
