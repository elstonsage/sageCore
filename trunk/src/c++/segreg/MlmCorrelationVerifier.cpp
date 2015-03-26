#include "segreg/MlmCorrelationVerifier.h"

#define rc residual_correlation_sub_model

#define FM rc::fm
#define MS rc::ms
#define MD rc::md
#define FS rc::fs
#define FD rc::fd
#define BB rc::bb
#define SS rc::ss
#define BS rc::bs

namespace SAGE
{
namespace SEGREG
{


MlmCorrelationVerifier::MlmCorrelationVerifier
    (const PedigreeDataSet::SubpedigreeCursor& speds,
     const rc&                                 resids,
     const binary_member_calculator&           bmc)
  : my_pairs_by_type(),
    my_bmc(bmc),
    my_resids(resids)
{
  // Determine which pairs we need (! fixed or != 0)
  bool use_pt[NUM_OF_CORRS];
  
  for(int i = 0; i < NUM_OF_CORRS; ++i)
  {
    use_pt[i] = !resids.is_correlation_fixed((rc::corr) i) ||
                resids.correlation((rc::corr) i) != 0;
  }

  // Build the pair list
  pair_class pc;

  for(PedigreeDataSet::SubpedigreeIterator
          subped = speds.first;
          subped != speds.second; ++subped)
  {
    for(FPED::FamilyConstIterator
            fam  = subped->family_begin();
            fam != subped->family_end(); ++fam)
    {
      const FPED::Member& mother = *fam->get_mother();
      const FPED::Member& father = *fam->get_father();
 
      // Add father-mother pair if appropriate
      if(use_pt[FM] && is_valid_pair(mother, father))
        add_pair(mother, father, FM);

      // For each child of the family
      for(FPED::OffspringConstIterator
              child = fam->offspring_begin();
              child != fam->offspring_end(); ++child)
      {
        bool child_is_female = child->is_female();

        // Call add routine on father-child
        pc = (child_is_female) ? FD : FS;

        if(use_pt[pc] && is_valid_pair(father, *child))
          add_pair(father, *child, pc);

        // Call add routine on mother-child
        pc = (child_is_female) ? MD : MS;

        if(use_pt[pc] && is_valid_pair(mother, *child))
          add_pair(mother, *child, pc);

        // For each later child of the family
        FPED::OffspringConstIterator next_child = child;
        
        for(++next_child ; next_child != fam->offspring_end(); ++next_child)
        {
          bool next_child_is_female = next_child->is_female();

          // Call add routine on child-child
          if(child_is_female)
          {
            pc = (next_child_is_female) ? SS : BS;
          }
          else
          {
            pc = (next_child_is_female) ? BS : BB;
          }

          if(use_pt[pc] && is_valid_pair(*next_child, *child))
            add_pair(*next_child, *child, pc);
        }
      }
    }
  }
}

/// Tests the correlations for bounds.  This is per note 7 of the SEGREG
/// formula document.
bool MlmCorrelationVerifier::current_correlation_estimates_are_valid() const
{
  // For each pair type
  for(size_t i = 0; i < 8; ++i)
  {
    double delta = my_resids.correlation((rc::corr) i);

    if(-1.0 < delta && delta < 1.0) continue;

    // For each pair in the list

    for(pair_list::const_iterator p  = my_pairs_by_type[i].begin();
                                  p != my_pairs_by_type[i].end(); ++p)
    {
      // For each genotype

      for(genotype_index g = index_AA; g != index_INVALID; ++g)
      {
        // Get theta_u(j) and theta_u(k)
        double tuj = my_bmc.get_expected_susc(*(p->first) , g);
        double tuk = my_bmc.get_expected_susc(*(p->second), g);

        // Get (e^theta)s
        double etuj = exp(tuj);
        double etuk = exp(tuk);

        // prod = (1+etuj) * (1+etuk)
        double prod = (1.0 + etuj) * (1.0 + etuk);

        if(delta > 1.0)
        {
          // mjk = 1.0 / max(etuj, etuk)
          double mjk = 1.0 / ((etuj < etuk) ? etuk : etuj);

          // Calculate upper bound on delta
          double ub = prod * mjk;

          // if ub < delta, we have a problem.
          if(ub < delta) return false;
        }
        else
        {
          double lb;
          
          // if tuj + tuk < 0
          if(tuj + tuk < 0)
            lb = -prod;
          else
            lb = -prod / exp(tuj + tuk);

          // if delta < lb, we have a problem
          if(delta < lb) return false;
        }
      }
    }
  }

  return true;
}

void MlmCorrelationVerifier::dump_valid_pairs() const
{
  for(size_t i = 0; i < 8; ++i)
  {
    if(my_pairs_by_type[i].size())
    {
      cout << "Pair Type " << i << endl;

      for(list<pair_type>::const_iterator pairs  = my_pairs_by_type[i].begin();
                                          pairs != my_pairs_by_type[i].end(); ++pairs)
      {
        cout << "    " << pairs->first->name() << "\t" << pairs->second->name() << endl;
      }
    }
  }
}

}}

#undef rc

#undef FM
#undef MS
#undef MD
#undef FS
#undef FD
#undef BB
#undef SS
#undef BS


