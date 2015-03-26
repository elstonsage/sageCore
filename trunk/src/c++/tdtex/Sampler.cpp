#include "error/internal_error.h"
#include "tdtex/Sampler.h"

namespace SAGE  {
namespace TDTEX {

//=============================================
//
//  Constructor
//
//=============================================
Sampler::Sampler(
  const Configuration          & config, 
  const RPED::MultiPedigree    & mp, 
        MPED::SexCode            sex,
        cerrorstream           & err) :

  my_config                (config),
  my_mp                    (mp),
  my_sex                   (sex),
  my_errors                (err),
  my_total_pedigrees       (mp.pedigree_count()),
  my_total_families        (0),
  my_total_pairs           (0),
  my_total_children        (0),
  my_informative_pedigrees (0),
  my_informative_families  (0),
  my_informative_pairs     (0),
  my_informative_children  (0),
  my_error_count           (0),
  my_informative           (false)
{ 
  // Create the correct transmission table:
  if(my_config.get_method() == Configuration::ALLELES)
  {
    my_table = TransmissionTablePtr(new AlleleTransmissionTable(get_gmodel()));
  }
  else // method == GENOTYPES
  {
    my_table = TransmissionTablePtr(new GenotypeTransmissionTable(get_gmodel()));
  }

  // Sample the pedigree data:
  for(RPED::PedigreeConstIterator ped = mp.pedigree_begin(); ped != mp.pedigree_end(); ++ped)
  {
    // Increment family count:
    my_total_families += ped->family_count();

    // Set up pedigree informativity:
    bool ped_informative = false;

    // Loop through all families and check them for informativity.
    for(RPED::FamilyConstIterator fam = ped->family_begin(); fam != ped->family_end(); ++fam)
    {
      ped_informative |= sample_family(*fam);
    }

    // If this pedigree is informative, increment the count:
    if(ped_informative)
      ++my_informative_pedigrees;

    // Logical-or the multipedigree informativity:
    my_informative |= ped_informative;
  }

//my_table->dump();
}

//=============================================
//
//  COPY Constructor
//
//=============================================
Sampler::Sampler(const Sampler & other) :
  my_config                (other.my_config),
  my_mp                    (other.my_mp),
  my_sex                   (other.my_sex),
  my_errors                (other.my_errors),
  my_total_pedigrees       (other.my_total_pedigrees),
  my_total_families        (other.my_total_families),
  my_total_pairs           (other.my_total_pairs),
  my_total_children        (other.my_total_children),
  my_informative_pedigrees (other.my_informative_pedigrees),
  my_informative_families  (other.my_informative_families),
  my_informative_pairs     (other.my_informative_pairs),
  my_informative_children  (other.my_informative_children),
  my_error_count           (other.my_error_count),
  my_informative           (other.my_informative)
{ 
  my_table = TransmissionTablePtr(other.my_table->clone());
}

//=============================================
//
//  generate_error(...)
//
//=============================================
void
Sampler::generate_error(
        size_t         marker,
  const RPED::Member & parent1,
  const RPED::Member & parent2,
  const RPED::Member & child)
{
  ++my_error_count;

  my_errors << priority(warning)
            << "Mendelian inheritance errors in pedigree '"
            << child.pedigree()->name() 
            << "' for parents '"

            << parent1 .name() << "' (" << ScoreCalculator::get_genotype_name(marker, parent1) << ") and '"
            << parent2 .name() << "' (" << ScoreCalculator::get_genotype_name(marker, parent2) << ") with child '"
            << child   .name() << "' (" << ScoreCalculator::get_genotype_name(marker, child)

            << ")." 
            << endl;
}

//=============================================
//
// is_affected
//
//=============================================
bool
Sampler::is_affected(const RPED::Member& ind, size_t trait) const
{
  double value = ind.pedigree()->info().trait(ind.index(), trait);

  return finite(value) && value > 0.0;
}

//=============================================
//
// is_affected_child(...)
//
//=============================================
bool
Sampler::is_affected_child(const RPED::Member& child) const
{
  return is_affected(child, my_config.get_trait());
}

//=============================================
//
//  is_affected_parent(...)
//
//=============================================
bool
Sampler::is_affected_parent(const RPED::Member& parent) const
{
  return is_affected(parent, my_config.get_parent_trait());
}

//=============================================
//
//  filter_parents(...)
//
//=============================================
void
Sampler::filter_parents(
  const RPED::Member     & parent1, 
  const RPED::Member     & parent2,
        TransmissionList & xmit) const
{
  // Set up local variables:
  bool                         parent_trait_set  = my_config.get_parent_trait() != (size_t)-1,
                               sex_selection_set = !MPED::is_sex_unknown(get_sex());
  TransmissionVector::iterator i                 = xmit.get_list().begin (),
                               e                 = xmit.get_list().end   ();

  while(i != e)
  {
    // Assume the transmission is valid to begin with...
    bool valid = true;
    
    // If source is not known:
    if(i->source == NULL)
    {
      // If we're supposed to subset on the basis of a parent's affectedness or sex, we can't
      // process this transmission because the source is unknown!
      if(parent_trait_set || sex_selection_set)
      {
        valid = false;
      }
    }
    // The source IS known:
    else
    {
      // Invalidate if the source does not come from one of the two given parents:
      valid &= (i->source == &parent1 || i->source == &parent2);
    
      // Invalidate if the source (parent) is not affected by the given parental trait, if specified:
      if(parent_trait_set)
        valid &= is_affected_parent(*(i->source));
    
      // Invalidate if the source (parent) is not the correct sex (if get_sex() is set):
      if(sex_selection_set)
        valid &= (i->source->get_effective_sex() == get_sex());
    }
    
    if(valid)
    {
      ++i;
    }
    else
    {
      xmit.get_list().erase(i++);
    }
  }
}

//=============================================
//
//  sample_family(...)
//
//=============================================
bool
Sampler::sample_family(const RPED::Family& fam)
{
  // If there are no offspring, skip this family:
  if(fam.offspring_count() == 0)
  {
    return false;
  }
    
  // Set up local variables:
  std::set<RPED::MemberConstPointer> skip_child;
  TransmissionList                   xmit;
  size_t                             informative_children = 0,
                                     informative_pairs    = 0;

  // Sample of parent/sib-pair quads
  if(my_config.get_max_sib_pairs() != 0)
  {
    // Select first sibling for the pair:
    for(RPED::OffspringConstIterator child1 = fam.offspring_begin(); child1 != fam.offspring_end(); ++child1)
    {
      // If we've processed more than maximum number of sib pairs from this
      // family, then we're done with this family.

      if(my_config.get_max_sib_pairs() != (size_t)-1 && informative_pairs >= my_config.get_max_sib_pairs())
        break;

      // If this child has already been added to the set of members
  
      if(skip_child.count(&*child1))
        continue;

      // If this child is not affected, then add it to the to-be-skipped list.
      // Then, go on to the next offspring.

      if(!is_affected_child(*child1))
      {
        skip_child.insert(&*child1);
        continue;
      }

      // Get a TransmissionList describing the transmission of the given
      // marker from two parents to a child.
      TransmissionList xmit1 = ScoreCalculator::score_child(my_config.get_marker(), *fam.parent1(), *fam.parent2(), *child1);

      // If there are any errors/inconsistencies in the TranmissionList,
      // report them:
      if(xmit1.get_error())
      {
        generate_error(my_config.get_marker(), *fam.parent1(), *fam.parent2(), *child1);
      }

      // If there was an error, or the TransmissionSize is empty, label
      // this child to-be-skipped and move on to the next offspring.
      if(xmit1.get_error() || !xmit1.get_list().size())
      {
        skip_child.insert(&*child1);
        continue;
      }

      // Now pick the second member of the sib pair:

      RPED::OffspringConstIterator child2 = child1;

      for(++child2; child2 != fam.offspring_end(); ++child2)
      {
        // If we've exceeded the maximum number of sib pairs, break out of this loop:
        if(informative_pairs >= my_config.get_max_sib_pairs())
          break;

        // If the first child is set to be skipped, then break out of the
        // 2nd-member search entirely:
        if(skip_child.count(&*child1))
          break;

        // If the second child is set to be skipped, advance to the next member:
        if(skip_child.count(&*child2))
          continue;

        // If the second child is not affected, set the child to be skipped
        // and advance to the next member:
        if(!is_affected_child(*child2))
        {
          skip_child.insert(&*child2);
          continue;
        }

        // Fetch a TransmissionList describing the transmission of the
        // indicated marker to the 2nd child:
        TransmissionList xmit2 = ScoreCalculator::score_child(my_config.get_marker(), *fam.parent1(), *fam.parent2(), *child2);

        // If there were any errors in the transmission, report them!
        if(xmit2.get_error())
          generate_error(my_config.get_marker(), *fam.parent1(), *fam.parent2(), *child2);

        // If there were any errors in transmission, or the size is null,
        // skip this child2:
        if(xmit2.get_error() || !xmit2.get_list().size())
        {
          skip_child.insert(&*child2);
          continue;
        }

        // Compare the TransmissionList's of each sib into a single TransmissionList:
        xmit = transmission_intersection(xmit1, xmit2);

        filter_parents(*fam.parent1(), *fam.parent2(), xmit);

        // Total pairs includes only pairs that consist of informative
        // individuals, but may not be informative as a pair
        ++my_total_pairs;

        if(xmit.get_error())
        {
          skip_child.insert(&*child1);
          skip_child.insert(&*child2);
          continue;
        }

        if(my_table->informative(xmit))
        {
          skip_child.insert(&*child1);
          skip_child.insert(&*child2);

          ++informative_pairs;
          ++my_informative_pairs;

          my_table->count_transmissions(xmit);
        }
      }

    } // End of family's offspring loop

  }

  // Sample over trios:
  // No need to add children to skip_child, since we do only a single iteration

  if(my_config.get_max_children() != 0)
  {
    for(RPED::OffspringConstIterator child = fam.offspring_begin(); child != fam.offspring_end(); ++child)
    {
      if(my_config.get_max_children() != (size_t)-1 && informative_children >= my_config.get_max_children())
        break;

      if(!is_affected_child(*child) || skip_child.count(&*child))
        continue;

      xmit = ScoreCalculator::score_child(my_config.get_marker(), *fam.parent1(), *fam.parent2(), *child);

      filter_parents(*fam.parent1(), *fam.parent2(), xmit);

      if(xmit.get_error())
        generate_error(my_config.get_marker(), *fam.parent1(), *fam.parent2(), *child);

      if(xmit.get_error() || !xmit.get_list().size())
        continue;

      // Total children are those that can be scored, but may not be
      // informative
      ++my_total_children;

      if(my_table->informative(xmit))
      {
        ++informative_children;
        ++my_informative_children;

        my_table->count_transmissions(xmit);
      }
    }
  }

  if(informative_children || informative_pairs)
  {
    ++my_informative_families;

    return true;
  }
  else
  {
    return false;
  }
}

} // End namespace TDTEX
} // End namespace SAGE
