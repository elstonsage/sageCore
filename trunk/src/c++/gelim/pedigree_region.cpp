#include <string>
#include "gelim/pedigree_region.h"
#include "gelim/ped_imodel_gen.h"

namespace SAGE
{

bool pedigree_region::build
    (const subpedigree& sped,
     const region&      r,
     bool               eliminate)
{
  reset();

  // Initialize primary data
  
  my_subpedigree = &sped;
  
  my_region = r;

  // Make sure of certain invariants
  
  if(!check_build_invariants()) return false;

  create_model_storage();

  // Create a pedigree_imodel_generator to build our markers,
  // and tell it whether to do genotype elimination

  pedigree_imodel_generator generator;

  generator.set_genotype_elimination(eliminate);

  for(size_type i = 0; i < my_markers.size(); ++i)
  {
    size_t m_index = my_subpedigree->multipedigree()->info().marker_find(my_region.locus(i).name());

    my_markers[i] = generator(*my_subpedigree, m_index);
    
    set_marker_status(i, generator);
  }

  return my_is_built = true;
}

bool pedigree_region::build
  (const subpedigree&     sped,
   const pedigree_region& pr,
   const vector<uint>&    pr_ids,
   bool                   eliminate)
{
  reset();

  // Initialize primary data
  
  my_subpedigree = &sped;
  
  my_region = pr.get_region();

  // Make sure of certain invariants
  
  if(!check_build_invariants()) return false;

  // Check the subpedigree versus the pr_ids.  This is a simple check to make
  // sure the sizes are correct.  It doesn't check the ids themselves for validity,
  // but we'll assume developers using this function will be smart and not
  // have values in pr_ids that are > r's subpedigree's member count.

  if(my_subpedigree->member_count() != pr_ids.size()) return false;

  create_model_storage();

  // Create a vector of phenotype ids.  This vector is just the subpedigree ids
  // vector + 1.
  
  vector<uint> phenotype_ids(pr_ids);
  
  // Add one to each element in phenotype_ids.
  
  std::transform(phenotype_ids.begin(), phenotype_ids.end(), phenotype_ids.begin(),
                 std::bind1st(std::plus<uint>(), 1));
  
  // Create a pedigree_imodel_generator to build our markers,
  // and tell it whether to do genotype elimination

  pedigree_imodel_generator generator;

  generator.set_genotype_elimination(eliminate);

  for(size_type i = 0; i < my_markers.size(); ++i)
  {
    size_t m_index = my_subpedigree->multipedigree()->info().marker_find(my_region.locus(i).name());

    my_markers[i] = generator(*my_subpedigree, m_index, pr[m_index], phenotype_ids);

    set_marker_status(i, generator);
  }

  return my_is_built = true;
}

void pedigree_region::create_model_storage()
{
  // Make sure there's space for the markers.

  my_markers            .resize(my_region.locus_count());
  my_model_consistencies.resize(my_region.locus_count(), true);
  my_model_informatives .resize(my_region.locus_count(), true);
}

bool pedigree_region::check_build_invariants()
{
  // Make sure of certain invariants
  
  if(!my_subpedigree->member_count())  return false;

  if(!my_region.valid())               return false;
  if(!my_region.locus_count())         return false;

  return true;
}

void pedigree_region::set_marker_status
  (size_t i,
   const pedigree_imodel_generator& generator)
{
  my_model_consistencies[i] = !generator.inconsistent();
  my_model_informatives [i] = generator.informative() && !generator.inconsistent();

  // If we're not being quiet, we check and report on errors

  if(!is_quiet())
  {
    // XXX Note:  These messages currently don't report on the pedigree
    // section.  This should be changed, since it may not be the whole
    // pedigree which is inconsistent or uninformative

//      cout << i << ' ' << my_region.locus(i).name() << endl; 

    string pedigree_name = my_subpedigree->pedigree()->name();
    string locus_name    = my_region.locus(i).name();

    if(!generator.informative())
    {
      my_errors << SAGE::priority(SAGE::warning)
                << "No marker phenotype data for pedigree '" 
                << pedigree_name << "' at marker '"
                << locus_name    << "'." << endl;
    }

    if(generator.inconsistent())
    {
      my_errors << SAGE::priority(SAGE::error) << "Pedigree '"
                << pedigree_name << "' is inconsistent at marker '"
                << locus_name    << "'.  No information will be "
                << "used at this marker for this pedigree." << endl;
    }
  }
}

}
