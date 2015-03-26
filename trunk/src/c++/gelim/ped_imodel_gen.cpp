#include "gelim/ped_imodel_gen.h"

#undef DEBUG

namespace SAGE
{

MLOCUS::inheritance_model
pedigree_imodel_generator::operator()(const subpedigree& sped, size_t marker) const
{
  // Set our bools initially to the expected (no error) state

  last_model_incon  = false;
  last_model_inform = true;

  // Find the original model:

  const FPED::Multipedigree::mpinfo_type& mpinfo = sped.multipedigree()->info();

  // If the marker isn't available, we have to return an empty model

  if( marker >= mpinfo.marker_count() )
    return empty_model(sped);

  // Otherwise, get the original marker

  const MLOCUS::inheritance_model& orig = mpinfo.marker_info(marker);

#ifdef DEBUG
  cout << "original model:" << endl;
  orig.print_info_sparse_matrix();
#endif

  // Create our vector of phenotypes
  
  vector<uint> phids(sped.member_count());

  for(uint ind = 0; ind < sped.member_count(); ++ind)
  {
    // Get the individual's index in the main subpedigree

    const FPED::FilteredMemberInfo& meminfo  = sped.member_index(ind).info();

    bool missing = meminfo.phenotype_missing(marker, orig);
    
    if(missing)
      phids[ind] = orig.get_missing_phenotype_id();
    else
      phids[ind] = meminfo.phenotype(marker);
  }

  // Use the other operator() for our work
  
  return operator()(sped, marker, orig, phids);
}

void
pedigree_imodel_generator::copy_missing_members
    (MLOCUS::inheritance_model&       model,
     const MLOCUS::inheritance_model& orig,
     const subpedigree&               sped,
     const vector<uint>&              pids) const
{
    for(uint ind = 0; ind < sped.member_count(); ++ind)
    {
      bool missing = !model.strict_phenotype(ind+1);

      if(missing)
      {
        model.copy_penetrance_sexed(orig, pids[ind], ind+1, sped.member_index(ind).get_effective_sex());
      }
    }
}

MLOCUS::inheritance_model
pedigree_imodel_generator::operator()
    (const subpedigree&               sped,
     size_t                           marker,
     const MLOCUS::inheritance_model& orig,
     const vector<uint>&              pids) const
{
  // Set our bools initially to the expected (no error) state

  last_model_incon  = false;
  last_model_inform = true;

  // If marker is passed the end of the set of markers in the multipedigree,
  // we have to return an empty model, since we can't verify state.

  if( marker >= sped.multipedigree()->info().marker_count() )
    return empty_model(sped);

#ifdef DEBUG
  cout << "original model:" << endl;
  orig.print_info_sparse_matrix();
#endif

  // First, make the initial model - copy the model

  MLOCUS::inheritance_model model(orig.gmodel());
  model.set_name(orig.name());

  // Now we add each of the individuals with their phenotype. 

  // Note that this could be done more efficiently by doing allele remapping
  // *before* we do all this, but for now, this is fine.

  for(uint ind = 0; ind < sped.member_count(); ++ind)
  {
#ifdef DEBUG
    uint  ped_index = sped.member_index(ind).index();
    uint sped_index = sped.member_index(ind).subindex();
  cout << ind  << ": sped_index = " << sped_index << ", index = " << ped_index
       << ", name = " << sped.member_index(ind).name() << endl;
#endif

    // Turn the ind's subpedigree index into a phenotype string - This string is 
    // a counter, so that individual phenotypes in the MLOCUS::inheritance_model
    // are the same order and index as the individual in the subpedigree + 1.
    // The +1 is necessary due to the missing phenotype being automatically
    // included for all phenotype models.

    string phen = long2str(ind, 6, 0, '0');

    // Is this individual missing?

    bool missing = sped.member_index(ind).info().phenotype_missing(marker, orig);

    model.add_phenotype(phen, !missing);

    if(!missing)
    {
#ifdef DEBUG
  cout << "not missing " << meminfo.phenotype(marker) << endl;
#endif

      model.copy_penetrance_sexed(orig, pids[ind], ind+1, sped.member_index(ind).get_effective_sex());
    }
    else
    {
#ifdef DEBUG
  cout << "missing" << endl;
#endif

      // Don't Do this yet!! It's expensive if there are lots of alleles to be
      // remapped
      // model.copy_penetrance(orig, orig.get_missing_phenotype_id(), ind+1);
    }
  }

  // Check for phenotypes that don't have all the same penetrance value for
  // every genotype.  If they're all uniform, we can't do much
  if(!model.penetrance_informative())
  {
    last_model_inform = false;
    last_model_incon  = false;

    // Check to see if we're remapping.  If we are, we can return the empty
    // model.  If we're not, we have to copy the missing phenotype information
    // to each individual with missing phenotype.
    if(prior_remap || post_remap)
      return empty_model(model, sped);
    else
    {
      copy_missing_members(model,orig,sped,pids);
      
      return model;
    }
  }

  // At this point we know that there is at least one phenotype where
  // every penetrance is not set to the same value.
  
  // We want to know if there's any genotype elimination or remapping even
  // possible.  This demands that there is at least one individual with both
  // a zero and non-zero penetrance.  If this isn't true, we can get out quickly
  if(!model.genotype_informative())
  {
    // Since we can't remap or eliminate, we simply copy our missing people and
    // get out.  Note that these 'missing' people might make the model genotype
    // informative, since their information might have missing genotypes, but
    // we assume, since they're declared missing, that that doesn't matter.  We
    // in fact assume that GE and such has been done.
    copy_missing_members(model,orig,sped,pids);

    last_model_inform = true;
    last_model_incon  = false;
    
    return model;
  }

  // At this point we have genotype_informative people.

  // We can try to do remapping and genotype elimination based upon this, 
  // but first must copy our missing people from the original model
  copy_missing_members(model,orig,sped,pids);

#ifdef DEBUG
  cout << "new model:" << endl;
  model.print_info_sparse_matrix();
#endif

  // Do pre-GE allele remapping

  if(prior_remap)
  {
    do_remap(model);

    if(!model.penetrance_informative())
    {
      last_model_inform = false;
      last_model_incon  = false;

      return empty_model(model, sped);
    }

#ifdef DEBUG
  cout << "new model: after prior_remap.." << endl;
  model.print_info_sparse_matrix();
#endif
  }
  
  // Do genotype elimination
  
  size_t state = 0;

  if(geno_elim)
  {
    state = do_gelim(sped, model, marker);

    if(state)
    {
      last_model_inform = true;
      last_model_incon  = true;

      return empty_model(model, sped);
    }

#ifdef DEBUG
  cout << "new model: after geno_elim.." << endl;
  model.print_info_sparse_matrix();
#endif
  }

  // Do post-GE remapping
  
  if(post_remap)
  {
    do_remap(model);

    // This if should never execute, because genotype elimination shouldn't
    // cause a marker to become uninformative when it was informative
    // before.
    if(!model.penetrance_informative())
    {
      last_model_inform = false;
      last_model_incon  = false;

      return empty_model(model, sped);
    }

#ifdef DEBUG
  cout << "new model: after post_remap.." << endl;
  model.print_info_sparse_matrix();
#endif
  }

  return model;
}

MLOCUS::inheritance_model pedigree_imodel_generator::empty_model(MLOCUS::inheritance_model model, const subpedigree& sped) const
{
  // Clear out all alleles if there's more than 1.
  if(model.allele_count() > 1)
  {
    model.mark_for_remap(model.allele_begin()->name());

    for(MLOCUS::allele_iterator a = ++model.allele_begin();
        a != model.allele_end(); ++a)
    {
      model.mark_for_remap(a->name());
    }
    model.remap();
  }

  // Make all phenotypes the missing phenotype.  This isn't a big deal, as
  // it is only one phenotype '~remap/~remap' (or X/X if there was only
  // one allele to begin with).  Note that this step is done regardless of 
  // what data might have been present.  It may not be necessary if there 
  // is only one allele, but it's a safety to remove any 'no valid genotypes'
  // cases for the members.

  size_t miss = model.get_missing_phenotype_id();

  for(size_t j = 1; j < model.phenotype_count()+1; ++j)
  {
    model.copy_penetrance_sexed(miss, j, sped.member_index(j-1).get_effective_sex(), true);
  }

  return model;
}

MLOCUS::inheritance_model pedigree_imodel_generator::empty_model(const subpedigree& sped) const
{
  last_model_inform = false;
  last_model_incon  = true;

  // Create empty alleles

  MLOCUS::inheritance_model model;

  model.add_allele("~remap", 1.0);

  // Make all phenotypes the missing phenotype.  This isn't a big deal, as
  // it is only one phenotype '~remap/~remap'.

  size_t miss = model.get_missing_phenotype_id();

  for(uint ind = 0; ind < sped.member_count(); ++ind)
  {
    // Turn the ind id into a phenotype string - This string is a counter,
    // so that individual phenotypes in the MLOCUS::inheritance_model are the same
    // order and index as the individual in the subpedigree + 1.  The
    // +1 is necessessary due to the missing phenotype being automatically
    // included for all phenotype models.

    string phen = long2str(ind, 6, 0, '0');

    model.add_phenotype(phen, false);

    model.copy_penetrance_sexed(miss, ind+1, sped.member_index(ind).get_effective_sex(), true);
  }

  return model;
}

size_t pedigree_imodel_generator::do_remap(MLOCUS::inheritance_model& model) const
{
  RemapSetType alleles;

  determine_allele_set(model, alleles);

  int count = 0;

  for(size_t i = 0; i < model.allele_count(); ++i)
  {
    if(alleles[i]) ++count;
  } 

  if(count < 2) return 0;

  int aid = 0;
  for(MLOCUS::allele_iterator a = model.allele_begin();
      a != model.allele_end(); ++a, ++aid)
  {
    if(alleles[aid])
    {
      model.mark_for_remap(a->name());
    
      a = model.allele_begin() + aid;
    }
  }

  model.remap();

  return count;
}

size_t pedigree_imodel_generator::
  do_gelim(const subpedigree& sped, MLOCUS::inheritance_model& model, size_t marker) const
{
  gelim.set_subpedigree(sped);

  return gelim.process(model, marker);
}

void pedigree_imodel_generator::determine_allele_set
    (const MLOCUS::inheritance_model& i, RemapSetType& bits) const
{
  // Initialize all the bits to be remapped.  We'll selectively mark them as we
  // detect they're not remappable
  for(size_t u = 0; u < i.allele_count(); ++u)
    bits[u] = true;

  // Iterate through phenotypes and find their penetrances for detecting alleles to
  // not remap
  MLOCUS::inheritance_model::phenotype_iterator ph = i.phenotype_begin();

  for( ; !bits.empty() && ph != i.phenotype_end(); ++ph)
  {
    // Only use strict phenotypes
    if(!i.strict_phenotype(*ph)) continue;

    // Check all phased genotypes

    MLOCUS::inheritance_model::phased_penetrance_iterator ppn = 
       i.phased_penetrance_begin(*ph);

    for( ; ppn != i.phased_penetrance_end(*ph); ++ppn)
    {
      if(*ppn != 0)
      {
        // determine the alleles
        MLOCUS::phased_genotype pgen = ppn.phased_geno();

        uint a1 = pgen.allele1().id();
        uint a2 = pgen.allele2().id();

        // Set them to null.
        bits[a1] = false;
        bits[a2] = false;
      }
    }

    // Check all unphased genotypes

    MLOCUS::inheritance_model::unphased_penetrance_iterator uppn = 
       i.unphased_penetrance_begin(*ph);

    for( ; uppn != i.unphased_penetrance_end(*ph); ++uppn)
    {
      if(*uppn != 0)
      {
        // determine the alleles
        MLOCUS::unphased_genotype upgen = uppn.unphased_geno();

        uint a1 = upgen.allele1().id();
        uint a2 = upgen.allele2().id();

        // Set them to null.
        bits[a1] = false;
        bits[a2] = false;
      }
    }
  }

#ifdef DEBUG
  cout << "Alleles Present: " << endl;

  for(size_t b = 0; b < i.allele_count(); ++b)
    if(!bits[b]) cout << i.get_allele(b).name() << "  ";

  cout << endl;

  cout << "Alleles not Present: " << endl;

  for(size_t b = 0; b < i.allele_count(); ++b)
    if(bits[b]) cout << i.get_allele(b).name() << "  ";

  cout << endl;

#endif

}

} // end namespace
