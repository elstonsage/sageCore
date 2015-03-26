#include "freq/Peeler.h"

// Uncomment the line below to enable peeler debugging:
// #define DEBUG_PEELER 1

// For debugging purposes ONLY!
int tab = 0;

#ifdef DEBUG_PEELER

  #define IN(msg)         for(int i=0;i<tab;++i) std::cout << "  "; std::cout << msg << std::endl; ++tab;
  #define OUT(msg) --tab; for(int i=0;i<tab;++i) std::cout << "  "; std::cout << msg << std::endl;
  
#else

  #define IN(msg) ;
  #define OUT(msg) ;
  
#endif

// End debugging stuff

namespace SAGE {
namespace FREQ {

//================================================================
//
//  Constructor
//
//================================================================
Peeler::Peeler(
  const MAXFUN::ParameterMgr   & mgr,
  const Sample::SpedInfo       & sped_info,
  const Sample                 & sample, 
        bool                     use_inbreeding,
        cerrorstream           & err)
  : 
  PeelerBase               (*(sped_info.sped)),
  my_errors                (err),
  my_sample                (sample),
  my_mgr                   (mgr),
  my_transm_prob_calc      (),
  my_geno_prob_calc        (sped_info.remapped_gmodel),
  my_allele_remapping      (sped_info.allele_remapping),
  my_remapped_allele_freqs (sped_info.allele_remapping.getRemappedAlleleCount()),
  my_lh                    (0.0),
  my_lh_count              (0)
{ 
  // Set inbreeding use:
  my_geno_prob_calc.setUseInbreeding(use_inbreeding);
  
  // Set up allele remapping:
  my_allele_remapping.setup(sample.getMarkerInfo().gmodel(), sped_info.remapped_gmodel);

  // Set up 'relevant' alleles bitfield:  
  my_relevant_alleles.resize(my_sample.getMarkerInfo().gmodel().allele_count(), false);

  for(FPED::MemberConstIterator ind = my_subpedigree.member_begin(); ind != my_subpedigree.member_end(); ++ind)
  {
    const GenotypeInfoVector & genotypes = my_sample.getGenotypeInfoVector(*ind);
    
    for(GenotypeInfoVector::const_iterator g = genotypes.begin(); g != genotypes.end(); ++g)
    {
      const MLOCUS::unphased_genotype genotype = g->genotype;
      
      std::vector<size_t> alleles(2);

      alleles.at(0) = my_allele_remapping.getUnmappedId(genotype.allele1().id()),
      alleles.at(1) = my_allele_remapping.getUnmappedId(genotype.allele2().id());
             
      for(size_t i = 0; i < alleles.size(); ++i)
      {
        if(alleles[i] == (size_t)-1) // If it's a remap allele!
        {
          for(std::set<size_t>::const_iterator j = my_allele_remapping.getRemapSet().begin(); j != my_allele_remapping.getRemapSet().end(); ++j)
          {
            my_relevant_alleles[*j] = true;
          }
        }
        else // It's not a remapped allele
        {
          my_relevant_alleles[alleles[i]] = true;
        }
      }
    }
  }

  // Setup cache:
  for(FPED::MemberConstIterator ind = my_subpedigree.member_begin(); ind != my_subpedigree.member_end(); ++ind)
    my_cache.get_individual_cache(*ind).setup(my_sample.getMarkerInfo().unphased_genotype_count(), *ind);

  // Populate transmission probability calc:
  my_transm_prob_calc.populateCache(sped_info.remapped_gmodel);

//  std::cout << "peeler constructor " << this << std::endl;
}

//================================================================
//
//  COPY Constructor
//
//================================================================
Peeler::Peeler(const Peeler & other) :
  PeelerBase               (other),
  my_errors                (other.my_errors),
  my_sample                (other.my_sample),
  my_mgr                   (other.my_mgr),
  my_transm_prob_calc      (other.my_transm_prob_calc),
  my_geno_prob_calc        (other.my_geno_prob_calc),
  my_allele_remapping      (other.my_allele_remapping),
  my_remapped_allele_freqs (other.my_remapped_allele_freqs),
  my_lh                    (other.my_lh),
  my_lh_count              (other.my_lh_count),
  my_relevant_alleles      (other.my_relevant_alleles)
{ 
//  std::cout << "peeler copy constructor from " << &other << " to " << this << std::endl;

}

//================================================================
//
//  Destructor
//
//================================================================
Peeler::~Peeler()
{ 
//  std::cout << "deallocating peeler " << this << " " << my_relevant_alleles.size() << std::endl;
}

//================================================================
//
//  calculateLikelihood()
//
//================================================================
log_double
Peeler::calculateLikelihood()
{
  IN("calculateLikelihood")

  // Update the allele remapping frequencies:
  bool changed_alleles = my_remapped_allele_freqs.update(my_mgr, my_allele_remapping);
  
  bool changed_inbreeding = my_geno_prob_calc.has_inbreeding_changed(my_mgr);

  // If we've already calculated the likelihood AT LEAST ONCE and relevant alleles
  // estimates have NOT changed, return the cached likelihood:
  if(my_lh_count && !changed_alleles && !changed_inbreeding)
  {
    return my_lh;
  }

  // Clear the cache:
  FPED::MemberConstIterator member_end_itr = my_subpedigree.member_end();

  for(FPED::MemberConstIterator ind = my_subpedigree.member_begin(); ind != member_end_itr; ++ind)
    my_cache.get_individual_cache(*ind).clearCachedVals();

  // Repopulate the genotype cache:
  my_geno_prob_calc.populateCache(my_mgr, my_remapped_allele_freqs);

  // Calculate the likelihood:

  my_lh = 0.0;

  const FPED::Member       & ind       = *my_subpedigree.member_begin();
  const GenotypeInfoVector & genotypes =  my_sample.getGenotypeInfoVector(ind);
  
  GenotypeInfoVector::const_iterator genotype_end_itr = genotypes.end();

  for(GenotypeInfoVector::const_iterator g = genotypes.begin(); g != genotype_end_itr; ++g)
  {
    // Grab anterior term:
    log_double anterior_term(anterior(ind, g->genotype));
    
    if(anterior_term.get_double() == 0.0)
      continue;
    
    // Grab posterior term:
    log_double posterior_term (1.0);
    
    for(FPED::MateConstIterator mate = my_sample.getMateBegin(ind); mate != my_sample.getMateEnd(ind); ++mate)
    {
      posterior_term *= posterior_with_mate(ind, mate->mate(), g->genotype);
    }
    
    // Multiply it all and add it to sped lh:
    my_lh += anterior_term * g->penetrance * posterior_term;
  }

  OUT("calculateLikelihood = " << my_lh << " ped name = " << get_subpedigree().pedigree()->name())

  ++my_lh_count;

  return my_lh;
}

//================================================================
//
//  internal_anterior_terminal(...)
//
//================================================================
const Peeler::result_type &
Peeler::internal_anterior_terminal(const member_type & ind, const data_type & g, result_type & iatg)
{
  IN("anterior terminal " << ind.name() << " g=" << g);
  
  iatg = my_geno_prob_calc.getProb(g.get_id());

  OUT("anterior terminal = " << iatg);

  return iatg;
}  

//================================================================
//
//  calculateAnteriorParent(...)          
//
//================================================================
log_double
Peeler::calculateAnteriorParent(const member_type & parent, const member_type & spouse_to_exclude, const GenotypeInfo & g)
{

  IN("anterior parent " << parent.name() << " exclude " << spouse_to_exclude.name() << " g=" << g.genotype);
  
  log_double result(0.0);

  if(g.penetrance)
  {
    log_double anterior_term = anterior(parent, g.genotype);
    
    if(anterior_term.get_double())
    {
      log_double post(1.0);
  
      for(FPED::MateConstIterator mate = my_sample.getMateBegin(parent); mate != my_sample.getMateEnd(parent); ++mate)
      {
        if(mate->mate().mpindex() != spouse_to_exclude.mpindex())
        {
          post *= posterior_with_mate(parent, mate->mate(), g.genotype);
        }
      }
      
      if(post.get_double())
      {
        result = anterior_term * g.penetrance * post;
      }
    }
  }

  OUT("anterior parent = " << result);

  return result;
}

//================================================================
//
//  calculateAnteriorChild(...)          
//
//================================================================
log_double
Peeler::calculateAnteriorChild(const member_type & ind, const data_type & g_m, const data_type & g_f)
{
  IN("anterior child " << ind.name() << " gm=" << g_m << " gf=" << g_f);

  log_double term(1.0);

  FPED::SiblingConstIterator sibling_end_itr = ind.sibling_end();
  
  for(FPED::SiblingConstIterator sib_itr = ind.sibling_begin(); sib_itr != sibling_end_itr; ++sib_itr)
  {
    log_double temp(0.0);

    const GenotypeInfoVector & genotypes = my_sample.getGenotypeInfoVector(*sib_itr);
    
    GenotypeInfoVector::const_iterator genotype_end_itr = genotypes.end();

    for(GenotypeInfoVector::const_iterator genotype_itr = genotypes.begin(); genotype_itr != genotype_end_itr; ++genotype_itr)
    {
      // Get penetrance term:
      if(!genotype_itr->penetrance)
        continue;
        
      // Get transmission prob:
      double transm_prob = my_transm_prob_calc.getProb(genotype_itr->genotype.get_id(),
                                                       g_m.get_id(), g_f.get_id());
      
      if(!transm_prob)
        continue;
        
      // Get posterior term:
      log_double post(1.0);
      
      for(FPED::MateConstIterator mate = my_sample.getMateBegin(*sib_itr); mate != my_sample.getMateEnd(*sib_itr); ++mate)
      {
        post *= posterior_with_mate(*sib_itr, mate->mate(), genotype_itr->genotype);
      }
      
      if(!post.get_double())
        continue;
      
      temp += transm_prob * genotype_itr->penetrance * post;
    }
    
    term *= temp;
  }
  
  OUT("anterior child = " << term);

  return term;
}

//================================================================
//
//  internal_anterior(...)          
//
//================================================================
const Peeler::result_type &
Peeler::internal_anterior(const member_type & ind, const data_type & g, result_type & ia)                        
{
  IN("anterior " << ind.name() << " g=" << g);
  
  ia = 0.0;

  const member_type    & p1            = *ind       . parent1           (),
                       & p2            = *ind       . parent2           ();
  const GenotypeInfoVector & p1_genotypes  =  my_sample . getGenotypeInfoVector (p1),
                       & p2_genotypes  =  my_sample . getGenotypeInfoVector (p2);

  GenotypeInfoVector::const_iterator p1_genotype_end_itr = p1_genotypes.end(),
                                 p2_genotype_end_itr = p2_genotypes.end();

  for(GenotypeInfoVector::const_iterator g1_itr = p1_genotypes.begin(); g1_itr != p1_genotype_end_itr; ++g1_itr)
  {
    // Parent1 term:
    log_double p1_anterior = calculateAnteriorParent(p1, p2, *g1_itr);
    
    if(!p1_anterior.get_double())
      continue;
    
    // Parent2 term:
    log_double p2_term(0.0);
    
    for(GenotypeInfoVector::const_iterator g2_itr = p2_genotypes.begin(); g2_itr != p2_genotype_end_itr; ++g2_itr)
    {
      // Transmission probability:
      double transmission_prob = my_transm_prob_calc.getProb(g.get_id(), g1_itr->genotype.get_id(), g2_itr->genotype.get_id());
      
      if(!transmission_prob)
        continue;
        
      // Anterior parent term:
      log_double anterior_parent = calculateAnteriorParent(p2, p1, *g2_itr);
      
      if(!anterior_parent.get_double())
        continue;
        
      // Anterior child term:
      log_double anterior_child = calculateAnteriorChild(ind, g1_itr->genotype, g2_itr->genotype);
      
      if(!anterior_child.get_double())
        continue;

      // Add it all up!
      p2_term += anterior_parent * transmission_prob * anterior_child;
    }
    
    if(!p2_term.get_double())
      continue;

    ia += p1_anterior * p2_term;
  }
  
  OUT("anterior = " << ia);

  return ia;
}

//================================================================
//  internal_posterior_with_mate(...)          
//================================================================
const Peeler::result_type & 
Peeler::internal_posterior_with_mate (const member_type& ind, const member_type& mate, const data_type & g_ind, result_type & ipem) 
{ 
  IN("posterior of ind " << ind.name() << " through mate=" << mate.name() << " gind=" << g_ind);
  
  ipem = 0.0;
  
  const GenotypeInfoVector & mate_genotypes = my_sample.getGenotypeInfoVector(mate);
  
  FPED::OffspringConstIterator offspring_begin_itr = ind.offspring_begin(mate),
                               offspring_end_itr   = ind.offspring_end();
                               
  GenotypeInfoVector::const_iterator mate_genotype_end_itr = mate_genotypes.end();

  // Begin mate genotype loop:
  for(GenotypeInfoVector::const_iterator g_mate_itr = mate_genotypes.begin(); g_mate_itr != mate_genotype_end_itr; ++g_mate_itr)
  {
    // Anterior term:
    log_double anterior_term = anterior(mate, g_mate_itr->genotype);

    if(!anterior_term.get_double())
      continue;

    // Calculate spouse term:
    
    log_double spouse_term(1.0);
    
    for(FPED::MateConstIterator mate_itr = my_sample.getMateBegin(mate); mate_itr != my_sample.getMateEnd(mate); ++mate_itr)
    {
      if(mate_itr->mate().mpindex() != ind.mpindex())
      {
        spouse_term *= posterior_with_mate(mate, mate_itr->mate(), g_mate_itr->genotype);
      }
    }

    if(!spouse_term.get_double())
      continue;

    // Calculate offspring term:
    
    log_double offspring_term(1.0);
    
    // Begin offspring loop:
    for(FPED::OffspringConstIterator offspring_itr = offspring_begin_itr; offspring_itr != offspring_end_itr; ++offspring_itr)
    {
      log_double temp(0.0);
      
      const GenotypeInfoVector & offspring_genotypes = my_sample.getGenotypeInfoVector(*offspring_itr);
      
      GenotypeInfoVector::const_iterator offspring_genotype_end_itr = offspring_genotypes.end();

      // Begin offspring genotype loop:
      for(GenotypeInfoVector::const_iterator g_offspring_itr = offspring_genotypes.begin(); g_offspring_itr != offspring_genotype_end_itr; ++g_offspring_itr)
      {
        double transmission_prob = my_transm_prob_calc.getProb(g_offspring_itr->genotype.get_id(), g_ind.get_id(), g_mate_itr->genotype.get_id());

        if(!transmission_prob)
          continue;
        
        log_double post(1.0);

        // Begin offspring mate loop:
        for(FPED::MateConstIterator offspring_mate_itr  = my_sample.getMateBegin (*offspring_itr); 
                                    offspring_mate_itr != my_sample.getMateEnd   (*offspring_itr); ++offspring_mate_itr)
        {
          post *= posterior_with_mate(*offspring_itr, offspring_mate_itr->mate(), g_offspring_itr->genotype);

        } // End of offspring mate loop
        
        temp += transmission_prob * g_offspring_itr->penetrance * post;

      } // End of offspring genotype loop

      offspring_term *= temp;

    } // End offspring loop

    if(!offspring_term.get_double())
      continue;

    // Add it all up!

    ipem += anterior_term * g_mate_itr->penetrance * spouse_term * offspring_term;

  } // End individual genotype loop

  OUT("posterior with mate = " << ipem);

  return ipem; 
}
    
//================================================================
// NOTE: None of the following functions are ever used.

//================================================================
//  internal_posterior_terminal(...)          
//================================================================
const Peeler::result_type & Peeler::internal_posterior_terminal(const member_type& ind, const data_type & g, result_type & ipt) { return ipt; }
    
//================================================================
//  internal_posterior_with_mate(...)          
//================================================================
const Peeler::result_type & Peeler::internal_posterior_except_mate (const member_type& ind, const member_type& mate, const data_type & g, result_type & ipwm) { return ipwm; }

//================================================================
//  internal_posterior(...)          
//================================================================
const Peeler::result_type & Peeler::internal_posterior(const member_type& ind, const data_type & g, result_type & ip) { return ip; }

} // End namespace FREQ
} // End namespace SAGE

