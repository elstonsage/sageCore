#include "freq/Estimator.h"

namespace SAGE {
namespace FREQ {

//============================================
//
//  runAnalysis(...)
//
//============================================
Results
Estimator::runAnalysis(const Configuration & config, const FPED::FilteredMultipedigree & mp, cerrorstream & errors)
{
  // Runtime output:
  std::cout << "Running next analysis:" << std::endl;
  
  // Create object to store results:
  Results results;
  
  results.config   =  config;
  results.analyses .  clear();
  results.analyses .  reserve(results.config.getMarkers().size());

  // Loop across markers for analysis:
  size_t m_count = 1;
  
  for(Configuration::MarkerVector::const_iterator marker_id = results.config.getMarkers().begin(); marker_id != results.config.getMarkers().end(); ++marker_id, ++m_count)
  {
    // Set up marker analysis:
    MarkerAnalysis marker_analysis;
    
    marker_analysis.marker_id = *marker_id;
    
    // Runtime output:
    std::cout << "  marker " << mp.info().marker_info(*marker_id).name() << " (" << m_count << " of " << results.config.getMarkers().size() << ")..." << std::endl;

    // Make sure there's at least 2 alleles in the model:
    if(mp.info().marker_info(*marker_id).gmodel().allele_count() < 2)
    {
      std::cout << "    There is insufficient data present at this marker for analysis. Skpping this marker..." << std::endl;
    }
    else
    {
      // Construct sample on current marker:
      Sample sample(mp, *marker_id);

      if(!sample.getValidSpedInfos().size() && !sample.getUnconnecteds().size())
      {
        // Report problem:
        errors << priority(error)
               << "Unable to construct a valid sample for this analysis."
               << " This may be due to mendelian inconsistencies and/or uninformative phenotypes for"
               << " the pedigree data. Skipping this marker..."
               << std::endl;
      }
      else
      {
        // Calculate allele frequencies:
        if(sample.getMarkerInfo().codominant())
        {
          std::cout << "    Calculating allele frequencies..." << std::endl;
        
          marker_analysis.allele_freqs = calculateAlleleFrequencies(results.config, sample);
        }
        else // We've got to set them to something!
        {
          size_t allele_count = sample.getMarkerInfo().gmodel().allele_count();
        
          marker_analysis.allele_freqs.founder_freqs . resize(allele_count, 0.0);
          marker_analysis.allele_freqs.all_freqs     . resize(allele_count, 0.0);
        }
      
        if(!results.config.skip_mle)
        {
          std::cout << "    Maximizing likelihood..." << std::endl;
        
          estimateLikelihood(results.config, sample, marker_analysis);
        }

        results.analyses.push_back(marker_analysis);

      } // End if-sample-is-valid

    } // End if-there's-at-least-2-alleles

  } // End of marker loop
  
  // Return results of analysis:
  return results;
}

//========================================================
//
//  calculateAlleleFrequencies(...)
//
//========================================================
AlleleFrequencies
Estimator::calculateAlleleFrequencies(const Configuration & config, const Sample & sample)
{
  // Set up local data storage:
  
  size_t              allele_count              = sample.getMarkerInfo().gmodel().allele_count(),
                      num_of_founder_alleles    = 0,
                      num_of_nonfounder_alleles = 0;
  std::vector<size_t> founder_allele_counts       (allele_count, 0),
                      nonfounder_allele_counts    (allele_count, 0);
  AlleleFrequencies   freqs;
     
  // Loop across singletons:
  for(Sample::MemberVector::const_iterator ind_itr = sample.getUnconnecteds().begin(); ind_itr != sample.getUnconnecteds().end(); ++ind_itr)
  {
    // If the phenotype is missing, skip this person:
    if(!sample.phenotypePresent(**ind_itr))
      continue;
      
    if(!sample.getGenotypeInfoVector(**ind_itr).size())
      continue;

    // Grab their (own and only!) genotype:
    const MLOCUS::unphased_genotype genotype = sample.getGenotypeInfoVector(**ind_itr)[0].genotype;
    
    // Increment allele count:
    num_of_founder_alleles += 2;
      
    // Increment the specific allele counts:
    ++founder_allele_counts[genotype.allele1().id()];
    ++founder_allele_counts[genotype.allele2().id()];
  }

  // Loop across speds:                 
  for(Sample::SpedInfoVector::const_iterator sped_itr  = sample.getValidSpedInfos().begin ();
                                             sped_itr != sample.getValidSpedInfos().end   (); ++sped_itr)
  {
    // Loop across inds:
    for(FPED::MemberConstIterator ind_itr = sped_itr->sped->member_begin(); ind_itr != sped_itr->sped->member_end(); ++ind_itr)
    {
      // If the phenotype is missing, skip this person:
      if(!sample.phenotypePresent(*ind_itr))
        continue;
      
      if(!sample.getGenotypeInfoVector(*ind_itr).size())
        continue;

      // Grab their (own and only!) genotype:
      const MLOCUS::unphased_genotype genotype = sample.getGenotypeInfoVector(*ind_itr)[0].genotype;
    
      // Extract allele ids:
      size_t allele1 = sped_itr->allele_remapping.getUnmappedId(genotype.allele1().id()),
             allele2 = sped_itr->allele_remapping.getUnmappedId(genotype.allele2().id());

      if(ind_itr->is_nonfounder()) // If the individual is a non-founder:
      {
        num_of_nonfounder_alleles += 2;
      
        ++nonfounder_allele_counts[allele1];
        ++nonfounder_allele_counts[allele2];
      }
      else // Individual is founder
      {
        num_of_founder_alleles += 2;
      
        ++founder_allele_counts[allele1];
        ++founder_allele_counts[allele2];
      }

    } // End individual loop
    
  } // End subpedigree loop
  
  // Calculate founder allele frequencies:
  freqs.founder_freqs.resize(allele_count, 0.0);
  
  if(num_of_founder_alleles)
    for(size_t i = 0; i < allele_count; ++i)
      freqs.founder_freqs[i] = (double)founder_allele_counts[i] / (double)num_of_founder_alleles;

  // Calculate non-founder allele frequencies:
  freqs.all_freqs.resize(allele_count);
  
  for(size_t i = 0; i < allele_count; ++i)
  {
    if(config.use_founder_weight && num_of_founder_alleles && num_of_nonfounder_alleles)
    {
      freqs.all_freqs[i] = (1 - config.founder_weight) * ((double)nonfounder_allele_counts [i] / (double)num_of_nonfounder_alleles) +
                                config.founder_weight  * ((double)founder_allele_counts    [i] / (double)num_of_founder_alleles);
    }
    else
    {
      freqs.all_freqs[i] = ((double)nonfounder_allele_counts[i] + (double)founder_allele_counts[i]) / 
                           ((double)num_of_founder_alleles      + (double)num_of_nonfounder_alleles);


    }
  }

  // Return the calculated frequencies:
  return freqs;
}

//========================================================
//
//  estimateLikelihood(...)
//
//========================================================
void
Estimator::estimateLikelihood(const Configuration & config, const Sample & sample, MarkerAnalysis & analysis)
{
  // Create ParameterMgr, MaxFunction, LikelihoodCalculator, SequenceCfg, and DebugCfg:
  MAXFUN::ParameterMgr mgr;
  MAXFUN::Function     func    (mgr);
  LikelihoodCalculator lh_calc (sample, mgr, false);
  MAXFUN::SequenceCfg  seq     (MAXFUN::SequenceCfg::DEFAULT_MAXIMIZATION, sample.getMarkerInfo().name());
  MAXFUN::DebugCfg     dbg     (config.maxfun_debug ? MAXFUN::DebugCfg::COMPLETE : MAXFUN::DebugCfg::NO_DEBUG_INFO);

  // Set function name:
  seq.setFunctionName(config.estimate_inbreeding ? "Log likelihood without inbreeding" : "Log likelihood");
  
  // Setup the components:
  setupComponents(config, sample, analysis.allele_freqs, mgr, func, lh_calc);
  
  // Runtime output:
  if(config.estimate_inbreeding)
    std::cout << "      Without inbreeding..." << std::endl;

  // Maximize the function and store the results:
  runMaximize(func, seq, dbg, analysis.maxfun_results);

  // If inbreeding estimation is turned on...
  if(config.estimate_inbreeding)
  {
    // Set the function name:
    seq.setFunctionName("Log likelihood with inbreeding");
    
    // Create ParameterMgr, MaxFunction, and LikelihoodCalculator:
    MAXFUN::ParameterMgr mgr2;
    MAXFUN::Function     func2    (mgr2);
    LikelihoodCalculator lh_calc2 (sample, mgr2, true);
  
    // Setup the components:
    setupComponents(config, sample, analysis.allele_freqs, mgr2, func2, lh_calc2);
  
    // Runtime output:
    std::cout << "      With inbreeding..." << std::endl;
    
    // Maximize and store results:
    runMaximize(func2, seq, dbg, analysis.maxfun_results_inbreeding);
  }
}

//========================================================
//
//  runMaximize(...)
//
//========================================================
void
Estimator::runMaximize(MAXFUN::Function & func, const MAXFUN::SequenceCfg & seq, const MAXFUN::DebugCfg & dbg, MAXFUN::Results & results)
{
  for(size_t i = 1; i <= 5; ++i)
  {
    double coeff = double(i) * 2.0;
    
    for(MAXFUN::ParameterIterator param  = func.getMgr().getParamBegin("_Alleles");
                                  param != func.getMgr().getParamEnd  ("_Alleles"); ++param)
    {
      param->setInitialEstimate(param->getInitialEstimate() * coeff);
    }
    
    results = MAXFUN::Maximizer::maximize(func, seq, dbg);
    
    if(results.getConverged())
      break;
  }
}

//========================================================
//
//  setupComponents(...)
//
//========================================================
void 
Estimator::setupComponents(
      const Configuration        & config,
      const Sample               & sample,
      const AlleleFrequencies    & init_freqs,
            MAXFUN::ParameterMgr & mgr,   
            MAXFUN::Function     & func,  
            LikelihoodCalculator & lhcalc)
{
  // Copy over initial allele frequencies, adjusting any that are 0.0:
  std::vector<double> adj_init_freqs(init_freqs.all_freqs);
  
  if(sample.getMarkerInfo().codominant())
  {
    double adj_sum = 0.0;
  
    for(size_t i = 0; i < adj_init_freqs.size(); ++i)
    {
      if(adj_init_freqs[i] == 0.0)
        adj_init_freqs[i] = 0.001;
      
      adj_sum += adj_init_freqs[i];
    }
  
    if(adj_sum != 1.0)
      for(size_t i = 0; i < adj_init_freqs.size(); ++i)
        adj_init_freqs[i] /= adj_sum;

  }

  // Add parameters:

  size_t allele_count = sample.getMarkerInfo().gmodel().allele_count();

  // Adjusted alleles:
  for(MLOCUS::allele_iterator al = sample.getMarkerInfo().allele_begin();
      al != sample.getMarkerInfo().allele_end(); ++al)
  {
    std::string name  = al->name();

    mgr.addParameter("Alleles",  name, MAXFUN::Parameter::DEPENDENT);
    mgr.getParameter("Alleles", name).setIncludePValue(false);
  }
  
  // Inbreeding-specific options:
  if(lhcalc.getEstimateInbreeding())
  {
    mgr.addParameter("Inbreeding coeff.", "Inbreeding coeff.", MAXFUN::Parameter::INDEPENDENT, 0.0, -1.0, 1.0);
    
    mgr.getParameter("Inbreeding coeff.", "Inbreeding coeff.").setIncludePValue(false);
  }

  // Unadjusted alleles:
  for(MLOCUS::allele_iterator al = sample.getMarkerInfo().allele_begin();
      al != sample.getMarkerInfo().allele_end(); ++al)
  {
    std::string name  = al->name();
    double      freq  = sample.getMarkerInfo().codominant() ? adj_init_freqs[al->id()] : 1.0 / (double)allele_count;

    mgr.addParameter("_Alleles", name, MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL, 2.0 * freq, 0.0);
    mgr.getParameter("_Alleles", name).setIncludeInOutput(false);
  }

  // Add steps to the function:
  func.addStep      (boost::bind(&LikelihoodCalculator::calculateDependents, boost::ref(lhcalc), _1), "Calculate dependents");

  if(lhcalc.getEstimateInbreeding())
  {
    func.addStep      (boost::bind(&LikelihoodCalculator::testInbredGenotypeFrequencies, boost::ref(lhcalc), _1), "Test genotype frequencies > 0");
  }

  func.setEvaluator (boost::bind(&LikelihoodCalculator::calculateLikelihood, boost::ref(lhcalc), _1));
}

} // End namespace FREQ
} // End namespace SAGE
