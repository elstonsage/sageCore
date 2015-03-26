#include "freq/LikelihoodCalculator.h"

namespace SAGE {
namespace FREQ {

//========================================
//
//  Constructor
//
//========================================
LikelihoodCalculator::LikelihoodCalculator(const Sample & sample, MAXFUN::ParameterMgr & mgr, bool use_inbreeding, const cerrorstream & err) :
  my_errors         (err),
  my_sample         (sample),
  my_use_inbreeding (use_inbreeding),
  my_singleton_calc (mgr, sample, use_inbreeding)
{
  // Set up peelers:
  my_peelers.clear();
  
  for(Sample::SpedInfoVector::const_iterator sped_itr  = my_sample.getValidSpedInfos().begin ();
                                             sped_itr != my_sample.getValidSpedInfos().end   (); ++sped_itr)
  {
    my_peelers[sped_itr->sped] = boost::shared_ptr<Peeler>(new Peeler(mgr, *sped_itr, my_sample, use_inbreeding));
  }                                                                 
}

//========================================
//
//  Copy Constructor
//
//========================================
LikelihoodCalculator::LikelihoodCalculator(const LikelihoodCalculator & other) :
  my_errors         (other.my_errors),
  my_sample         (other.my_sample),
  my_peelers        (other.my_peelers),
  my_use_inbreeding (other.my_use_inbreeding),
  my_singleton_calc (other.my_singleton_calc)
{ }

//========================================
//
//  Destructor
//
//========================================
LikelihoodCalculator::~LikelihoodCalculator()
{ }

//========================================
//
//  calculateDependents
//
//========================================
int 
LikelihoodCalculator::calculateDependents(MAXFUN::ParameterMgr & mgr) const
{
  // Add up all the independent allele freqs, then set the last one to 1.0 minus the sum:
  size_t allele_count = my_sample.getMarkerInfo().gmodel().allele_count();
  double freq_sum     = 0.0;

  for(size_t i = 0; i < allele_count; ++i)
  {
    freq_sum += mgr.getParameter("_Alleles", i).getCurrentEstimate();
  }
                 
  for(size_t i = 0; i < allele_count; ++i)
  {
    double unadj_est = mgr.getParameter("_Alleles", i).getCurrentEstimate();
           
    mgr.getParameter("Alleles", i).setCurrentEstimate(unadj_est / freq_sum);
  }

  // Return success:
  return 0;
}

int
LikelihoodCalculator::testInbredGenotypeFrequencies(MAXFUN::ParameterMgr& mgr) const
{
  if(!getEstimateInbreeding()) return 0;
  
  size_t allele_count = my_sample.getMarkerInfo().gmodel().allele_count();

  double inbreeding = mgr.getParameter("Inbreeding coeff.", "Inbreeding coeff.").getCurrentEstimate();
                 
  for(size_t i = 0; i < allele_count; ++i)
  {
    double a1_freq = mgr.getParameter("Alleles", i).getCurrentEstimate();

    if(GenotypeProbCalc::getInbredHomozygousFrequency(a1_freq, inbreeding) < 0.0) 
    {
      return 1;
    }
    
    for(size_t j = i+1; j < allele_count; ++j)
    {
      double a2_freq = mgr.getParameter("Alleles", i).getCurrentEstimate();
      
      if(GenotypeProbCalc::getInbredHeterozygousFrequency(a1_freq, a2_freq, inbreeding) < 0.0) return 2;
    }
  }

  // Return success:
  return 0;
}

//========================================
//
//  calculateLikelihood
//
//========================================
double 
LikelihoodCalculator::calculateLikelihood(MAXFUN::ParameterMgr & mgr) const
{
  // Calculate likelihood:
  log_double lh(1.0);
  
  // For every subpedigree:
  for(PeelerMap::iterator itr = my_peelers.begin(); itr != my_peelers.end(); ++itr)
    lh *= itr->second->calculateLikelihood();
  
  // For all the singletons:
  lh *= my_singleton_calc.calculateLikelihood();

  return lh.get_log();
}
    
} // End namespace FREQ
} // End namespace SAGE

