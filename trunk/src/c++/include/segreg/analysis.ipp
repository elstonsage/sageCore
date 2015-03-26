#ifndef SEGREG_ANALYSIS_H
#include "segreg/analysis.h"
#endif

namespace SAGE   {
namespace SEGREG {

//=========================
// primary_analysis inlines
//=========================

inline 
primary_analysis::primary_analysis
  (APP::Output_Streams& o)
  : out(o),
    errors(out.errors()),
    my_quality(true),
    my_models(), 
    my_results()
{ }

inline
primary_analysis::~primary_analysis()
{ }

inline
const primary_analysis_results& primary_analysis::get_results() const
{
  return my_results;
}

inline
bool primary_analysis::is_quality() const
{
  return my_quality;
}

inline
void primary_analysis::set_quality(bool b)
{
  my_quality = b;
}

//=================================
// primary_analysis_results inlines
//=================================

inline 
primary_analysis_results::primary_analysis_results()
  : my_valid(false),
    my_ped_data(), 
    my_model(), 
    my_maxfun_results(),
    my_subpedigree_count(0),
    my_unconnected_count(0),
    my_type_probs(),
    my_pfunc_probs(),
    my_mlm_resid_corr()
{ }

//lint -e{1554}
inline 
primary_analysis_results::primary_analysis_results
  (const primary_analysis_results& p)
  : my_valid            (p.my_valid),
    my_ped_data         (p.my_ped_data), 
    my_model            (p.my_model), 
    my_maxfun_results   (p.my_maxfun_results),
    my_subpedigree_count(p.my_subpedigree_count),
    my_unconnected_count(p.my_unconnected_count),
    my_type_probs       (p.my_type_probs),
    my_pfunc_probs      (p.my_pfunc_probs),
    my_mlm_resid_corr   (p.my_mlm_resid_corr)
{ }

inline 
primary_analysis_results::~primary_analysis_results()
{ }

inline 
primary_analysis_results&
    primary_analysis_results::operator=(const primary_analysis_results& p)
{
  my_valid             = p.my_valid;

  //lint -e{1555}
  my_ped_data          = p.my_ped_data; 
  my_model             = p.my_model; 
  my_maxfun_results    = p.my_maxfun_results;
  my_subpedigree_count = p.my_subpedigree_count;
  my_unconnected_count = p.my_unconnected_count;
  my_type_probs        = p.my_type_probs;
  my_pfunc_probs       = p.my_pfunc_probs;
  my_mlm_resid_corr    = p.my_mlm_resid_corr;

  return *this;
}

inline
bool primary_analysis_results::is_valid() const
{
  return my_valid;
}

inline
const model& primary_analysis_results::get_final_model() const
{
  return my_model;
}

inline
const MAXFUN::Results& primary_analysis_results::get_maxfun_results() const
{
  return my_maxfun_results;
}

inline
uint primary_analysis_results::get_subpedigree_count() const
{
  return my_subpedigree_count;
}

inline
uint primary_analysis_results::get_unconnected_count() const
{
  return my_unconnected_count;
}

inline
const pg::post_geno_map& primary_analysis_results::get_type_probs() const
{
  return my_type_probs;
}

inline
const pf::pen_func_map& primary_analysis_results::get_pfunc_probs() const
{
  return my_pfunc_probs;
}

inline
const MlmResidCorrelationCalculator& primary_analysis_results::get_mlm_resid_corr() const
{
  return my_mlm_resid_corr;
}

}
}
