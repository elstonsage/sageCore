#ifndef PEDIGREE_ANALYSIS_SAMPLE_H
#include "mlod/PedigreeAnalysisSample.ipp"
#endif

namespace SAGE {
namespace MLOD {

inline
PedigreeAnalysisSample::~PedigreeAnalysisSample()
{ }

inline
size_t PedigreeAnalysisSample::get_original_member_count      () const
{
  return my_original_rped.member_count();
}

inline
size_t PedigreeAnalysisSample::get_original_subpedigree_count      () const
{
  return my_original_sped_count;
}

inline
size_t PedigreeAnalysisSample::get_original_pedigree_count    () const
{
  return my_original_rped.pedigree_count();
}

inline
size_t PedigreeAnalysisSample::get_member_count      () const
{
  return my_fped.member_count();
}

inline
size_t PedigreeAnalysisSample::get_subpedigree_count () const
{
  return my_filtered_sped_count;
}

inline
size_t PedigreeAnalysisSample::get_pedigree_count    () const
{
  return my_fped.pedigree_count();
}

inline
size_t PedigreeAnalysisSample::get_largest_bit_count () const
{
  return my_largest_bit_count;
}

inline
bool PedigreeAnalysisSample::is_member_included(const RPED::Member& mem)
{
  MemberInclusionType t = get_member_type(mem);
  
  return t == MIT_DIRECTLY_INFORMATIVE ||
         t == MIT_STRUCTURALLY_INFORMATIVE;
}

inline
PedigreeAnalysisSample::MemberInclusionType 
  PedigreeAnalysisSample::get_member_type(const RPED::Member& mem)
{
  return my_member_types[mem.mpindex()];
}

inline
PedigreeAnalysisSample::MemberInclusionType PedigreeAnalysisSample::get_member_type(const FPED::Member& mem)
{
  return my_member_types[mem.info().get_source_member()->mpindex()];
}

inline
const RPED::MultiPedigree& PedigreeAnalysisSample::get_original_multipedigree() const
{
  return my_original_rped;
}

inline
const FPED::Multipedigree& PedigreeAnalysisSample::get_filtered_multipedigree() const
{
  return my_fped;
}
inline
PedigreeAnalysisSample::SpedIterator PedigreeAnalysisSample::get_subpedigree_begin() const
{
  return my_spedigrees.begin();
}

inline
PedigreeAnalysisSample::SpedIterator PedigreeAnalysisSample::get_subpedigree_end()   const
{
  return my_spedigrees.end();
}
    
inline
const FPED::Subpedigree& PedigreeAnalysisSample::get_subpedigree(size_t spid) const
{
  return *my_spedigrees[spid];
}

}
}

