#ifndef __MPOINT_LIKE_H
#include "lvec/mpoint_like.h"
#endif

namespace SAGE
{

// ================
// Inline Functions
// ================

inline void mpoint_likelihood_data::request_resource(ref_type r)
{
  ++ref_count[r];

  if(ref_count[r] == 1)
    my_build = false;

  if(r == MP_COMBINED && !compute[r])
    ++ref_count[0];
}

inline void mpoint_likelihood_data::release_resource(ref_type r)
{
  --ref_count[r];

  if(r == MP_COMBINED && !compute[r])
    --ref_count[0];
}

inline bool mpoint_likelihood_data::set_markers
  (const pedigree_region& pr)
{
  my_valid = false;

  my_markers = pr;

  my_region = pr.get_region();

  if(compute[2]) ref_count[0] += ref_count[2];

  compute[0] = compute[1] = compute[2] = false;

  delete dots;

  dots = NULL;
  
  return true;
}

inline void mpoint_likelihood_data::set_meiosis_map(const meiosis_map& _mm)
{
    my_valid = false;
    my_mm = _mm;

    if(compute[2])
      ref_count[0] += ref_count[2];

    compute[0] = compute[1] = compute[2] = false;
}

inline bool mpoint_likelihood_data::instantiate_dot_formatter()
{
  if(my_verbose)
  {
    if(dots) delete dots;

    dots = new text_dot_formatter(cout);
    
    dots->set_prefix_width(43);
    
    dots->set_trigger_count(my_region.locus_count());
  }

  return my_verbose;
}

inline const meiosis_map& 
    mpoint_likelihood_data::get_meiosis_map() const
{
  return my_mm;
}

inline const mpoint_likelihood_data::lvector&
    mpoint_likelihood_data::single_point_vector(long v)
{
  if(!valid() ||!built()
              || (compute[2] && !ref_count[0])
              || (!compute[2] && ref_count[0] <= ref_count[2])
              || v < 0 || v > lvector_count())
    return bad_lvector;

  if(!compute[0]) compute_single_point(single_point);

  return single_point[v].l;
}

inline const mpoint_likelihood_data::lvector&
    mpoint_likelihood_data::right_sided_vector (long v)
{
  if(!valid() || !built() || !ref_count[1] || v < 0 || v > lvector_count())
    return bad_lvector;

  if(!compute[1]) compute_separate();

  return right[v].l;
}

inline const mpoint_likelihood_data::lvector&
    mpoint_likelihood_data::left_sided_vector  (long v)
{
  if(!valid() || !built() || !ref_count[1] || v < 0 || v > lvector_count())
    return bad_lvector;

  if(!compute[1]) compute_separate();
  
  return left[v].l;
}

inline const mpoint_likelihood_data::lvector&
    mpoint_likelihood_data::multi_point_vector (long v)
{
  if(!built() || !ref_count[2] || v < 0 || v > lvector_count())
    return bad_lvector;

  if(!compute[2]) 
  {
    compute_combined();
  
    ref_count[0] -= ref_count[2];
  }
  
  return multi_point[v].l;
}

inline long mpoint_likelihood_data::lvector_count() const
{
  return my_region.locus_count();
}

inline bool mpoint_likelihood_data::built() const
{
  return my_build;
} 

}
