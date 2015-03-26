#ifndef IBD_ANALYSIS_H
#define IBD_ANALYSIS_H

//==========================================================================
//  File:      ibd_analysis.h
//
//  Author:
//
//  History:   Initial implementation.
//             Updated to new libraries.                         yjs Jul. 03 
//
//  Notes:
//
//  Copyright (c) 2003 R.C. Elston
//    All Rights Reserved
//==========================================================================

#include "ibd/ibd.h"

namespace SAGE {

class ibd_analysis
{
  public:

    ibd_analysis();
    
    ~ibd_analysis() { }

    bool compute(const lvector&, const meiosis_map&, bool i_state = false);

    bool valid() const;

    double prob_share(mem_pointer, mem_pointer, long n) const;
    double prob_share(long, long,               long n) const;

    double sib_prob_share(mem_pointer, mem_pointer, long n) const;
    double sib_prob_share(long, long,               long n) const;

    void   get_lvec_probability(vector<double>& lvec_prob) const;
    void   get_pair_ibd_state(mem_pointer p1, mem_pointer p2, vector<size_t>& pair_ibd_state) const;

    const a_marker_ibd_state& get_ibd_state()  const;

  private:

    size_t find_id(mem_pointer) const;

    bool   is_sib(const meiosis_map&, long, long) const;
    bool   is_hsib(const meiosis_map&, long, long) const;
    bool   is_maternal_hsib(const meiosis_map&, long, long) const;
    bool   is_paternal_hsib(const meiosis_map&, long, long) const;
    bool   is_brother_brother(const meiosis_map&, long, long) const;

    double& my_prob_share(long, long, long n);
    double& my_sib_prob_share(long, long, long n);

    // Data objects

    meiosis_map        my_meiosis_map;

    vector<double>     my_ibd_sharing;
    vector<double>     my_sib_ibd_sharing;

    a_marker_ibd_state my_ibd_state;

    bool               my_valid;
};

#include "ibd/ibd_analysis.ipp"

} // end of namespace SAGE

#endif
