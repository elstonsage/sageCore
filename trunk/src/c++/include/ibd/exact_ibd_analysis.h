#ifndef RELTEST_EXACT_IBD_H
#define RELTEST_EXACT_IBD_H

//==========================================================================
//  File:       exact_ibd_analysis.h
//
//  Author:     Qing Sun
//
//  History:    initial coding of the file                                99
//              Updated to new libraries                         yjs Jul. 03
//
//  Notes:      Exact multipoint IBD generator
//
//  Copyright (c) 1999 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "ibd/ibd_analysis.h"

namespace SAGE {

class exact_ibd_analysis
{
  public:

    class exact_ibd_params
    {
      public:
        bool   generate_ibd_intervals;
        size_t ibd_count;
    };
    
    exact_ibd_analysis(cerrorstream& e = sage_cerr,
                       ostream& output = std::cout,
                       bool verbose    = false);

    ~exact_ibd_analysis();

    cerrorstream get_errors()  const;
    cerrorstream set_errors(cerrorstream& s);

    bool         built()       const;
    bool         valid()       const;

    IBD*         ibd_adaptor() const;

    bool         build(long max_markers, long max_bits,
                       bool multi_point  = true,
                       bool single_point = false,
                       bool intervals    = true);

    bool         set_pedigree(const meiosis_map&, const pedigree_region&);

    bool         build_ibds(bool use_intervals);

    size_t       add_pair(mem_pointer i1, mem_pointer i2, pair_type pt);

    bool         compute(const string& title, bool sp = false,
                         bool intervals = true, bool i_state = false);

  private:

    typedef Likelihood_Vector            lvector;
    typedef mpoint_likelihood_data       ldata;
    
    void compute_single_point(const lvector&, long, bool);
    void compute_two_point(region_type&, long, long, long, bool);

    // Member
    //
    exact_ibd_params            my_params;

    lvector                     temp1;
    lvector                     temp2;
    
    ibd_analysis                my_ibd_analysis;

    ldata                       my_ldata;

    meiosis_map                 my_meiosis_map;

    region_type                 my_region;
    pedigree_region             my_ped_region;

    IBD*                        my_ibds;

    text_dot_formatter*         my_dots;
    
    bool                        my_built;
    bool                        my_valid;
    bool                        my_verbose;

    ostream&                    my_output;

    cerrorstream                errors;
};

#include "ibd/exact_ibd_analysis.ipp"

} // end of namespace SAGE

#endif
