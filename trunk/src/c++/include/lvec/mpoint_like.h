#ifndef __MPOINT_LIKE_H
#define __MPOINT_LIKE_H

//
// Multipoint Likelihood Vector Data
//
// Author: Geoff Wedig (wedig@darwin.cwru.edu)
//
// History: 0.01 gcw Initial Interface Development Jun 15 1998
//
// Copyright (c) 1998 R. C. Elston

#include "LSF/LSF.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "rped/genome_description.h"
#include "util/dots.h"
#include "gelim/pedigree_region.h"
#include "lvec/lvector.h"

namespace SAGE
{

class mpoint_likelihood_data
{
  public:

    typedef RPED::genome_description::region_type region;
    typedef Likelihood_Vector               lvector;
            
    enum ref_type { SINGLE_POINT = 0, MP_SEPARATE = 1, MP_COMBINED = 2 };

    mpoint_likelihood_data(SAGE::cerrorstream& = SAGE::sage_cerr,
                           bool verbose = true);

    ~mpoint_likelihood_data();

    // Constructor variables

    bool set_markers(const pedigree_region&);
    void set_meiosis_map(const meiosis_map&);

    const meiosis_map& get_meiosis_map() const;

    // Construction/Destruction

    // Construct the memory space for the vectors.  This requires the region
    // to be set
    bool build(size_t max_markers, size_t max_size);
    bool built() const;
    
    // Test to see if the object is currently valid (pedigree and mm are
    // good and the same, etc)
    bool valid();

    // Error setup

    SAGE::cerrorstream get_errors() const { return err; }
    SAGE::cerrorstream set_errors(SAGE::cerrorstream& s) { return err = s; }

    // Referencing

    void request_resource(ref_type);
    
    void release_resource(ref_type);

    // Sizing and obtaining of the likelihood vectors

    const lvector& single_point_vector(long v);
    const lvector& right_sided_vector (long v);
    const lvector& left_sided_vector  (long v);
    const lvector& multi_point_vector (long v);

    // If valid, returns # of vectors.  If not valid, returns 0
    long lvector_count() const;

    void dump_times(ostream&) const;

  protected:

    typedef clock_t clock_type;

    struct lvector_profile
    {
      lvector_profile() : self_time(0), inherited_time(0), copy_time(0), operation_time(0) { }
      lvector_profile(const lvector& _l) : l(_l), self_time(0), inherited_time(0), copy_time(0), operation_time(0) { }
      lvector_profile(const lvector_profile& p)
      {
        *this = p;
      }

      void operator=(const lvector_profile& p)
      {
        l = p.l;
        self_time = p.self_time;
        inherited_time = p.inherited_time;
        copy_time = p.copy_time;
        operation_time = p.operation_time;
      }

      lvector l;

      clock_type self_time;
      clock_type inherited_time;
      clock_type copy_time;
      clock_type operation_time;
    };

    bool instantiate_dot_formatter();

    // Allocate memory for the types
    bool build_single_point();
    bool build_separate();
    bool build_combined();

    // Each of the later functions call the earlier if needed, storing them as
    // temporaries if unnecessary.
    void compute_single_point(vector<lvector_profile>&);
    void compute_separate    ();
    void compute_combined    ();

    SAGE::cerrorstream err;

    size_t             my_size, my_mcount;

    meiosis_map        my_mm;

    pedigree_region    my_markers;

    region             my_region;
    
    vector<lvector_profile> single_point;
    vector<lvector_profile> right;
    vector<lvector_profile> left;
    vector<lvector_profile> multi_point;

    lvector bad_lvector;
    lvector temp;
    
    text_dot_formatter* dots;

    long ref_count[3];
    bool compute[3];

    bool my_build;
    bool my_valid;
    bool my_verbose;
};

}

#include "lvec/mpoint_like.ipp"

#endif
