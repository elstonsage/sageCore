#ifndef  GENIBD_SIM_PAIR_H
#define  GENIBD_SIM_PAIR_H

//==========================================================================
//  File:    sim_pair.h
//
//  Author:  Qing Sun
//           Yeunjoo Song
//
//  History: Version 0.90
//           Maternal & paternal bit split for sib pair done.   - yjs Jun 02
//           bool sib to size_t my_pair_class to add hsib.      - yjs Oct 02
//            my_pair_class = 0 default
//                            1 sib
//                            2 hsib
//           Took out simulation pair structure from old MP_mcmc.h
//           & updated to new libraries.                        - yjs Sep 04 
//
//  Notes:   This header defines the basic relative pair data structures for
//           ibd simulation method.
//
//  Copyright (c) 1998 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "genibd/definitions.h"

namespace SAGE
{

namespace GENIBD
{

class sim_relative_pair 
{
  public:
    
    sim_relative_pair(size_t m = 0);
    sim_relative_pair(size_t m,
                      fmember_const_pointer m1, fmember_const_pointer m2,
                      pair_type pt);

    ~sim_relative_pair();

    sim_relative_pair& operator=(const sim_relative_pair& r);

    fmember_const_pointer get_first_ind()  const;
    fmember_const_pointer get_second_ind() const;

    bool   is_sib()             const;
    bool   is_hsib()            const;
    bool   is_brother_brother() const;
    bool   is_sister_sister()   const;

    void   set_last_sharing    (size_t m, size_t i);
    void   set_last_sib_sharing(size_t m, size_t i);

    size_t get_last_sharing    (size_t m)  const;
    size_t get_last_sib_sharing(size_t m)  const;

    void   set_values(size_t m, double f0, double f2, size_t steps);
    void   set_values(size_t m, double f0, double f1mp, double f2, size_t steps);

    double get_value    (size_t m, size_t i) const;

    void   increment_value(size_t m, size_t increment = 1, bool x_linked = false);

  protected:

    fmember_const_pointer      my_member_one;
    fmember_const_pointer      my_member_two;

    pair_type                  my_pair_type;

    struct marker_data
    {
      size_t f0;
      size_t f2;
      size_t last_sharing;
      size_t total;

      // Added for maternal & paternal bit split. - yjs Jul. 2002
      //
      size_t f1m;
      size_t f1p;
      size_t last_sib_sharing;

      marker_data()
      {
        f0 = f2 = f1m = f1p = total = 0;
        last_sharing = last_sib_sharing = (size_t) -1;
      }
    };

    vector<marker_data> my_data;
};

#include "genibd/sim_pair.ipp"

} // end of namespace GENIBD

} // end of namespace SAGE

#endif
