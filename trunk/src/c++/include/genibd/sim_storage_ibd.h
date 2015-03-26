#ifndef GENIBD_SIM_STORAGE_IBD_H
#define GENIBD_SIM_STORAGE_IBD_H

//==========================================================================
//  File:    sim_storage_ibd.h
//
//  Author:  Qing Sun
//           Yeunjoo Song (song@darwin.cwru.edu)                                 
//
//  History: Qing Sun - Version 0.90
//           yjs Jun 02 - Maternal & paternal bit split for sib pair done.
//           yjs Sep 04 - Old sim_IBD, sim_sib_IBD, and Pairs structures
//                        were combined & updated to make it more consistent
//                        with basic_storage_ibd structure for exact &
//                        single method.
//
//  Notes:   This header defines the conatiners of pair data structure
//           and analysis result for ibd simulation program.
//
//  Copyright (c) 1997 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "genibd/sim_pair.h"

namespace SAGE
{

namespace GENIBD
{

class sim_storage_ibd : public IBD
{
  public:

    sim_storage_ibd(const meiosis_map&         mmap,
                    const mcmc_parameters&     par,
                    vector<sim_relative_pair>& pairs);

    virtual ~sim_storage_ibd();

    virtual void            build();
    virtual bool            built() const;

    virtual bool            has_pedigree();

    virtual size_t          pair_count() const;

    virtual const id_pair   get_pair(size_t i) const;
    virtual id_pair         get_pair(size_t i);

    virtual const id_pair   get_pair(const std::string &ped,
                                     const std::string &i1,
                                     const std::string &i2,
                                     error_t &e) const;

    virtual id_pair         get_pair(const std::string &ped,
                                     const std::string &i1,
                                     const std::string &i2,
                                     error_t &e);

    virtual bool            get_pair(size_t i, std::string &ped,
                                     std::string &i1,
                                     std::string &i2) const;
            
    virtual bool            use_pair(size_t i) const;
    virtual bool            use_pair(const mem_pointer i1, const mem_pointer i2) const;

    virtual bool            valid_pair(size_t i) const;
    virtual bool            invalidate_pair(size_t i) const;

    virtual size_t          add_pair(mem_pointer i1, mem_pointer i2, pair_type pt);
    virtual size_t          pair_index(const mem_pointer i1, const mem_pointer i2) const;

    virtual bool set_ibd(size_t i, size_t m, double f0, double f2);
    virtual bool set_ibd(size_t i, const std::vector<double> &f0,
                                   const std::vector<double> &f2);

    virtual bool set_ibd(size_t i, size_t m, double f0, double f1, double f2);
    virtual bool set_ibd(size_t i, const std::vector<double> &f0,
                                   const std::vector<double> &f1mp,
                                   const std::vector<double> &f2);

    virtual bool get_ibd(size_t i, size_t m, double &f0, double &f2) const;
    virtual bool get_ibd(size_t i, std::vector<double> &f0,
                                   std::vector<double> &f2) const;

    virtual bool get_ibd(size_t i, size_t m, double &f0, double &f1mp, double &f2) const;
    virtual bool get_ibd(size_t i, std::vector<double> &f0,
                                   std::vector<double> &f1mp,
                                   std::vector<double> &f2) const;

    virtual const sped_pointer get_subped(size_t i) const;
    virtual sped_pointer       get_subped(size_t i);

    virtual bool set_ibd_state(size_t m, const a_marker_ibd_state& i_state);
    virtual bool get_ibd_state(size_t m,       a_marker_ibd_state& i_state) const;

    virtual bool set_ibd_state(const sped_pointer sp, const ibd_state_info& i_info);
    virtual bool get_ibd_state(const sped_pointer sp, ibd_state_info& i_info) const;

  protected:

    friend class sim_storage_ibd_comparison;

    bool test_pedigree(const std::string& = string()) const;

    const meiosis_map&          my_pedigree;

    const mcmc_parameters&      my_params;
    vector<sim_relative_pair>&  my_pair_ibds;

    bool                        my_built;
};

// ============================
//  sim_storage_ibd_comparison
// ============================

class sim_storage_ibd_comparison
{
  public:

    bool operator()(const sim_relative_pair& pair1,
                    const sim_relative_pair& pair2) const
    {
      if(pair1.get_first_ind()->name()  <  pair2.get_first_ind()->name())    return true;

      if(pair1.get_first_ind()->name()  == pair2.get_first_ind()->name()  &&
         pair1.get_second_ind()->name() <  pair2.get_second_ind()->name())   return true;

      return false;
    }
};

#include "genibd/sim_storage_ibd.ipp"

} // end of namespace GENIBD

} // end of namespace SAGE

#endif
