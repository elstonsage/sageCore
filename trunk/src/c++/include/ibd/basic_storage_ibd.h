#ifndef BASIC_STORAGE_IBD_H
#define BASIC_STORAGE_IBD_H

#include "ibd/ibd.h"

namespace SAGE {

// ================
// basic_storage_ibd
// ================

class basic_storage_ibd : public SAGE::IBD
{
  public:

    basic_storage_ibd(const meiosis_map& mmap);
    basic_storage_ibd(const meiosis_map& mmap, 
                      const region_type& region,
                      bool               use_intervals);

    virtual ~basic_storage_ibd();

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
    virtual bool get_ibd_state(const sped_pointer sp,       ibd_state_info& i_info) const;

  protected:

    friend class basic_storage_ibd_comparison;

    bool test_pedigree(const std::string& = string()) const;

    meiosis_map                    my_pedigree;

    vector<ibd_pair_info>          my_pairs;
    vector<ibd_probability_info>   my_ibd_probs;

    markers_ibd_state              my_ibd_states;

    bool                        my_built;
};

#include "ibd/basic_storage_ibd.ipp"

} // end of namespace SAGE

#endif
