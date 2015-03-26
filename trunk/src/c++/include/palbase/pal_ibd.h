#ifndef REL_PAIR_IBD_H
#define REL_PAIR_IBD_H

#include "palbase/relative_pairs.h"

namespace SAGE    {
namespace PALBASE {

class pair_analysis_ibd : public IBD
{
  public:

    pair_analysis_ibd( relative_pairs &p );

    virtual ~pair_analysis_ibd();

    virtual void            build();
    virtual bool            built() const;

    virtual bool            has_pedigree();

    virtual size_t          add_marker(const string& name, double dist, gmodel_type mt);

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

    virtual bool set_ibd(size_t i, size_t m, double f0, double f1mp, double f2);
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
    virtual bool get_ibd_state(size_t m, a_marker_ibd_state& i_state) const;

    virtual bool set_ibd_state(const sped_pointer sp, const ibd_state_info& i_info);
    virtual bool get_ibd_state(const sped_pointer sp, ibd_state_info& i_info) const;

  protected:

    relative_pairs*  my_pairs;
};

#include "pal_ibd.ipp"

} // end of namespace PALBASE
} // end of namespace SAGE

#endif
