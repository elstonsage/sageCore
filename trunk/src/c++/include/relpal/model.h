#ifndef RELPAL_MODEL_H
#define RELPAL_MODEL_H

//=============================================================================
// File:    model.h
//
// Author:  Yeunjoo Song
//
// History: Initial implementation                                  yjs Mar. 07
//
// Notes:   This file contains definition for following data structures.
//            class regression_model
//
// Copyright (c) 2007   R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "relpal/params.h"

namespace SAGE {
namespace RELPAL {

class regression_model
{
  public:

    friend class relpal_analysis;

    regression_model();

    void    clear();
    void    invalidate();

    bool    valid()                    const;

    size_t  get_trait_count()          const;
    size_t  get_ind_parameter_count()  const;
    size_t  get_ped_parameter_count()  const;

    string get_name()                  const;

    const analysis_options&        get_analysis_options() const;
    const pvalue_options&          get_pvalue_options()   const;

    dependent_variable&            get_trait(size_t t);
    independent_variable&          get_ind_parameter(size_t i);
    independent_variable&          get_ped_parameter(size_t i);

    const dependent_variable&      get_trait(size_t t)         const;
    const independent_variable&    get_ind_parameter(size_t i) const;
    const independent_variable&    get_ped_parameter(size_t i) const;

    trait_iterator                 trait_begin();
    trait_iterator                 trait_end();

    parameter_iterator             ind_parameter_begin();
    parameter_iterator             ind_parameter_end();

    parameter_iterator             ped_parameter_begin();
    parameter_iterator             ped_parameter_end();

    trait_const_iterator           trait_begin()          const;
    trait_const_iterator           trait_end()            const;

    parameter_const_iterator       ind_parameter_begin()  const;
    parameter_const_iterator       ind_parameter_end()    const;

    parameter_const_iterator       ped_parameter_begin()  const;
    parameter_const_iterator       ped_parameter_end()    const;

    void      clear_traits();

    tib_value add_trait(const dependent_variable& t);

    void      clear_ind_parameters();
    void      clear_ped_parameters();
    void      clear_markers();

    pib_value add_ind_parameter(const independent_variable& p);
    pib_value add_ped_parameter(const independent_variable& p);

    pib_value add_ind_covariate(const covariate_type& c);
    pib_value add_ped_covariate(const covariate_type& c);

    void      add_intercept(size_t t = 0);
    void      add_random_err_variance(size_t t1 = 0, size_t t2 = 0);
    void      add_common_env_variance(size_t t1 = 0, size_t t2 = 0);
    void      add_polygenic_variance(size_t t1 = 0, size_t t2 = 0);

    void      set_analysis_options(const analysis_options& op);
    void      set_pvalue_options(const pvalue_options& op);

    void      set_name(const string& s);

    bool      operator==(const regression_model& r) const;
    bool      operator!=(const regression_model& r) const;

    void      dump_model(const relative_pairs& relpairs, ostream &out) const;

  protected:

    void      validate();

  private:

    string            my_name;

    trait_vector      my_traits;

    parameter_vector  my_ind_parameters;
    parameter_vector  my_ped_parameters;

    analysis_options  my_analysis_opt;

    pvalue_options    my_pvalue_opt;

    bool              my_valid;
};

#include "relpal/model.ipp"

} // end of namespace RELPAL
} // end of namespace SAGE

#endif
