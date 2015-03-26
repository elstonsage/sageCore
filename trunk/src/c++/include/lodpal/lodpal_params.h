#ifndef LODPAL_PARAMS_H
#define LODPAL_PARAMS_H

//****************************************************************************
//* File:      lodpal_params.h                                               *
//*                                                                          *
//* Author:    Kevin Jacobs & Yeunjoo Song                                   *
//*                                                                          *
//* History:   Version 0.0 Initial implementation                kbj         *
//*            1.0 one-parameter model added.                    yjs Nov. 00 *
//*            1.1 covariate added.                              yjs Nov. 00 *
//*            1.2 bad_sib_pair storage & func added.            yjs Mar. 01 *
//*            1.3 rel_pair_map storage & func added.            yjs Mar. 01 * 
//*            1.3 multipoint, singlepoint seperated.            yjs Apr. 01 *
//*            1.4 dsp, re-parameterization added.               yjs Apr. 01 *
//*            1.5 evaluate & update_bound for dsp added.        yjs May. 01 *
//*            1.6 diagnostic option added.                      yjs Jun. 01 *
//*            1.7 parameters seperated from arp.h.              yjs Jul. 01 *
//*            2.0 weight added.                                 yjs Jan. 02 *
//*            3.0 x-linkage added.                              yjs May. 02 *
//*            4.0 parent-of-origin added.                       yjs Feb. 03 *
//*            contrast option added.                            yjs Feb. 06 *
//*                                                                          *
//* Notes:     This header file defines ARPTest class.                       *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "lodpal/definitions.h"

namespace SAGE   {
namespace LODPAL {

class parameter_estimate
{
  public:
    parameter_estimate();
    parameter_estimate(double d);

    void            clear();

    void            fix_value();
    void            fix_value(double v);
    void            set_value(double v);
    void            release_value();
    void            set_initial_value(double d);
    void            set_first_derivative(double d);
    void            set_stderr(double d);
    
    double          initial_value()      const;
    double          value()              const;
    double          first_derivative()   const;
    bool            fixed()              const;
    double          get_stderr()             const;
      
  private:
    double          my_initial_value;
    double          my_value;
    double          my_first_derivative;
    double          my_stderr;
    bool            my_fixed;
};

struct trait_parameter
{
  enum trait_type { conaff, condisc, noconunaff, contrast };

  explicit trait_parameter(size_t t = (size_t)-1, trait_type tt = conaff, double cp = 0.0);

  string name(const RelativePairs& pairs) const;

  bool operator==(const trait_parameter& t) const;
  bool operator!=(const trait_parameter& t) const;
  bool operator< (const trait_parameter& t) const;

  size_t         trait;
  double         cutpoint;
  trait_type     pair_select;
  SampleInfo     info;
  bool           valid;
};

struct autosomal_model_type
{
  enum model_type      { one_parameter, two_parameter };
  enum constraint_type { constrained, unconstrained };
  enum poo_fixed_type  { none, maternal, paternal };

  explicit autosomal_model_type(model_type      m    = one_parameter,
                                constraint_type type = constrained,
                                poo_fixed_type  po_f = none,
                                bool poo             = false,
                                double al = -std::numeric_limits<double>::infinity());

  string name()                                       const;
  string poo_name()                                   const;

  bool   operator==(const autosomal_model_type& m)    const;
  bool   operator!=(const autosomal_model_type& m)    const;

  model_type            model;
  constraint_type       constraint;
  poo_fixed_type        fixed;

  bool                  parent_of_origin;
  double                alpha;
};

struct x_linkage_model_type
{
  explicit x_linkage_model_type(bool lambda1 = true,
                                bool lambda2 = true,
                                double al = 2.634);

  string name()                                       const;

  bool   operator==(const x_linkage_model_type& m)    const;
  bool   operator!=(const x_linkage_model_type& m)    const;

  bool                lambda1_equal;
  bool                lambda2_fixed;
  double              alpha;
};

struct covariate_type
{
  enum cov_ad { dsp, mean, minimum, prop };
  enum cov_op { none, sum, diff, prod, both, avg, single, pair };

  covariate_type();

  explicit covariate_type(size_t c, cov_ad ad = mean, cov_op op = sum,
                          double ad_val = std::numeric_limits<double>::quiet_NaN(),
                          double powr = 1.0,
                          double init_d1 = 0.1,
                          double init_d2 = 0.1);

  string name             (const RelativePairs& pairs) const;
  string effect_name      ()                           const;
  string short_effect_name()                           const;

  bool   operator==(const covariate_type& c)           const;
  bool   operator!=(const covariate_type& c)           const;

  size_t              covariate;

  double              power;
  cov_op              operation;
  cov_ad              adjust;
  double              adjust_value;

  double              variance;
  double              variance_y;

  SampleInfo          info;
  SampleInfo          info_y;

  SampleInfo          info_unweighted;
  SampleInfo          info_y_unweighted;
  
  // params for autosomal models.
  parameter_estimate  delta1;
  parameter_estimate  delta2;

  // params for x-linked models.
  parameter_estimate  delta1mm;
  parameter_estimate  delta1mf;
  parameter_estimate  delta1ff;
  parameter_estimate  delta2ff;

  // params for parent-of-origin models.
  parameter_estimate  delta1m;
  parameter_estimate  delta1p;

  bool                valid;
};

struct marker_type
{
  enum inheritance_type { autosomal, x_linked, y_linked };
  
  explicit marker_type(size_t m = (size_t)-1,
                       inheritance_type in = autosomal,
                       double init_b1 = 0.1,
                       double init_b2 = 0.1);

  string name             (const RelativePairs& pairs) const;
  string effect_name      ()                           const;
  string short_effect_name()                           const;

  bool   operator==(const marker_type& m)              const;
  bool   operator!=(const marker_type& m)              const;

  size_t              marker;

  inheritance_type    inheritance;

  SampleInfo          info;

  // params for autosomal models.
  parameter_estimate  beta1;
  parameter_estimate  beta2;

  parameter_estimate  lambda1;
  parameter_estimate  lambda2;

  // params for x-linked models.
  parameter_estimate  beta1mm;
  parameter_estimate  beta1mf;
  parameter_estimate  beta1ff;
  parameter_estimate  beta2ff;

  parameter_estimate  lambda1mm;
  parameter_estimate  lambda1mf;
  parameter_estimate  lambda1ff;
  parameter_estimate  lambda2ff;

  // params for parent-of-origin models.
  parameter_estimate  beta1m;
  parameter_estimate  beta1p;

  parameter_estimate  lambda1m;
  parameter_estimate  lambda1p;

  bool                valid;
};

struct independent_variable
{
  independent_variable();
  
  string name       (const RelativePairs& pairs)           const;
  string effect_name()                                     const;
  
  string covariate_name       (const RelativePairs& pairs) const;
  string covariate_effect_name()                           const;
  string marker_name          (const RelativePairs& pairs) const;
  string marker_effect_name   ()                           const;
  
  bool   operator==(const independent_variable& c)         const;
  bool   operator!=(const independent_variable& c)         const;
  
  void clear();
  
  typedef  vector<marker_type>       marker_vector;
  typedef  vector<covariate_type>    covariate_vector;
  
  marker_vector                      markers;
  covariate_vector                   covariates;
  
  bool                               valid;
};

struct lodpal_weight_type
{
  enum weight_op { none, single, pair };

  explicit lodpal_weight_type(size_t w = (size_t)-1, weight_op = single);

  string name(const RelativePairs& pairs) const;

  bool operator==(const lodpal_weight_type& w) const;
  bool operator!=(const lodpal_weight_type& w) const;

  size_t         weight;

  weight_op      operation;

  SampleInfo     info;
  bool           valid;
};

class lodpal_parameters
{
  public:
    friend class lodpal_pairs;
    friend class ARP_base_analysis;
    friend class ARP_one_analysis;
    friend class ARP_two_analysis;
    friend class ARP_x_one_analysis;
    friend class ARP_x_two_analysis;
    friend class ARP_po_one_analysis;
    friend class ARP_po_two_analysis;

    typedef std::vector<trait_parameter>           trait_vector;
    typedef trait_vector::iterator                 trait_iterator;
    typedef trait_vector::const_iterator           trait_const_iterator;

    typedef independent_variable::marker_vector    marker_vector;
    typedef marker_vector::iterator                marker_iterator;
    typedef marker_vector::const_iterator          marker_const_iterator;

    typedef independent_variable::covariate_vector covariate_vector;
    typedef covariate_vector::iterator             covariate_iterator;
    typedef covariate_vector::const_iterator       covariate_const_iterator;


    lodpal_parameters();

    void    clear();

    bool    valid()                                      const;
    void    invalidate();

    bool    skip_uninformative_pairs()                   const;
    bool    use_pair_cache()                             const;
    bool    turn_off_default()                           const;
    bool    multipoint()                                 const;
    bool    print_lambda()                               const;
    bool    sib_pairs_only()                             const;
    bool    autosomal_marker_exist()                     const;
    bool    x_linked_marker_exist()                      const;
    bool    use_mm_pair()                                const;
    bool    use_mf_pair()                                const;
    bool    use_ff_pair()                                const;
    size_t  used_pair_type()                             const;
    size_t  diagnostic_marker()                          const;
    size_t  max_marker_name_size()                       const;
    
    void    set_skip_uninformative_pairs(bool s);
    void    set_use_pair_cache(bool c);
    void    set_turn_off_default(bool d);
    void    set_multipoint(bool m);
    void    set_print_lambda(bool l);
    void    set_sib_pairs_only(bool s);
    void    set_autosomal_marker_exist(bool x);
    void    set_x_linked_marker_exist(bool x);
    void    set_x_linkage_pair_type(bool mm, bool mf, bool ff);
    void    set_diagnostic_marker(size_t dm);
    void    set_max_marker_name_size(size_t s);

    void    dump_parameters(const RelativePairs &pairs)  const;
    
    size_t  trait_count()                                const;
    size_t  subset_count()                               const;
    size_t  marker_count()                               const;
    size_t  covariate_count()                            const;
    size_t  parameter_count()                            const;

  protected:
    marker_iterator                marker_begin();
    marker_iterator                marker_end();

    covariate_iterator             covariate_begin();
    covariate_iterator             covariate_end();

  public:
    marker_const_iterator          marker_begin()                     const;
    marker_const_iterator          marker_end()                       const;

    covariate_const_iterator       covariate_begin()                  const;
    covariate_const_iterator       covariate_end()                    const;

    const trait_parameter&         trait_parameters(size_t t)         const;
    const trait_parameter&         subset_parameters(size_t t)        const;
    const marker_type&             marker_parameters(size_t m)        const;
    const covariate_type&          covariate_parameters(size_t c)     const;

    const independent_variable&    parameters()                       const;
    const lodpal_weight_type&      weight_parameter()                 const;
    const autosomal_model_type&    autosomal_model()                  const;
    const x_linkage_model_type&    x_linkage_model()                  const;

    trait_parameter&               trait_parameters(size_t t);
    trait_parameter&               subset_parameters(size_t t);
    marker_type&                   marker_parameters(size_t m);
    covariate_type&                covariate_parameters(size_t c);

    independent_variable&          parameters();
    lodpal_weight_type&            weight_parameter();
    autosomal_model_type&          autosomal_model();
    x_linkage_model_type&          x_linkage_model();
    
    typedef std::pair<trait_iterator, bool>      tib_value;
    typedef std::pair<marker_iterator, bool>     mib_value;
    typedef std::pair<covariate_iterator, bool>  cib_value;

    void      clear_subsets();
    tib_value add_subset(const trait_parameter& t);
    tib_value set_subset(const trait_parameter& t);
    tib_value add_subset(size_t t);
    tib_value set_subset(size_t t);

    void      clear_traits();
    tib_value add_trait(const trait_parameter& t);
    tib_value set_trait(const trait_parameter& t);
    tib_value add_trait(size_t t, trait_parameter::trait_type tt = trait_parameter::conaff, double cp = 0.0);
    tib_value set_trait(size_t t, trait_parameter::trait_type tt = trait_parameter::conaff, double cp = 0.0);

    void      clear_parameters();
    void      clear_markers();
    void      clear_covariates();
    mib_value add_parameter(const marker_type& m);
    cib_value add_parameter(const covariate_type& c);

    mib_value add_marker(size_t m,
                         marker_type::inheritance_type in = marker_type::autosomal,
                         double b1 = 0.1,
                         double b2 = 0.1);
                         
    cib_value add_covariate(size_t c,
                            covariate_type::cov_ad ad = covariate_type::mean,
                            covariate_type::cov_op op = covariate_type::sum,
                            double ad_val = std::numeric_limits<double>::quiet_NaN(),
                            double power = 1.0,
                            double d1 = 0.1,
                            double d2 = 0.1);

    void      clear_weight();
    bool      add_weight(const lodpal_weight_type& w);
    bool      set_weight(const lodpal_weight_type& w);
    bool      add_weight(size_t w, lodpal_weight_type::weight_op op = lodpal_weight_type::pair);
    bool      set_weight(size_t w, lodpal_weight_type::weight_op op = lodpal_weight_type::pair);

    void      clear_autosomal_model();
    bool      add_autosomal_model(const autosomal_model_type& m);
    bool      set_autosomal_model(const autosomal_model_type& m);
    bool      add_autosomal_model(autosomal_model_type::model_type      m = autosomal_model_type::one_parameter,
                                  autosomal_model_type::constraint_type cont = autosomal_model_type::constrained,
                                  autosomal_model_type::poo_fixed_type  fixed = autosomal_model_type::none,
                                  bool poo_test = false,
                                  double al = 2.634);
    bool      set_autosomal_model(autosomal_model_type::model_type      m = autosomal_model_type::one_parameter,
                                  autosomal_model_type::constraint_type cont = autosomal_model_type::constrained,
                                  autosomal_model_type::poo_fixed_type  fixed = autosomal_model_type::none,
                                  bool poo_test = false,
                                  double al = 2.634);

    void      clear_x_linkage_model();
    bool      add_x_linkage_model(const x_linkage_model_type& m);
    bool      set_x_linkage_model(const x_linkage_model_type& m);
    bool      add_x_linkage_model(bool lambda1 = true, bool lambda2 = true, double al = 2.634);
    bool      set_x_linkage_model(bool lambda1 = true, bool lambda2 = true, double al = 2.634);

  private:
    bool                      my_valid;
    bool                      my_skip_uninformative_pairs;
    bool                      my_use_pair_cache;
    bool                      my_multipoint;
    bool                      my_turn_off_default;
    bool                      my_print_lambda;
    bool                      my_sib_pairs_only;
    bool                      my_autosomal_marker_exist;
    bool                      my_x_linked_marker_exist;
    bool                      my_use_mm_pair;
    bool                      my_use_mf_pair;
    bool                      my_use_ff_pair;

    size_t                    my_total_used_pair_type;
    size_t                    my_diagnostic_marker;
    size_t                    my_max_marker_name_size;
    
    trait_vector              my_traits;
    trait_vector              my_subsets;
    independent_variable      my_parameters;
    lodpal_weight_type        my_weight;

    autosomal_model_type      my_autosomal_model;
    x_linkage_model_type      my_x_linkage_model;
};


#include "lodpal/lodpal_params.ipp"

} // end of namespace LODPAL
} // end of namespace SAGE

#endif
