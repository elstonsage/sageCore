#ifndef SEGREG_CALCULATOR_H
#define SEGREG_CALCULATOR_H
//====================================================================================
// File:        segreg_calculator.h
//
// Purpose:     Calculation for segregation likelihood. this class creates calculator
//              objects based on the models and tests valid subpedigrees by using
//              multipedigree utilities.
//
// Author:      Kai He
//
// History:     07/26/2001 initialed
//
// Copyright (c) 2001 R. C. Elston
//===================================================================================


#include "segreg/LikelihoodElements.h"
#include "segreg/regressive_peeler.h"
#include "segreg/FPMM_peeler.h"
#include "segreg/mlm_peeler.h"
#include "segreg/MlmCorrelationVerifier.h"
#include "segreg/mlm_resid_corr.h"
#include "segreg/SL_calculator.h"
#include "segreg/RegUnconnectedLikelihood.h"
#include "segreg/PedigreeDataSet.h"
#include "segreg/segreg_datatypes.h"
#include "segreg/pen_function.h"
#include "mped/sp.h"
#include "mped/mp.h"
#include "mped/mp_utilities.h"
#include "maxfun/maxfun.h"
#include "numerics/kahan.h"
#include "numerics/log_double.h"
#include "numerics/cephes.h"
#include "numerics/functions.h"
#include <cmath>

namespace SAGE
{
namespace SEGREG
{

/// The segreg_calculator calculates the likelihood of the pedigree data under a model.

/// Likelihood for the pedigree:
/// \f[In \ Regressive \ Models \ and \ MLM(multivariate \ logistic \ model):\f]
/// \f[L(P)=\sum_{u_i}L(P, \ u_i)\f]
/// \f[In \ FPMM:\f]
/// \f[L(P)=\sum_{u_i}\sum_{v_i}L(P, \ u_i, \ v_i)\f]
///
/// Posterior probability of Genotype Ui for individual i, given all the pedigree data.
/// \f[In \ Regressive \ Models \ and \ MLM:\f]
/// \f[P(u_i|P)=\frac{L(P, \ u_i)}{L(P)} \quad for \ i = 1,2,...,n\f]
/// \f[In \ FPMM:\f]
/// \f[P(u_i|P)=\frac{\sum_{v_i}L(P, \ u_i, \ v_i)}{L(P)} \quad for \ i=1,2,...,n\f]
///
/// Joint Likelihood for the Pedigree and the genotype Ui:
/// \f[In \ Regressive \ Model:\f]
/// \f[L(P,\ u_i)=ant_i(u_i)\prod_{s \in S_i}pos_{is}(u_i)\f]
/// \f[In \ FPMM \ with \ genotype \ Ui \ and \ polygenic \ number \ Vi:\f]
/// \f[L(P, \ u_i, \ v_i)=ant_i(u_i,v_i)\prod_{s \in S_i}pos_{is}(u_i, v_i)\f]
///
/// Likelihood of regressive model for unconnected individual.
/// \f[L(y_1,y_2,...,y_k)=\prod_{j=1}^{k}\sum_{u_j}\psi_j(u_j)P(y_j|u_j)\f]
///
/// Likelihood of MLM model for unconnected individual.
/// \f[L(y_1,y_2,...,y_k)=\prod_{j=1}^{k}\sum_{u_j}\psi_j(u_j)frac{e^{\theta_u(i)y_i}}{1 + e^{\theta_u(i)}}\f]
///
/// Likelihood of MLM model.
/// \f[L(P,\ u_i)=sum_{\Lambda_0}{L(F_0)[\Pi_{\Omega_1}\sum_{\Lambda_{1j}}\frac{L(F_{1j})}{L(1j)}
///  [...\Pi_{\Omega_{kj}}\sum_{\Lambda_{kj}}\frac{L(F_{kj})}{L(kj)}]...]}\f]

class segreg_calculator : public MaxFunction
{
  public:

    // Constructors  and destructors

    segreg_calculator  (const PedigreeDataSet&,
                        const model&); 
    ~segreg_calculator ();

    // Create segreg_calculator components

    void set_continuous_penalty_component (double c);

    double calculate_continuous_penalty() const;

    log_double calculate_prevalence_penalty() const;

    // Get the count of subpedigrees and unconnecteds used in the analysis

    uint get_subpedigree_count() const;
    uint get_unconnected_count() const;


    // Likelihood calculation including both connected and unconnected
    // members in all pedigrees
 
    log_double calculate();

    // Penetrance calculations for all individuals
    
    void calculate_post_geno_probs(pg::post_geno_map&);

    void calculate_pen_func_probs(pf::pen_func_map&);

    // Calculat mlm residual aasociations and correlations
        
    void calculate_mlm_resid_corr(MlmResidCorrelationCalculator&) const;    

  protected:

    typedef FPED::SubpedigreeConstIterator   subpedigree_iterator;
    typedef FPED::SubpedigreeConstPointer    subpedigree_pointer;
    typedef FPED::Member                     member_type;
    typedef FPED::MemberConstPointer         member_const_pointer;
    typedef FPED::MemberConstIterator        member_const_iterator;
    typedef FPED::Family                     family_type;
    typedef FPED::FamilyConstPointer         family_const_pointer;

    void setup_components             ();

    void setup_regressive_components  ();
    void setup_mlm_components         ();
    void setup_fpmm_components        ();
    
    bool using_ascertainment() const;

    //MaxFunction' methods from maxfun.h

    virtual double evaluate          (parameter_vector&);
    virtual int    update_bounds     (parameter_vector&); 
            double internal_evaluate (parameter_vector&);

    //----------------------------------------------------------------------
    // Calculation based on model type
    
    log_double calc_likelihood ();

    //----------------------------------------------------------------------
    // Regressive model Pedigree likelihood calculation

    void calc_connected_regressive   ( );
    void calc_connected_FPMM         ( );
    void calc_connected_MLM          ( );

    log_double calc_connected_regressive   (const FPED::Subpedigree& ps, bool asc = false);
    log_double calc_connected_FPMM         (const FPED::Subpedigree& ps, bool asc = false);
    log_double calc_connected_MLM          (const FPED::Subpedigree& ps, bool asc = false);

    void calc_unconnected_regressive ( );
    void calc_unconnected_FPMM       ( );
    void calc_unconnected_MLM        ( );

    log_double calc_unconnected_regressive (const FPED::Member&, bool asc = false);
    log_double calc_unconnected_FPMM       (const FPED::Member&, bool asc = false);
    log_double calc_unconnected_MLM        (const FPED::Member&, bool asc = false);

    //----------------------------------------------------------------------

    void calc_connected_pgeno_regressive (const FPED::Subpedigree&, pg::post_geno_map&);
    void calc_connected_pgeno_FPMM       (const FPED::Subpedigree&, pg::post_geno_map&);
    void calc_connected_pgeno_MLM        (const FPED::Subpedigree&, pg::post_geno_map&);

    void calc_unconnected_pgeno_regressive (const FPED::Member&, pg::post_geno_map&);
    void calc_unconnected_pgeno_FPMM       (const FPED::Member&, pg::post_geno_map&);
    void calc_unconnected_pgeno_MLM        (const FPED::Member&, pg::post_geno_map&);

    void calc_member_pgeno_regressive  (pg::post_geno_map& p, regressive_peeler& peeler,
                                        const member_type& ind);
    void calc_member_pgeno_FPMM        (pg::post_geno_map& p, const member_type& ind);
    void calc_member_pgeno_MLM         (pg::post_geno_map& p, const member_type& ind);

    //----------------------------------------------------------------------
    // Penetrance functions for continous data

    void calc_connected_pen_func   (const FPED::Subpedigree&, pf::pen_func_map&);

    void calc_unconnected_pen_func (const FPED::Member&, pf::pen_func_map&);

    void calc_nonfounder_pen_func  (const FPED::Member&, pf::pen_func_map&);
    void calc_founder_pen_func     (const FPED::Member&, pf::pen_func_map&);

    //----------------------------------------------------------------------
    // Calculates subpedigree likelihood given a member and their genotype
    
    log_double calc_subped_given_member_regressive(regressive_peeler*  rpl,
                                                   const FPED::Member& indi,
                                                   const TypeDescription::State& genotype);

    log_double calc_subped_given_member_MLM       (mlm_peeler*         mpl, 
                                                   const FPED::Member& i,  
                                                   const TypeDescription::State& genotype);

    log_double calc_subped_given_member_FPMM      (FPMM_peeler*        fpl,
                                                   const FPED::Member& indi,
                                                   genetic_info        genoinfo);
		

    /// Initializes likelihood1 and likelihood2 to 0.0 at the beginning of a calculation.

    void clear_likelihood();

    /// Accumulates pedigree likelihoods into likelihood1 and likelihood2. 
    /// Returns true if the likelihood remains finite.
    bool accumulate_likelihood(log_double likelihood);

    /// Combines likelihood1 and likelihood2 (as required) into a final likelihood1.

    bool finalize_likelihood();

    
    //----------------------------------------------------------------------
    // Data members
    //----------------------------------------------------------------------

    const PedigreeDataSet&            my_ped_data;
    const model                    &  md;

    LikelihoodElements             *  my_like_elts;
    LikelihoodElements             *  my_asc_like_elts;
    FPMM_SL                        *  fpmmsl;
    FPMM_SL                        *  asc_fpmmsl;
    mlm_peeler                     *  mpl;

    binary_member_calculator       *  bmc;
    binary_member_calculator       *  abmc;

    std::auto_ptr<MlmCorrelationVerifier>  my_mlm_corr_verifier;    
    //MlmResidCorrelationCalculator       *  mlm_resid_corr;

    double                            last_likelihood;

    double                            c_penalty_value;

    log_double                        likelihood1;
    log_double                        likelihood2;

    static vector<double> reg_ant_vec; // due to JA for likelihood rescaling
    static vector<double> reg_pos_vec; // see above
    static vector<double> mlm_ant_vec; // due to JA for likelihood rescaling
    static vector<double> mlm_pos_vec; // see above
    static vector<double> fpmm_ant_vec; // due to JA for likelihood rescaling
    static vector<double> fpmm_pos_vec; // see above

};

}
}

#include "segreg/segreg_calculator.ipp"

#endif
