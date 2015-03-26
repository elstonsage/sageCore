#ifndef MLM_RESID_CORR_H
#define MLM_RESID_CORR_H
//======================================================================================
//  File:       mlm_resid_corr.h
//
//  Author:     Kai He
//
//  History:    09.2003 
//
//
//  Copyright (c) 2001 R. C. Elston
//  All rights reserved
//======================================================================================

#include "segreg/member_calculator.h"
#include "segreg/PedigreeDataSet.h"
#include "segreg/model.h"
#include "boost/array.hpp"
#include "numerics/sinfo.h"

namespace SAGE
{
namespace SEGREG
{

const int NUM_CORR      = 8;
const int NUM_CLASSES   = 4;
const int NUM_GENOTYPES = 3;

/// \brief Calculates Residual Correlations based upon the Residual Associations in the model.
///
/// This class implements SAGE4.4 internal Documentation Appendix G.
///
/// It performs its calculations based upon the original covariates as provided in the
/// multipedigree, not the centered, etc. values in the binary member calculator.
/// This is according to Yi and Dr. Elston.
class MlmResidCorrelationCalculator
{
  friend class primary_analysis_results;
  
  public:
    typedef residual_correlation_sub_model                      resid_model;
    typedef pair<string, double>                                dopair;
    typedef residual_correlation_sub_model::corr                corr;
     
    /// \name Object Management
    //@{
     
    /// Default Constructor
    ///
    MlmResidCorrelationCalculator ();
    
    /// Constructor
    ///
    /// This version calculates our correlations based upon a set of subpedigrees
    /// and a model (using the calculate() function).
    ///
    /// \param b binary_member_calculator used for testing individual validity, primarily
    /// \param l All subpedigrees we'll use for our calculation
    /// \param m The model containing the residual associations and other parameters
    MlmResidCorrelationCalculator (const binary_member_calculator& b,
                            const PedigreeDataSet::SubpedigreeCursor& l, 
                            const model& m);
    
    /// Copy Constructor
    ///
    MlmResidCorrelationCalculator (const MlmResidCorrelationCalculator&);
    
    /// Copy Operator
    ///
    MlmResidCorrelationCalculator& operator=(const MlmResidCorrelationCalculator&);
                                      
    /// Destructor
    ///
    ~MlmResidCorrelationCalculator() {}
    //@}
    
    /// Calculates, for a given set of subpedigrees, four correlations: the
    /// father-mother, father-offspring, mother-offspring, and sib-sib.
    /// In cases where there are covariates, the father-offspring and mother-offspring
    /// are calculated separately, even if there is a single parental association,
    /// which may result in slightly different correlations.
    ///
    /// \param b binary_member_calculator used for testing individual validity, primarily
    /// \param l All subpedigrees we'll use for our calculation
    /// \param m The model containing the residual associations and other parameters
    void calculate(const binary_member_calculator& b,
                   const PedigreeDataSet::SubpedigreeCursor& l, 
                   const model& m);
    
    /// For testing and output checking, dumps the correlations and calculations
    /// to the screen.
    ///
    /// \param m Model used for determining names of variables and such.
    void   dump_data(const model& m);
    
    /// Returns the correlations for the particular pair type.  Note that
    /// children are not sorted by sex, so father-son = father-daughter,
    /// mother-son = mother-daughter and
    /// brother-brother = brother-sister = sister-sister.  This may change in
    /// the future.
    ///
    /// \param c The correlation desired
    double get_corr_ij  (resid_model::corr c);
    
    /// Returns the correlations for the particular pair type for the particular
    /// type (AA, AB, or BB).  Note that
    /// children are not sorted by sex, so father-son = father-daughter,
    /// mother-son = mother-daughter and
    /// brother-brother = brother-sister = sister-sister.  This may change in
    /// the future.
    ///
    /// \param c     The correlation desired
    /// \param gtype The genotype of the correlation desired.
    double get_corr_ij_u(resid_model::corr c,size_t gtype);

  private:

    /// \brief stores information about each covariate while performing our calculations
    ///
    /// This struct is used internally by the calculate() function to store information
    /// about each covariate.
    struct class_cov
    {
      size_t index;
      string name;
      double coefficient;
      bool   interaction;
      double i_taus[NUM_GENOTYPES];
      SampleInfo mo_sinfo;
      SampleInfo fa_sinfo;
      SampleInfo pa_sinfo;
      SampleInfo ch_sinfo;
    };
    
    typedef std::list<const FPED::Member*> member_list;

    /// Classifies each member in the subpedigree list into the father, mother
    /// and/or child lists.  Note that a member can belong to the child list
    /// and one of the parental lists simultaneously.
    void classify_members(const binary_member_calculator& b,
                          const PedigreeDataSet::SubpedigreeCursor& spl);

    /// Calculate the \f$theta_u(i)\f$ for each member class (father, mother, child)
    ///
    /// \param m The model
    void calculate_theta_u_i(const model& m);

    /// Calculate the \f$theta(i)\f$ for each member class (father, mother, child)
    ///
    /// \param m The model
    void calculate_theta_i  (const model& m);

    /// Calculate the correlations given the model
    ///
    /// \param m The model holding the residual correlations and type frequencies
    void   calculate_corrs(const model&);
    
    /// Calculates the correlations for the specific relative pair type
    ///
    /// \param c        The pair type to calculate
    /// \param delta_ij The residual association for the pair type
    /// \param theta_i  The theta (see appendix) for members of type i
    /// \param theta_j  The theta (see appendix) for members of type j
    void   calculate_ij_corrs(resid_model::corr c, double delta_ij, size_t theta_i, size_t thetat_j);

    /// Calculates the correlation, given the association (delta) and the
    /// means for members groups i and j.  See Appendix G.  Works for both
    /// general and for type-specific correlations.
    ///
    /// \param delta_ij The residual association for the pair type
    /// \param theta_i  The theta (see appendix) for members of type i
    /// \param theta_j  The theta (see appendix) for members of type j
    double   calculate_corr(double delta_ij, double theta_i, double theta_j) const;
    
    /// \name Evaluation parameters
    //@{
    /// Storage of the covariate mean sums for each individaul class
    ///
    boost::array<double, NUM_CLASSES> my_means;          
    /// Storage of the \f$theta(i)\f$ for each individual class
    ///
    boost::array<double, NUM_CLASSES> my_theta_i;
    /// Storage of the \f$theta_u(i)\f$ for each genotype for each individual class
    ///
    boost::array<boost::array<double, NUM_GENOTYPES>, NUM_CLASSES> my_theta_u_i;

    /// Storage of the covariate interaction mean sums for each genotype for each individual class
    ///
    boost::array<boost::array<double, NUM_GENOTYPES>, NUM_CLASSES> my_tau_means;
    //@}
    
    /// \name Member Lists
    //@{
    /// List of members who are female and a parent
    ///
    member_list my_mother_list;

    /// List of members who are male and a parent
    ///
    member_list my_father_list;
    
    /// List of members who are a parent.  Is the join of the mother and father lists
    ///
    member_list my_parent_list;

    /// List of members who have parents in the data set
    ///
    member_list my_child_list;
    //@}
    
    /// Correlation parameters
    //@{
    /// Storage of the calculated correlations for each pair type
    ///
    boost::array<double, NUM_CORR>                              my_corr_ij;
    /// Storage of the calculated correlations for each genotype for each pair type
    ///
    boost::array<boost::array<double, NUM_GENOTYPES>, NUM_CORR> my_corr_ij_u;     
    //@}
};

}
}

#endif


