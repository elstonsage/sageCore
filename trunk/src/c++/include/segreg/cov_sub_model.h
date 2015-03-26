#ifndef SEGREG_COV_SUB_MODEL_H
#define SEGREG_COV_SUB_MODEL_H
//============================================================================
// File:      cov_sub_model.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   7/9/01 - created.                            djb
//                                                                          
// Notes:     defines the mean covariate sub-models for SEGREG.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include "sub_model_base.h"
#include "type_sub_model.h"
#include "fped/fped.h"
#include "rped/rped.h"

namespace SAGE
{

namespace SEGREG
{

/// @name covariate sub-model constants
//@{
extern const std::string  MEAN_COV_NAME;
extern const std::string  VAR_COV_NAME;
extern const std::string  SUSCEPT_COV_NAME;
extern const std::string  COMPOSITE_TRAIT_NAME;
extern const double       COV_DEFAULT_VALUE;
extern const bool         COV_DEFAULT_FIXED;
extern const bool         COV_DEFAULT_INTERACTION;
extern const double       INTERACTION_INIT_VALUE;
//@}

// Covariates know the model exists to permit friend access, but really don't 
//need to know anything else about it.
class model;

// We don't want default constructor for the covariate models, so let's turn off
// the lint warnings for all the covariate models:
//lint -esym(1712, *CovariateSubmodel, CompositeTraitSubmodel)

/// \brief base class for all covariate related sub models
///
/// The CovariateSubmodel is the base class upon which all covariate
/// related submodels are built.  This includes the mean, variance, susceptibility
/// and composite trait.
///
/// Covariates have an associated 'mean', a double value which is actually controlled
/// externally and should represent the mean value of the covariate. XXX This
/// should be FIXED. XXX (sampling will do it)
class CovariateSubmodel : public SegregSubmodel
{
  public:
    friend class model;
  
    /// Enumeration describing the types of covariates
    enum CovariateTypeEnum
    {
        ct_MEAN, ///< Mean Covariate
        ct_VAR,  ///< Variance Covariate
        ct_SUSC, ///< Susceptibility Covariate
        ct_COMP  ///< Composite Trait Covariate
    };

    /// \brief Storage for a covariate.
    ///
    /// The covariate struct stores all the information needed to include
    /// a particualr covariate in a SEGREG analysis.
    struct covariate
    { 
      covariate();
      covariate(const std::string& trait,
                size_t             index,
                double             coeff,
                bool               inter,
                bool               fixed);

      std::string  trait_name;              ///< Name of trait in FPED::Multipedigree
      size_t       trait_index;             ///< Index of trait in FPED::Multipedigree
      double       coefficient;             ///< Starting/Current coefficient
      bool         has_interaction;         ///< Indicates an interaction
      bool         is_fixed;                ///< Indicates a fixed value
      double       i_taus[NUM_OF_TYPES];    ///< When interactions are used, coefficients
                                            ///< for the interaction terms.
      
      size_t       maxfun_index;            ///< \internal 
                                            ///< Index of first maxfun parameter that applies to
                                            ///< this covariate.
    };
  
    /// \name Basic Object Information
    //@{
    virtual std::string  name() const = 0;

    virtual std::string option_description() const;

    /// Returns the name of the sub-block for the current type
    ///
    virtual std::string get_subblock_name() const = 0;

    /// Returns the name of the sub-block for the given covariate type
    static std::string get_cov_type_subblock_name(CovariateTypeEnum type);

    /// Dump the current status to the ostream
    ///
    /// \param out output stream
    void  dump(std::ostream& out) const;

    /// Pretty prints a covariate table to the ostream
    /// 
    /// \param out    output stream
    /// \param indent Number of characters to indent the table
    void  print_covariate_table(std::ostream& out, size_t indent = 0) const;
    //@}

    /// \name Covariate access
    ///
    //@{

    /// Returns the number of covariates in the model
    ///
    size_t get_covariate_count() const;

    /// Returns the currently stored mean of the covariate at i
    ///
    /// \param i The covariate's index
    double get_covariate_mean(size_t i) const;
    /// Sets the covariate's mean for the covariate at i to d.
    ///
    /// \param i the covariate index
    /// \param d the value to set the mean to
    void   set_covariate_mean(size_t i, double d);
    
    /// Calculate the covariate means based on the RefMultiPedigree's data
    ///
    void   calculate_covariate_means(const FPED::Multipedigree&);

    /// Calculate the mean of the covariate at i based upon the RefMultiPedigree's
    /// data
    void   calculate_covariate_mean(const FPED::Multipedigree&, size_t i);

    /// Get a reference to a vector of the current covariates with their values.
    ///
    const std::vector<covariate>&  covariates() const;

    /// Returns \c true if there is a covariate with the same name as the covariate
    /// given, \c false otherwise
    bool  has_covariate(const covariate& cov) const;
    
    /// Returns \c true if there is a covariate with the name given, 
    /// \c false otherwise.
    bool  has_covariate(const std::string& cov) const;

    /// Return the index of the covariate with the name given 
    size_t covariate_index(const std::string&) const;
    //@}

  protected:
  
    /// \name Object Management
    //@{
    CovariateSubmodel(const genotype_specific_mean_susc_sub_model* m_ptr,
                          cerrorstream& errors);
    CovariateSubmodel(const CovariateSubmodel& other);
    CovariateSubmodel&  operator=(const CovariateSubmodel& other);

    virtual ~CovariateSubmodel();
    //@}

    /// Adds a covariate given the values specified.  Note
    /// that using this version of the function when the
    /// covariate type is ct_COMP will result in a compile-time failure.
    ///
    /// \param mp           The multipedigree (for looking up data
    /// \param trait_name   The name of the trait
    /// \param coefficient  Initial value of the covariate coefficient
    /// \param interaction  Covariate has interactions
    /// \param fixed        Is the coefficient fixed?
    /// \param type_missing Interactions are only permitted when there
    ///                     are two or three main types (or the main type is
    ///                     not known)
    bool base_add_covariate
          (const RPED::MultiPedigree*    mp,
           const string&                 trait_name,
           double                        coefficient,
           bool                          interaction,
           bool                          fixed, 
           bool                          type_missing);

  private:

    /// \name MAXFUN related functions
    ///
    /// These functions are used by MAXFUN (or other MAXFUN related functions)
    /// to set this submodel up for maximization
    //@{
    virtual int  update();
    /// This function builds the my_parameters structure for inclusion into the
    /// maxfunapi. We don't do this until right before the maximization because
    /// the interactions change the number of maximization parameters according
    /// to the number of means, which we have no control over. Thus, we delay
    /// the building of the my_parameters until the last possible instant.
    ///
    /// \return 0 Always Successful.
    virtual int finalizeConfiguration();
    
    /// Adds interaction parameters to maxfun's list of parameters when
    /// required.
    void initialize_interactions(vector<covariate>::const_iterator c_iter);
    //@}
    
    std::vector<covariate>::iterator  find(const covariate& cov);

    void  insert_covariate(const covariate& cov);                 

    /// \name Data
    //@{ 
    mutable std::vector<covariate>  my_covariates;
    std::vector<double>     my_covariate_means;
    
    const genotype_specific_mean_susc_sub_model*  my_m_ptr;
    //@}
};

/// \brief Submodel for covariate storage
///
/// The CovariateSubmodel stores covariates of a given type.
/// Types include mean, variance, susceptibility
/// and composite trait.
///
/// The CovariateSubmodel is templatized on the type of covariates it contains
/// so that certain features do not have to be done multiple times
template <CovariateSubmodel::CovariateTypeEnum COV_TYPE>
class TypedCovariateSubmodel : public CovariateSubmodel
{
  public:
    friend class model;
  
    /// \name Object Management
    //@{
    TypedCovariateSubmodel(const genotype_specific_mean_susc_sub_model* m_ptr,
                          cerrorstream& errors);
    TypedCovariateSubmodel(const TypedCovariateSubmodel& other);
    TypedCovariateSubmodel&  operator=(const TypedCovariateSubmodel& other);

    virtual ~TypedCovariateSubmodel();
    //@}
    
    /// \name Basic Object Information
    //@{
    /// Returns the 'name' of the covariate model type.
    ///
    virtual std::string  name() const;

    /// Returns the name of the sub-block for the current type
    ///
    virtual std::string get_subblock_name() const;

    /// \name Covariate access
    ///
    //@{

    /// Adds a covariate given the values specified.  Note
    /// that using this version of the function when the
    /// covariate type is ct_COMP will result in a compile-time failure.
    ///
    /// \param mp           The multipedigree (for looking up data
    /// \param trait_name   The name of the trait
    /// \param coefficient  Initial value of the covariate coefficient
    /// \param interaction  Covariate has interactions
    /// \param fixed        Is the coefficient fixed?
    /// \param type_missing Interactions are only permitted when there
    ///                     are two or three main types (or the main type is
    ///                     not known)
    bool add_covariate
          (const RPED::MultiPedigree*    mp,
           const string&                 trait_name,
           double                        coefficient,
           bool                          interaction,
           bool                          fixed, 
           bool                          type_missing);

    /// Adds a covariate given the values specified.  Note
    /// that using this version of the function when the
    /// covariate type is NOT ct_COMP will result in a compile-time failure.
    ///
    /// \param mp           The multipedigree (for looking up data
    /// \param trait_name   The name of the trait
    /// \param coefficient  Initial value of the covariate coefficient
    /// \param fixed        Is the coefficient fixed?
    /// \param type_missing Interactions are only permitted when there
    ///                     are two or three main types (or the main type is
    ///                     not known)
    bool add_covariate
          (const RPED::MultiPedigree*    mp,
           const string&                 trait_name,
           double                        coefficient,
           bool                          fixed, 
           bool                          type_missing);

  private:
  
};

typedef TypedCovariateSubmodel<CovariateSubmodel::ct_MEAN> MeanCovariateSubmodel;
typedef TypedCovariateSubmodel<CovariateSubmodel::ct_VAR>  VarianceCovariateSubmodel;
typedef TypedCovariateSubmodel<CovariateSubmodel::ct_SUSC> SusceptibilityCovariateSubmodel;
typedef TypedCovariateSubmodel<CovariateSubmodel::ct_COMP> CompositeTraitSubmodel;

/// Returns a list of trait names which appear in both CovariateSubmodel
///
/// \param sm1 The first  submodel
/// \param sm2 The second submodel
std::list<std::string> get_shared_covariates(const CovariateSubmodel& sm1,
                                             const CovariateSubmodel& sm2);

/// Compares the trait lists of two CovariateSubmodel objects.  Returns \c true
/// if there are no traits shared in the two lists, \c false otherwise
///
/// \param sm1 The first  submodel
/// \param sm2 The second submodel
bool are_covariates_exclusive(const CovariateSubmodel& sm1,
                              const CovariateSubmodel& sm2);


} 
} 
#include "segreg/cov_sub_model.ipp"

#endif

