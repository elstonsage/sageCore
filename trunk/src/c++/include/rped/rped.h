#ifndef RPED_NEW_H
#define RPED_NEW_H

//============================================================================
//  Reference Pedigree structure -- Definition of General Reference pedigree
//                                  storage classes
//
//  Author: Kevin Jacobs (jacobs@darwin.cwru.edu)
//
//  History:   0.1  kbj Initial implementation                   May 07 98
//                  yjs Added a new class PhenotypeReaderInfo    Jan 2003
//                       & new read_..() function in RefMpedInfo.
//                  djb Added a check to PedigreeSort().          8/1/3
//
//  Copyright (c) 1998  R.C. Elston
//    All Rights Reserved
//============================================================================

#include <string>
#include <cmath>
#include <vector>
#include <set>
#include <iostream>
#include "numerics/isnan.h"
#include "globals/config.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "error/errorbuf.h"
#include "error/internal_error.h"
#include "fortran/Tokenizer.h"
#include "fortran/Token_func.h"
#include "mped/mp.h"
#include "mped/sp.h"
#include "mlocus/imodel.h"

namespace SAGE {
namespace RPED {

template <class GI=MPED::no_info, class FI=MPED::no_info, class SI=MPED::no_info, class PI=MPED::no_info, class MI=MPED::no_info>
class PedData
{
public:
    typedef GI       indinfo_type;
    typedef FI       faminfo_type;
    typedef SI       subinfo_type;
    typedef PI       pedinfo_type;
    typedef MI       mpinfo_type;

    PedData();

    PedData(const MPED::multipedigree_base& m);

    size_t     pedigree_info_count()           const;
    size_t  subpedigree_info_count(size_t p)   const;
    size_t       family_info_count(size_t p)   const;
    size_t   individual_info_count(size_t p)   const;

    size_t  subpedigree_info_count(const MPED::pedigree_base& p)   const;
    size_t       family_info_count(const MPED::pedigree_base& p)   const;
    size_t   individual_info_count(const MPED::pedigree_base& p)   const;

    const mpinfo_type&   multipedigree_info()                   const;
    const pedinfo_type&       pedigree_info(size_t p)           const;
    const subinfo_type&    subpedigree_info(size_t p, size_t s) const;
    const subinfo_type&         family_info(size_t p, size_t f) const;
    const indinfo_type&     individual_info(size_t p, size_t i) const;

    mpinfo_type&   multipedigree_info();
    pedinfo_type&       pedigree_info(size_t p);
    subinfo_type&    subpedigree_info(size_t p, size_t s);
    subinfo_type&         family_info(size_t p, size_t f);
    indinfo_type&     individual_info(size_t p, size_t i);

    const pedinfo_type&        pedigree_info(const MPED::pedigree_base& p)               const;
    const subinfo_type&     subpedigree_info(const MPED::pedigree_base& p, size_t s)     const;
    const subinfo_type&     subpedigree_info(const MPED::subpedigree_base& s)            const;
    const faminfo_type&          family_info(const MPED::pedigree_base& p, size_t f)     const;
    const faminfo_type&          family_info(const MPED::subpedigree_base& s, size_t f)  const;
    const faminfo_type&          family_info(const MPED::family_base& f)                 const;
    const faminfo_type&          family_info(const MPED::member_base& m)                 const;
    const indinfo_type&      individual_info(const MPED::pedigree_base& p, size_t i)     const;
    const indinfo_type&      individual_info(const MPED::subpedigree_base& s, size_t i)  const;
    const indinfo_type&      individual_info(const MPED::member_base& m)                 const;

    pedinfo_type&        pedigree_info(const MPED::pedigree_base& p);
    subinfo_type&     subpedigree_info(const MPED::pedigree_base& p, size_t s);
    subinfo_type&     subpedigree_info(const MPED::subpedigree_base& s);
    faminfo_type&          family_info(const MPED::pedigree_base& p, size_t f);
    faminfo_type&          family_info(const MPED::subpedigree_base& s, size_t f);
    faminfo_type&          family_info(const MPED::family_base& f);
    faminfo_type&          family_info(const MPED::member_base& m);
    indinfo_type&      individual_info(const MPED::pedigree_base& p, size_t i);
    indinfo_type&      individual_info(const MPED::subpedigree_base& s, size_t i);
    indinfo_type&      individual_info(const MPED::member_base& m);

    mpinfo_type&                  info();
    pedinfo_type&                 info(const MPED::pedigree_base& p);
    subinfo_type&                 info(const MPED::subpedigree_base& s);
    faminfo_type&                 info(const MPED::family_base& f);
    indinfo_type&                 info(const MPED::member_base& m);

    const mpinfo_type&            info()                                          const;
    const pedinfo_type&           info(const MPED::pedigree_base& p)              const;
    const subinfo_type&           info(const MPED::subpedigree_base& s)           const;
    const faminfo_type&           info(const MPED::family_base& f)                const;
    const indinfo_type&           info(const MPED::member_base& m)                const;

    pedinfo_type&           operator[](const MPED::pedigree_base& p);
    subinfo_type&           operator[](const MPED::subpedigree_base& s);
    faminfo_type&           operator[](const MPED::family_base& f);
    indinfo_type&           operator[](const MPED::member_base& m);

    const pedinfo_type&     operator[](const MPED::pedigree_base& p)              const;
    const subinfo_type&     operator[](const MPED::subpedigree_base& s)           const;
    const faminfo_type&     operator[](const MPED::family_base& f)                const;
    const indinfo_type&     operator[](const MPED::member_base& m)                const;

    void update_storage(const MPED::multipedigree_base& m);
    void update_storage(const MPED::pedigree_base& p);
    void update_storage(const MPED::subpedigree_base& s);
    void update_storage(const MPED::family_base& f);
    void update_storage(const MPED::member_base& m);

private:
    typedef std::vector<indinfo_type>    indinfo_storage_type;
    typedef std::vector<faminfo_type>    faminfo_storage_type;
    typedef std::vector<subinfo_type>    subinfo_storage_type;

    struct ped_storage
    {
        pedinfo_type                     ped_info;
        subinfo_storage_type             sub_data;
        faminfo_storage_type             fam_data;
        indinfo_storage_type             ind_data;
    };

    struct mped_storage
    {
        mpinfo_type                      mped_info;
        std::vector<ped_storage>         ped_data;
    };

    mped_storage                         my_storage;
};

/** \brief Describes a single trait
  *
  * \par Introduction
  *
  * RefTraitInfo describes a single trait, and is subsequently used within RefMPedInfo.
  *
  * \par Basic attributes
  *
  * Every trait has the name three basic traits: name, type, and usage.
  *
  * The first (name) is simple enough: a name associated with the trait.
  *
  * Trait type (specified by the enum trait_t) describes the type of data that the trait represents.
  * Most of the time, traits will be either continuous, discrete, or binary. There are two other types (unknown
  * and invalid) available in case such distinctions are needed.
  *
  * Trait usage (specified by the enum trait_use) is useful for categorizing the intended analytical use of
  * the trait. A SAGE program may be designed to iterate through all traits of type trait_variate, for instance.
  *
  * \par Missingness
  *
  * Since individuals may very well be uninformative for some traits, the RefTraitInfo must be configured to
  * handle missingness. Furthermore, since trait data is read in
  * as a string value (but is then converted to a numeric value), the missing code system must handle both strings
  * and numbers.
  *
  * In setting up your RefTraitInfo instance, make sure to set both the string-type and numeric-type missing codes.
  * This is done with set_numeric_missing_code() and set_string_missing_code().
  *
  * \par Binary traits
  *
  * If your RefTraitInfo's type is binary_trait, then you'll need to set up the RefTraitInfo for reading
  * in binary values.
  *
  * Please note that the RefTraitInfo is designed to organize binary values conceptually as affected/unaffected;
  * although it's entirely reasonable to understand binary values as true/false, the RefTraitInfo was designed with
  * the analysis of disease in mind. Consequently, the internal vocabulary reflects the affected/unaffected concept.
  *
  * You'll need to set the affected/unaffected code (as indicated by the source data), both in numeric and string
  * form. This is accomplish with the functions set_string_affected_code(), set_string_unaffected_code(), 
  * set_numeric_affected_code(), set_numeric_unaffected_code().
  *
  * You can also set a numeric threshold value to determine affectedness. Values greater than the indicated threshold
  * will be interpreted as \b affected, while values less-than-or-equal-to the threshold will be interpreted as
  * \b unaffected.
  *
  */
class RefTraitInfo
{
public:
  /// Identifies trait category
  enum trait_t   { unknown_trait,    /*!< Trait is of unknown type.                     */
                   invalid_trait,    /*!< Trait is an invalid type.                     */
                   continuous_trait, /*!< Trait is continuous (-INF to +INF, eg. 0.55). */
                   binary_trait,     /*!< Trait is binary (0/1, true/false, etc.).      */
                   discrete_trait,   /*!< Trait is discrete (-INF to +INF, eg. 3).      */
                   categorical_trait /*!< Trait is categorical ("red", "blue", etc.).   */
                 }; 

  /// Identifies the intended use of the trait
  enum trait_use { unknown_use,     /*!< Trait's use is unknown.                       */
                   trait_variate,   /*!< Trait is intended for use as a primary trait. */
                   trait_covariate  /*!< Trait is intended for use as a covariate.     */
                 };

  /// A vector of string values; each value is a valid categorical value. Its index is 
  /// the actual value that is stored for an individual.
  typedef std::vector<std::string> CategoryVector;

  /// @name Constructors
  //@{

    ///
    /// Constructor.
    RefTraitInfo();

    ///
    /// Constructor
    /// \param trait_name Name of the trait
    /// \param trait_type Type of the trait
    /// \param usage Intended use of the trait
    explicit RefTraitInfo(string    trait_name, 
                          trait_t   trait_type = continuous_trait,
                          trait_use usage      = unknown_use);

  //@}

  /// @name Basic attributes
  //@{

    ///
    /// Sets the trait name.
    /// \param n The name of the trait
    void set_name(const string &n);

    ///
    /// Sets the trait alias name.
    /// \param n The alias name of the trait
    void set_alias_name(const string &n);

    ///
    /// Returns the trait name.
    const string & name() const;

    ///
    /// Returns the trait alias name.
    const string & alias_name() const;

    ///
    /// Sets the trait type.
    /// \param typ The type of the trait
    void set_type(trait_t typ);

    ///
    /// Returns the trait type.
    trait_t type() const;

    ///
    /// Sets the trait usage.
    /// \param tru The usage of the trait
    void set_usage(trait_use tru);

    ///
    /// Returns the trait usage.
    trait_use usage() const;
  
  //@}

  /// @name Missingness
  //@{

    ///
    /// Sets the trait's missing code as a numeric value.
    /// \param nm The missing code
    void set_numeric_missing_code(double nm);

    ///
    /// Returns the trait's missing code.
    double numeric_missing_code() const;

    ///
    /// Sets the trait's missing code as a string.
    /// \param sm The missing code
    void set_string_missing_code(const string &sm);

    ///
    /// Returns the trait's missing code.
    const string &string_missing_code() const;

  //@}

  /// @name Binary trait functions
  //@{

    ///
    /// Sets the string-type affectedness code for reading in binary traits.
    /// \param ac Affectedness code
    void set_string_affected_code(const string &ac);

    ///
    /// Returns the string-type affectedness code for reading in binary traits.
    const string &string_affected_code() const;

    ///
    /// Sets the string-type unaffectedness code for reading in binary traits.
    /// \param uc Unaffectedness code
    void set_string_unaffected_code(const string &uc);

    ///
    /// Returns the string-type unaffectedness code for reading in binary traits.
    const string &string_unaffected_code() const;

    ///
    /// Sets the numeric-type affectedness code for reading in binary traits.
    /// \param ac Affectedness code
    void set_numeric_affected_code(double ac);

    ///
    /// Returns the numeric-type affectedness code for reading in binary traits.
    double numeric_affected_code() const;

    ///
    /// Sets the numeric-type unaffectedness code for reading in binary traits.
    /// \param uc Unaffectedness code
    void set_numeric_unaffected_code(double uc);

    ///
    /// Returns the numeric-type unaffectedness code for reading in binary traits.
    double numeric_unaffected_code() const;

    ///
    /// Sets the numeric threshold value for reading in binary traits (see binary trait section of detailed description).
    /// \param th Threshold value
    void set_threshold(double th);

    ///
    /// Returns the numeric threshold value for reading in binary traits.
    double threshold() const;

  //@}
  
  /// @name Category functions
  //@{
  
    ///
    /// Returns the categories read in for this trait.
    const CategoryVector & get_categories() const { return my_categories; }

    ///
    /// Returns the categories read in for this trait.
    CategoryVector & get_categories() { return my_categories; }
    
    ///
    /// Adds the named category (if it's not already in the list, and returns the id of that category.
    size_t add_category(const std::string & v) { if(get_id(v) == (size_t)-1) { my_categories.push_back(v); } return get_id(v); }

    ///
    /// Looks up the category name and, if it has been read in, returns the corresponding id.
    /// Otherwise, returns (size_t)-1.
    size_t get_id(const std::string & category) const;  
    
    ///
    /// Returns whether or not the list of categories can be added to.
    bool get_lockout() const { return my_lockout; }

    ///
    /// Sets the lockout status (see get_lockout() ).
    void set_lockout(bool b) { my_lockout = b; }

  //@}

private:
  void init();

  // Basic info:
  string    my_name;
  string    my_alias_name;
  trait_t   my_type;
  trait_use my_use;

  double  my_threshold;

  // Data members for missing value stuff:
  double  my_numeric_missing_code;
  string  my_string_missing_code;

  // Data members for binary trait stuff:
  double  my_numeric_affected_code;
  double  my_numeric_unaffected_code;
  string  my_string_affected_code;
  string  my_string_unaffected_code;

  // Data members for categorical trait stuff:
  CategoryVector my_categories;
  bool           my_lockout;
};

/** \brief Stores information about traits whose type is not addressed in rped
  *
  * In the event that pedigree data has traits that cannot be classified according as
  * phenotypes, traits, covariates, markers, alleles, or trait markers, the RefStringInfo
  * type is available. Although the data read in for any particular RefStringInfo cannot
  * be used directly by an analysis (since the nature of its data is unspecifiable), rped
  * still allows the user to note that the trait does exist.
  *
  * This feature is still used, though, in a program's output. If the program lists values
  * on an individual-by-individual basis, it may be helpful to include unclassifiable traits
  * in such output. For this reason, the RefStringInfo class is made available.
  *
  */
class RefStringInfo
{
public:

  /// @name Constructors
  //@{

    ///
    /// Constructor.
    RefStringInfo();

    ///
    /// Constructor.
    /// \param name The name of the string-type trait
    explicit RefStringInfo(string name);

  //@}

  /// @name Name
  //@{

    ///
    /// Sets the name of the trait.
    /// \param n The new name
    void set_name(const string &n);

    ///
    /// Returns the name of the trait.
    const string &name()  const;
  
  //@}

  /// @name Missingness
  //@{

    ///
    /// Sets the missing code (for reading in the trait value from pedigree data).
    /// \param sm The missing code
    void set_string_missing_code(const string &sm);

    ///
    /// Returns the missing code (for reading in the trait value from pedigree data).
    const string &string_missing_code() const;
  
  //@}

private:
  void init();

  string  my_name;
  string  my_string_missing_code;
};

/** \brief Please see MLOCUS::inheritance_model for more information.
  *
  * Please see mlocus::inheritance_model for more information.
  */
typedef MLOCUS::inheritance_model RefMarkerInfo;

/** \brief Storage for information read from the marker block
  *
  * \par Introduction
  *
  * The PhenotypeReaderInfo specifies options for reading in marker/allele information from pedigree data.
  * It stores information specified in the \c marker block of a parameter file (see user manual for more information).
  *
  * \par Source data syntax
  *
  * There are a few syntax-related functions available for specifying how an allele is to be read in from pedigree 
  * data. Assuming marker data takes the form "[allele_name][delimeter][allele_name]", you can use the
  * syntax-related functions to specify the nature of that delimeter, as well as the allele missing code.
  *
  * The syntax functions (delimeter and missing code) apply generally to all markers read from the source data. If, 
  * however, a specific marker entry in the pedigree block specifies a delimeter or missing code attribute, that
  * attribute \b overrides such attributes specified in the marker block.
  *
  * \par Allele frequency adjustment
  *
  * For programs that use/estimate allele frequencies, it may be analytically useful to adjust the allele frequency
  * calculated from the source data. That is, a given set of source data may in fact have zero occurences of allele
  * 'A'; the population estimate, however, may be non-zero. Furthermore, in order to successfully arrive at a
  * population estimate, the initial estimate may need to be non-zero (for mathematical reasons).
  *
  * Consequently, it is possible to specify allele adjustment options in the marker block of a parameter file.
  * This is accomplished with the set_allele_adjustment() function, as well as its helper functions 
  * (set_min_allele_freq(), set_max_allele_freq()). Please see PhenotypeReaderInfo::allele_adj for more information
  * on the type of adjustments available.
  */
class PhenotypeReaderInfo
{
public:

  /** \brief Allele frequency adjustment option
    *
    * Allele frequencies given by the source data can be overridden in a number of ways.
    * This enum describes those override methods.
    */
  enum allele_adj { none,   /*!< No frequency adjustment                                                    */
                    equal,  /*!< All allele frequencies set to be equal                                     */
                    comp,   /*!< For any frequency i, i = 1.0 - i (the complement of i)                     */
                    min,    /*!< For any frequency i such that i < minimum frequency, i = minimum frequency */
                    max,    /*!< For any frequency i such that i > maximum frequency, i = maximum frequency */
                    min_max /*!< Equivalent to min \b and max options together                              */
                  };
  
  /// @name Constructor
  //@{

    ///
    /// Constructor.
    PhenotypeReaderInfo();

  //@}

  /// @name Source data syntax
  //@{

    ///
    /// When marker data is read in, it may take the form "allele_name/allele_name". The delimiter
    /// (in this case, the forward slash) may be any single character. That character is set in the
    /// PhenotypeReaderInfo with this function.
    /// \param s The allele delimeter character
    void set_allele_delimiter(char s);

    ///
    /// Returns the allele delimeter character.
    char get_allele_delimiter() const;

    ///
    /// The missing code for an allele is stored as a string. This function sets the missing code.
    /// \param s The missing code for the allele
    void set_allele_missing(const string s);

    ///
    /// Returns the allele's missing code.
    const string get_allele_missing() const;

  //@}

  /// @name Allele frequency adjustment
  //@{

    ///
    /// Sets the allele frequency adjustment option.
    /// \param a The allele frequency adjustment option
    void set_allele_adjustment(allele_adj a = none);

    ///
    /// Returns the allele frequency adjustment option.
    allele_adj get_allele_adjustment() const;

    ///
    /// Sets the minimum allele frequency (in conjunction with frequency adjustment options PhenotypeReaderInfo::min 
    /// and PhenotypeReaderInfo::min_max)
    /// \param m The minimum allele frequency
    void set_min_allele_freq(double m = 0.0);

    ///
    /// Returns the minimum allele frequency.
    double get_min_allele_freq() const;

    ///
    /// Sets the maximum allele frequency (in conjunction with frequency adjustment options PhenotypeReaderInfo::max
    /// and PhenotypeReaderInfo::min_max)
    /// \param m The maximum allele frequency
    void set_max_allele_freq(double m = 1.0);

    ///
    /// Returns the maximum allele frequency.
    double get_max_allele_freq() const;

  //@}

    void set_covariate_moi   (const string cm = "");
    void set_covariate_allele(const string ca = "");
    void set_allow_hemizygote(bool b = false);

    const string get_covariate_moi()    const;
    const string get_covariate_allele() const;
    bool         get_allow_hemizygote() const;

private:

  allele_adj my_allele_adjustment;
  double     my_minimum_allele_freq;
  double     my_maximum_allele_freq;

  char       my_allele_delimiter;
  string     my_allele_missing;

  string     my_cov_moi;
  string     my_cov_allele;
  bool       my_allow_hemizygote;
};

/** \brief Storage for trait names, missing codes and other meta-information for RPED::RefMultiPedigree::multipedigree_type
  *
  * \par Introduction
  *
  * Recall that mped's derived classes are templatized; to turn them into usable objects we need to
  * create a template specialization of those classes. In rped, the mped derived classes are
  * specialized with the RefMPedInfo (multipedigree info) and RefPedInfo (pedigree info) classes.
  *
  * RefMPedInfo's purpose is to keep track of source data options (such as formatting, delimeters, etc.).
  * In addition, it also makes certain basic marker data statistics available.
  *
  * To use RefMPedInfo, the first thing you need to do is set up the pedigree structure configuration options.
  * These options are specified via the "Basic data importation options" section of this class.
  *
  * Then, you should scan through the pedigree block of the parameter file, and create traits/markers as necessary.
  * This is done via the three "Adding / editing / removing [object]" sections.
  *
  * \par Basic data importation options
  *
  * There are four options that you need to set for reading in pedigree data:
  * 
  * o Individual missing code - see RefMPedInfo::set_individual_missing_code
  * 
  * o Sex code for male - see set_sex_code_male()
  * 
  * o Sex code for female - see set_sex_code_female()
  * 
  * o Sex code for unknown sex - see set_sex_code_unknown()
  * 
  * \par Traits
  *
  * RefMPedInfo tracks trait information, such as trait name, type (binary, continuous, etc.), and so on. The
  * individual trait values themselves are stored in the RefTraitInfo object (see SAGE::RPED::RefTraitInfo).
  *
  * Before reading in actual pedigree data, you should create entries for all traits that you plan to store.
  * You can do this with the add_trait(), add_continuous_trait(), add_binary_trait(), add_continuous_covariate(),
  * and add_binary_covariate() functions.
  *
  * In scanning the pedigree block of the parameter file, you should create a trait entry for all indicated
  * traits/covariates, including the various trait-specific options. For instance, if the user has the following 
  * pedigree block:
  *
  * \code
  * pedigree
  * {
  *   ...
  *   covariate = ENV_01, missing = -1
  *   ...
  * }
  * \endcode
  *
  * Your code should do something like this:
  *
  * \code
  * SAGE::RPED::RefMultiPedigree my_RMP();
  * int id = my_RMP.info().add_trait("ENV_01", SAGE::RPED::RefTraitInfo::continuous_trait, SAGE::RPED::RefTraitInfo::trait_covariate);
  *
  * my_RMP.info().trait_info(id).set_numeric_missing_code(-1);
  * \endcode
  *
  * \par Markers
  *
  * You should also create marker entries for all markers specified in the parameter file. This is done via
  * the functions in the "Adding / editing / removing markers" section.
  *
  * For instance, if the pedigree block has a marker declared like this:
  *
  * \code
  * pedigree
  * {
  *   ...
  *   marker = LOC_02, delimiter = /, name = LOC_02, missing = -1
  *   ...
  * }
  * \endcode
  *
  * Your code should do something like this:
  *
  * \code
  * SAGE::RPED::RefMultiPedigree my_RMP();
  * int id = my_RMP.info().add_marker("LOC_02");
  *
  * my_RMP.info().marker_info(id)...
  * \endcode
  *
  * Also, please note that there is a group of functions designed to read in marker block parameters as well. 
  * If you encounter a marker block in parsing the parameter file, you can pass it to the 
  * RefMPedInfo::read_pheno_reader_info function. Although individual marker listings (in the pedigree block)
  * override the options in the marker block, any markers lacking such specific options will use those specified
  * in the marker block.
  *
  * \par String-type fields
  *
  * The user can also specify (in the pedigree block) fields whose type is not addressed by SAGE. In this
  * case, the data is indeed read in for each user, even though it cannot be used in an analysis. It can,
  * however, still be used in the program's output (in a table of individual field values, for instance).
  *
  * In scanning the pedigree block of the parameter file, you should create an entry for any string-type
  * field the user specifies. For instance, if the pedigree block looks like this:
  *
  * \code
  * pedigree
  * {
  *   ...
  *   string = FOO
  *   ...
  * }
  * \endcode
  *
  * Your code should do something like this:
  *
  * \code
  * SAGE::RPED::RefMultiPedigree my_RMP();
  *
  * my_RMP.info().add_string_field("FOO");
  * \endcode
  *
  */
class RefMPedInfo
{
public:
 
  /// @name Constructor
  //@{

    ///
    /// Constructor.
    RefMPedInfo();

  //@}

  /// @name Basic data importation options
  //@{

    ///
    /// Sets the missing code for an individual (from the pedigree data).
    /// \param imc Individual missing code
    void set_individual_missing_code(const string &imc);

    ///
    /// Returns the missing code for an individual (from the pedigree data).
    const string &individual_missing_code() const;

    ///
    /// Sets the sex code for male (from the pedigree data).
    /// \param c Sex code for male
    void set_sex_code_male(const string &c);

    ///
    /// Returns the sex code for male (from the pedigree data).
    const string &sex_code_male() const;

    ///
    /// Sets the sex code for female (from the pedigree data).
    /// \param c Sex code for female
    void set_sex_code_female(const string &c);

    ///
    /// Returns the sex code for female (from the pedigree data).
    const string &sex_code_female() const;

    ///
    /// Sets the sex code for unknown sex (from the pedigree data).
    /// \param c Sex code for unknown sex
    void set_sex_code_unknown(const string &c);

    ///
    /// Returns the sex code for unknown sex (from the pedigree data).
    const string &sex_code_unknown() const;

  //@}

  /// @name Adding / editing / removing traits
  //@{

    ///
    /// Adds an entry for a new trait.
    /// \param trait_name The name of the trait
    /// \param trait_type The type of the trait (see RefTraitInfo::trait_t)
    /// \param use The usage purpose of the trait (see RefTraitInfo::trait_use)
    /// \returns The id number of the newly created trait
    size_t add_trait(const string&                 trait_name, 
                           RefTraitInfo::trait_t   trait_type = RefTraitInfo::continuous_trait,
                           RefTraitInfo::trait_use use        = RefTraitInfo::unknown_use);

    ///
    /// Adds an entry for a new continuous trait.
    /// \param trait_name The name of the trait
    /// \param use The usage purpose of the trait (see RefTraitInfo::trait_use)
    /// \returns The id number of the newly created trait
    size_t add_continuous_trait(const string&                 trait_name, 
                                      RefTraitInfo::trait_use use = RefTraitInfo::unknown_use);

    ///
    /// Adds an entry for a new binary trait.
    /// \param trait_name The name of the trait
    /// \param use The usage purpose of the trait (see RefTraitInfo::trait_use)
    /// \returns The id number of the newly created trait
    size_t add_binary_trait(const string &trait_name,
                   RefTraitInfo::trait_use use = RefTraitInfo::unknown_use);

    ///
    /// Adds an entry for a new continuous covariate.
    /// \param covariate_name The name of the trait
    /// \returns The id number of the newly created trait
    size_t add_continuous_covariate(const string &covariate_name);

    ///
    /// Adds an entry for a new binary covariate.
    /// \param covariate_name The name of the trait
    /// \returns The id number of the newly created trait
    size_t add_binary_covariate(const string &covariate_name);

    ///
    /// Returns the non-const RefTraitInfo instance associated with the indicated trait.
    /// \param t The trait id
    RefTraitInfo & trait_info(size_t t);

    ///
    /// Removes the trait entry for the most recently added trait.
    void remove_last_trait();

    ///
    /// Removes the trait entry for the indicated trait.
    /// \param t_id The id number of the trait.
    void remove_trait_info(size_t t_id);

  //@}

  /// @name Adding / editing / removing markers
  //@{

    ///
    /// Adds an entry for the indicated marker.
    /// \param marker_name The name of the new marker
    /// \returns The id number of the newly created marker
    size_t add_marker(const string &marker_name);

    ///
    /// Adds an entry for the indicated marker.
    /// \param marker_name The name of the new marker
    /// \param model The MLOCUS::inheritance_model associated with this marker
    /// \returns The id number of the newly created marker
    size_t add_marker(const string &marker_name, const MLOCUS::inheritance_model& model);
  
    ///
    /// Returns the non-const MLOCUS::inheritance_model_map instance associated with this RefMPedinfo.
    MLOCUS::inheritance_model_map & markers();

    ///
    /// Returns the non-const RefMarkerInfo instance associated with the indicated marker.
    /// \param m The id number of the marker in question
    RefMarkerInfo & marker_info(size_t m);

    ///
    /// Returns the non-const RefMarkerInfo instance associated with the indicated marker.
    /// \param m The name of the marker in question
    RefMarkerInfo & marker_info(const string& m);

    ///
    /// Removes the entry for the indicated marker.
    /// \param m_id The id number of the marker
    void remove_marker_info(size_t m_id);

    ///
    /// Removes the entry for the indicated marker.
    /// \param m_name The name of the marker
    void remove_marker_info(string m_name);

  //@}

  /// @name Marker block functions
  //@{

    ///
    /// Reads in marker block options from the indicated parameter set.
    /// \param params The data object containing the marker block configuration
    /// \param errors The error stream to which error messages will be directed
    void read_pheno_reader_info(const LSFBase* params, cerrorstream &errors);

    ///
    /// Returns the non-const PhenotypeReaderInfo instance associated with this RefMPedInfo.
    PhenotypeReaderInfo& get_pheno_reader_info();

    ///
    /// Returns the const PhenotypeReaderInfo instance associated with this RefMPedInfo.
    const PhenotypeReaderInfo& get_pheno_reader_info() const;

  //@}  

  /// @name Adding / editing / removing string-type fields (unanalysable)
  //@{

    ///
    /// Adds an entry for the indicated string-type field.
    /// \param name The name of the string-type field
    size_t add_string_field(const string &name);

    ///
    /// Returns the non-const RefStringInfo instance associated with the indicated string-type field.
    /// \param s The id number of the string-type field
    RefStringInfo &string_info(size_t s);

  //@}

  /// @name Querying trait status
  //@{

    ///
    /// Returns the const RefTraitInfo instance associated with the indicated trait.
    /// \param t The trait id
    const RefTraitInfo &trait_info(size_t t) const;

    ///
    /// Returns the number of traits created.
    size_t trait_count() const;

    ///
    /// Returns \c true if the trait exists, \c false if it does not exist.
    /// \param trait_name The name of the trait in question
    bool trait_exists(const string& trait_name) const;

    ///
    /// Returns the type of the trait in question.
    /// \param trait_name The name of the trait in question
    RefTraitInfo::trait_t get_trait_type(const string& trait_name) const;

    ///
    /// Returns the id number of the indicated trait.
    /// \param name The name of the trait in question
    /// \retval -1 Trait was not found
    /// \retval >=0 The id number of the trait
    size_t trait_find(const std::string &name) const;

  //@}

  /// @name Querying marker status
  //@{

    ///
    /// Returns the number of markers created.
    size_t marker_count() const;

    ///
    /// Returns the const MLOCUS::inheritance_model_map associated with this RefMPedInfo.
    const MLOCUS::inheritance_model_map& markers() const;

    ///
    /// Returns the const RefMarkerInfo instance associated with the indicated marker.
    /// \param m The id number of the marker
    const RefMarkerInfo &marker_info(size_t m) const;

    ///
    /// Returns \c true if the marker exists, \c false if it does not exist.
    /// \param marker_name The name of the marker in question
    bool marker_exists(const string& marker_name) const;

    ///
    /// Returns the id number of the indicated trait.
    /// \param name The name of the trait in question
    /// \retval -1 Trait was not found
    /// \retval >=0 The id number of the trait
    size_t marker_find(const std::string &name) const;

  //@}

  /// @name Querying string-type field status
  //@{

    ///
    /// Returns the number of stirng-type fields in this RefMPedInfo.
    size_t string_count() const;

    ///
    /// Returns the const RefStringInfo instance associated with the indicated string-type field.
    /// \param s The name of the string-type field in question
    const RefStringInfo & string_info(size_t s) const;

    ///
    /// Returns the id number of the string-type field.
    /// \retval -1 Field was not found
    /// \retval >=0 The id number of the field
    size_t string_find(const std::string &name) const;

  //@}

private:
  string my_individual_missing_code;
  string my_sex_code_male;
  string my_sex_code_female;
  string my_sex_code_unknown;

  typedef vector<RefTraitInfo>     trait_vector;
  typedef vector<RefStringInfo>    string_vector;  
  typedef MLOCUS::inheritance_model_map    marker_map;

  trait_vector         my_trait_info;
  string_vector        my_string_info;
  marker_map           my_marker_info;
  PhenotypeReaderInfo  my_pheno_reader_info;  
};

/** \brief Storage individual trait values
  *
  * \par Introduction
  *
  * The object is designed to be attached to a SAGE::MPED::pedigree. It is attached to the object through
  * RPED::RefMultiPedigree as a template parameter.
  *
  * Once you have created an instance of this object (for association with a SAGE::RPED::RefPedigree), you
  * should do the following:
  *
  * \c 1 Build it (see Initial build)
  *
  * \c 2 Set the correct number of markers/traits/string-fields (see Adding/removing trait/string/phenotype fields)
  *
  * \c 3 Add the individual field values (see Setting individual trait/string/phenotype values)
  *
  * Later, when you want to fetch individual field values, you should use the functions from the 
  * "Fetching total number of trait/string/marker/member values" section.
  */
class RefPedInfo
{
public:

  /// @name Constructor
  //@{

    ///
    /// Constructor.
    RefPedInfo();

  //@}

  /// @name Initial build
  //@{

    ///
    /// Once you have added all the members, you should invoke this function to finalize the internal data storage.
    void build(MPED::pedigree_base &ped);

  //@}

  /// @name Adding & removing trait / string / phenotype fields
  //@{

    ///
    /// Sets the number of trait fields.
    /// \param traits The number of trait fields
    void resize_traits(size_t traits);

    ///
    /// Sets the number of marker fields.
    /// \param markers The number of markers
    /// \param mped_info The const reference of RefMPedInfo instance
    void resize_markers(size_t markers, const RefMPedInfo& mped_info);

    ///
    /// Sets the number of string-type fields.
    /// \param s The number of s
    void resize_strings(size_t s);

    ///
    /// Removes the entry for the indicated trait field.
    /// \param t_id The id number of the trait field in question
    void remove_trait(size_t t_id);

    ///
    /// Removes the entry for the indicated marker field.
    /// \param m_id The id number of the marker field in question
    void remove_marker(size_t m_id);

  //@}

  /// @name Setting individual trait / string / phenotype values
  //@{

    ///
    /// Sets the trait value for an individual (according to the source data).
    /// \param i The index of the individual
    /// \param t The number of the trait field
    /// \param d The trait value
    /// \retval true Value set successfully
    /// \retval false Value \b not set successfully (either \c i or \c t is out of range)
    bool set_trait(size_t i, size_t t, double d);

    ///
    /// Sets the trait value for an individual (according to the source data).
    ///
    /// Please note that since the trait value being set is in string form, it must
    /// be compared against configuration parameters for that trait (such as missing code
    /// and numeric threshold).
    /// \param i The index of the individual
    /// \param t The number of the trait field
    /// \param value The trait value (in string form)
    /// \param trait_info The RefTraitInfo instance associated with the trait being added
    /// \retval 0 Trait set ok
    /// \retval 1 Trait set ok, but missing
    /// \retval 2 Bad trait value, assumed missing
    /// \retval 3 Invalid invididual or trait id
    /// \retval 4 Trait not set due to trait settings
    int set_trait(size_t i, 
                  size_t t, 
                  const std::string  & value, 
                  RefTraitInfo & trait_info);

    ///
    /// Sets the phenotype (marker) value for an individual (according to the source data).
    /// \param i The index of the individual
    /// \param m The number of the marker field
    /// \param p The phenotype value (see RefMarkerInfo for more information on how to obtain phenotype values)
    /// \retval true Value set successfully
    /// \retval false Value \b not set successfully (either \c i or \c t is out of range)
    bool set_phenotype(size_t i, 
                       size_t m, 
                       uint   p);

    ///
    /// Sets the phenotype (marker) value for an individual (according to the source data).
    /// \param i The index of the individual
    /// \param m The number of the marker field
    /// \param allele1 The first allele
    /// \param allele2 The first allele
    /// \param marker The RefMarkerInfo instance associated with this marker
    /// \retval 0 Marker set ok
    /// \retval 1 Marker set ok, but missing
    /// \retval 2 Bad marker value, assumed missing
    /// \retval 3 Invalid invididual or marker id    
    int set_phenotype(size_t i, 
                      size_t m, 
                      const std::string & allele1,
                      const std::string & allele2,
                      RefMarkerInfo & marker);

    ///
    /// For a marker that is X- or Y-linked, use this function to set the phenotypic value. It will
    /// be checked against the X/Y linkage options in the RefMarkerInfo.
    /// \param ind_id The index of the individual
    /// \param ind_sex The sex of the individual
    /// \param marker_id The number of the marker field
    /// \param allele1 The first allele
    /// \param allele2 The first allele
    /// \param marker The RefMarkerInfo instance associated with this marker
    /// \retval 0 Marker set ok
    /// \retval 1 Marker set ok, but missing
    /// \retval 2 Bad marker value, assumed missing
    /// \retval 3 Invalid invididual or marker id    
    int set_phenotype(size_t ind_id, 
                      MPED::SexCode ind_sex, 
                      size_t marker_id, 
                      const std::string & allele1,
                      const std::string & allele2,
                      RefMarkerInfo & marker);

    ///
    /// Sets the string-type field value for an individual.
    /// \param i The index of the individual
    /// \param s The id number of the string-type field
    /// \param val The string-type value
    /// \retval true Value was set successfully
    /// \retval false Value was \b not set successfully
    bool set_string(size_t i, 
                    size_t s, 
                    const string& val);

  //@}

  /// @name Re-organizing members
  //@{

    ///
    /// Swaps all entries for the two indicated members.
    /// \param m1 The index of the first member
    /// \param m2 The index of the second member
    void swap_members(size_t m1, size_t m2);

  //@}

  /// @name Fetching total number of trait / string / marker / member values
  //@{

    ///
    /// Returns the number of trait fields present in this object.
    size_t trait_count() const;

    ///
    /// Returns the number of string-type fields present in this object.
    size_t string_count() const;

    ///
    /// Returns the number of marker fields present in this object.
    size_t marker_count() const;

    ///
    /// Returns the number of members present in this object.
    size_t member_count() const;

  //@}

  /// @name Fetching individual trait / string / phenotype values
  //@{

    ///
    /// Returns an individual's trait value.
    /// \param i The index of the individual
    /// \param t The id number of the trait field
    double trait(size_t i, size_t t) const;

    ///
    /// Returns \c true if the trait value is missing, \c false if it is present.
    /// \param i The index of the individual
    /// \param t The id number of the trait field
    bool trait_missing(size_t i, size_t t) const;

    ///
    /// Returns an individual's phenotypic value.
    /// \param i The index of the individual
    /// \param m The id number of the marker field
    uint phenotype(size_t i, size_t m) const;

    ///
    /// Returns \c true if the phenotypic value is missing, \c false if it is present.
    /// \param i The index of the individual
    /// \param m The id number of the marker field
    /// \param mi The RefMarkerInfo instance associated with this marker
    bool phenotype_missing(size_t i, size_t m, const RefMarkerInfo& mi) const;

    ///
    /// Returns an individual's string-type value.
    /// \param i The index of the individual
    /// \param s The id number of the string-type field
    string get_string(size_t i, size_t s) const;

  //@}

private:
    int set_autosomal_phenotype(size_t i,  size_t m, const std::string & allele1,
                                const std::string & allele2, RefMarkerInfo & marker);
    int set_x_linked_phenotype(size_t i,
                               MPED::SexCode ind_sex, 
                               size_t m, const std::string & allele1,
                               const std::string & allele2, RefMarkerInfo & marker);
    int set_y_linked_phenotype(size_t i,
                               MPED::SexCode ind_sex, 
                               size_t m, const std::string & allele1,
                               const std::string & allele2, RefMarkerInfo & marker);
                      
  typedef vector<double> dvector;
  typedef vector<string> svector;
  typedef vector<uint>   uvector;

  vector<dvector> my_traits;
  vector<svector> my_strings;
  vector<uvector> my_markers;
  size_t my_member_count;
  size_t my_marker_count;
  size_t my_trait_count;
  size_t my_string_count;
};

/** \brief MPED::multipedigree instantiated with RefPedInfo and RefMPedInfo as template parameters
  *
  * MPED::multipedigree specialized on RefPedInfo and RefMPedInfo as template parameters.
  */
typedef MPED::multipedigree <MPED::no_info, MPED::no_info, MPED::no_info, RefPedInfo, RefMPedInfo> RefMultiPedigree;

/** @name Typedef-ed from RefMultiPedigree
  *
  * All of the following typedefs come from RefMultiPedigree.
  */
//@{

  typedef RefMultiPedigree::pedigree_type     RefPedigree;
  typedef RefMultiPedigree::subpedigree_type  RefSubpedigree;
  typedef RefMultiPedigree::family_type       RefFamily;
  typedef RefMultiPedigree::member_type       RefMember;

//@}

/** \brief MPED::multipedigree instantiated with RefPedInfo and RefMPedInfo as template parameters
  *
  * MPED::multipedigree specialized on RefPedInfo and RefMPedInfo as template parameters.
  */
typedef MPED::multipedigree <MPED::no_info, MPED::no_info, MPED::no_info, RefPedInfo, RefMPedInfo> MultiPedigree;

// Old-style typedefs
typedef MultiPedigree::member_const_pointer        member_const_pointer;  

/** @name Typedef-ed from MultiPedigree.
  *
  * All of the following typedefs come from MultiPedigree.
  */
//@{

  typedef MultiPedigree::pedigree_type     Pedigree;
  typedef MultiPedigree::subpedigree_type  Subpedigree;
  typedef MultiPedigree::family_type       Family;
  typedef MultiPedigree::member_type       Member;

  typedef MultiPedigree::member_pointer        MemberPointer;
  typedef MultiPedigree::family_pointer        FamilyPointer;
  typedef MultiPedigree::subpedigree_pointer   SubpedigreePointer;
  typedef MultiPedigree::pedigree_pointer      PedigreePointer;
  typedef MultiPedigree::multipedigree_pointer MultipedigreePointer;

  typedef MultiPedigree::member_const_pointer        MemberConstPointer;
  typedef MultiPedigree::family_const_pointer        FamilyConstPointer;
  typedef MultiPedigree::subpedigree_const_pointer   SubpedigreeConstPointer;
  typedef MultiPedigree::pedigree_const_pointer      PedigreeConstPointer;
  typedef MultiPedigree::multipedigree_const_pointer MultipedigreeConstPointer;

  typedef MultiPedigree::family_iterator        FamilyIterator;
  typedef MultiPedigree::mate_iterator          MateIterator;
  typedef MultiPedigree::member_iterator        MemberIterator;
  typedef MultiPedigree::offspring_iterator     OffspringIterator;
  typedef MultiPedigree::parent_iterator        ParentIterator;
  typedef MultiPedigree::pedigree_iterator      PedigreeIterator;
  typedef MultiPedigree::progeny_iterator       ProgenyIterator;
  typedef MultiPedigree::sibling_iterator       SiblingIterator;
  typedef MultiPedigree::subpedigree_iterator   SubpedigreeIterator;

  typedef MultiPedigree::family_const_iterator        FamilyConstIterator;
  typedef MultiPedigree::mate_const_iterator          MateConstIterator;
  typedef MultiPedigree::member_const_iterator        MemberConstIterator;
  typedef MultiPedigree::offspring_const_iterator     OffspringConstIterator;
  typedef MultiPedigree::parent_const_iterator        ParentConstIterator;
  typedef MultiPedigree::pedigree_const_iterator      PedigreeConstIterator;
  typedef MultiPedigree::progeny_const_iterator       ProgenyConstIterator;
  typedef MultiPedigree::sibling_const_iterator       SiblingConstIterator;
  typedef MultiPedigree::subpedigree_const_iterator   SubpedigreeConstIterator;

//@}

// Sort a pedigree into topological order based on lineage
void OldPedigreeSort(RefPedigree &p);

void PedigreeSort(RefPedigree &p);

MemberConstPointer ancestor_error(std::set<size_t>& ancestors, const MemberConstPointer ind);

bool CheckPedigreeSort(RefPedigree &p);

} // End namespace RPED
} // End namespace SAGE

#include "rped/rped.ipp"

#endif
