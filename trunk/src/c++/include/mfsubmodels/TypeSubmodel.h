#ifndef MFSUBMODELS_TYPESUBMODEL_H
#define MFSUBMODELS_TYPESUBMODEL_H

#include <math.h>
#include "maxfunapi/maxfunapi.h"
#include "error/internal_error.h"
#include "mfsubmodels/TypeConstraint.h"
#include "globals/SAGEConstants.h"

namespace SAGE        {
namespace MFSUBMODELS {

const int DUMP_PRECISION = 12;
const int NUM_OF_TYPES   = 3;

/** \brief Base class for all genotype specific sub_models.
  * Please note that this class can never be directly instantiated, as it is pure virtual.
  *
  * \par Purpose
  *
  * \par Parsing logic
  *
  * Although parsing is technically handled by the TypeParser object, much of the underlying logic
  * is handled in this object. In particular, the isInitialSetupConsistent() function takes care
  * of examining the user's initial settings (from the type_foo sub-block) and figuring out what
  * to do with them. Please see the documentation for that function for more information.
  *
  * \par Special options for two genotypes
  *
  * If the ConstraintType is set to TWO (not TWO_DOM or TWO_REC, simply TWO), then
  * the submodel will have to assume that the TWO option is implicitly dominant
  * or recessive. This is because it is meaningless to specify two genotypic foo-types
  * without accounting for the third genotype (in a diallelic model).
  * 
  * Therefore, there exists the functions isTwoDomByDefault() and setTwoDomByDefault().
  *
  * If the ConstraintType is in fact set to TWO, then the submodel will query the
  * isTwoDomByDefault() function to see how it should interpret the TWO option. You can
  * change such default behavior with the setTwoDomByDefault() function.
  *
  */
class TypeSpecificSubmodel : public MAXFUN::Submodel
{
  public:

  /// @name Friend classes
  //@{
  
    friend class TypeSpecificParser;

  //@}
  
  /// @name Enumerations
  //@{
  
    ///
    /// Indicates which category of content this submodel will be describing.
    enum CategoryEnum
    {
      MEAN           = 0, /// For submodels describing means
      VARIANCE       = 1, /// For submodels describing variances
      SUSCEPTIBILITY = 2  /// For submodels describing susceptibilities
    };

  //@}

  /// @name Constructors
  //@{

    ///
    /// This function is protected so it can only be used by a derived type.
    /// \param c The category of content this submodel will be describing
    /// \param errors The errorstream to use
    TypeSpecificSubmodel(CategoryEnum c, cerrorstream& errors = sage_cerr);

    ///
    /// This function is protected so it can only be used by a derived type.
    TypeSpecificSubmodel(const TypeSpecificSubmodel& other);

    ///
    /// This function is protected so it can only be used by a derived type.
    TypeSpecificSubmodel& operator=(const TypeSpecificSubmodel& other);

  //@}

  /// @name Category-specific functions
  //@{
  
    ///
    /// Returns the category of this object.
    CategoryEnum getCategory() const;

    ///
    /// Returns the brief name of the type of data stored by class ("mean", "var", "suscept").
    string getBriefName() const;

    ///
    /// Returns the name of the type of data stored by class ("mean", "variance", "susceptibility")
    string getSingularName() const;
    
    ///
    /// Returns the name of the type of data stored by class ("means", "variances", "susceptibilities")
    string getPluralName() const;
  
    ///
    /// Returns \c true if the given ConstraintTypeEnum is allowed, \c false otherwise.
    /// \param e The ConstraintTypeEnum to evaluate
    bool isConstraintTypeAllowed(TypeConstraint::OptionEnum e) const;

    ///
    /// If your analysis uses both MEAN and VARIANCE submodels, you must make sure that the settings assigned
    /// in both submodels are compatible with each other. In particular, you need to make sure that
    /// the number of variances does not exceed the number of means. This function carries out such an evaluation.
    /// \param mean_submodel The mean submodel
    /// \param var_submodel The variance submodel
    /// \retval 0 The submodels are compatible
    /// \retval 1 The submodels are incompatible because the number of variances exceeds the number of means
    static int areMeanAndVarianceCompatible(const TypeSpecificSubmodel & mean_submodel, const TypeSpecificSubmodel & var_submodel);

  //@}

  /// @name Querying configuration
  //@{
  
    ///
    /// Returns the TypeConstraintEnum for this object.
    const TypeConstraint & getTypeConstraint() const;

    ///
    /// Returns a string describing both the genotype constraint and the type of data stored by the class.
    /// For instance: "three means, decreasing"
    string convertConfigurationToString() const;  

    ///
    /// Returns the current value corresponding to the given index.
    /// \param gt The index in question (0, 1, or 2)
    double getValue(size_t gt) const;

    ///
    /// Compares the TypeConstraint against the initial values for the types, and indicates if
    /// there are any conflicts.
    ///
    /// First of all, if the cateogry of this submodel is VARIANCE, and any given value is less than zero,
    /// an error will be produced.
    ///
    /// To better understand this function, here is a comprehensive list of all the situations examined:
    ///
    /// \par Any OPTION
    ///
    /// If no values given, then parsing succeeds.
    ///
    /// \par OPTION == ONE
    ///
    /// If at least two different values are given, then an error is produced.
    /// 
    /// If all given values are identical, then parsing succeeds.
    ///
    /// \par OPTION == TWO / TWO_DOM / TWO_REC
    ///
    /// If only one value is given, an error is produced.
    ///
    /// \par OPTION == TWO
    ///
    /// If AA and AB are given, and AA == AB, then an error is produced.
    ///
    /// If AA and AB are given, and AA != AB, then the option is interpreted as TWO_REC, and BB = AB (with warning #1)
    ///
    /// If AA and BB are given, and AA == BB, then the option is interpreted as getDefaultTwoTreatment(), and AB = AA.
    ///
    /// If AA and BB are given, and AA != BB, then the option is interpreted as getDefaultTwoTreatment(), and AB = the appropriate value.
    ///
    /// If AB and BB are given, and AB == BB, then an error is produced.
    ///
    /// If AB and BB are given, and AB != BB, then the option is interpreted as TWO_DOM, and AA = AB (with warning)
    ///
    /// If all values are given, and all values are unique, then an error is produced.
    ///
    /// If all values are given, and all values are the same, then the option is interpreted as getDefaultTwoTreatment().
    ///
    /// If all values are given, AA == AB, AA != BB, then the option is interpreted as TWO_DOM.
    /// 
    /// If all values are given, AA == BB, AA != AB, then an error is produced.
    /// 
    /// If all values are given, AB == BB, AA != AB, then the option is interpreted as TWO_REC.
    ///
    /// \par OPTION == TWO_DOM
    ///
    /// If AA and AB are given, and AA != AB, then an error is produced.
    //
    /// If AA and AB are given, and BB is not set, then an error is produced.
    ///
    /// If AA or AB is given, and BB is given, then AA or AB (whichever is missing) is assigned to the appropriate value.
    ///
    /// \par OPTION == TWO_REC
    ///
    /// If AB and BB are set, and AB != BB, then an error is produced.
    ///
    /// If AB and BB are given, and AA is not set, then an error is produced.
    ///
    /// If AB or BB is given, and AA is given, then AB or BB (whichever is missing) is assigned to the appropriate value.
    ///
    /// \par OPTION == THREE
    ///
    /// If only one or two values are given, an error is produced.
    ///
    /// \par OPTION == THREE_ADD
    ///
    /// If three values are given, and they do not meet the THREE_ADD requirements, an error is produced.
    ///
    /// \par OPTION == THREE_INC
    ///
    /// If three values are given, and they do not meet the THREE_INC requirements, an error is produced.
    ///
    /// \par OPTION == THREE_DEC
    ///
    /// If three values are given, and they do not meet the THREE_DEC requirements, an error is produced.
    ///        
    /// \par Evaluating fixedness
    ///
    /// The next step in evaluating the consistency of the initial settings is making sure there are no
    /// conflicting fixedness specifiers.
    ///
    /// \par OPTION = any
    ///
    /// If no fixed statuses are given, then everything's fine.
    ///
    /// \par OPTION = ONE
    ///
    /// If all given fixed statuses = X, then any remaining unassigned fix statuses are set to X.
    ///
    /// If more than one fixed statuses are given and they are different, an error is produced.
    /// 
    /// \par OPTION = TWO (interpreted as TWO_DOM) or TWO_DOM
    ///
    /// If fixed status is given for AA and AB, and said statuses are different, an error is produced.
    ///
    /// If fixed status X is given for AA, then the fixed status of AB is set to X.
    ///
    /// If fixed status X is given for AB, then the fixed status of AA is set to X.
    ///
    /// \par OPTION = TWO (interpreted as TWO_REC) or TWO_REC
    ///
    /// If fixed status is given for AB and BB, and said statuses are different, an error is produced.
    ///
    /// If fixed status X is given for AB, then the fixed status of BB is set to X.
    ///
    /// If fixed status X is given for BB, then the fixed status of AB is set to X.
    ///
    /// \param warning_code This value is set by isInitialSetupConsistent() to indicate any warning conditions encountered.
    /// Possible values include:
    /// 0 - No warning condition
    /// 1 - OPTION=TWO, AA and BB are given but not equal; interpreting as default behavior given by getDefaultTwoTreatment()
    /// 2 - OPTION=TWO, AB and BB are given but not equal; interpreting as TWO_DOM
    ///
    /// \retval 0 Everything is ok
    /// \retval 1 OPTION=ONE, but at least two different values given
    /// \retval 2 OPTION=TWO_DOM, but different values given for AA and AB
    /// \retval 3 OPTION=TWO_DOM, valid value given for dominant type (AA/AB), but no value given for recessive type (BB)
    /// \retval 4 OPTION=TWO_REC, but different values given for AB and BB
    /// \retval 5 OPTION=TWO_REC, valid value given for dominant type (AB/BB), but no value given for recessive type (AA)
    /// \retval 6 OPTION=TWO/TWO_DOM/TWO_REC, but only one value specified
    /// \retval 7 OPTION=THREE_ADD, three values given, but do not meet relational requirement
    /// \retval 8 OPTION=THREE_DEC, three values given, but do not meet relational requirement
    /// \retval 9 OPTION=THREE_INC, three values given, but do not meet relational requirement
    /// \retval 10 OPTION=THREE_ADD/THREE_DEC/THREE_INC, at least one value given but not all three
    /// \retval 11 OPTION=TWO, but three different values given
    /// \retval 12 OPTION=TWO, AA equals AB, but BB not given; not sure how to interpret this
    /// \retval 13 OPTION=TWO, AB equals BB, but AA not given; not sure how to interpret this
    /// \retval 14 OPTION=TWO, AA equals BB, and AB != AA or BB; not sure how to interpret this
    /// \retval 100 OPTION=ONE, at least two different fixedness statuses given
    /// \retval 101 OPTION=TWO (interpreted as TWO_DOM) or TWO_DOM, conflicting fixed status given for AA and AB
    /// \retval 102 OPTION=TWO (interpreted as TWO_REC) or TWO_REC, conflicting fixed status given for AB and BB
    /// \retval 200 Variance sub-block indicated, value less than zero given
    int isInitialSetupConsistent(int & warning_code) const;

  //@}
    
  /// \brief Setting configuration
  //@{
  
    ///
    /// Sets this object's TypeConstraint to the given value.
    /// Please note that this option should be set \b prior to setting the values
    /// for individual genotypes (setTypeValues() ). Those set functions
    /// depend on the TypeConstraint.
    /// \param constraint The new TypeConstraint value
    /// \retval true The type was successfully set
    /// \retval false The type was not successfully set, because the derived submodel does not allow it
    bool setTypeConstraint(const TypeConstraint & constraint);
    
    ///
    /// Sets this object's TypeConstraint to the given string value.
    /// \param t String for the constraint type, such as "ONE", "THREE_ADD", etc. (see TypeConstraint::OptionEnum)
    /// \retval true The type was successfully set
    /// \retval false The type was not successfully set, because the derived submodel does not allow it
    bool setTypeConstraint(const string & t);

    ///
    /// Sets the options for the given genotype.
    /// \param i The index in question (see genotype_index)
    /// \param initial_val The initial value for the genotype
    /// \param fixedness_set Boolean indicating whether or not the user has set the fixedness
    /// \param fixedness If the user did indeed set a fixedness, indicates whether or not it is fixed
    /// \retval true If the new options satisfy the relational requirements given by getTypeConstraint()
    /// \retval false If the the new options conflict with the relational requirements given by getTypeConstraint()
    bool setTypeValues(size_t i, double initial_val, bool fixedness_set, bool fixedness);

  //@}

  /// @name Special option for two (geno)types
  //@{

    ///
    /// If the user specified OPTION=TWO, but did not specify any genotypic values, this function
    /// returns how the option will be interpreted.
    static TypeConstraint::OptionEnum getDefaultTwoTreatment();
    
    ///
    /// Sets the default interpretive behavior for ConstraintType = TWO when no genotypic values
    /// are specified by the user. The parameter in question can \b only be TypeConstraint::TWO_DOM
    /// or TypeConstraint::TWO_REC. 
    /// \param e Indicates default behavior (see above)
    /// \retval true Value set successfully
    /// \retval false Value not set successfully (because \c e was neither TWO_DOM nor TWO_REC)
    static bool setDefaultTwoTreatment(TypeConstraint::OptionEnum e);
    
    ///
    /// If the user specified OPTION=TWO, but did not specify any genotypic values, then the parser
    /// will have figured out how to interpret this option. You can find out what the parser decided
    /// to do with OPTION=TWO by using this function.
    TypeConstraint::OptionEnum treatThisTwoAs() const;

    ///
    /// Sets the interpretive behavoir for this specific submodel. The parser should use this function
    /// once it has determined how to treat OPTION=TWO.
    ///
    /// Note: Why is this function const? Because it is called by isInitialSetupConsistent(), which is
    /// const.
    /// \retval true Value set successfully
    /// \retval false Value not set successfully (because \c e was neither TWO_DOM nor TWO_REC)
    bool setTreatThisTwoAs(TypeConstraint::OptionEnum e) const;

  //@}

  /// @name Blah

    /// Tests the submodel for completeness, which is defined as all
    /// parameters specified as finite values.
    bool is_complete() const;

  //@}

  /// @name Debugging
  //@{
    
    /// \internal
    /// Dumps the contents of this object.
    void dump(std::ostream& out) const;

  //@}

  protected:

  /// @name Required virtual interface
  //@{
  
    virtual int update ();
    
  //@}

  /// \brief Maximization Initialization
  ///
  //@{
      
    ///
    /// This function synchronizes the my_parameters vector with what's in
    /// my_types and my_types fixed.  Necessary for inclusion in maxfun
    /// or other times where that must be synched.
    virtual int finalizeConfiguration();

    ///
    /// Initialize my_parameters if TypeConstraint = ONE
    void initializeOne();

    ///
    /// Initialize my_parameters for TypeConstraint = TWO_DOM
    /// (A dominant over B (\f$t_AB = t_AA\f$) )
    void initializeTwoDom();

    ///
    /// Initialize my_parameters for TypeConstraint = TWO_REC
    /// (A recessive under B (\f$t_AB = t_BB\f$) )
    void initializeTwoRec();

    ///
    /// Initialize my_parameters for TypeConstraint = THREE
    void initializeThree();

    ///
    /// Initialize my_parameters for TypeConstraint = THREE_ADD
    /// (\f$t_AB = (t_AA + t_BB) / 2\f$)
    void initializeThreeAdd();

  //@}

  /// @name Option specific synchronization functions (sync w. Maxfun).
  //@{

    ///
    /// Synchronizes current parameter estimates with Maxfun for ConstraintType = ONE.
    /// \retval 0 Synchronization successful
    /// \retval !=0 Synchronization not successful
    int synchronizeOne();

    ///
    /// Synchronizes current parameter estimates with Maxfun for ConstraintType = TWO.
    /// \retval 0 Synchronization successful
    /// \retval !=0 Synchronization not successful
    int synchronizeTwo();

    ///
    /// Synchronizes current parameter estimates with Maxfun for ConstraintType = THREE.
    /// \retval 0 Synchronization successful
    /// \retval !=0 Synchronization not successful
    int synchronizeThree();

    ///
    /// Synchronizes current parameter estimates with Maxfun for ConstraintType = TWO_DOM.
    /// \retval 0 Synchronization successful
    /// \retval !=0 Synchronization not successful
    int synchronizeTwoDom();

    ///
    /// Synchronizes current parameter estimates with Maxfun for ConstraintType = TWO_REC.
    /// \retval 0 Synchronization successful
    /// \retval !=0 Synchronization not successful
    int synchronizeTwoRec();

    ///
    /// Synchronizes current parameter estimates with Maxfun for ConstraintType = THREE_ADD.
    /// \retval 0 Synchronization successful
    /// \retval !=0 Synchronization not successful
    int synchronizeThreeAdd();

    ///
    /// Synchronizes current parameter estimates with Maxfun for ConstraintType = THREE_DEC.
    /// \retval 0 Synchronization successful
    /// \retval !=0 Synchronization not successful
    int synchronizeThreeDec();

    ///
    /// Synchronizes current parameter estimates with Maxfun for ConstraintType = THREE_INC.
    /// \retval 0 Synchronization successful
    /// \retval !=0 Synchronization not successful
    int synchronizeThreeInc();

  //@}    

  /// Constraint Checking
  //@{

    ///
    /// Returns \c true if the relational constraints of THREE_ADD are met, \c false otherwise.
    bool isThreeAddConstraintMet() const;

    ///
    /// Returns \c true if the relational constraints of THREE_DEC are met, \c false otherwise.
    bool isThreeDecConstraintMet() const;

    ///
    /// Returns \c true if the relational constraints of THREE_INC are met, \c false otherwise.
    bool isThreeIncConstraintMet() const;

  //@}
        

  private:

  /// @name Data members
  //@{
  
    CategoryEnum   my_category;
    TypeConstraint my_type_constraint;

    mutable double my_types     [NUM_OF_TYPES];
    mutable bool   my_fixed     [NUM_OF_TYPES];
    mutable bool   my_fixed_set [NUM_OF_TYPES]; // Indicates whether or not the user specified a fixity...

    ///
    /// If the user specifies OPTION=TWO, then at parse time we have to figure out whether to
    /// treat it as TWO_DOM or TWO_REC. This value will be set at parse time.
    mutable TypeConstraint::OptionEnum my_treat_two_as;

  //@}
  
  /// @name Static data
  //@{
    
    ///
    /// Brief names for each category ("mean", "var", "suscept")
    static const string s_brief_names[NUM_OF_TYPES];

    ///
    /// Singular complete names for each category ("mean", "variance", "susceptibility")
    static const string s_singular_names[NUM_OF_TYPES];

    ///
    /// Plural complete names for each category ("means", "variances", "susceptibilitys")
    static const string s_plural_names[NUM_OF_TYPES];

    ///
    /// Indicates whether a particular constraint type is allowed on a particular category
    static const bool s_type_constraint_allowed [NUM_OF_TYPES][NUM_OF_OPTIONS];
    
    ///
    /// If the user specifies option=TWO, but doesn't specify any genotypic values, then
    /// the parser will have to assume this means either TWO_DOM or TWO_REC. The
    /// static value s_default_two_behavior should be set to either TWO_DOM or TWO_REC,
    /// depending on what kind of default behavior we want.
    static TypeConstraint::OptionEnum s_default_two_treatment;

    ///
    /// Default value for a type mean
    static const double s_mean_default_value;
    
    ///
    /// Default fixedness for a type mean
    static const bool s_mean_default_fixed;

    ///
    /// Default value fro a type variance
    static const double s_var_default_value;
    
    ///
    /// Default fixedness for a type variance
    static const bool s_var_default_fixed;
    
    ///
    /// Minimum epsilon value (distance from zero)
    static const double s_var_epsilon;

    ///
    /// Default lower bound for a type variance
    static const double s_var_lower_bound;

  //@}
      
};

/** \brief Base class for all genotype specific sub_models.
  * Please note that this class can never be directly instantiated, as it is pure virtual.
  *
  * \par Purpose
  *
  * \par Parsing logic
  *
  * Although parsing is technically handled by the TypeParser object, much of the underlying logic
  * is handled in this object. In particular, the isInitialSetupConsistent() function takes care
  * of examining the user's initial settings (from the type_foo sub-block) and figuring out what
  * to do with them. Please see the documentation for that function for more information.
  *
  * \par Special options for two genotypes
  *
  * If the ConstraintType is set to TWO (not TWO_DOM or TWO_REC, simply TWO), then
  * the submodel will have to assume that the TWO option is implicitly dominant
  * or recessive. This is because it is meaningless to specify two genotypic foo-types
  * without accounting for the third genotype (in a diallelic model).
  * 
  * Therefore, there exists the functions isTwoDomByDefault() and setTwoDomByDefault().
  *
  * If the ConstraintType is in fact set to TWO, then the submodel will query the
  * isTwoDomByDefault() function to see how it should interpret the TWO option. You can
  * change such default behavior with the setTwoDomByDefault() function.
  *
  */
class NewTypeSpecificSubmodel : public MAXFUN::NewSubmodel
{
  public:

  /// @name Friend classes
  //@{
  
    friend class MAXFUN::ParameterMgr;
    friend class NewTypeSpecificParser;

  //@}
  
  /// @name Enumerations & typedefs
  //@{
  
    ///
    /// Indicates which category of content this submodel will be describing.
    enum CategoryEnum
    {
      MEAN           = 0, /// For submodels describing means
      VARIANCE       = 1, /// For submodels describing variances
      SUSCEPTIBILITY = 2  /// For submodels describing susceptibilities
    };

  //@}

  /// @name Required virtual interface
  //@{
  
    virtual int update ();

    virtual MAXFUN::NewSubmodelShPtr clone();
            
  //@}

  /// @name Category-specific functions
  //@{
  
    ///
    /// Returns the category of this object.
    CategoryEnum getCategory() const;

    ///
    /// Returns the brief name of the type of data stored by class ("mean", "var", "suscept").
    string getBriefName() const;

    ///
    /// Returns the name of the type of data stored by class ("mean", "variance", "susceptibility")
    string getSingularName() const;
    
    ///
    /// Returns the name of the type of data stored by class ("means", "variances", "susceptibilities")
    string getPluralName() const;
  
    ///
    /// Returns \c true if the given ConstraintTypeEnum is allowed, \c false otherwise.
    /// \param e The ConstraintTypeEnum to evaluate
    bool isConstraintTypeAllowed(TypeConstraint::OptionEnum e) const;

    ///
    /// If your analysis uses both MEAN and VARIANCE submodels, you must make sure that the settings assigned
    /// in both submodels are compatible with each other. In particular, you need to make sure that
    /// the number of variances does not exceed the number of means. This function carries out such an evaluation.
    /// \param mean_submodel The mean submodel
    /// \param var_submodel The variance submodel
    /// \retval 0 The submodels are compatible
    /// \retval 1 The submodels are incompatible because the number of variances exceeds the number of means
    static int areMeanAndVarianceCompatible(const NewTypeSpecificSubmodel & mean_submodel, const NewTypeSpecificSubmodel & var_submodel);

  //@}

  /// @name Querying configuration
  //@{
  
    ///
    /// Returns the TypeConstraintEnum for this object.
    const TypeConstraint & getTypeConstraint() const;

    ///
    /// Returns a string describing both the genotype constraint and the type of data stored by the class.
    /// For instance: "three means, decreasing"
    string convertConfigurationToString() const;  

    ///
    /// Returns the current value corresponding to the given index.
    /// \param gt The index in question (0, 1, or 2)
    double getValue(size_t gt) const;

    ///
    /// Compares the TypeConstraint against the initial values for the types, and indicates if
    /// there are any conflicts.
    ///
    /// First of all, if the cateogry of this submodel is VARIANCE, and any given value is less than zero,
    /// an error will be produced.
    ///
    /// To better understand this function, here is a comprehensive list of all the situations examined:
    ///
    /// \par Any OPTION
    ///
    /// If no values given, then parsing succeeds.
    ///
    /// \par OPTION == ONE
    ///
    /// If at least two different values are given, then an error is produced.
    /// 
    /// If all given values are identical, then parsing succeeds.
    ///
    /// \par OPTION == TWO / TWO_DOM / TWO_REC
    ///
    /// If only one value is given, an error is produced.
    ///
    /// \par OPTION == TWO
    ///
    /// If AA and AB are given, and AA == AB, then an error is produced.
    ///
    /// If AA and AB are given, and AA != AB, then the option is interpreted as TWO_REC, and BB = AB (with warning #1)
    ///
    /// If AA and BB are given, and AA == BB, then the option is interpreted as getDefaultTwoTreatment(), and AB = AA.
    ///
    /// If AA and BB are given, and AA != BB, then the option is interpreted as getDefaultTwoTreatment(), and AB = the appropriate value.
    ///
    /// If AB and BB are given, and AB == BB, then an error is produced.
    ///
    /// If AB and BB are given, and AB != BB, then the option is interpreted as TWO_DOM, and AA = AB (with warning)
    ///
    /// If all values are given, and all values are unique, then an error is produced.
    ///
    /// If all values are given, and all values are the same, then the option is interpreted as getDefaultTwoTreatment().
    ///
    /// If all values are given, AA == AB, AA != BB, then the option is interpreted as TWO_DOM.
    /// 
    /// If all values are given, AA == BB, AA != AB, then an error is produced.
    /// 
    /// If all values are given, AB == BB, AA != AB, then the option is interpreted as TWO_REC.
    ///
    /// \par OPTION == TWO_DOM
    ///
    /// If AA and AB are given, and AA != AB, then an error is produced.
    //
    /// If AA and AB are given, and BB is not set, then an error is produced.
    ///
    /// If AA or AB is given, and BB is given, then AA or AB (whichever is missing) is assigned to the appropriate value.
    ///
    /// \par OPTION == TWO_REC
    ///
    /// If AB and BB are set, and AB != BB, then an error is produced.
    ///
    /// If AB and BB are given, and AA is not set, then an error is produced.
    ///
    /// If AB or BB is given, and AA is given, then AB or BB (whichever is missing) is assigned to the appropriate value.
    ///
    /// \par OPTION == THREE
    ///
    /// If only one or two values are given, an error is produced.
    ///
    /// \par OPTION == THREE_ADD
    ///
    /// If three values are given, and they do not meet the THREE_ADD requirements, an error is produced.
    ///
    /// \par OPTION == THREE_INC
    ///
    /// If three values are given, and they do not meet the THREE_INC requirements, an error is produced.
    ///
    /// \par OPTION == THREE_DEC
    ///
    /// If three values are given, and they do not meet the THREE_DEC requirements, an error is produced.
    ///        
    /// \par Evaluating fixedness
    ///
    /// The next step in evaluating the consistency of the initial settings is making sure there are no
    /// conflicting fixedness specifiers.
    ///
    /// \par OPTION = any
    ///
    /// If no fixed statuses are given, then everything's fine.
    ///
    /// \par OPTION = ONE
    ///
    /// If all given fixed statuses = X, then any remaining unassigned fix statuses are set to X.
    ///
    /// If more than one fixed statuses are given and they are different, an error is produced.
    /// 
    /// \par OPTION = TWO (interpreted as TWO_DOM) or TWO_DOM
    ///
    /// If fixed status is given for AA and AB, and said statuses are different, an error is produced.
    ///
    /// If fixed status X is given for AA, then the fixed status of AB is set to X.
    ///
    /// If fixed status X is given for AB, then the fixed status of AA is set to X.
    ///
    /// \par OPTION = TWO (interpreted as TWO_REC) or TWO_REC
    ///
    /// If fixed status is given for AB and BB, and said statuses are different, an error is produced.
    ///
    /// If fixed status X is given for AB, then the fixed status of BB is set to X.
    ///
    /// If fixed status X is given for BB, then the fixed status of AB is set to X.
    ///
    /// \param warning_code This value is set by isInitialSetupConsistent() to indicate any warning conditions encountered.
    /// Possible values include:
    /// 0 - No warning condition
    /// 1 - OPTION=TWO, AA and BB are given but not equal; interpreting as default behavior given by getDefaultTwoTreatment()
    /// 2 - OPTION=TWO, AB and BB are given but not equal; interpreting as TWO_DOM
    ///
    /// \retval 0 Everything is ok
    /// \retval 1 OPTION=ONE, but at least two different values given
    /// \retval 2 OPTION=TWO_DOM, but different values given for AA and AB
    /// \retval 3 OPTION=TWO_DOM, valid value given for dominant type (AA/AB), but no value given for recessive type (BB)
    /// \retval 4 OPTION=TWO_REC, but different values given for AB and BB
    /// \retval 5 OPTION=TWO_REC, valid value given for dominant type (AB/BB), but no value given for recessive type (AA)
    /// \retval 6 OPTION=TWO/TWO_DOM/TWO_REC, but only one value specified
    /// \retval 7 OPTION=THREE_ADD, three values given, but do not meet relational requirement
    /// \retval 8 OPTION=THREE_DEC, three values given, but do not meet relational requirement
    /// \retval 9 OPTION=THREE_INC, three values given, but do not meet relational requirement
    /// \retval 10 OPTION=THREE_ADD/THREE_DEC/THREE_INC, at least one value given but not all three
    /// \retval 11 OPTION=TWO, but three different values given
    /// \retval 12 OPTION=TWO, AA equals AB, but BB not given; not sure how to interpret this
    /// \retval 13 OPTION=TWO, AB equals BB, but AA not given; not sure how to interpret this
    /// \retval 14 OPTION=TWO, AA equals BB, and AB != AA or BB; not sure how to interpret this
    /// \retval 100 OPTION=ONE, at least two different fixedness statuses given
    /// \retval 101 OPTION=TWO (interpreted as TWO_DOM) or TWO_DOM, conflicting fixed status given for AA and AB
    /// \retval 102 OPTION=TWO (interpreted as TWO_REC) or TWO_REC, conflicting fixed status given for AB and BB
    /// \retval 200 Variance sub-block indicated, value less than zero given
    int isInitialSetupConsistent(int & warning_code) const;

  //@}
    
  /// \brief Setting configuration
  //@{
  
    ///
    /// Sets the category for this submodel (see CategoryEnum).
    bool setCategory(CategoryEnum c);

    ///
    /// Sets this object's TypeConstraint to the given value.
    /// Please note that this option should be set \b prior to setting the values
    /// for individual genotypes (setTypeValues() ). Those set functions
    /// depend on the TypeConstraint.
    /// \param constraint The new TypeConstraint value
    /// \retval true The type was successfully set
    /// \retval false The type was not successfully set, because the derived submodel does not allow it
    bool setTypeConstraint(const TypeConstraint & constraint);
    
    ///
    /// Sets this object's TypeConstraint to the given string value.
    /// \param t String for the constraint type, such as "ONE", "THREE_ADD", etc. (see TypeConstraint::OptionEnum)
    /// \retval true The type was successfully set
    /// \retval false The type was not successfully set, because the derived submodel does not allow it
    bool setTypeConstraint(const string & t);

    ///
    /// Sets the options for the given genotype.
    /// \param i The index in question (see genotype_index)
    /// \param initial_val The initial value for the genotype
    /// \param fixedness_set Boolean indicating whether or not the user has set the fixedness
    /// \param fixedness If the user did indeed set a fixedness, indicates whether or not it is fixed
    /// \retval true If the new options satisfy the relational requirements given by getTypeConstraint()
    /// \retval false If the the new options conflict with the relational requirements given by getTypeConstraint()
    bool setTypeValues(size_t i, double initial_val, bool fixedness_set, bool fixedness);

  //@}

  /// @name Special option for two (geno)types
  //@{

    ///
    /// If the user specified OPTION=TWO, but did not specify any genotypic values, this function
    /// returns how the option will be interpreted.
    static TypeConstraint::OptionEnum getDefaultTwoTreatment();
    
    ///
    /// Sets the default interpretive behavior for ConstraintType = TWO when no genotypic values
    /// are specified by the user. The parameter in question can \b only be TypeConstraint::TWO_DOM
    /// or TypeConstraint::TWO_REC. 
    /// \param e Indicates default behavior (see above)
    /// \retval true Value set successfully
    /// \retval false Value not set successfully (because \c e was neither TWO_DOM nor TWO_REC)
    static bool setDefaultTwoTreatment(TypeConstraint::OptionEnum e);
    
    ///
    /// If the user specified OPTION=TWO, but did not specify any genotypic values, then the parser
    /// will have figured out how to interpret this option. You can find out what the parser decided
    /// to do with OPTION=TWO by using this function.
    TypeConstraint::OptionEnum treatThisTwoAs() const;

    ///
    /// Sets the interpretive behavoir for this specific submodel. The parser should use this function
    /// once it has determined how to treat OPTION=TWO.
    ///
    /// Note: Why is this function const? Because it is called by isInitialSetupConsistent(), which is
    /// const.
    /// \retval true Value set successfully
    /// \retval false Value not set successfully (because \c e was neither TWO_DOM nor TWO_REC)
    bool setTreatThisTwoAs(TypeConstraint::OptionEnum e) const;

  //@}

  /// @name Blah

    /// Tests the submodel for completeness, which is defined as all
    /// parameters specified as finite values.
    bool is_complete() const;

  //@}

  /// @name Debugging
  //@{
    
    /// \internal
    /// Dumps the contents of this object.
    void dump(std::ostream& out) const;

  //@}

protected:

  /// @name Constructors
  //@{

    ///
    /// This function is protected so it can only be used by a derived type.
    /// \param c The category of content this submodel will be describing
    /// \param errors The errorstream to use
    NewTypeSpecificSubmodel(CategoryEnum c = MEAN, cerrorstream& errors = sage_cerr);

    ///
    /// This function is protected so it can only be used by a derived type.
    NewTypeSpecificSubmodel(const NewTypeSpecificSubmodel& other);

    ///
    /// This function is protected so it can only be used by a derived type.
    NewTypeSpecificSubmodel& operator=(const NewTypeSpecificSubmodel& other);

  //@}


  /// \brief Maximization Initialization
  ///
  //@{
      
    ///
    /// This function synchronizes the my_parameters vector with what's in
    /// my_types and my_types fixed.  Necessary for inclusion in maxfun
    /// or other times where that must be synched.
    virtual int finalizeConfiguration();

    ///
    /// Initialize my_parameters if TypeConstraint = ONE
    void initializeOne();

    ///
    /// Initialize my_parameters for TypeConstraint = TWO_DOM
    /// (A dominant over B (\f$t_AB = t_AA\f$) )
    void initializeTwoDom();

    ///
    /// Initialize my_parameters for TypeConstraint = TWO_REC
    /// (A recessive under B (\f$t_AB = t_BB\f$) )
    void initializeTwoRec();

    ///
    /// Initialize my_parameters for TypeConstraint = THREE
    void initializeThree();

    ///
    /// Initialize my_parameters for TypeConstraint = THREE_ADD
    /// (\f$t_AB = (t_AA + t_BB) / 2\f$)
    void initializeThreeAdd();

  //@}

  /// @name Option specific synchronization functions (sync w. Maxfun).
  //@{

    ///
    /// Synchronizes current parameter estimates with Maxfun for ConstraintType = ONE.
    /// \retval 0 Synchronization successful
    /// \retval !=0 Synchronization not successful
    int synchronizeOne();

    ///
    /// Synchronizes current parameter estimates with Maxfun for ConstraintType = TWO.
    /// \retval 0 Synchronization successful
    /// \retval !=0 Synchronization not successful
    int synchronizeTwo();

    ///
    /// Synchronizes current parameter estimates with Maxfun for ConstraintType = THREE.
    /// \retval 0 Synchronization successful
    /// \retval !=0 Synchronization not successful
    int synchronizeThree();

    ///
    /// Synchronizes current parameter estimates with Maxfun for ConstraintType = TWO_DOM.
    /// \retval 0 Synchronization successful
    /// \retval !=0 Synchronization not successful
    int synchronizeTwoDom();

    ///
    /// Synchronizes current parameter estimates with Maxfun for ConstraintType = TWO_REC.
    /// \retval 0 Synchronization successful
    /// \retval !=0 Synchronization not successful
    int synchronizeTwoRec();

    ///
    /// Synchronizes current parameter estimates with Maxfun for ConstraintType = THREE_ADD.
    /// \retval 0 Synchronization successful
    /// \retval !=0 Synchronization not successful
    int synchronizeThreeAdd();

    ///
    /// Synchronizes current parameter estimates with Maxfun for ConstraintType = THREE_DEC.
    /// \retval 0 Synchronization successful
    /// \retval !=0 Synchronization not successful
    int synchronizeThreeDec();

    ///
    /// Synchronizes current parameter estimates with Maxfun for ConstraintType = THREE_INC.
    /// \retval 0 Synchronization successful
    /// \retval !=0 Synchronization not successful
    int synchronizeThreeInc();

  //@}    

  /// Constraint Checking
  //@{

    ///
    /// Returns \c true if the relational constraints of THREE_ADD are met, \c false otherwise.
    bool isThreeAddConstraintMet() const;

    ///
    /// Returns \c true if the relational constraints of THREE_DEC are met, \c false otherwise.
    bool isThreeDecConstraintMet() const;

    ///
    /// Returns \c true if the relational constraints of THREE_INC are met, \c false otherwise.
    bool isThreeIncConstraintMet() const;

  //@}
        
private:

  /// @name Data members
  //@{
  
    CategoryEnum   my_category;
    TypeConstraint my_type_constraint;

    mutable double my_types     [NUM_OF_TYPES];
    mutable bool   my_fixed     [NUM_OF_TYPES];
    mutable bool   my_fixed_set [NUM_OF_TYPES]; // Indicates whether or not the user specified a fixity...

    ///
    /// If the user specifies OPTION=TWO, then at parse time we have to figure out whether to
    /// treat it as TWO_DOM or TWO_REC. This value will be set at parse time.
    mutable TypeConstraint::OptionEnum my_treat_two_as;

  //@}
  
  /// @name Static data
  //@{
    
    ///
    /// Brief names for each category ("mean", "var", "suscept")
    static const string s_brief_names[NUM_OF_TYPES];

    ///
    /// Singular complete names for each category ("mean", "variance", "susceptibility")
    static const string s_singular_names[NUM_OF_TYPES];

    ///
    /// Plural complete names for each category ("means", "variances", "susceptibilitys")
    static const string s_plural_names[NUM_OF_TYPES];

    ///
    /// Indicates whether a particular constraint type is allowed on a particular category
    static const bool s_type_constraint_allowed [NUM_OF_TYPES][NUM_OF_OPTIONS];
    
    ///
    /// If the user specifies option=TWO, but doesn't specify any genotypic values, then
    /// the parser will have to assume this means either TWO_DOM or TWO_REC. The
    /// static value s_default_two_behavior should be set to either TWO_DOM or TWO_REC,
    /// depending on what kind of default behavior we want.
    static TypeConstraint::OptionEnum s_default_two_treatment;

    ///
    /// Default value for a type mean
    static const double s_mean_default_value;
    
    ///
    /// Default fixedness for a type mean
    static const bool s_mean_default_fixed;

    ///
    /// Default value fro a type variance
    static const double s_var_default_value;
    
    ///
    /// Default fixedness for a type variance
    static const bool s_var_default_fixed;
    
    ///
    /// Minimum epsilon value (distance from zero)
    static const double s_var_epsilon;

    ///
    /// Default lower bound for a type variance
    static const double s_var_lower_bound;

  //@}
      
};

extern const MAXFUN::SMType<MFSUBMODELS::NewTypeSpecificSubmodel> new_type_specific_submodel;

} // End namespace MFSUBMODELS
} // End namespace SAGE

#include "mfsubmodels/TypeSubmodel.ipp"

#endif
