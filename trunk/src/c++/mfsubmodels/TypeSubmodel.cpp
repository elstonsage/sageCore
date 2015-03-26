#include "util/StringUtils.h"
#include "globals/SAGEConstants.h"
#include "mfsubmodels/TypeSubmodel.h"
#include "output/Output.h"
#include "output/ViewPrettyPrint.h"

namespace SAGE        {
namespace MFSUBMODELS {

// Extern variables instantiated here:

const std::string  TYPE_MEAN_NAME        = "Type means";
const std::string  TYPE_SUSCEPT_NAME     = "Type susceptibilities";
const std::string  TYPE_VAR_NAME         = "Type variances";

// Static variables instantiated here:

const double TypeSpecificSubmodel::s_mean_default_value = numeric_limits<double>::quiet_NaN();  // 0;
const bool   TypeSpecificSubmodel::s_mean_default_fixed = false;
const double TypeSpecificSubmodel::s_var_default_value  = numeric_limits<double>::quiet_NaN();  // 1;
const bool   TypeSpecificSubmodel::s_var_default_fixed  = false;
const double TypeSpecificSubmodel::s_var_epsilon        = 0.00001;
const double TypeSpecificSubmodel::s_var_lower_bound    = 0;

const string TypeSpecificSubmodel::s_brief_names    [NUM_OF_TYPES] = { "mean",  "var",       "suscept"};
const string TypeSpecificSubmodel::s_singular_names [NUM_OF_TYPES] = { "mean",  "variance",  "susceptibility"};
const string TypeSpecificSubmodel::s_plural_names   [NUM_OF_TYPES] = { "means", "variances", "susceptibilities"};

const bool TypeSpecificSubmodel::s_type_constraint_allowed[NUM_OF_TYPES][NUM_OF_OPTIONS] = { 
  { true, true, true, true, true, true, true,  true  },
  { true, true, true, true, true, true, false, false },
  { true, true, true, true, true, true, true,  true  } };

TypeConstraint::OptionEnum TypeSpecificSubmodel::s_default_two_treatment = TypeConstraint::TWO_DOM;

//=============================================
//
//  Constructor
//
//=============================================
TypeSpecificSubmodel::TypeSpecificSubmodel(CategoryEnum c, cerrorstream & errors)  
    : 
    Submodel           (errors),
    my_category        (c),
    my_type_constraint (TypeConstraint::ONE)
{
  my_types[0] = SAGE::QNAN;
  my_types[1] = SAGE::QNAN;
  my_types[2] = SAGE::QNAN;

  my_fixed[0] = false;
  my_fixed[1] = false;
  my_fixed[2] = false;

  my_fixed_set[0] = false;
  my_fixed_set[1] = false;
  my_fixed_set[2] = false;
}

//===============================================
//
//  Copy constructor
//
//===============================================
TypeSpecificSubmodel::TypeSpecificSubmodel(const TypeSpecificSubmodel& other)
    : MAXFUN::Submodel   (other),
      my_category        (other.my_category),
      my_type_constraint (other.my_type_constraint),
      my_treat_two_as    (other.my_treat_two_as)
{
  for(int i = 0; i < NUM_OF_TYPES; ++i)
  {
    my_types     [i] = other.my_types     [i];
    my_fixed     [i] = other.my_fixed     [i];
    my_fixed_set [i] = other.my_fixed_set [i];
  }
}

//===============================================
//
//  operator=(...)
//
//===============================================
TypeSpecificSubmodel&
TypeSpecificSubmodel::operator=(const TypeSpecificSubmodel& other)
{
  if(this != &other)
  {
    MAXFUN::Submodel::operator=(other);
    
    my_category        = other.my_category;
    my_type_constraint = other.my_type_constraint;
    my_treat_two_as    = other.my_treat_two_as;
    
    for(int i = 0; i < NUM_OF_TYPES; ++i)
    {
      my_types     [i] = other.my_types     [i];
      my_fixed     [i] = other.my_fixed     [i];
      my_fixed_set [i] = other.my_fixed_set [i];
    }
  }
  
  return *this;
}

//==============================================
//
//  dump()
//
//==============================================
void
TypeSpecificSubmodel::dump(std::ostream& out) const
{
  OUTPUT::Table t(getSingularName());
  
  t << (OUTPUT::TableRow() << ("type_" + getBriefName()))
    << (OUTPUT::TableRow() << "option" 
                           << getTypeConstraint().toString() 
                           << (getTypeConstraint().getOption() == TypeConstraint::TWO ? "# interpreted as " + TypeConstraint(treatThisTwoAs()).toString() : ""))
    << (OUTPUT::TableRow() << getBriefName() << "=AA, val=" << my_types[0] << ", fixed=" << my_fixed[0])
    << (OUTPUT::TableRow() << getBriefName() << "=AB, val=" << my_types[1] << ", fixed=" << my_fixed[1])
    << (OUTPUT::TableRow() << getBriefName() << "=BB, val=" << my_types[2] << ", fixed=" << my_fixed[2]);
    
  out << t;
}

//===============================================================
//
//  convertConfigurationToString()
//
//===============================================================
string 
TypeSpecificSubmodel::convertConfigurationToString() const
{
  string type;                        

  type = my_type_constraint.getTypeCountAsString() + " " +
         (my_type_constraint.getTypeCount() == 1 ? getSingularName() : getPluralName()) + 
         (my_type_constraint.hasTypeRelation() ? ", " + my_type_constraint.getTypeRelationAsString() : "");
         
  if(!type.size()) 
    SAGE_internal_error();

  return type;
}

//===============================================================
//
//  setTypeValues(...)
//
//===============================================================
bool
TypeSpecificSubmodel::setTypeValues(
  size_t i, 
  double initial_val, 
  bool   fixed_set,
  bool   fixedness)
{
  // Set the value and the fixity:
  my_types     [i] = initial_val;
  my_fixed     [i] = fixedness;
  my_fixed_set [i] = fixed_set;

  // Return success!
  return true;
}


//==============================================================
//
//  isInitialSetupConsistent()
//
//==============================================================
int 
TypeSpecificSubmodel::isInitialSetupConsistent(int & warning_code) const
{
  // Grab necessary values for analysis:

  double AA                  = my_types[0],
         AB                  = my_types[1],
         BB                  = my_types[2];
  bool   AA_set              = !SAGE::isnan(AA),
         AB_set              = !SAGE::isnan(AB),
         BB_set              = !SAGE::isnan(BB),

         no_values_set       = !AA_set && !AB_set && !BB_set,

         one_value_set       = ( AA_set && !AB_set && !BB_set) ||
                               (!AA_set &&  AB_set && !BB_set) ||
                               (!AA_set && !AB_set &&  BB_set),

         two_values_set      = ( AA_set &&  AB_set && !BB_set) ||
                               ( AA_set && !AB_set &&  BB_set) ||
                               (!AA_set &&  AB_set &&  BB_set),
 
         three_values_set    = ( AA_set &&  AB_set &&  BB_set),
 
         at_least_two_values_set = two_values_set || three_values_set,

         all_set_values_same = one_value_set 

                               ||

                               (( AA_set &&  AB_set && !BB_set) && (AA == AB) ||
                                ( AA_set && !AB_set &&  BB_set) && (AA == BB) ||
                                (!AA_set &&  AB_set &&  BB_set) && (AB == BB)) 

                               ||

                               (three_values_set && AA == AB && AA == BB),

         all_set_values_unique = three_values_set && (AA != AB) && (AB != BB) && (AA != BB);

  double any_value = AA_set ? AA : (AB_set ? AB : (BB_set ? BB : QNAN));

  //===================================================================
  // If all preset values are QNAN, then it's all ok, because the 
  // analysis will figure out valid initial values later on.
  //===================================================================

  if(no_values_set)
  {
    return 0;
  }

  //===================================================================
  // If this is a VARIANCE submodel, and any value is <= 0, return an 
  // error!
  //===================================================================

  else if((getCategory() == VARIANCE) && ((AA_set && AA <= s_var_epsilon + s_var_lower_bound) || 
                                          (AB_set && AB <= s_var_epsilon + s_var_lower_bound) || 
                                          (BB_set && BB <= s_var_epsilon + s_var_lower_bound)))
  {
    return 200;
  }

  //===================================================================
  // If option is ONE...
  //===================================================================

  else if(getTypeConstraint().getTypeCount() == 1)
  {
    if(all_set_values_same)
    {
      my_types[0] = my_types[1] = my_types[2] = any_value;
    }
    else
    {
      return 1;
    }
  }

  //===================================================================
  // If option is TWO...
  //===================================================================

  else if(getTypeConstraint().getTypeCount() == 2)
  {
    if(at_least_two_values_set)
    {
      if(getTypeConstraint().getOption() == TypeConstraint::TWO)
      {
        if(two_values_set)
        {
          if(AA_set && AB_set && AA == AB)
          {
            return 12;
          }
          else if(AA_set && AB_set && AA != AB)
          {
            setTreatThisTwoAs(TypeConstraint::TWO_REC);

            my_types[2] = my_types[1];
          }
          else if(AA_set && BB_set && AA == BB)
          {
            setTreatThisTwoAs(getDefaultTwoTreatment());

            my_types[1] = my_types[0];
          }
          else if(AA_set && BB_set && AA != BB)
          {
            setTreatThisTwoAs(getDefaultTwoTreatment());

            my_types[1] = treatThisTwoAs() == TypeConstraint::TWO_DOM ? my_types[0] : my_types[1];

            warning_code = 1;
          }
          else if(AB_set && BB_set && AB == BB)
          {
            return 13;
          }
          else if(AB_set && BB_set && AB != BB)
          {
            setTreatThisTwoAs(TypeConstraint::TWO_DOM);

            my_types[0] = my_types[1];

            warning_code = 2;
          }
        }
        else if(three_values_set && all_set_values_same)
        {
          setTreatThisTwoAs(getDefaultTwoTreatment());
        }
        else if(three_values_set && all_set_values_unique)
        {
          return 11;
        }
        else if(three_values_set) // We can assume there are only two unique values set
        {
          if(AA == AB)
          {
            setTreatThisTwoAs(TypeConstraint::TWO_DOM);
          }
          else if(AA == BB)
          {
            return 14;
          }
          else if(AB == BB)
          {
            setTreatThisTwoAs(TypeConstraint::TWO_REC);
          }
        }
      }
      else if(getTypeConstraint().getOption() == TypeConstraint::TWO_DOM)
      {
        if(AA_set && AB_set && (AA != AB))
        {
          return 2;
        }
        else if(!BB_set)
        {
          return 3;
        }
        else
        {
          AA_set ? my_types[1] = my_types[0] : my_types[0] = my_types[1];
        }
      }
      else // option == TWO_REC
      {
        if(AB_set && BB_set && (AB != BB))
        {
          return 4;
        }
        else if(!AA_set)
        {
          return 5;
        }
        else
        {
          AB_set ? my_types[2] = my_types[1] : my_types[1] = my_types[2];
        }
      }
    }
    else // !at_least_two_values_set
    {
      return 6;
    }
  }

  //====================================================================
  // If option is THREE...
  //====================================================================

  else if(getTypeConstraint().getTypeCount() == 3)
  {
    if(three_values_set)
    {
      if((getTypeConstraint().getOption() == TypeConstraint::THREE_ADD) && !isThreeAddConstraintMet())
      {
        return 7;
      }
      else if((getTypeConstraint().getOption() == TypeConstraint::THREE_DEC) && !isThreeDecConstraintMet())
      {
        return 8;
      }
      else if((getTypeConstraint().getOption() == TypeConstraint::THREE_INC) && !isThreeIncConstraintMet())
      {
        return 9;
      }
    }
    else // !all_three_values_set
    {
      return 10;
    }
  }

  // Evaluate fixedness stuff:

  bool AA_fixed = my_fixed[0],
       AB_fixed = my_fixed[1],
       BB_fixed = my_fixed[2],

       AA_fixed_set = my_fixed_set[0],
       AB_fixed_set = my_fixed_set[1],
       BB_fixed_set = my_fixed_set[2],

       one_fixed_set = ( AA_fixed_set && !AB_fixed_set && !BB_fixed_set) ||
                       (!AA_fixed_set &&  AB_fixed_set && !BB_fixed_set) ||
                       (!AA_fixed_set && !AB_fixed_set &&  BB_fixed_set),

       two_fixed_set = ( AA_fixed_set &&  AB_fixed_set && !BB_fixed_set) ||
                       ( AA_fixed_set && !AB_fixed_set &&  BB_fixed_set) ||
                       (!AA_fixed_set &&  AB_fixed_set &&  BB_fixed_set),

       three_fixed_set = AA_fixed_set && AB_fixed_set && BB_fixed_set,

       at_least_one_fixed_set = one_fixed_set || two_fixed_set || three_fixed_set,

       any_fixed = AA_fixed_set ? AA_fixed : (AB_fixed_set ? AB_fixed : BB_fixed),

       all_set_fixed_same = one_fixed_set

                             ||

                             (( AA_fixed_set &&  AB_fixed_set && !BB_fixed_set) && (AA_fixed == AB_fixed) ||
                              ( AA_fixed_set && !AB_fixed_set &&  BB_fixed_set) && (AA_fixed == BB_fixed) ||
                              (!AA_fixed_set &&  AB_fixed_set &&  BB_fixed_set) && (AB_fixed == BB_fixed))

                             ||

                             (three_values_set && AA_fixed == AB_fixed && AA_fixed == BB_fixed);

  if(!at_least_one_fixed_set)
  {
    ;
  }
  else if(getTypeConstraint().getTypeCount() == 1)
  {
    if(all_set_fixed_same)
    {
      my_fixed[0] = my_fixed[1] = my_fixed[2] = any_fixed;
    }
    else
    {
      return 100;
    }
  }
  else if(getTypeConstraint().getTypeCount() == 2)
  {
    if(getTypeConstraint().getOption() == TypeConstraint::TWO_DOM || 
       (getTypeConstraint().getOption() == TypeConstraint::TWO && treatThisTwoAs() == TypeConstraint::TWO_DOM))
    {
      if(AA_fixed_set && AB_fixed_set && (AA_fixed != AB_fixed))
      {
        return 101;
      }
      else if(AA_fixed_set && !AB_fixed_set)
      {
        my_fixed[1] = my_fixed[0];
      }
      else if(!AA_fixed_set && AB_fixed_set)
      {
        my_fixed[0] = my_fixed[1];
      }
    }
    else if(getTypeConstraint().getOption() == TypeConstraint::TWO_REC || 
           (getTypeConstraint().getOption() == TypeConstraint::TWO && treatThisTwoAs() == TypeConstraint::TWO_REC))
    {
      if(AB_fixed_set && BB_fixed_set && (AB_fixed != BB_fixed))
      {
        return 102;
      }
      else if(AB_fixed_set && !BB_fixed_set)
      {
        my_fixed[2] = my_fixed[1];
      }
      else if(!AB_fixed_set && BB_fixed_set)
      {
        my_fixed[1] = my_fixed[2];
      }
    }
  }

  // Return success:

  return 0;
}

//==============================================================
//
//  areMeanAndVarianceCompatible(...)
//
//==============================================================
int 
TypeSpecificSubmodel::areMeanAndVarianceCompatible(const TypeSpecificSubmodel & mean_submodel, const TypeSpecificSubmodel & var_submodel)
{
  if(var_submodel.getTypeConstraint().getTypeCount() > mean_submodel.getTypeConstraint().getTypeCount())
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

//===============================================================
//
//  finalizeConfiguration()
//
//===============================================================
int
TypeSpecificSubmodel::finalizeConfiguration()
{
  switch(my_type_constraint.getOption())
  {
    case TypeConstraint::ONE:       initializeOne      ();  break;
    case TypeConstraint::TWO:       treatThisTwoAs() == TypeConstraint::TWO_DOM ?  
                                    initializeTwoDom   () : 
                                    initializeTwoRec   ();  break;
    case TypeConstraint::THREE:     initializeThree    ();  break;
    case TypeConstraint::TWO_DOM:   initializeTwoDom   ();  break;
    case TypeConstraint::TWO_REC:   initializeTwoRec   ();  break;
    case TypeConstraint::THREE_ADD: initializeThreeAdd ();  break;
    case TypeConstraint::THREE_DEC:
    case TypeConstraint::THREE_INC:

      // For threeDec and threeInc, they have the same values as the basic
      // three model, but have a functional relationship to one another.

      initializeThree();
      
      // Redefine any non-fixed as independent functional.
      
      my_parameters[0].initial_type = my_fixed[0] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL;
      my_parameters[1].initial_type = my_fixed[1] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL;
      my_parameters[2].initial_type = my_fixed[2] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL;

      break;
  }

  return 0;
}

//===============================================================
//
//  initializeOne()
//
//===============================================================
void
TypeSpecificSubmodel::initializeOne()
{
  my_parameters.resize(1);

  my_parameters[0] = MAXFUN::ParameterInput(
      toUpper(getPluralName()),
      getBriefName(),
      my_fixed[0] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
      my_types[0],
      my_category == VARIANCE ? s_var_lower_bound + s_var_epsilon : NEGATIVE_INF,
      POSITIVE_INF);
}

//===============================================================
//
//  initializeTwoDom()
//
//===============================================================
void
TypeSpecificSubmodel::initializeTwoDom()
{
  my_parameters.resize(2);

  my_parameters[0] = MAXFUN::ParameterInput
    ( toUpper(getPluralName()),
      getBriefName() + "_AA_AB",
      my_fixed[0] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
      my_types[0],
      my_category == VARIANCE ? s_var_lower_bound + s_var_epsilon : NEGATIVE_INF,
      POSITIVE_INF);

  my_parameters[1] = MAXFUN::ParameterInput
    ( toUpper(getPluralName()),
      getBriefName() + "_BB",
      my_fixed[2] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
      my_types[2],
      my_category == VARIANCE ? s_var_lower_bound + s_var_epsilon : NEGATIVE_INF,
      POSITIVE_INF);
}

//===============================================================
//
//  initializeTwoRec()
//
//===============================================================
void
TypeSpecificSubmodel::initializeTwoRec()
{ 
  my_parameters.resize(2);

  my_parameters[0] = MAXFUN::ParameterInput
    ( toUpper(getPluralName()),
      getBriefName() + "_AA",
      my_fixed[0] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
      my_types[0],
      my_category == VARIANCE ? s_var_lower_bound + s_var_epsilon : NEGATIVE_INF,
      POSITIVE_INF);

  my_parameters[1] = MAXFUN::ParameterInput
    ( toUpper(getPluralName()),
      getBriefName() + "_AB_BB",
      my_fixed[2] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
      my_types[2],
      my_category == VARIANCE ? s_var_lower_bound + s_var_epsilon : NEGATIVE_INF,
      POSITIVE_INF);
}


//===============================================================
//
//  initializeThreeAdd()
//
//===============================================================
void
TypeSpecificSubmodel::initializeThreeAdd()
{
  my_parameters.resize(3);

  my_parameters[0] = MAXFUN::ParameterInput
    ( toUpper(getPluralName()),
      getBriefName() + "_AA",
      my_fixed[0] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL,
      my_types[0],
      my_category == VARIANCE ? s_var_lower_bound + s_var_epsilon : NEGATIVE_INF,
      POSITIVE_INF);

  my_parameters[1] = MAXFUN::ParameterInput
    ( toUpper(getPluralName()),
      getBriefName() + "_AB",
      my_fixed[1] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::DEPENDENT,
      my_types[1],
      my_category == VARIANCE ? s_var_lower_bound + s_var_epsilon : NEGATIVE_INF,
      POSITIVE_INF);

  my_parameters[2] = MAXFUN::ParameterInput
    ( toUpper(getPluralName()),
      getBriefName() + "_BB",
      my_fixed[2] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL,
      my_types[2],
      my_category == VARIANCE ? s_var_lower_bound + s_var_epsilon : NEGATIVE_INF,
      POSITIVE_INF);
}

//===============================================================
//
//  initializeThree()
//
//===============================================================
void
TypeSpecificSubmodel::initializeThree()
{
  my_parameters.resize(3);

  my_parameters[0] = MAXFUN::ParameterInput
    ( toUpper(getPluralName()),
      getBriefName() + "_AA",
      my_fixed[0] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
      my_types[0],
      my_category == VARIANCE ? s_var_lower_bound + s_var_epsilon : NEGATIVE_INF,
      POSITIVE_INF);

  my_parameters[1] = MAXFUN::ParameterInput
    ( toUpper(getPluralName()),
      getBriefName() + "_AB",
      my_fixed[1] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
      my_types[1],
      my_category == VARIANCE ? s_var_lower_bound + s_var_epsilon : NEGATIVE_INF,
      POSITIVE_INF);

  my_parameters[2] = MAXFUN::ParameterInput
    ( toUpper(getPluralName()),
      getBriefName() + "_BB",
      my_fixed[2] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
      my_types[2],
      my_category == VARIANCE ? s_var_lower_bound + s_var_epsilon : NEGATIVE_INF,
      POSITIVE_INF);
}

//=================================================================
//
//  update()
//
//=================================================================
int
TypeSpecificSubmodel::update()
{
  int return_value = 1;

  switch(my_type_constraint.getOption())
  {
    case TypeConstraint::ONE       : return_value = synchronizeOne      (); break;
    case TypeConstraint::TWO       : return_value = synchronizeTwo      (); break;
    case TypeConstraint::THREE     : return_value = synchronizeThree    (); break;
    case TypeConstraint::TWO_DOM   : return_value = synchronizeTwoDom   (); break;
    case TypeConstraint::TWO_REC   : return_value = synchronizeTwoRec   (); break;
    case TypeConstraint::THREE_ADD : return_value = synchronizeThreeAdd (); break;
    case TypeConstraint::THREE_DEC : return_value = synchronizeThreeDec (); break;
    case TypeConstraint::THREE_INC : return_value = synchronizeThreeInc (); break;
    default                        : SAGE_internal_error();
  }
   
  return return_value;
}

//===============================================================
//
//  synchronizeOne()
//
//===============================================================
int
TypeSpecificSubmodel::synchronizeOne()
{
  my_types[0] = getParam(0);
  my_types[1] = getParam(0);
  my_types[2] = getParam(0);

  return 0;
}

//===============================================================
//
//  synchronizeTwo()
//
//===============================================================
int  
TypeSpecificSubmodel::synchronizeTwo()
{
  return getDefaultTwoTreatment() == TypeConstraint::TWO_DOM ? synchronizeTwoDom() : synchronizeTwoRec();
}
 
//===============================================================
//
//  synchronizeThree()
//
//===============================================================
int  
TypeSpecificSubmodel::synchronizeThree()
{
  my_types[0] = getParam(0);
  my_types[1] = getParam(1);
  my_types[2] = getParam(2);

  return 0;
}

//===============================================================
//
//  synchronizeTwoDom()
//
//===============================================================
int  
TypeSpecificSubmodel::synchronizeTwoDom()
{
  my_types[0] = getParam(0);
  my_types[1] = getParam(0);
  my_types[2] = getParam(1);

  return 0;
}
 
//===============================================================
//
//  synchronizeTwoRec()
//
//===============================================================
int  
TypeSpecificSubmodel::synchronizeTwoRec()
{
  my_types[0] = getParam(0);
  my_types[1] = getParam(1);
  my_types[2] = getParam(1);

  return 0;
}

//===============================================================
//
//  synchronizeThreeAdd()
//
//===============================================================
int  
TypeSpecificSubmodel::synchronizeThreeAdd()
{
  my_types[0] = getParam(0);
  my_types[2] = getParam(1);

  // For this option one parameter is a function of the other two.

  my_types[1] = (my_types[0] + my_types[2]) / 2;
  
  // Update maxfun w. newly calculated value of 'AB' parameter.

  getParam(1) = my_types[1];
  
  return 0;
}
 
//==============================================================
//
//  isThreeAddConstraintMet()
//
//==============================================================
bool
TypeSpecificSubmodel::isThreeAddConstraintMet() const
{
  return my_types[1] == (my_types[0] + my_types[2]) / 2;
}

// - Values of genotype specific means must be in decreasing order. 
//   Comparisons involving a QNAN always return true.
//
bool
TypeSpecificSubmodel::isThreeDecConstraintMet() const
{
  bool  AA_or_AB_nan = SAGE::isnan(my_types[0]) || SAGE::isnan(my_types[1]);
  bool  AB_or_BB_nan = SAGE::isnan(my_types[1]) || SAGE::isnan(my_types[2]);
  bool  AA_or_BB_nan = SAGE::isnan(my_types[0]) || SAGE::isnan(my_types[2]);
  
  bool  AA_gte_AB = AA_or_AB_nan || my_types[0] >= my_types[1];
  bool  AB_gte_BB = AB_or_BB_nan || my_types[1] >= my_types[2];
  bool  AA_gte_BB = AA_or_BB_nan || my_types[0] >= my_types[2];

  if(AA_gte_AB && AB_gte_BB && AA_gte_BB)
  {
    return true;
  }
  else
  {
    return false;
  }
}

// - Values of genotype specific means must be in increasing order.
//   Comparisons involving a QNAN always return true.
//
bool
TypeSpecificSubmodel::isThreeIncConstraintMet() const
{
  bool  AA_or_AB_nan = SAGE::isnan(my_types[0]) || SAGE::isnan(my_types[1]);
  bool  AB_or_BB_nan = SAGE::isnan(my_types[1]) || SAGE::isnan(my_types[2]);
  bool  AA_or_BB_nan = SAGE::isnan(my_types[0]) || SAGE::isnan(my_types[2]);
  
  bool  AA_lte_AB = AA_or_AB_nan || my_types[0] <= my_types[1];
  bool  AB_lte_BB = AB_or_BB_nan || my_types[1] <= my_types[2];
  bool  AA_lte_BB = AA_or_BB_nan || my_types[0] <= my_types[2];

  if(AA_lte_AB && AB_lte_BB && AA_lte_BB)
  {
    return true;
  }
  else
  {
    return false;
  }
}

//===============================================================
//
//  synchronizeThreeDec()
//
//===============================================================
int TypeSpecificSubmodel::synchronizeThreeDec()
{
  my_types[0] = getParam(0);
  my_types[1] = getParam(1);
  my_types[2] = getParam(2);

  if(isThreeDecConstraintMet())
  {
    return 0;
  }
  else
  {
    // - sub-model will now be in an inconsistent state, but this will be 
    //   corrected on the next successful synchronization.  In the interim
    //   maxfun will not call evaluate() so we do not have a problem 
    //   (per discussion w. gcw).
    //
    return 1;    
  }
}
 
//==================================================================
//
//  synchronizeThreeInc()
//
//==================================================================
int TypeSpecificSubmodel::synchronizeThreeInc()
{
  my_types[0] = getParam(0);
  my_types[1] = getParam(1);
  my_types[2] = getParam(2);

  if(isThreeIncConstraintMet())
  {
    return 0;
  }
  else
  {
    return 1;   
  }
}





//==================================================================================



const MAXFUN::SMType<MFSUBMODELS::NewTypeSpecificSubmodel> new_type_specific_submodel = MAXFUN::SMType<MFSUBMODELS::NewTypeSpecificSubmodel>();

// Static variables instantiated here:

const double NewTypeSpecificSubmodel::s_mean_default_value = numeric_limits<double>::quiet_NaN();  // 0;
const bool   NewTypeSpecificSubmodel::s_mean_default_fixed = false;
const double NewTypeSpecificSubmodel::s_var_default_value  = numeric_limits<double>::quiet_NaN();  // 1;
const bool   NewTypeSpecificSubmodel::s_var_default_fixed  = false;
const double NewTypeSpecificSubmodel::s_var_epsilon        = 0.00001;
const double NewTypeSpecificSubmodel::s_var_lower_bound    = 0;

const string NewTypeSpecificSubmodel::s_brief_names    [NUM_OF_TYPES] = { "mean",  "var",       "suscept"};
const string NewTypeSpecificSubmodel::s_singular_names [NUM_OF_TYPES] = { "mean",  "variance",  "susceptibility"};
const string NewTypeSpecificSubmodel::s_plural_names   [NUM_OF_TYPES] = { "means", "variances", "susceptibilities"};

const bool NewTypeSpecificSubmodel::s_type_constraint_allowed[NUM_OF_TYPES][NUM_OF_OPTIONS] = { 
  { true, true, true, true, true, true, true,  true  },
  { true, true, true, true, true, true, false, false },
  { true, true, true, true, true, true, true,  true  } };

TypeConstraint::OptionEnum NewTypeSpecificSubmodel::s_default_two_treatment = TypeConstraint::TWO_DOM;

//============================================================================
//  clone()
//============================================================================
MAXFUN::NewSubmodelShPtr
NewTypeSpecificSubmodel::clone()
{
  return MAXFUN::NewSubmodelShPtr(new NewTypeSpecificSubmodel(*this));
}

//=============================================
//
//  Constructor
//
//=============================================
NewTypeSpecificSubmodel::NewTypeSpecificSubmodel(CategoryEnum c, cerrorstream & errors)  
    : 
    NewSubmodel        (s_singular_names[c], errors),
    my_type_constraint (TypeConstraint::ONE)
{
  setCategory(c);

  my_types[0] = QNAN;
  my_types[1] = QNAN;
  my_types[2] = QNAN;

  my_fixed[0] = false;
  my_fixed[1] = false;
  my_fixed[2] = false;

  my_fixed_set[0] = false;
  my_fixed_set[1] = false;
  my_fixed_set[2] = false;
}

//===============================================
//
//  Copy constructor
//
//===============================================
NewTypeSpecificSubmodel::NewTypeSpecificSubmodel(const NewTypeSpecificSubmodel& other)
    : MAXFUN::NewSubmodel (other),
      my_category         (other.my_category),
      my_type_constraint  (other.my_type_constraint),
      my_treat_two_as     (other.my_treat_two_as)
{
  for(int i = 0; i < NUM_OF_TYPES; ++i)
  {
    my_types     [i] = other.my_types     [i];
    my_fixed     [i] = other.my_fixed     [i];
    my_fixed_set [i] = other.my_fixed_set [i];
  }
}

//===============================================
//
//  operator=(...)
//
//===============================================
NewTypeSpecificSubmodel&
NewTypeSpecificSubmodel::operator=(const NewTypeSpecificSubmodel& other)
{
  if(this != &other)
  {
    MAXFUN::NewSubmodel::operator=(other);
    
    my_category        = other.my_category;
    my_type_constraint = other.my_type_constraint;
    my_treat_two_as    = other.my_treat_two_as;
    
    for(int i = 0; i < NUM_OF_TYPES; ++i)
    {
      my_types     [i] = other.my_types     [i];
      my_fixed     [i] = other.my_fixed     [i];
      my_fixed_set [i] = other.my_fixed_set [i];
    }
  }
  
  return *this;
}

//==============================================
//
//  dump()
//
//==============================================
void
NewTypeSpecificSubmodel::dump(std::ostream& out) const
{
  OUTPUT::Table t(getSingularName());
  
  t << (OUTPUT::TableRow() << ("type_" + getBriefName()))
    << (OUTPUT::TableRow() << "option" 
                           << getTypeConstraint().toString() 
                           << (getTypeConstraint().getOption() == TypeConstraint::TWO ? "# interpreted as " + TypeConstraint(treatThisTwoAs()).toString() : ""))
    << (OUTPUT::TableRow() << getBriefName() << "=AA, val=" << my_types[0] << ", fixed=" << my_fixed[0])
    << (OUTPUT::TableRow() << getBriefName() << "=AB, val=" << my_types[1] << ", fixed=" << my_fixed[1])
    << (OUTPUT::TableRow() << getBriefName() << "=BB, val=" << my_types[2] << ", fixed=" << my_fixed[2]);
    
  out << t;
}

//===============================================================
//
//  convertConfigurationToString()
//
//===============================================================
string 
NewTypeSpecificSubmodel::convertConfigurationToString() const
{
  string type;                        

  type = my_type_constraint.getTypeCountAsString() + " " +
         (my_type_constraint.getTypeCount() == 1 ? getSingularName() : getPluralName()) + 
         (my_type_constraint.hasTypeRelation() ? ", " + my_type_constraint.getTypeRelationAsString() : "");
         
  if(!type.size()) 
    SAGE_internal_error();

  return type;
}

//===============================================================
//
//  setTypeValues(...)
//
//===============================================================
bool
NewTypeSpecificSubmodel::setTypeValues(
  size_t i, 
  double initial_val, 
  bool   fixed_set,
  bool   fixedness)
{
  // Set the value and the fixity:
  my_types     [i] = initial_val;
  my_fixed     [i] = fixedness;
  my_fixed_set [i] = fixed_set;

  // Return success!
  return true;
}


//==============================================================
//
//  isInitialSetupConsistent()
//
//==============================================================
int 
NewTypeSpecificSubmodel::isInitialSetupConsistent(int & warning_code) const
{
  // Grab necessary values for analysis:

  double AA                  = my_types[0],
         AB                  = my_types[1],
         BB                  = my_types[2];
  bool   AA_set              = !SAGE::isnan(AA),
         AB_set              = !SAGE::isnan(AB),
         BB_set              = !SAGE::isnan(BB),

         no_values_set       = !AA_set && !AB_set && !BB_set,

         one_value_set       = ( AA_set && !AB_set && !BB_set) ||
                               (!AA_set &&  AB_set && !BB_set) ||
                               (!AA_set && !AB_set &&  BB_set),

         two_values_set      = ( AA_set &&  AB_set && !BB_set) ||
                               ( AA_set && !AB_set &&  BB_set) ||
                               (!AA_set &&  AB_set &&  BB_set),
 
         three_values_set    = ( AA_set &&  AB_set &&  BB_set),
 
         at_least_two_values_set = two_values_set || three_values_set,

         all_set_values_same = one_value_set 

                               ||

                               (( AA_set &&  AB_set && !BB_set) && (AA == AB) ||
                                ( AA_set && !AB_set &&  BB_set) && (AA == BB) ||
                                (!AA_set &&  AB_set &&  BB_set) && (AB == BB)) 

                               ||

                               (three_values_set && AA == AB && AA == BB),

         all_set_values_unique = three_values_set && (AA != AB) && (AB != BB) && (AA != BB);

  double any_value = AA_set ? AA : (AB_set ? AB : (BB_set ? BB : QNAN));

  //===================================================================
  // If all preset values are QNAN, then it's all ok, because the 
  // analysis will figure out valid initial values later on.
  //===================================================================

  if(no_values_set)
  {
    return 0;
  }

  //===================================================================
  // If this is a VARIANCE submodel, and any value is <= 0, return an 
  // error!
  //===================================================================

  else if((getCategory() == VARIANCE) && ((AA_set && AA <= s_var_epsilon + s_var_lower_bound) || 
                                          (AB_set && AB <= s_var_epsilon + s_var_lower_bound) || 
                                          (BB_set && BB <= s_var_epsilon + s_var_lower_bound)))
  {
    return 200;
  }

  //===================================================================
  // If option is ONE...
  //===================================================================

  else if(getTypeConstraint().getTypeCount() == 1)
  {
    if(all_set_values_same)
    {
      my_types[0] = my_types[1] = my_types[2] = any_value;
    }
    else
    {
      return 1;
    }
  }

  //===================================================================
  // If option is TWO...
  //===================================================================

  else if(getTypeConstraint().getTypeCount() == 2)
  {
    if(at_least_two_values_set)
    {
      if(getTypeConstraint().getOption() == TypeConstraint::TWO)
      {
        if(two_values_set)
        {
          if(AA_set && AB_set && AA == AB)
          {
            return 12;
          }
          else if(AA_set && AB_set && AA != AB)
          {
            setTreatThisTwoAs(TypeConstraint::TWO_REC);

            my_types[2] = my_types[1];
          }
          else if(AA_set && BB_set && AA == BB)
          {
            setTreatThisTwoAs(getDefaultTwoTreatment());

            my_types[1] = my_types[0];
          }
          else if(AA_set && BB_set && AA != BB)
          {
            setTreatThisTwoAs(getDefaultTwoTreatment());

            my_types[1] = treatThisTwoAs() == TypeConstraint::TWO_DOM ? my_types[0] : my_types[1];

            warning_code = 1;
          }
          else if(AB_set && BB_set && AB == BB)
          {
            return 13;
          }
          else if(AB_set && BB_set && AB != BB)
          {
            setTreatThisTwoAs(TypeConstraint::TWO_DOM);

            my_types[0] = my_types[1];

            warning_code = 2;
          }
        }
        else if(three_values_set && all_set_values_same)
        {
          setTreatThisTwoAs(getDefaultTwoTreatment());
        }
        else if(three_values_set && all_set_values_unique)
        {
          return 11;
        }
        else if(three_values_set) // We can assume there are only two unique values set
        {
          if(AA == AB)
          {
            setTreatThisTwoAs(TypeConstraint::TWO_DOM);
          }
          else if(AA == BB)
          {
            return 14;
          }
          else if(AB == BB)
          {
            setTreatThisTwoAs(TypeConstraint::TWO_REC);
          }
        }
      }
      else if(getTypeConstraint().getOption() == TypeConstraint::TWO_DOM)
      {
        if(AA_set && AB_set && (AA != AB))
        {
          return 2;
        }
        else if(!BB_set)
        {
          return 3;
        }
        else
        {
          AA_set ? my_types[1] = my_types[0] : my_types[0] = my_types[1];
        }
      }
      else // option == TWO_REC
      {
        if(AB_set && BB_set && (AB != BB))
        {
          return 4;
        }
        else if(!AA_set)
        {
          return 5;
        }
        else
        {
          AB_set ? my_types[2] = my_types[1] : my_types[1] = my_types[2];
        }
      }
    }
    else // !at_least_two_values_set
    {
      return 6;
    }
  }

  //====================================================================
  // If option is THREE...
  //====================================================================

  else if(getTypeConstraint().getTypeCount() == 3)
  {
    if(three_values_set)
    {
      if((getTypeConstraint().getOption() == TypeConstraint::THREE_ADD) && !isThreeAddConstraintMet())
      {
        return 7;
      }
      else if((getTypeConstraint().getOption() == TypeConstraint::THREE_DEC) && !isThreeDecConstraintMet())
      {
        return 8;
      }
      else if((getTypeConstraint().getOption() == TypeConstraint::THREE_INC) && !isThreeIncConstraintMet())
      {
        return 9;
      }
    }
    else // !all_three_values_set
    {
      return 10;
    }
  }

  // Evaluate fixedness stuff:

  bool AA_fixed = my_fixed[0],
       AB_fixed = my_fixed[1],
       BB_fixed = my_fixed[2],

       AA_fixed_set = my_fixed_set[0],
       AB_fixed_set = my_fixed_set[1],
       BB_fixed_set = my_fixed_set[2],

       one_fixed_set = ( AA_fixed_set && !AB_fixed_set && !BB_fixed_set) ||
                       (!AA_fixed_set &&  AB_fixed_set && !BB_fixed_set) ||
                       (!AA_fixed_set && !AB_fixed_set &&  BB_fixed_set),

       two_fixed_set = ( AA_fixed_set &&  AB_fixed_set && !BB_fixed_set) ||
                       ( AA_fixed_set && !AB_fixed_set &&  BB_fixed_set) ||
                       (!AA_fixed_set &&  AB_fixed_set &&  BB_fixed_set),

       three_fixed_set = AA_fixed_set && AB_fixed_set && BB_fixed_set,

       at_least_one_fixed_set = one_fixed_set || two_fixed_set || three_fixed_set,

       any_fixed = AA_fixed_set ? AA_fixed : (AB_fixed_set ? AB_fixed : BB_fixed),

       all_set_fixed_same = one_fixed_set

                             ||

                             (( AA_fixed_set &&  AB_fixed_set && !BB_fixed_set) && (AA_fixed == AB_fixed) ||
                              ( AA_fixed_set && !AB_fixed_set &&  BB_fixed_set) && (AA_fixed == BB_fixed) ||
                              (!AA_fixed_set &&  AB_fixed_set &&  BB_fixed_set) && (AB_fixed == BB_fixed))

                             ||

                             (three_values_set && AA_fixed == AB_fixed && AA_fixed == BB_fixed);

  if(!at_least_one_fixed_set)
  {
    ;
  }
  else if(getTypeConstraint().getTypeCount() == 1)
  {
    if(all_set_fixed_same)
    {
      my_fixed[0] = my_fixed[1] = my_fixed[2] = any_fixed;
    }
    else
    {
      return 100;
    }
  }
  else if(getTypeConstraint().getTypeCount() == 2)
  {
    if(getTypeConstraint().getOption() == TypeConstraint::TWO_DOM || 
       (getTypeConstraint().getOption() == TypeConstraint::TWO && treatThisTwoAs() == TypeConstraint::TWO_DOM))
    {
      if(AA_fixed_set && AB_fixed_set && (AA_fixed != AB_fixed))
      {
        return 101;
      }
      else if(AA_fixed_set && !AB_fixed_set)
      {
        my_fixed[1] = my_fixed[0];
      }
      else if(!AA_fixed_set && AB_fixed_set)
      {
        my_fixed[0] = my_fixed[1];
      }
    }
    else if(getTypeConstraint().getOption() == TypeConstraint::TWO_REC || 
       (getTypeConstraint().getOption() == TypeConstraint::TWO && treatThisTwoAs() == TypeConstraint::TWO_REC))
    {
      if(AB_fixed_set && BB_fixed_set && (AB_fixed != BB_fixed))
      {
        return 102;
      }
      else if(AB_fixed_set && !BB_fixed_set)
      {
        my_fixed[2] = my_fixed[1];
      }
      else if(!AB_fixed_set && BB_fixed_set)
      {
        my_fixed[1] = my_fixed[2];
      }
    }
  }

  // Return success:

  return 0;
}

//==============================================================
//
//  areMeanAndVarianceCompatible(...)
//
//==============================================================
int 
NewTypeSpecificSubmodel::areMeanAndVarianceCompatible(const NewTypeSpecificSubmodel & mean_submodel, const NewTypeSpecificSubmodel & var_submodel)
{
  if(var_submodel.getTypeConstraint().getTypeCount() > mean_submodel.getTypeConstraint().getTypeCount())
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

//===============================================================
//
//  finalizeConfiguration()
//
//===============================================================
int
NewTypeSpecificSubmodel::finalizeConfiguration()
{
  switch(my_type_constraint.getOption())
  {
    case TypeConstraint::ONE:       initializeOne      ();  break;
    case TypeConstraint::TWO:       treatThisTwoAs() == TypeConstraint::TWO_DOM ?  
                                    initializeTwoDom   () : 
                                    initializeTwoRec   ();  break;
    case TypeConstraint::THREE:     initializeThree    ();  break;
    case TypeConstraint::TWO_DOM:   initializeTwoDom   ();  break;
    case TypeConstraint::TWO_REC:   initializeTwoRec   ();  break;
    case TypeConstraint::THREE_ADD: initializeThreeAdd ();  break;
    case TypeConstraint::THREE_DEC:
    case TypeConstraint::THREE_INC:

      // For threeDec and threeInc, they have the same values as the basic
      // three model, but have a functional relationship to one another.

      initializeThree();
      
      // Redefine any non-fixed as independent functional.
      
      my_parameters[0].initial_type = my_fixed[0] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL;
      my_parameters[1].initial_type = my_fixed[1] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL;
      my_parameters[2].initial_type = my_fixed[2] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL;

      break;
  }

  return 0;
}

//===============================================================
//
//  initializeOne()
//
//===============================================================
void
NewTypeSpecificSubmodel::initializeOne()
{
  my_parameters.resize(1);

  my_parameters[0] = MAXFUN::ParameterInput(
      toUpper(getPluralName()),
      getBriefName(),
      my_fixed[0] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
      my_types[0],
      my_category == VARIANCE ? s_var_lower_bound + s_var_epsilon : NEGATIVE_INF,
      POSITIVE_INF);
}

//===============================================================
//
//  initializeTwoDom()
//
//===============================================================
void
NewTypeSpecificSubmodel::initializeTwoDom()
{
  my_parameters.resize(2);

  my_parameters[0] = MAXFUN::ParameterInput
    ( toUpper(getPluralName()),
      getBriefName() + "_AA_AB",
      my_fixed[0] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
      my_types[0],
      my_category == VARIANCE ? s_var_lower_bound + s_var_epsilon : NEGATIVE_INF,
      POSITIVE_INF);

  my_parameters[1] = MAXFUN::ParameterInput
    ( toUpper(getPluralName()),
      getBriefName() + "_BB",
      my_fixed[2] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
      my_types[2],
      my_category == VARIANCE ? s_var_lower_bound + s_var_epsilon : NEGATIVE_INF,
      POSITIVE_INF);
}

//===============================================================
//
//  initializeTwoRec()
//
//===============================================================
void
NewTypeSpecificSubmodel::initializeTwoRec()
{ 
  my_parameters.resize(2);

  my_parameters[0] = MAXFUN::ParameterInput
    ( toUpper(getPluralName()),
      getBriefName() + "_AA",
      my_fixed[0] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
      my_types[0],
      my_category == VARIANCE ? s_var_lower_bound + s_var_epsilon : NEGATIVE_INF,
      POSITIVE_INF);

  my_parameters[1] = MAXFUN::ParameterInput
    ( toUpper(getPluralName()),
      getBriefName() + "_AB_BB",
      my_fixed[2] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
      my_types[2],
      my_category == VARIANCE ? s_var_lower_bound + s_var_epsilon : NEGATIVE_INF,
      POSITIVE_INF);
}


//===============================================================
//
//  initializeThreeAdd()
//
//===============================================================
void
NewTypeSpecificSubmodel::initializeThreeAdd()
{
  my_parameters.resize(3);

  my_parameters[0] = MAXFUN::ParameterInput
    ( toUpper(getPluralName()),
      getBriefName() + "_AA",
      my_fixed[0] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL,
      my_types[0],
      my_category == VARIANCE ? s_var_lower_bound + s_var_epsilon : NEGATIVE_INF,
      POSITIVE_INF);

  my_parameters[1] = MAXFUN::ParameterInput
    ( toUpper(getPluralName()),
      getBriefName() + "_AB",
      my_fixed[1] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::DEPENDENT,
      my_types[1],
      my_category == VARIANCE ? s_var_lower_bound + s_var_epsilon : NEGATIVE_INF,
      POSITIVE_INF);

  my_parameters[2] = MAXFUN::ParameterInput
    ( toUpper(getPluralName()),
      getBriefName() + "_BB",
      my_fixed[2] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL,
      my_types[2],
      my_category == VARIANCE ? s_var_lower_bound + s_var_epsilon : NEGATIVE_INF,
      POSITIVE_INF);
}

//===============================================================
//
//  initializeThree()
//
//===============================================================
void
NewTypeSpecificSubmodel::initializeThree()
{
  my_parameters.resize(3);

  my_parameters[0] = MAXFUN::ParameterInput
    ( toUpper(getPluralName()),
      getBriefName() + "_AA",
      my_fixed[0] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
      my_types[0],
      my_category == VARIANCE ? s_var_lower_bound + s_var_epsilon : NEGATIVE_INF,
      POSITIVE_INF);

  my_parameters[1] = MAXFUN::ParameterInput
    ( toUpper(getPluralName()),
      getBriefName() + "_AB",
      my_fixed[1] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
      my_types[1],
      my_category == VARIANCE ? s_var_lower_bound + s_var_epsilon : NEGATIVE_INF,
      POSITIVE_INF);

  my_parameters[2] = MAXFUN::ParameterInput
    ( toUpper(getPluralName()),
      getBriefName() + "_BB",
      my_fixed[2] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
      my_types[2],
      my_category == VARIANCE ? s_var_lower_bound + s_var_epsilon : NEGATIVE_INF,
      POSITIVE_INF);
}

//=================================================================
//
//  update()
//
//=================================================================
int
NewTypeSpecificSubmodel::update()
{
  int return_value = 1;

  switch(my_type_constraint.getOption())
  {
    case TypeConstraint::ONE       : return_value = synchronizeOne      (); break;
    case TypeConstraint::TWO       : return_value = synchronizeTwo      (); break;
    case TypeConstraint::THREE     : return_value = synchronizeThree    (); break;
    case TypeConstraint::TWO_DOM   : return_value = synchronizeTwoDom   (); break;
    case TypeConstraint::TWO_REC   : return_value = synchronizeTwoRec   (); break;
    case TypeConstraint::THREE_ADD : return_value = synchronizeThreeAdd (); break;
    case TypeConstraint::THREE_DEC : return_value = synchronizeThreeDec (); break;
    case TypeConstraint::THREE_INC : return_value = synchronizeThreeInc (); break;
    default                        : SAGE_internal_error();
  }
   
  return return_value;
}

//===============================================================
//
//  synchronizeOne()
//
//===============================================================
int
NewTypeSpecificSubmodel::synchronizeOne()
{
  my_types[0] = getParam(0);
  my_types[1] = getParam(0);
  my_types[2] = getParam(0);

  return 0;
}

//===============================================================
//
//  synchronizeTwo()
//
//===============================================================
int  
NewTypeSpecificSubmodel::synchronizeTwo()
{
  return getDefaultTwoTreatment() == TypeConstraint::TWO_DOM ? synchronizeTwoDom() : synchronizeTwoRec();
}
 
//===============================================================
//
//  synchronizeThree()
//
//===============================================================
int  
NewTypeSpecificSubmodel::synchronizeThree()
{
  my_types[0] = getParam(0);
  my_types[1] = getParam(1);
  my_types[2] = getParam(2);

  return 0;
}

//===============================================================
//
//  synchronizeTwoDom()
//
//===============================================================
int  
NewTypeSpecificSubmodel::synchronizeTwoDom()
{
  my_types[0] = getParam(0);
  my_types[1] = getParam(0);
  my_types[2] = getParam(1);

  return 0;
}
 
//===============================================================
//
//  synchronizeTwoRec()
//
//===============================================================
int  
NewTypeSpecificSubmodel::synchronizeTwoRec()
{
  my_types[0] = getParam(0);
  my_types[1] = getParam(1);
  my_types[2] = getParam(1);

  return 0;
}

//===============================================================
//
//  synchronizeThreeAdd()
//
//===============================================================
int  
NewTypeSpecificSubmodel::synchronizeThreeAdd()
{
  my_types[0] = getParam(0);
  my_types[2] = getParam(1);

  // For this option one parameter is a function of the other two.

  my_types[1] = (my_types[0] + my_types[2]) / 2;
  
  // Update maxfun w. newly calculated value of 'AB' parameter.

  getParam(1) = my_types[1];
  
  return 0;
}
 
//==============================================================
//
//  isThreeAddConstraintMet()
//
//==============================================================
bool
NewTypeSpecificSubmodel::isThreeAddConstraintMet() const
{
  return my_types[1] == (my_types[0] + my_types[2]) / 2;
}

// - Values of genotype specific means must be in decreasing order. 
//   Comparisons involving a QNAN always return true.
//
bool
NewTypeSpecificSubmodel::isThreeDecConstraintMet() const
{
  bool  AA_or_AB_nan = SAGE::isnan(my_types[0]) || SAGE::isnan(my_types[1]);
  bool  AB_or_BB_nan = SAGE::isnan(my_types[1]) || SAGE::isnan(my_types[2]);
  bool  AA_or_BB_nan = SAGE::isnan(my_types[0]) || SAGE::isnan(my_types[2]);
  
  bool  AA_gte_AB = AA_or_AB_nan || my_types[0] >= my_types[1];
  bool  AB_gte_BB = AB_or_BB_nan || my_types[1] >= my_types[2];
  bool  AA_gte_BB = AA_or_BB_nan || my_types[0] >= my_types[2];

  if(AA_gte_AB && AB_gte_BB && AA_gte_BB)
  {
    return true;
  }
  else
  {
    return false;
  }
}

// - Values of genotype specific means must be in increasing order.
//   Comparisons involving a QNAN always return true.
//
bool
NewTypeSpecificSubmodel::isThreeIncConstraintMet() const
{
  bool  AA_or_AB_nan = SAGE::isnan(my_types[0]) || SAGE::isnan(my_types[1]);
  bool  AB_or_BB_nan = SAGE::isnan(my_types[1]) || SAGE::isnan(my_types[2]);
  bool  AA_or_BB_nan = SAGE::isnan(my_types[0]) || SAGE::isnan(my_types[2]);
  
  bool  AA_lte_AB = AA_or_AB_nan || my_types[0] <= my_types[1];
  bool  AB_lte_BB = AB_or_BB_nan || my_types[1] <= my_types[2];
  bool  AA_lte_BB = AA_or_BB_nan || my_types[0] <= my_types[2];

  if(AA_lte_AB && AB_lte_BB && AA_lte_BB)
  {
    return true;
  }
  else
  {
    return false;
  }
}

//===============================================================
//
//  synchronizeThreeDec()
//
//===============================================================
int NewTypeSpecificSubmodel::synchronizeThreeDec()
{
  my_types[0] = getParam(0);
  my_types[1] = getParam(1);
  my_types[2] = getParam(2);

  if(isThreeDecConstraintMet())
  {
    return 0;
  }
  else
  {
    // - sub-model will now be in an inconsistent state, but this will be 
    //   corrected on the next successful synchronization.  In the interim
    //   maxfun will not call evaluate() so we do not have a problem 
    //   (per discussion w. gcw).
    //
    return 1;    
  }
}
 
//==================================================================
//
//  synchronizeThreeInc()
//
//==================================================================
int NewTypeSpecificSubmodel::synchronizeThreeInc()
{
  my_types[0] = getParam(0);
  my_types[1] = getParam(1);
  my_types[2] = getParam(2);

  if(isThreeIncConstraintMet())
  {
    return 0;
  }
  else
  {
    return 1;   
  }
}












} // End namespace MFSUBMODELS
} // End namespace SAGE
