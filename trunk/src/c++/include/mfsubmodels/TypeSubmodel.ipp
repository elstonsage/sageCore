#ifndef MFSUBMODELS_TYPESUBMODEL_H
#include "TypeSubmodel.h"
#endif

namespace SAGE        {
namespace MFSUBMODELS {

inline string TypeSpecificSubmodel::getBriefName    ()         const { return s_brief_names    [getCategory()];                                  }
inline string TypeSpecificSubmodel::getSingularName ()         const { return s_singular_names [getCategory()];                                  }
inline string TypeSpecificSubmodel::getPluralName   ()         const { return s_plural_names   [getCategory()];                                  }
inline double TypeSpecificSubmodel::getValue        (size_t i) const { return my_types[i];                                                       }
inline bool   TypeSpecificSubmodel::is_complete     ()         const { return finite(getValue(0)) && finite(getValue(1)) && finite(getValue(2)); }

inline bool TypeSpecificSubmodel::isConstraintTypeAllowed (TypeConstraint::OptionEnum e) const { return s_type_constraint_allowed[getCategory()][e]; }

inline       TypeSpecificSubmodel::CategoryEnum TypeSpecificSubmodel::getCategory            () const { return my_category;             }
inline       TypeConstraint::OptionEnum         TypeSpecificSubmodel::getDefaultTwoTreatment ()       { return s_default_two_treatment; }
inline       TypeConstraint::OptionEnum         TypeSpecificSubmodel::treatThisTwoAs         () const { return my_treat_two_as;         }
inline const TypeConstraint&                    TypeSpecificSubmodel::getTypeConstraint      () const { return my_type_constraint;      }



//===============================================================
//  setDefaultTwoTreatment(...)
//===============================================================
inline bool TypeSpecificSubmodel::setDefaultTwoTreatment(TypeConstraint::OptionEnum e)
{
  if(e != TypeConstraint::TWO_DOM && e != TypeConstraint::TWO_REC)
    return false;

  s_default_two_treatment = e;

  return true;
}

//===============================================================
//  setTreatThisTwoAs(...)
//===============================================================
inline bool TypeSpecificSubmodel::setTreatThisTwoAs(TypeConstraint::OptionEnum e) const
{
  if(e != TypeConstraint::TWO_DOM && e != TypeConstraint::TWO_REC)
    return false;

  my_treat_two_as = e;

  return true;
}

//==============================================================
//  setTypeConstraint(...) #1
//==============================================================
inline bool TypeSpecificSubmodel::setTypeConstraint(const TypeConstraint & constraint)
{
  if(isConstraintTypeAllowed(constraint.getOption()) == true)
  {
    my_type_constraint = constraint;

    return true;
  }
  else
  {
    return false;
  }
}

//==============================================================
//  setTypeConstraint(...) #2
//==============================================================
inline bool TypeSpecificSubmodel::setTypeConstraint(const string & t)
{
       if(t == "ONE")       return setTypeConstraint(TypeConstraint::ONE);
  else if(t == "TWO")       return setTypeConstraint(TypeConstraint::TWO);
  else if(t == "THREE")     return setTypeConstraint(TypeConstraint::THREE);
  else if(t == "TWO_DOM")   return setTypeConstraint(TypeConstraint::TWO_DOM);
  else if(t == "TWO_REC")   return setTypeConstraint(TypeConstraint::TWO_REC);
  else if(t == "THREE_ADD") return setTypeConstraint(TypeConstraint::THREE_ADD);
  else if(t == "THREE_DEC") return setTypeConstraint(TypeConstraint::THREE_DEC);
  else if(t == "THREE_INC") return setTypeConstraint(TypeConstraint::THREE_INC);
  else                      return false;
}


//========================================================================================================






inline string NewTypeSpecificSubmodel::getBriefName    ()         const { return s_brief_names    [getCategory()];                                  }
inline string NewTypeSpecificSubmodel::getSingularName ()         const { return s_singular_names [getCategory()];                                  }
inline string NewTypeSpecificSubmodel::getPluralName   ()         const { return s_plural_names   [getCategory()];                                  }
inline double NewTypeSpecificSubmodel::getValue        (size_t i) const { return my_types[i];                                                       }
inline bool   NewTypeSpecificSubmodel::is_complete     ()         const { return finite(getValue(0)) && finite(getValue(1)) && finite(getValue(2)); }

inline bool NewTypeSpecificSubmodel::isConstraintTypeAllowed (TypeConstraint::OptionEnum e) const { return s_type_constraint_allowed[getCategory()][e]; }

inline       NewTypeSpecificSubmodel::CategoryEnum NewTypeSpecificSubmodel::getCategory            () const { return my_category;             }
inline       TypeConstraint::OptionEnum         NewTypeSpecificSubmodel::getDefaultTwoTreatment ()       { return s_default_two_treatment; }
inline       TypeConstraint::OptionEnum         NewTypeSpecificSubmodel::treatThisTwoAs         () const { return my_treat_two_as;         }
inline const TypeConstraint&                    NewTypeSpecificSubmodel::getTypeConstraint      () const { return my_type_constraint;      }

inline bool NewTypeSpecificSubmodel::setCategory(CategoryEnum c) 
{ 
  my_category = c; 

  setName(s_singular_names[my_category]);

  return true; 
}

//===============================================================
//  setDefaultTwoTreatment(...)
//===============================================================
inline bool NewTypeSpecificSubmodel::setDefaultTwoTreatment(TypeConstraint::OptionEnum e)
{
  if(e != TypeConstraint::TWO_DOM && e != TypeConstraint::TWO_REC)
    return false;

  s_default_two_treatment = e;

  return true;
}

//===============================================================
//  setTreatThisTwoAs(...)
//===============================================================
inline bool NewTypeSpecificSubmodel::setTreatThisTwoAs(TypeConstraint::OptionEnum e) const
{
  if(e != TypeConstraint::TWO_DOM && e != TypeConstraint::TWO_REC)
    return false;

  my_treat_two_as = e;

  return true;
}

//==============================================================
//  setTypeConstraint(...) #1
//==============================================================
inline bool NewTypeSpecificSubmodel::setTypeConstraint(const TypeConstraint & constraint)
{
  if(isConstraintTypeAllowed(constraint.getOption()) == true)
  {
    my_type_constraint = constraint;

    return true;
  }
  else
  {
    return false;
  }
}

//==============================================================
//  setTypeConstraint(...) #2
//==============================================================
inline bool NewTypeSpecificSubmodel::setTypeConstraint(const string & t)
{
       if(t == "ONE")       return setTypeConstraint(TypeConstraint::ONE);
  else if(t == "TWO")       return setTypeConstraint(TypeConstraint::TWO);
  else if(t == "THREE")     return setTypeConstraint(TypeConstraint::THREE);
  else if(t == "TWO_DOM")   return setTypeConstraint(TypeConstraint::TWO_DOM);
  else if(t == "TWO_REC")   return setTypeConstraint(TypeConstraint::TWO_REC);
  else if(t == "THREE_ADD") return setTypeConstraint(TypeConstraint::THREE_ADD);
  else if(t == "THREE_DEC") return setTypeConstraint(TypeConstraint::THREE_DEC);
  else if(t == "THREE_INC") return setTypeConstraint(TypeConstraint::THREE_INC);
  else                      return false;
}








// - For debugging.
//
inline std::ostream&  
operator<<(std::ostream& out, const model_input& mi)
{
  out << "value : " << mi.value << "\n"
      << std::boolalpha
      << "fixed:  " << mi.fixed << std::endl;
  
  return out;
}

} // End namespace MFSUBMODELS
} // End namespace SAGE
