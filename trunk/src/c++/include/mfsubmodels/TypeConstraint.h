#ifndef MFSUBMODELS_TYPECONSTRAINT_H
#define MFSUBMODELS_TYPECONSTRAINT_H

namespace SAGE        {
namespace MFSUBMODELS {

const int NUM_OF_OPTIONS = 8;

/** \brief Describes a type-based constraint
  *
  * Although this class is primarily intended to describe constraints between \b genotypes,
  * it does not have to be used for that purpose.
  */
class TypeConstraint
{
  public:

  ///
  /// Describes a set of constraints on genotypes.
  /// Note that some options are only available for certain sub_models. 
  /// Specifically, three_dec and three_inc are only valid for mean and
  /// susceptibility models.
  ///
  /// Please note that NUM_OF_OPTIONS \b must be set to the number of values
  /// in this enum.
  enum OptionEnum
  {
    ONE = 0,       ///< One mean/var/susc.
    TWO,           ///< Two means/vars/suscs
    THREE,         ///< Three means/vars/suscs
    TWO_DOM,       ///< Two means/vars/suscs with uAB = uAA
    TWO_REC,       ///< Two means/vars/suscs with uAB = uBB
    THREE_ADD,     ///< Three means/vars/suscs with uAB = (uAA + uBB) / 2
    THREE_DEC,     ///< Three means/vars/suscs with uAA > uAB > uBB
    THREE_INC      ///< Three means/vars/suscs with uAA < uAB < uBB
  };

  /// @name Constructors & operator
  //@{
  
    ///
    /// Basic constructor. Value defaults to TypeConstraint::ONE.
    TypeConstraint();
    
    ///
    /// Another constructor. Initializes its ConstraintType to the given parameter.
    /// \param c The initial ConstraintType.
    TypeConstraint(OptionEnum c);
    
    ///
    /// Copy constructor.
    /// \param other The object to copy
    TypeConstraint(const TypeConstraint & other);
  
    ///
    /// Assignment operator.
    /// \param other The object to copy
    TypeConstraint& operator=(const TypeConstraint & other);
  
  //@}

  /// @name Basic functionality
  //@{
  
    ///
    /// Returns this object's ConstraintType (see OptionEnum).
    const OptionEnum & getOption() const;
    
    ///
    /// Sets this object's Option (see OptionEnum).
    /// \param c The new Option
    void setOption(const OptionEnum & c);

    ///
    /// Sets this object's Option from a string.
    /// Valid values are: "ONE", "TWO", "THREE", "TWO_DOM", "TWO_REC",
    /// "THREE_ADD", "THREE_DEC", "THREE_INC".
    /// \param t The string from which to set the Option
    /// \retval true The Option was set successfully
    /// \retval false The Option was not set successfully
    bool setOption(const string & t);

    ///
    /// Given a TypeConstraint, reports the number of genotypes.
    /// \param option The TypeConstraintEunm in question
    size_t getTypeCount() const;

    ///
    /// Given a TypeConstraint, returns \c true if there is a preset relation between the types
    /// (for two_dom, two_rec, three_add, three_dec, and three_inc), \false otherwise.
    /// \param option The TypeConstraintEunm in question
    bool hasTypeRelation() const;

  //@}

  /// @name String extraction
  //@{
  
    ///
    /// Converts this object's Option to a string.
    /// (Ie: "TWO_REC", "THREE_INC")
    /// \param option The TypeConstraintEunm in question
    string toString() const;

    ///
    /// Reports the number of genotypes as a string.
    /// (Ie: "one", "two", or "three"
    /// \param option The TypeConstraintEnum in question
    string getTypeCountAsString() const;
    
    ///
    /// Given a TypeConstraint, reports the relationship between those genotypes as a string.
    /// (Ie: "decreasing", "increasing")
    /// \param option The TypeConstraintEunm in question
    string getTypeRelationAsString() const;

  //@}

  private:

    OptionEnum my_option;
};

//============================
//  INLINE FUNCTIONS
//============================

inline
TypeConstraint::TypeConstraint() 
{
  my_option = ONE;
}
    
inline
TypeConstraint::TypeConstraint(OptionEnum c) 
{
  my_option = c;
}
    
inline
TypeConstraint::TypeConstraint(const TypeConstraint & other) 
{
  my_option = other.my_option;
}
  
inline
TypeConstraint& 
TypeConstraint::operator=(const TypeConstraint & other) 
{
  my_option = other.my_option;
  return *this;
}
  
inline const TypeConstraint::OptionEnum & 
TypeConstraint::getOption() const 
{ 
  return my_option; 
}

inline void 
TypeConstraint::setOption(const OptionEnum & c) 
{ 
  my_option = c; 
}

inline bool
TypeConstraint::setOption(const string & opt) 
{ 
       if(opt == "ONE")       { my_option = ONE;       return true; }
  else if(opt == "TWO")       { my_option = TWO;       return true; }
  else if(opt == "THREE")     { my_option = THREE;     return true; }
  else if(opt == "TWO_DOM")   { my_option = TWO_DOM;   return true; }
  else if(opt == "TWO_REC")   { my_option = TWO_REC;   return true; }
  else if(opt == "THREE_ADD") { my_option = THREE_ADD; return true; }
  else if(opt == "THREE_DEC") { my_option = THREE_DEC; return true; }
  else if(opt == "THREE_INC") { my_option = THREE_INC; return true; }
  else return false;
}

inline string 
TypeConstraint::toString() const
{
  switch(my_option)
  {
    case ONE       : return "one";      
    case TWO       : return "two";      
    case THREE     : return "three";    
    case TWO_DOM   : return "two_dom";  
    case TWO_REC   : return "two_rec";  
    case THREE_ADD : return "three_add";
    case THREE_DEC : return "three_dec";
    case THREE_INC : return "three_inc";
    default        : return "";
  }                                          
}

inline size_t 
TypeConstraint::getTypeCount() const 
{
  switch(my_option)
  {
    case ONE       : return 1;          
    case TWO       : return 2;
    case THREE     : return 3;
    case TWO_DOM   : return 2;
    case TWO_REC   : return 2;
    case THREE_ADD : return 3;
    case THREE_DEC : return 3;
    case THREE_INC : return 3;
    default        : return 0;
  }
}

inline string 
TypeConstraint::getTypeCountAsString() const 
{
  size_t cnt = getTypeCount();
  
       if(cnt == 1) return "one";
  else if(cnt == 2) return "two";
  else if(cnt == 3) return "three";
  else return "";
}
    
inline bool 
TypeConstraint::hasTypeRelation() const 
{
  switch(my_option)
  {
    case ONE       : return false;
    case TWO       : return false;
    case THREE     : return false;
    case TWO_DOM   : return true; 
    case TWO_REC   : return true; 
    case THREE_ADD : return true; 
    case THREE_DEC : return true; 
    case THREE_INC : return true; 
    default        : return false;
  }
}

inline string 
TypeConstraint::getTypeRelationAsString() const 
{
  switch(my_option)
  {
    case ONE       : return "";
    case TWO       : return "";
    case THREE     : return "";
    case TWO_DOM   : return "A dominant";
    case TWO_REC   : return "A recessive";
    case THREE_ADD : return "additive";   
    case THREE_DEC : return "decreasing"; 
    case THREE_INC : return "increasing";
    default        : return "";
  }
}

} // End namespace MFSUBMODELS
} // End namespace SAGE

#endif
