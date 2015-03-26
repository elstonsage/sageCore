#ifndef FUNC_EXPRESSION
#define FUNC_EXPRESSION

#include <vector>
#include <string>
#include <iostream>

namespace SAGE {
namespace FUNC {

class Expression
{
  public:
  
    typedef std::vector<std::string> VariableNames;

  /// @name Constructors
  //@{

    ///
    /// Constructor.
    Expression();
    
    ///
    /// Copy constructor.
    Expression(const Expression &);
    
    /// 
    /// Assignment operator.
    Expression& operator= (const Expression &);
    
  //@}

  /// @name Misc
  //@{

    ///
    /// Sets the python expression (ie: "x * y").
    /// Returns success / failure.
    bool setExpression(const std::string & expression);
    
    ///
    /// Returns the list of variables that will have to be set for the function to evaluate.
    const VariableNames & getVariableNames() const;
    
    ///
    /// Sets the named variable to the given value.
    /// Returns \c false if the given name is not in the VariableNames (given by getVariableNames() ).
    bool setVariable(const std::string & name, double value);
    
    ///
    /// Evaluates the expression.
    double evaluate() const;
    
    ///
    /// Debugging.
    void dump() const;

  //@}

  private:
  
    typedef std::vector<double> VariableValues;

    std::string    my_expression;
    VariableNames  my_variable_names;
    VariableValues my_variable_values;

};


} // End namespace FUNC
} // End namespace SAGE


#endif
