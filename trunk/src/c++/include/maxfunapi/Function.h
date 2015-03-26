#ifndef MAXFUN_FUNCTION_SCRIPT_H
#define MAXFUN_FUNCTION_SCRIPT_H

#include "globals/SAGEConstants.h"
#include "boost/shared_ptr.hpp"
#include "boost/lambda/lambda.hpp"
#include "functors/functors.h"
#include "maxfun/maxfun.h"
#include "maxfunapi/ParameterMgr.h"

namespace SAGE   {
namespace MAXFUN {

/// \brief Represents a function for maxfun to maximize.
///
/// A "Function", conceptually speaking, is a process that takes in input and produces output.
/// In this instance, the SAGE::MAXFUN::Function is customized for the purposes of function maximization.
///
/// \par Creating a Function
///
/// Since function maximization presupposes the presence (and organization) of parameters, a SAGE::MAXFUN::Function
/// must be bound to a SAGE::MAXFUN::ParameterMgr. This ParameterMgr must be independently instantiated before
/// the Function itself can be created. Once the ParameterMgr has been created, the Function can in turn be created
/// (and bound at the point of construction to the ParameterMgr).
///
/// For example:
///
/// \code
///
/// SAGE::MAXFUN::ParameterMgr m;
/// SAGE::MAXFUN::Function f(m);
///
/// \endcode
///
/// \par Function maximization
///
/// Function maximization is an iterative process: The maxfun API, having been given a ParameterMgr and a Function,
/// sets the "current" parameters estimates in the ParameterMgr to a set of values. It then queries the Function
/// to find out if those parameter estimates are acceptable (see "the Sequence" below). If the estimates are
/// acceptable, the API queries the Function to get the Function's final value (see "the final value" below). 
/// If the estimates are not acceptable, the API reverts to the last known "good" estimates are tries alternate possibilities.
///
/// \par Calculating a Function
///
/// The calculation of a Function has two distinct steps: the execution of a sequence of one or more steps, 
/// and the execution of a final, separate step that returns the function's value.
///
/// 1. The Sequence
///
/// A Function's internal "sequence" of steps consists of 1 or more functors that are executed one after the other.
/// Each functor is required to take as input a non-const ParameterMgr reference (the current parameter estimates,
/// that is) and return an integer indicating the exit status of that step. For the exit status, 0 indicates success,
/// non-zero indicates failure. When the sequence executes, it evalutes the exit status of each step. If any step
/// indicates failure, the sequence aborts and the maximization process is notified.
///
/// 2. The final value
///
/// The actual function value is given by a separate functor executed after the sequence is successfully run.
/// It must be a functor that takes a non-const ParameterMgr& and returns a double (the final function value).
///
class Function : public SAGE::MaxFunction
{
private:

  ///
  /// A functor that will assign the value of a parameter according to the double
  /// returned by a given functor.
  ///
  /// Please note that, in order to allow full use of lambda expressions, you cannot
  /// directly instantiate this class. Instead, use the makeParamSetter() function to
  /// do so. 
  class DependentParamSetter
  {
    public:
  
      typedef SAGE::FUNCTORS::make_function<double, boost::mpl::vector<ParameterMgr&> >::type FunctorType;
      
      DependentParamSetter();
    
      DependentParamSetter(int id, const FunctorType& f) : target_id(id), getVal(f) { }
    
      int operator()(ParameterMgr& mgr) const
      {
        mgr.getParameter(target_id).setCurrentEstimate(getVal(mgr));
        
        return 0;
      }
    
    private:
  
      int         target_id;
      FunctorType getVal;
  };

public:

  /// @name Typedefs
  //@{
  
    /// The sequence type.
    typedef SAGE::FUNCTORS::Sequence<int, boost::mpl::vector<ParameterMgr&> > SequenceType;

    /// The required type for the evaluator functor.
    typedef SAGE::FUNCTORS::make_function<double, boost::mpl::vector<ParameterMgr&> >::type EvaluatorFunctor;
    
  //@}

  /// @name Constructors
  //@{
  
    ///
    /// Constructor.
    /// \param mgr The ParameterMgr to which this Function is bound.
    explicit Function(MAXFUN::ParameterMgr& mgr) : my_mgr(mgr) 
    { 
      addStep(SAGE::FUNCTORS::ReturnConstValue<int, SAGE::MAXFUN::Function::SequenceType::args> (0), "[Default step]");
    }
    
    ///
    /// Copy constructor.
    Function(const Function & other) : my_mgr(other.my_mgr), my_seq(other.my_seq), my_eval(other.my_eval) { }
  
  //@}
    
  /// @name Required virtual interface
  //@{

    ///
    /// Required for the maxfunapi to work. Ignore this function!
    virtual double evaluate(vector<double> & params) { return my_eval(my_mgr); }

    ///
    /// Required for the maxfunapi to work. Ignore this function!
    virtual int update_bounds(vector<double>& params) { getMgr().update(params); return my_seq(my_mgr); }

  //@}

  /// @name Adding sequence steps / setting evaluation step:
  //@{
  
    ///
    /// Adds a step to the "update" process. Each step must be a functor that takes a ParameterMgr& and returns an error
    /// code (and int). The int returned by the step will be interpreted as an error if it is nonzero.
    void addStep(const SequenceType::FunctorType& step, const std::string& name)
    {
      my_seq.addStep(step, boost::lambda::_1 == 0, name);
    }
    
    ///
    /// Adds a step to the "update" process. Each step must be a functor that takes a ParameterMgr& and returns an error
    /// code (and int). The int returned by the step will be interpreted as an error if it is nonzero.
    template<typename F>
    void addDependentParamStep(const string& group_name, const string& param_name, const F& f)
    {
      int  id   = my_mgr.getParamID(group_name, param_name);
      
      string name = "Calculate " + group_name + ":" + param_name;
      addStep(DependentParamSetter(id, f), name);
    }
    
    template<typename F>
    void addDependentParamStep(size_t param_idx, const F& f)
    {
      Parameter&  param = my_mgr.getParameter(param_idx);
      
      string name = "Calculate " + param.getGroupName() + ":" + param.getName();
      addStep(DependentParamSetter(param_idx, f), name);
    }    
    
    ///
    /// Sets the "evaluator" function (the function that will return the Function's value). The evaluator must take
    /// a ParameterMgr& as input and return a double.
    void setEvaluator(const EvaluatorFunctor& eval)
    {
      my_eval = eval;
    }
  
  //@}

  /// @name Getting ParameterMgr:
  //@{
  
    ///
    /// Returns a const reference to the ParameterMgr to which this Function is bound.
    const ParameterMgr& getMgr() const { return my_mgr; }
    
    ///
    /// Returns a non-const reference to the ParameterMgr to which this Function is bound.
    ParameterMgr& getMgr() { return my_mgr; }
    
  //@}
  
  private:

  /// @name Disallowed
  //@{
  
    Function& operator= (const Function & other);
    
  //@}

  /// @name Data members
  //@{
  
    ParameterMgr&      my_mgr;
    SequenceType       my_seq;
    EvaluatorFunctor   my_eval;
 
  //@}
};

#define MF_PARAM(param_id) \
boost::lambda::bind<double>(&SAGE::MAXFUN::ParameterMgr::getEst, boost::lambda::_1, param_id)

#define MF_CONSTANT(x) \
boost::lambda::constant(x)      


} // End namespace MAXFUN
} // End namespace SAGE

#endif
