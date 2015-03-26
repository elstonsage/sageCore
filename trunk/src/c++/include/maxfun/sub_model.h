#ifndef SUB_MODEL_H
#define SUB_MODEL_H
//============================================================================
// File:      sub_model.h
//                                                                          
// Author:    Geoff Wedig (wedig@darwin.cwru.edu)
//                                                                          
// History:   0.1 gcw Initial Implementation                    Apr 2001   
//                djb reformatted, added  misc. global 
//                    declarations/definitions                  Jun 2001
//            1.0 gcw Moved out of SEGREG                       Aug 2002
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

/** @file
 *            The purpose of this file is to create an easy method of
 *            inclusion of different sub-models into maxfun.
 *
 *
 */



#include <math.h>
#include <assert.h>
#include <cmath>
#include <list>
#include <vector>
#include <functional>
#include <algorithm>
#include <limits>
#include <ostream>
#include <sstream>
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "maxfun/maxfun.h"
#include "LSF/parse_ops.h"
#include "globals/SAGEConstants.h"

using std::string;
using std::ostream;

namespace SAGE
{

/// Externs that are often used with maxfun parameters
//@{

extern const double  POSITIVE_INF;
extern const double  NEGATIVE_INF;
extern const double  DEFAULT_LOWER_BOUND;
extern const double  DEFAULT_UPPER_BOUND;

//@}

// Declared so the sub_model_sequencer knows it's coming.
class  sub_model;

///       pstatus details the various options for a parameter in maxfun
/**
 *        Note that pstatus only includes the input status options.
 *        Output status is controlled elsewhere. (Q: Should be moved to maxfun.h?
 */
enum pstatus { indep_func = 1, indep_non_func = 2, dependent = 3, fixed = 4 };

/// struct for specifying a parameter to maxfun.
///
/// It includes the value, bounds, status flag and a text string
/// identifier of the parameter.

struct maxfun_parameter  
{
  maxfun_parameter(string   ident    = "",
                   pstatus  status   = fixed,
                   double   val      =   0,
                   double   l_bound  = DEFAULT_LOWER_BOUND,
                   double   u_bound  = DEFAULT_UPPER_BOUND );
                   
  maxfun_parameter(const maxfun_parameter& other);
  maxfun_parameter&  operator=(const maxfun_parameter& other);
                   
  double   value;
  double   lower_bound;
  double   upper_bound;
  pstatus  status;
  string   identifier;
};

/// struct often used for input of variables into a sub_model.
///
/// The model_input struct is used to input variables into a sub_model

struct model_input
{
  model_input(double i = QNAN, bool f = false);
  model_input(const model_input& other);
  model_input&  operator=(const model_input& other);

  double  value;
  bool    fixed;
};

/**
 *            The sub_model_sequencer ontrols the sequencing and updating of
 *            sub_models that are given to maxfun.  For each sub_model
 *            given, it adds the sub_model's parameters into maxfun.
 *
 *            It creates a new function to wrapper the Maxfunction it is given. This
 *            new function handles the updating and control of the sub_models. When
 *            update_bounds() is called, each sub_model in its sub_model list is called using
 *            synchronize().  If any fail, it returns that fail state.  Otherwise, it
 *            calls update_bounds on it's internal function, which may perform
 *            additional checks.
 *
 *            The evaluate function simply passes through to the internal Maxfunction
 *            object.  All sub_models should be up to date at this point.
 */
class sub_model_sequencer
{
  protected:

    // - A reference to a sub_model with the location of the parameters.
    //
    struct sub_model_reference
    {
      sub_model_reference(sub_model*, size_t);

      sub_model* mod;
      size_t pbegin;
    };

    typedef std::list<sub_model_reference> sub_model_list;

  public:

    class sub_model_sequencer_function : public MaxFunction
    {
      public:
        inline sub_model_sequencer_function(MaxFunction& mf, sub_model_list&);
      
        virtual  inline ~sub_model_sequencer_function();

      protected:
        virtual inline double  evaluate(parameter_vector& theta);
        virtual inline int     update_bounds(parameter_vector& theta);

        // Data members.
        MaxFunction&     my_maxfunction;
        sub_model_list&  my_sub_models;
    };

  public:
    sub_model_sequencer(MaxFunction& mf);

    ~sub_model_sequencer();

          Maxfun&  maxfun();
    const Maxfun&  maxfun() const;

    void  add_sub_model(sub_model* m);

  private:
  
    // Data members.
    sub_model_list                my_sub_models;
    sub_model_sequencer_function  my_function;
    Maxfun                        my_maxfun;
};


/**           base class to all sub_models.
 *            It is used to assign, sequence, and access the maxfun
 *            estimated model parameters without the need to access them
 *            directly.
 *                                                                          
 *            In essence, a 'sub-model' is a class that contains a set of
 *            parameters for maxfun to use.  There is only a minimal public
 *            interface to the sub-model; most of the interface is supplied
 *            in the derived class.  The sub-model provides the initial
 *            estimate, bounds, and status type for each parameter to a
 *            sub-model sequencer.  The sub-model sequencer then inserts
 *            this into the list of parameters on its maxfun object.
 *
 *            In addition to the construction phase, the sub-model sequencer
 *            wrappers the function calls maxfun makes to MaxFunction.  This
 *            wrapper copies out the data from the maxfun parameter list
 *            back into the sub-models, which should do any dependency
 *            checking.  The evaluate() function can then use the sub-model
 *            instead of the raw vector supplied to it from maxfun.  This
 *            allows a nearly transparent use of the sub-models, without the
 *            need for manual sequencing, which is prone to errors.
 *
 *            In addition, the modular framework of the sub-model makes it
 *            quite easy to maintain. With proper use (each sub-model
 *            dealing with only a few parameters), it becomes nearly trivial
 *            to take out one sub-model and replace it with another, without
 *            disturbing code that does not rely on that sub-model. With the
 *            raw vector, this is not the case.  Knowing which parameter is
 *            which index is of prime importance.  The sub-model interface
 *            takes care of that detail so that the programmer doesn't have
 *            to.
 */
class sub_model
{
  public:
    friend ostream&  operator<<(ostream& out, const sub_model& sm);
  
    inline sub_model(cerrorstream& errors = sage_cerr);
    inline sub_model(const sub_model& other);
    inline sub_model&  operator=(const sub_model& other);
  
    // - Virtual destructor is required due to virtual functions.
    //
    virtual inline ~sub_model();

    // - This public interface is a simple one and derived classes should add to
    //   it.  Primarily, this is for access to the 'underlying' sub-model.
    //
    virtual inline string  option_description() const; 
    virtual inline string  name() const;
    
    inline bool    in_use() const;
    inline size_t  parameter_count() const; 
    inline string  parameter_identifier(size_t i) const;

    virtual inline double  parameter(size_t i) const; 
    
    virtual  inline const std::vector<maxfun_parameter>&   get_parameters();

  protected:
    friend class sub_model_sequencer;
    friend class sub_model_sequencer::sub_model_sequencer_function;
    typedef std::vector<double>         parameter_vector;
    typedef parameter_vector::iterator  parameter_iterator;

    // - synchronize is called by update_bounds.  It should do
    //   any work necessary to complete the public interface of the derived
    //   class.  Should return 0 if all ok, positive number if there is
    //   a problem.
    //
    //   Note that any dependent variables, or any variables that are modified
    //   from what is given by maxfun need to be copied into both my_parameters
    //   and the appropriate maxfun variable.
    //
    //   Default version copies the parameters from maxfun into my_parameters
    //   and returns 0;
    //
    virtual inline int synchronize(parameter_iterator start);

    // Data members.
    bool  my_in_use;                          // Set by the sub_model_sequencer.

    mutable cerrorstream  my_errors; //< Stream for reporting errors

    mutable std::vector<maxfun_parameter>  my_parameters;
};

//lint -esym(534, sub_model::synchronize)

/// Debugging functions
//@{
ostream&  operator<<(ostream& out, const maxfun_parameter& max_par);

string    pstatus_2_string(pstatus status);
string    value_phrase(const maxfun_parameter& mp);
//@}
 
#include "maxfun/sub_model.ipp"

}

#endif


