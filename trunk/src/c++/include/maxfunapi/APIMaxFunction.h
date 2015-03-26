#ifndef APIMAXFUNCTION_H
#define APIMAXFUNCTION_H

#include <time.h>
#include "maxfun/maxfun.h"
#include "maxfunapi/Datatypes.h"
#include "maxfunapi/DebugCfg.h"
#include "maxfunapi/ParameterMgr.h"

namespace SAGE   {
namespace MAXFUN {

//======================================================================
//
//  class APIMaxFunction
//
//======================================================================
class APIMaxFunction : public SAGE::MaxFunction
{
  public:
    //==================================================================
    // Constructor:
    //==================================================================

    APIMaxFunction(SAGE::MaxFunction &, ParameterMgr &, const DebugCfg &);
    
    APIMaxFunction(const APIMaxFunction &);

    //==================================================================
    // Public utility functions:
    //==================================================================

    virtual double evaluate      (vector<double> & params);
    virtual int    update_bounds (vector<double> & params);

    //==================================================================
    // Public accessors:
    //==================================================================

    ParameterMgr & getParameterMgr() const;

    ///
    /// Tells this object to which OUTPUT::Table it should direct runtime output.
    void setTable(OUTPUT::Table *);

    
  private:
    APIMaxFunction & operator= (const APIMaxFunction &);

    //==================================================================
    // Data members:
    //==================================================================

    SAGE::MaxFunction&  my_max_function;
    ParameterMgr&       my_parameter_mgr;
    const DebugCfg&     my_debug_cfg;
    OUTPUT::Table*      my_table;
    double lastvalue; // due to JA
    void output_iteration_end_results();
};

}} // End namespace

#endif
