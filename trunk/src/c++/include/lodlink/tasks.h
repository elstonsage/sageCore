#ifndef LODLINK_TASKS_H
#define LODLINK_TASKS_H
//============================================================================
// File:      tasks.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   10/18/2 - created.                                   djb
//                                                                          
// Notes:     base classes for program tasks.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include <ostream>
#include <string>
#include <vector>
#include <utility>
#include <functional>
#include <algorithm>
#include "boost/smart_ptr.hpp"
#include "numerics/log_double.h"
#include "numerics/cephes.h"
#include "fped/fped.h"
#include "maxfun/maxfun.h"
#include "maxfunapi/maxfunapi.h"
#include "lodlink/instructions.h"
#include "lodlink/mle_sub_model.h"
#include "lodlink/output.h"
#include "lodlink/results.h"
#include "lodlink/likelihood.h"
#include "lodlink/max_opt.h"

using std::ostream;
using std::string;
using std::vector;
using std::pair;

namespace SAGE
{

namespace LODLINK
{

typedef boost::shared_ptr<task_result>  result_ptr;
typedef boost::shared_ptr<non_ss_alt_result>  non_ss_alt_result_ptr;
typedef boost::shared_ptr<ss_alt_result>  ss_alt_result_ptr;

bool  non_ss_has_marker(const non_ss_alt_result_ptr& ptr, const string& marker);

template<typename Arg1>
struct has_marker : public binary_function<Arg1, string, bool>
{
  bool  operator()(Arg1 ptr, string marker) const
  {
    return  ptr->marker == marker;
  }
};

template<class T, class R> void  do_task_calculations(T& task);

//----------------------------------------------------------------------------
//  Class:    task
//                                                                          
//  Purpose:  base class for a family of classes representing actions
//            specified by the user in a LODLINK analysis block.
//                                                                          
//----------------------------------------------------------------------------
//
class task
{
  public:
    typedef SAGE::FPED::FilteredMultipedigree::pedigree_const_iterator     pedigree_const_iterator;
    typedef SAGE::FPED::FilteredMultipedigree::subpedigree_type            subpedigree;
    typedef SAGE::FPED::FilteredMultipedigree::subpedigree_const_iterator  subpedigree_const_iterator;
    typedef SAGE::FPED::FilteredMultipedigree::member_const_iterator       member_const_iterator;
  
    task(cerrorstream& errors, const FPED::FilteredMultipedigree& mped, const instructions& instr);
    virtual ~task();
  
    virtual void  announce_start() const = 0;
    virtual void  calculate() = 0;
    void  write(ostream& summary, ostream& detail) const;

  protected:
    MAXFUN::Results  maximize(mle_sub_model& mle, size_t trait_index, size_t marker_index) const;    
    MAXFUN::Results  maximize(const group& g, mle_sub_model& mle, size_t trait_index, size_t marker_index);

    virtual void  write_summary(ostream& out) const = 0;
    virtual void  write_detail(ostream& out) const = 0;
    
    bool  likelihood_finite(mle_sub_model& mle, size_t trait, size_t marker, const string& test);
    bool  likelihood_finite(const string& group_name, const group& g, mle_sub_model& mle, size_t trait, 
                                                      size_t marker, const string& test);
    string  marker_name(size_t marker_index);
    
  protected:
    cerrorstream&  my_errors;
    bool  completed;        
    const FPED::FilteredMultipedigree&  my_mped;
    const instructions&  my_instructions;
    vector<result_ptr>  my_results;
};

enum sf_type { sf_LINKAGE, sf_HOMOGENEITY };    // Smith/Faraway type

//----------------------------------------------------------------------------
//  Class:    non_ss_smiths_faraways_test
//                                                                          
//  Purpose:  represents Smith's test for linkage homogeneity and Faraway's 
//            test for linkage under heterogeneity w a single recombination
//            fraction.
//            
//                                                                          
//----------------------------------------------------------------------------
//
class non_ss_smiths_faraways_test : public task
{
  public:
    non_ss_smiths_faraways_test(cerrorstream& errors, const FPED::FilteredMultipedigree& mped, 
                                const instructions& instr, sf_type t);
                                
    void  announce_start() const;                                
    void  calculate();
    static void  clear();

  private:
    void  calculate_alt(size_t trait_index, size_t marker_index, non_ss_smiths_faraways_result& result);
    void  calculate_null(size_t trait_index, size_t marker_index, non_ss_smiths_faraways_result& result);
    void  calculate_posteriors(size_t trait_index, size_t marker_index,
                                  non_ss_smiths_faraways_result& result);

    void  write_summary(ostream& out) const;
    void  write_detail(ostream& out) const;
    void  write_vc_matrix(ostream& out) const;
                                      
    sf_type  my_type;
    
    static bool  alt_calculated;             // altenative hypothesis calculations for
                                             // Smith's model performed for this analysis?
    static vector<non_ss_alt_result_ptr>  alt_results;
};

//----------------------------------------------------------------------------
//  Class:    ss_smiths_faraways_test
//                                                                          
//  Purpose:  represents Smith's test for linkage homogeneity and Faraway's 
//            test for linkage under heterogeneity w two recombination
//            fractions.
//            
//                                                                          
//----------------------------------------------------------------------------
//
class ss_smiths_faraways_test : public task
{
  public:
    ss_smiths_faraways_test(cerrorstream& errors, const FPED::FilteredMultipedigree& mped, 
                            const instructions& instr, sf_type t);

    void  announce_start() const;                                
    void  calculate();
    
    static void  clear();

  private:
    void  calculate_alt(size_t trait_index, size_t marker_index, ss_smiths_faraways_result& result);
    void  calculate_null(size_t trait_index, size_t marker_index, ss_smiths_faraways_result& result);
    void  calculate_posteriors(size_t trait_index, size_t marker_index,
                                      ss_smiths_faraways_result& result);
                                      
    void  write_summary(ostream& out) const;
    void  write_detail(ostream& out) const;
    void  write_vc_matrix(ostream& out) const;
                                      
    sf_type  my_type;
    
    static bool  alt_calculated;              // altenative hypothesis calculations for
                                              // Smith's model performed for this analysis?    
    static vector<ss_alt_result_ptr>  alt_results;
};

#include "lodlink/tasks.ipp"
}
}

#endif

