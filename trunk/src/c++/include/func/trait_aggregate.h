#ifndef TRAIT_AGGREGATE_H
#define TRAIT_AGGREGATE_H
//============================================================================
// File:      trait_aggregate.h
//                                                                          
// Author:    Kai He
//                                                                          
// History:   Created 12/11/02
//                                                                          
// Notes:     
//                                                                          
// Copyright (c) 2000 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include <string>
#include <map>
#include <vector>
#include <list>
#include <cctype>
#include <limits>
#include <cassert>
#include <setjmp.h>
#include <signal.h>
#include <unistd.h>
#include <iomanip>

#include "rped/rped.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "numerics/sinfo.h"
#include "numerics/histogram.h"
#include "numerics/mt.h"
#include "numerics/functions.h"
#include "Function.h"
#include "FunctionParser.h"

namespace SAGE {
namespace FUNC {

class trait_aggregate
{
  public:
    enum function_string { isMissing, 
    		           aggregate_mean, 
    		           aggregate_var, 
    		           aggregate_mean_var  };

    typedef std::map<string,double>           map;
    typedef std::map<string,double>::iterator miter;

    struct record
    {
      string trait_name;
      size_t trait_no;
      map    mem_value;
    };

    typedef std::string                   string;
    typedef std::list<char>               list;
    typedef std::list<char>::iterator	  list_iterator;
    typedef std::list<string>             list_string;
    typedef std::list<string>::iterator   list_string_iterator;  
    typedef std::vector<double>           trait_value_vector;
    typedef std::vector<record>           v_record;
    typedef std::vector<record>::iterator v_record_iterator;
    
    trait_aggregate           ()                    { }    
    trait_aggregate           (cerrorstream errors):my_errors(errors) { } 
    trait_aggregate           (const trait_aggregate& );
    trait_aggregate& operator=(const trait_aggregate& );
    
    void   dump             ();
    void   dump             (SampleInfo,string);
    void   print_out        ();
    size_t get_trait_no	    () const { return my_trait_no;   	}
    string get_trait_name   () const { return my_trait_name; 	}	    
 
    size_t get_trait_no     () { return my_trait_no;      }
    string get_trait_name   () { return my_trait_name;    }
                                
    trait_value_vector& get_traits() { return my_traits;	}
                                                            
    string parse_expression   (RPED::RefMultiPedigree&, string);
    string parse_aggregation  (RPED::RefMultiPedigree&, string);
    string parse_merge        (string);
    string parse_expression2  (RPED::RefMultiPedigree&, string);
    bool   is_in_class        (RPED::RefMultiPedigree&, RPED::RefPedInfo&, string, size_t);   
    void   process_aggregation(RPED::RefMultiPedigree&, string, string);    			     
    void   create_adjusted_trait(RPED::RefMultiPedigree&, string);
    void   create_adjusted_trait(RPED::RefMultiPedigree&, const string&, string);
    void   create_adjusted_trait(RPED::RefMultiPedigree&, const FunctionParser::TraitData&);
    void   merge_class_trait    (RPED::RefMultiPedigree&, const FunctionParser::TraitData&);
    void   create_class_trait   (RPED::RefMultiPedigree&, const FunctionParser::TraitData&);
        
    void   write_error_msg(const string&);
    void   write_test_msg (const string&); //test parsing aggregation functions
                
  private:
  
    cerrorstream        my_errors;
    size_t              my_trait_no;
    string              my_trait_name;
    trait_value_vector  my_traits;
    SampleInfo          my_sinfo;
    record              my_record;
    v_record            my_v_record;
            
};

} // End namespace FUNC
} // End namespace SAGE

#endif
