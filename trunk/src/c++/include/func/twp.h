#ifndef TWP_H
#define TWP_H
//============================================================================
// File:      twp.h (trimming-winsorization process)
//
// Author:    Kai He
//
//
// History:   12/2003 Initial version
//
// Notes:
//
// Copyright (c) 2002 R.C. Elston
// All Rights Reserved
//============================================================================


#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <limits>
#include <list>
#include <map>
#include <cctype>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <utility>

#include "error/errorstream.h"
#include "error/errormanip.h"
#include "numerics/cephes.h"
#include "numerics/constants.h"
#include "numerics/sinfo.h"
#include "numerics/mt.h"
#include "numerics/functions.h"

#include "mlocus/genotype.h"
#include "mped/mp.h"
#include "mped/sp.h"
#include "mped/mp_utilities.h"
#include "rped/rped.h"
#include "func/FunctionParser.h"

namespace SAGE {
namespace FUNC {

class trim_winsor_process
{
  public:
    typedef RPED::RefMultiPedigree::pedigree_iterator	 	pedigree_iterator;
    typedef RPED::RefMultiPedigree::subpedigree_iterator 	subpedigree_iterator;   
    typedef RPED::RefPedigree::family_iterator 		family_iterator;
    typedef RPED::RefPedigree::member_iterator 		member_iterator;
    typedef list<string>				parameter_list;
    typedef const FunctionParser::TraitData&  		parser_data;
                        
    struct Person
    {
      size_t		position;
      double            trait;
      pedigree_iterator ped;
      member_iterator   mem;
    };

    typedef multimap<double,Person>	trait_map;
    typedef trait_map::iterator		trait_map_iterator;
    typedef pair<double,Person>		pair;
    
    typedef multimap<size_t,Person>     pos_map;
    typedef pos_map::iterator           pos_map_iterator;
                
    
    trim_winsor_process();
    trim_winsor_process(cerrorstream& err = sage_cerr);
    trim_winsor_process(const trim_winsor_process& );
    trim_winsor_process& operator =(const trim_winsor_process& );
    ~trim_winsor_process();

    size_t check_parenthese(string);
    string deparenthese    (string);
    size_t add_trait_number(RPED::RefMPedInfo&, parser_data);

    void parse_twp_expression (RPED::RefMultiPedigree&, string);
    void create_twp(RPED::RefMultiPedigree&, parser_data);
    void create_twp_trait(RPED::RefMultiPedigree&, parser_data);

    void build_data           (RPED::RefMultiPedigree&, parser_data);
    void sort_positioning_data();
    void trim_data            ();
    void winsor_data          ();
    
    string get_process(){ return my_process;}
    string get_trait()  { return my_trait;}
    double get_gamma()  { return my_gamma;}
    size_t get_members(){ return my_people.size();}
    size_t get_missing_members(){ return my_missing_members;}
    
    /// for error output
    
    enum error_message { invalid_trait, duplicate_trait, invalid_gamma };
    void write_error_message(error_message, const string);

    /// for testing and debugging
        
    void dump_members();
    void dump_sort_members();
    void dump_pos_members();
    
  protected:

  private:
  
    vector<Person>	my_people;
    vector<Person>      my_sort_people;
    trait_map		my_trait_map;
    pos_map		my_pos_map;
    size_t		my_missing_members;

    cerrorstream	my_errors;
    string		my_process;
    string		my_trait;
    double		my_gamma;    

};    

} // End namespace FUNC
} // End namespace SAGE
    
#endif
