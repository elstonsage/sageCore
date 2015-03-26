#ifndef TAI_H
#define TAI_H
//============================================================================
// File:      tai.h (Transmitted Allele Status)
//
// Author:    Kai He
//
//
// History:   10/2003 Initial version
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

/** \class transmitted_allele_indicator
  * \brief Class that calculates the transmitted allele indicator variable.
  *
  */
class transmitted_allele_indicator
{
public:

    typedef RPED::RefMultiPedigree::pedigree_iterator    pedigree_iterator;
    typedef RPED::RefMultiPedigree::subpedigree_iterator subpedigree_iterator;   
    typedef RPED::RefPedigree::family_iterator           family_iterator;
    typedef RPED::RefMultiPedigree::member_const_pointer member_const_pointer;
    typedef RPED::RefPedigree::member_iterator           member_iterator;
    typedef RPED::RefPedigree::offspring_iterator        offspring_iterator;
    typedef list<string>                                 parameter_list;
                    

  /// @name Constructors & operators
  //@{

    ///
    /// Constructor
    /// \param err The errorstream to which error messages will be sent
    explicit transmitted_allele_indicator (cerrorstream& err = sage_cerr);

    /// 
    /// Copy constructor
    /// \param other Object that will be copied
    transmitted_allele_indicator (const transmitted_allele_indicator & other);

    ///
    /// Assignment operator
    /// \param other Object that will be copied
    transmitted_allele_indicator& operator= (const transmitted_allele_indicator & other);

    ///
    /// Destructor
    ~transmitted_allele_indicator();

  //@}

  /// @name Creating the TAI status variable
  //@{

    ///
    /// Creates a TAI trait within the indicated RefMultiPedigree.
    /// \param mp The RefMultiPedigree instance in which to create the TAI status variable
    /// \param data Parser data indicating on which marker/allele to base the TAI status calculation
    void createTaiStatusTrait(RPED::RefMultiPedigree & mp, const FunctionParser::TraitData & data);
    
  //@}

  private:

  /// @name Parsing functions
  //@{

    void parseTaiExpression (string expression);

    bool verifyParameters(const RPED::RefMultiPedigree & mp);

  //@}
  
  /// @name RefMultiPedigree value-populating functions
  //@{
  
    void populateRefMultiPedigree(RPED::RefMultiPedigree& mp, const FunctionParser::TraitData & pd);

    double getTaiStatus(const RPED::RefMember & mem);

  //@}

    bool is_homozygous(string allele1_name, string allele2_name, string comparison_allele_name);
    bool is_heterozygous(string, string, string);
    bool missing_allele(string);
    bool allele_in_genotype(string, string, string);

    size_t check_parenthese(string);
    string deparenthese    (string);
    size_t add_trait_number(RPED::RefMPedInfo&, const FunctionParser::TraitData &);

    string get_allele_1(string);
    string get_allele_2(string);
    string get_marker();
    string get_allele();
    
    cerrorstream	my_errors;
    string		my_marker_name;
    string		my_allele_name;

    enum TaiType { TAI, UTAI }; // Whether this is to be a TAI calculation or a UTAI calculation
    
    TaiType my_type;
};    
    
//========================
//  INLINE FUNCTIONS
//========================

inline string transmitted_allele_indicator::get_marker() { return my_marker_name; }
inline string transmitted_allele_indicator::get_allele() { return my_allele_name; }
        
} // End namespace FUNC
} // End namespace SAGE

#endif
