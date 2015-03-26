#ifndef MP_UTILITIES_H
#define MP_UTILITIES_H

#include "mped/sp.h"
#include "mped/mp.h"

namespace SAGE {
namespace MPED {

/** \brief Provides a number of static utility functions for the multipedigree library
  *
  *              The mp_utilities file defines subroutine's functions which
  *		will be used in SAGE programs.
  *		
  * 		nuclear_family_count: the number of nuclear families to which
  * 		a member belongs. For unconnected member, this number is zero.
  *		For connected member, it includes both families for which the
  *		member is a child and those for which it is a parent.
  *
  *		connector: a member who belongs to (at least 2) multiple
  *		families, at the most time, it is a child in a family and
  *		is a parent in other family or families.
  *
  *		terminal family: a nuclear family with 1 or fewer connectors.
  *		
  *
  */
class mp_utilities
{
  public:
  typedef pedigree_base::subpedigree_const_iterator subpedigree_const_iterator;
  typedef pedigree_base::subpedigree_const_pointer  subpedigree_const_pointer;
  typedef pedigree_base::pedigree_const_iterator    pedigree_const_iterator;
  typedef pedigree_base::pedigree_const_pointer     pedigree_const_pointer;
  typedef pedigree_base::family_const_iterator      family_const_iterator;
  typedef pedigree_base::family_const_pointer       family_const_pointer;
  typedef pedigree_base::member_const_iterator      member_const_iterator;
  typedef pedigree_base::member_const_pointer       member_const_pointer;
  typedef pedigree_base::mate_const_iterator        mate_const_iterator;
  typedef pedigree_base::subpedigree_iterator 	    subpedigree_iterator;
  typedef pedigree_base::subpedigree_pointer        subpedigree_pointer;
  typedef pedigree_base::family_iterator            family_iterator;
  typedef pedigree_base::family_pointer       	    family_pointer;
  typedef pedigree_base::member_iterator            member_iterator;
  typedef pedigree_base::member_pointer             member_pointer;
  typedef pedigree_base::mate_iterator              mate_iterator;
  typedef pedigree_base::offspring_const_iterator   offspring_const_iterator;
  typedef pedigree_base::offspring_iterator         offspring_iterator;

  /// @name Graph-related structure
  //@{

    /** @ingroup UtilityFunctions
      *
      * Indicates whether or not the subpedigree has loops.
      * \param subpedigree The subpedigree in question
      * \retval true The subpedigree has loops
      * \retval false The subpedigree does not have loops
      */
    static bool has_loops( const subpedigree_base& subpedigree);

    /** @ingroup UtilityFunctions
      *
      * Indicates whether or not the subpedigree has loops.
      * \param subpedigree The subpedigree in question
      * \retval true The subpedigree does not have loops
      * \retval false The subpedigree has loops
      */
    static bool  no_loops( const subpedigree_base& subpedigree);

    /** @ingroup UtilityFunctions
      *
      * Indicates whether or not the subpedigree has mating chains.
      * \param subpedigree The subpedigree in question
      * \retval true The subpedigree has mating chains
      * \retval false The subpedigree does not have mating chains
      */
    static bool has_chains(const subpedigree_base& subpedigree);

    /** @ingroup UtilityFunctions
      *
      * Indicates whether or not the subpedigree has mating chains.
      * \param subpedigree The subpedigree in question
      * \retval true The subpedigree does not have mating chains
      * \retval false The subpedigree has mating chains
      */
    static bool  no_chains(const subpedigree_base& subpedigree);

  /// @name Statistics

    /** @ingroup UtilityFunctions
      *
      * Returns the size of the largest mating cluster in the subpedigree.
      * \param subpedigree The subpedigree in question
      */
    static int max_cluster_size(const subpedigree_base& subpedigree);

    /// \ingroup UtilityFunctions
    ///
    /// Returns the size of the largest mating cluster in the subpedigree to which
    /// the member belongs.  Returns 0 if the member is unconnected.
    ///
    /// \param mem A member of the subpedigree in question
//    static int max_cluster_size(const member_base& mem);

    /** @ingroup UtilityFunctions
      *
      * Hoopla!
      */
    static int getMaxSibshipSize(const multipedigree_base &);

  //@}

  /// @name Nuclear family counts
  //@{

    /** @ingroup UtilityFunctions
      *
      * Returns the number of nuclear families of which the indicated individual is a member.
      */
    static int nuclear_family_count(member_const_pointer  );

    /** @ingroup UtilityFunctions
      *
      * Hoopla!
      */
    static int nuclear_family_count(member_const_iterator );

    /** @ingroup UtilityFunctions
      *
      * Hoopla!
      */
    static int nuclear_family_count(member_pointer  );

    /** @ingroup UtilityFunctions
      *
      * Hoopla!
      */
    static int nuclear_family_count(member_iterator );
 
  //@}

  /// @name Connector counts
  //@{

    /** @ingroup UtilityFunctions
      *
      * Returns the number of connecting individuals in the indicated family.
      */
    static int connector_count(family_const_pointer  );

    /** @ingroup UtilityFunctions
      *
      * Hoopla!
      */
    static int connector_count(family_const_iterator );

    /** @ingroup UtilityFunctions
      *
      * Hoopla!
      */
    static int connector_count(family_pointer  );

    /** @ingroup UtilityFunctions
      *
      * Hoopla!
      */
    static int connector_count(family_iterator );

  //@}

  /// @name Individual status
  //@{

    /** @ingroup UtilityFunctions
      *
      * Returns whether or not the indicated member belongs to at least two nuclear families.
      */
    static bool is_connector(member_const_iterator );

    /** @ingroup UtilityFunctions
      *
      * Hoopla!
      */
    static bool is_connector(member_const_pointer  );

    /** @ingroup UtilityFunctions
      *
      * Hoopla!
      */
    static bool is_connector(member_iterator );

    /** @ingroup UtilityFunctions
      *
      * Hoopla!
      */
    static bool is_connector(member_pointer  );

    /** @ingroup UtilityFunctions
      *
      * Hoopla!
      */
    static bool is_connector(const member_base&  );

    /** @ingroup UtilityFunctions
      *
      * Returns whether or not an individual is a founder
      */
    static bool is_founder(member_const_iterator );

    /** @ingroup UtilityFunctions
      *
      * Hoopla!
      */
    static bool is_founder(member_iterator );

    /** @ingroup UtilityFunctions
      *
      * Hoopla!
      */
    static bool is_founder(member_const_pointer );

    /** @ingroup UtilityFunctions
      *
      * Hoopla!
      */
    static bool is_founder(member_pointer );

    /** @ingroup UtilityFunctions
      *
      * Hoopla!
      */
    static bool is_founder(const member_base&  );

    /** @ingroup UtilityFunctions
      *
      * Returns whether or not the indicated individual is in the indicated family.
      */
    static bool in_family(member_const_pointer,  family_const_pointer  ); 

    /** @ingroup UtilityFunctions
      *
      * Hoopla!
      */
    static bool in_family(member_const_iterator, family_const_iterator );

    /** @ingroup UtilityFunctions
      *
      * Hoopla!
      */
    static bool in_family(member_pointer,  family_pointer  );

    /** @ingroup UtilityFunctions
      *
      * Hoopla!
      */
    static bool in_family(member_iterator, family_iterator );

  /// @name Family status
  //@{

    /** @ingroup UtilityFunctions
      *
      * Returns whether or not family has one and only one connector.
      */
    static bool is_terminal_family(family_const_iterator );

    /** @ingroup UtilityFunctions
      *
      * Hoopla!
      */
    static bool is_terminal_family(family_const_pointer  );

    /** @ingroup UtilityFunctions
      *
      * Hoopla!
      */
    static bool is_terminal_family(family_iterator );

    /** @ingroup UtilityFunctions
      *
      * Hoopla!
      */
    static bool is_terminal_family(family_pointer  );

  //@}
};

} // End namespace MPED
} // End namespace SAGE

#include "mped/mp_utilities.ipp"

#endif
