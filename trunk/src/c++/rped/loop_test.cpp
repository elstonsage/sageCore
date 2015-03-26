#include <iostream>
#include <iomanip>
#include <cassert>
#include <fstream>
#include <iomanip>
#include "rped/loop.h"
#include "mped/mp.h"
#include "rped/rpfile.h"
#include "error/errorstream.h"

namespace SAGE {
namespace RPED {

typedef RefPedigree ped_t;
typedef ped_t::family_const_iterator      f_iter;
typedef ped_t::mate_const_iterator        a_iter;
typedef ped_t::member_const_iterator      m_iter;
typedef ped_t::offspring_const_iterator   o_iter;
typedef ped_t::parent_const_iterator      p_iter;
typedef ped_t::sibling_const_iterator     s_iter;
typedef ped_t::subpedigree_const_iterator u_iter;

bool    
load_data(const string &filename, string format, RefMultiPedigree& P)
{
    RefFortranPedigreeFile ped_reader(SAGE::sage_cerr);

    ped_reader.set_format(format);

    ped_reader.add_pedigree_id();
    ped_reader.add_individual_id();
    ped_reader.add_sex_field();
    ped_reader.add_parent_id();
    ped_reader.add_parent_id();

    if( !ped_reader.input(P, filename) )
    {
      cerr << "Error reading pedigree data" << endl;
      return FALSE;
    }
    return TRUE;
}

} // End namespace RPED
} // End namespace SAGE

int
main(int argc, char* argv[])
{
    if (argc != 3)
    {
        std::cerr << "usage: " << argv[0] << " filename format" << std::endl;
        return -1;
    }

    SAGE::RPED::RefMultiPedigree P;
    load_data(argv[1], argv[2], P);

    SAGE::RPED::LoopChecker  lck;

    SAGE::RPED::RefMultiPedigree::pedigree_iterator i;
    for( i = P.pedigree_begin(); i != P.pedigree_end(); ++i)
    {
      cout << endl << "Checking pedigree " << i->name() << ':';

      if(!lck.check_pedigree( &*i ))
        cout << "   Error checking pedigree:  bad pedigree." ;
      else if(lck.loops())
        cout << endl << "   Loops detected:"
             << "  " << lck.marriage_loops() << " marriage loops."
             << "  " << lck.non_marriage_loops() << " non-marriage loops.";
      else
        cout << " normal pedigree.";

      cout << endl;
    }

    return 0;
}

