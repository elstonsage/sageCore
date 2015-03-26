// Copyright(c) 1998 RC Elston
#if _MSC_VER
    #include <mpfwd.h>
    #pragma hdrstop
#endif

#include <fstream>
#include <iomanip>
#include "mped/mp.h"
#include "rped/rped.h"
#include "rped/rpfile.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "output/Output.h"

namespace SAGE {
namespace RPED {

typedef RefPedigree ped_t;

typedef ped_t::family_const_iterator f_iter;
typedef ped_t::mate_iterator        a_iter;
typedef ped_t::member_iterator      m_iter;
typedef ped_t::offspring_const_iterator   o_iter;
typedef ped_t::parent_iterator      p_iter;
typedef ped_t::sibling_iterator     s_iter;
typedef ped_t::subpedigree_iterator u_iter;

bool    
load_data(const string &filename, string format, string outfile,
          RefMultiPedigree& P)
{
    cerrorstream errors;
    errors << prefix("[%P] ") << line_width(80);

    RefFortranPedigreeFile ped_reader(errors);

    ped_reader.set_format(format);

    ped_reader.add_pedigree_id();
    ped_reader.add_individual_id();
    ped_reader.add_sex_field();
    ped_reader.add_parent_id();
    ped_reader.add_parent_id();

    P.info().set_sex_code_male("1");
    P.info().set_sex_code_female("2");
    P.info().set_sex_code_unknown(" ");
//  P.info().set_individual_missing_code("0");

    if( !ped_reader.input(P, filename) )
    {
      cerr << "Error reading pedigree data" << endl;
      return FALSE;
    }

    ped_reader.output(P, outfile);
    return TRUE;
}

void
print_summary(const ped_t& P)
{
  SAGE::OUTPUT::Table t1("Pedigree statistics");
  
  t1 << (SAGE::OUTPUT::TableRow() << "Members"          << P.member_count())
     << (SAGE::OUTPUT::TableRow() << "Nuclear families" << P.family_count());
  
  std::cout << t1;

  SAGE::OUTPUT::Table t2("Nuclear families");
  
  f_iter  ff = P.family_begin();
  f_iter  fl = P.family_end();

  for (;  ff != fl;  ++ff)
  {
    t2 << (SAGE::OUTPUT::TableRow() << "Parents" << ff->name1() << ff->name2());
      
    o_iter  of = ff->offspring_begin();
    o_iter  ol = ff->offspring_end();

    SAGE::OUTPUT::TableRow r;
      
    r << "Children";

    for (;  of != ol;  ++of)
    {
      r << of->name();
    }
        
    t2 << r;
  }
    
  std::cout << t2;
}

void check_order(RefPedigree &p)
{
  RefPedigree::member_iterator i = p.member_begin();
  for( ; i != p.member_end(); ++i)
  {
    int p1,p2,mp;
    p1=p2=-1;
    if(i->parent1()) p1 = i->parent1()->index();
    if(i->parent2()) p2 = i->parent2()->index();
    mp = std::max(p1,p2);
    if( mp > -1 && i->index() < (size_t)mp)
      cout << "Ordering violated for individual '" << i->name() << "'."
           << "(" << p1 << "," << p2 << ") -> " << i->info() << std::endl;
  }
}

} // End namespace RPED
} // End namespace SAGE

int
main(int argc, char* argv[])
{
    if (argc != 4)
    {
        std::cerr << "usage: " << argv[0] << " filename format <ped file>" << std::endl;
        return -1;
    }

    SAGE::RPED::RefMultiPedigree P;
    load_data(argv[1], argv[2], argv[3], P);

    SAGE::RPED::RefMultiPedigree::pedigree_iterator i;
    for( i = P.pedigree_begin(); i != P.pedigree_end(); ++i)
    {
      print_summary(*i);
      check_order(*i);
    }

    for( i = P.pedigree_begin(); i != P.pedigree_end(); ++i)
      cout << i->name() << " = " << i->index() << endl;

    return 0;
}

