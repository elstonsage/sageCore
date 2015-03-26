#include <iostream>
#include <iomanip>
#include <cassert>
#include <fstream>
#include <iomanip>
#include "mped/mp.h"
#include "rped/rpfile.h"
#include "fped/fped.h"
#include "error/errorstream.h"

using namespace std;
using namespace SAGE;

template <class MTYPE>
void print_member_row(ostream& o, const MTYPE& m, bool newline = true)
{
  o << setw(5) << right << m.name() << "   ";
  
  if(m.parent1())
  {
    o << setw(5) << right << m.parent1()->name() << "   "
      << setw(5) << right << m.parent2()->name() << "   ";
  }
  else
  {
    o << " --Founder--    ";
  }
  
  if(m.pedigree()->info().trait_missing(m.index(), 0))
    o << " - ";
  else if(m.pedigree()->info().trait(m.index(), 0) == 0)
    o << " N ";
  else
    o << " Y ";
  
  if(newline)
    o << endl;
}

void print_member_row(ostream& o, const FPED::Member& m)
{
  print_member_row(o, m, false);

  o << "    ";
  
  // Get the source individual
  const MPED::member_base* src = m.info().get_source_member();
  
  if(src)
  {
    o << setw(5) << right << src->name()
      << endl;
  }
  else
  {
    o << "---" << endl;
  }
}

void print_header_rows       (ostream& o, const RPED::RefPedigree& p)
{
  o << "  MEM    Par1    Par2   Inf" << endl;
  o << "=====   =====   =====   ===" << endl;
}

void print_header_rows       (ostream& o, const FPED::Pedigree& p)
{
  o << "  MEM    Par1    Par2   Inf   Source" << endl;
  o << "=====   =====   =====   ===   ======" << endl;
}

template <class Iterator>
void print_member_rows(ostream& o, Iterator b, Iterator e)
{
  for( ; b != e; ++b)
  {
    print_member_row(o, *b);
  }
}

template <class FTYPE>
void print_family_table(ostream& o, const FTYPE& f)
{
  o << "Family of Pedigree: " << f.pedigree()->name() << endl << endl;
  
  print_header_rows(o, *f.pedigree());
  
  print_member_rows(o, f.parent_begin(),    f.parent_end());
  print_member_rows(o, f.offspring_begin(), f.offspring_end());
  
  o << endl;
}

template <class STYPE>
void print_subpedigree_table(ostream& o, const STYPE& s)
{
  o << "Subpedigree of Pedigree: " << s.name() << endl << endl;
  
  print_header_rows(o, s.pedigree());
  
  print_member_rows(o, s.member_begin(), s.member_end());
  
  o << endl;
}

template <class PTYPE>
void print_pedigree_table(ostream& o, const PTYPE& p)
{
  o << "Pedigree: " << p.name() << endl << endl;
  
  print_header_rows(o, p);
  
  print_member_rows(o, p.member_begin(), p.member_end());
  
  o << endl;
}

template <class MPTYPE>
void print_multipedigree(ostream& o, const MPTYPE& m)
{
  
  o << "Multipedigree: " << endl << endl;

  for(typename MPTYPE::pedigree_const_iterator p = m.pedigree_begin(); p != m.pedigree_end(); ++p)
    print_pedigree_table(o, *p);
}

template <class MTYPE>
class member_informative : public unary_function<MTYPE,bool>
{
  public:
    bool operator() (const MTYPE& m) const
    {
      bool inf = !m.pedigree()->info().trait_missing(m.index(), 0) &&
                  m.pedigree()->info().trait(m.index(), 0);
      
      return inf;
    }
};

template<class MPTYPE>
void test_filter_all          (const MPTYPE& P)
{
    // Run Tests using RefMultiPedigree
    
    test_filter_multipedigree(P);
    
    // Run tests using the RefPedigree
    
    for(typename MPTYPE::pedigree_const_iterator i = P.pedigree_begin();
        i != P.pedigree_end(); ++i)
    {
      test_filter_pedigree(*i);
    }

    // Run family tests using the RefFamily
    
    for(typename MPTYPE::pedigree_const_iterator i = P.pedigree_begin();
        i != P.pedigree_end(); ++i)
    {
      for(typename MPTYPE::family_const_iterator j = i->family_begin();
          j != i->family_end(); ++j)
      {
        test_filter_family(*j);
      }
    }
    
    // Create a nuclear family based upon the first and last members of the 
    // filrst pedigree in the multipedigree
    
    const typename MPTYPE::pedigree_type& ped = P.pedigree_index(0);
    
    const typename MPTYPE::member_type& ind1 = ped.member_index(0);
    const typename MPTYPE::member_type& ind2 = ped.member_index(ped.member_count()-1);
    
    FPED::FilteredMultipedigree f(P);
    
    FPED::filter_to_sib_pair(f, ind1, ind2);

    f.construct();
    
    print_multipedigree(cout, f);  
}

template<class MPTYPE>
void test_filter_multipedigree(const MPTYPE& m)
{
  cout << "Multipedigree Filtration Test" << endl
       << "=======================================" << endl << endl;
  
  cout << "Original ";
     
  print_multipedigree(cout, m);

  // Test1: Create a filtered_multipedigree and insert all members from source object

  FPED::FilteredMultipedigree f1(m);

  cout << "=======================================" << endl << endl;

  cout << "Filtering all members...";

  FPED::MPFilterer::add_multipedigree(f1, m);

  f1.construct();

  cout << "Done." << endl << endl;

  // Do some basic checks -- We have included everyone, so everything should match
  // perfectly.
  
  assert(f1.pedigree_count() == m.pedigree_count());
  assert(f1.member_count()   == m.member_count());
  
  // Dump it to screen
  
  print_multipedigree(cout, f1);  

  cout << "=======================================" << endl << endl;

  // Test2: Create a filtered_multipedigree with only those members who are informative
  //        included

  FPED::FilteredMultipedigree f2(m);

  cout << "Filtering informative members...";

  FPED::FilterResults results = FPED::MPFilterer::
      add_multipedigree_filtered_by_members(f2, m, member_informative<typename MPTYPE::member_type>());

  f2.construct();

  cout << "Done." << endl << endl;

  // Dump it to screen
  
  print_multipedigree(cout, f2);  

  assert(results.get_included_member_count() + results.get_excluded_member_count() == m.member_count());

  cout << "=======================================" << endl << endl;
  
  // Test3: Create a filtered_multipedigree with members who are informative
  //        and those who are structurally informative included

  FPED::FilteredMultipedigree f3(m);

  cout << "Filtering structurally informative members...";

  results = FPED::MPFilterer::add_multipedigree_filtered_by_members
              (f3, m, FPED::is_inf_within_sped(member_informative<typename MPTYPE::member_type>()));

  f3.construct();

  cout << "Done." << endl << endl;

  // Dump it to screen
  
  print_multipedigree(cout, f3);  

  assert(results.get_included_member_count() + results.get_excluded_member_count() == m.member_count());

  cout << "=======================================" << endl << endl;
}

// Simply return the RefMultiPedigree.  This is done to support the upcasting function
// below
const RPED::RefMultiPedigree& get_mped(const RPED::RefMultiPedigree& r)
{
  return r;
}

// Do an upcast of the multipedigree type to the filtered multipedigree to make
// our tests work. This shouldn't be needed in most code, as the
// base type is ok in most instances.
const FPED::FilteredMultipedigree& get_mped(const FPED::FilteredMultipedigree::multipedigree_type& r)
{
  return ((const FPED::FilteredMultipedigree&) r);
}

template<class PTYPE>
void test_filter_pedigree(const PTYPE& p)
{
  cout << "Pedigree Filtration Test" << endl
       << "=======================================" << endl << endl;
  
  cout << "Original ";
     
  print_pedigree_table(cout, p);

  // Test1: Create a filtered_multipedigree and insert all members from source object

  FPED::FilteredMultipedigree f1(get_mped(*p.multipedigree()));

  cout << "=======================================" << endl << endl;

  cout << "Filtering all members...";

  FPED::MPFilterer::add_pedigree(f1, p);

  f1.construct();

  cout << "Done." << endl << endl;

  // Do some basic checks -- We have included everyone, so everything should match
  // perfectly.
  
  assert(f1.pedigree_count() == 1);
  assert(f1.member_count()   == p.member_count());
  
  // Dump it to screen
  
  print_multipedigree(cout, f1);  

  cout << "=======================================" << endl << endl;

  // Test2: Create a filtered_multipedigree with only those members who are informative
  //        included

  FPED::FilteredMultipedigree f2(get_mped(*p.multipedigree()));

  cout << "Filtering informative members...";

  FPED::FilterResults results = FPED::MPFilterer::
      add_pedigree_filtered_by_members(f2, p, member_informative<typename PTYPE::member_type>());

  f2.construct();

  cout << "Done." << endl << endl;

  // Dump it to screen
  
  print_multipedigree(cout, f2);  

  assert(results.get_included_member_count() + results.get_excluded_member_count() == p.member_count());

  cout << "=======================================" << endl << endl;
  
  // Test3: Create a filtered_multipedigree with members who are informative
  //        and those who are structurally informative included

  FPED::FilteredMultipedigree f3(get_mped(*p.multipedigree()));

  cout << "Filtering structurally informative members...";

  results = FPED::MPFilterer::add_pedigree_filtered_by_members
              (f3, p,FPED::is_inf_within_sped(member_informative<typename PTYPE::member_type>()));

  f3.construct();

  cout << "Done." << endl << endl;

  // Dump it to screen
  
  print_multipedigree(cout, f3);  

  assert(results.get_included_member_count() + results.get_excluded_member_count() == p.member_count());

  cout << "=======================================" << endl << endl;
}

template<class FTYPE>
void test_filter_family(const FTYPE& f)
{
  cout << "Family Filtration Test" << endl
       << "=======================================" << endl << endl;
  
  cout << "Original ";
     
  print_family_table(cout, f);

  // Test1: Create a filtered_multipedigree and insert all members from source object

  FPED::FilteredMultipedigree f1(get_mped(*f.multipedigree()));

  cout << "=======================================" << endl << endl;

  cout << "Filtering all members...";

  FPED::MPFilterer::add_family(f1, f);

  f1.construct();

  cout << "Done." << endl << endl;

  // Do some basic checks -- We have included everyone, so everything should match
  // perfectly.
  
  assert(f1.pedigree_count() == 1);
  assert(f1.member_count()   == f.offspring_count() + 2);
  
  // Dump it to screen
  
  print_multipedigree(cout, f1);  

  cout << "=======================================" << endl << endl;

  // Test2: Create a filtered_multipedigree with only those members who are informative
  //        included

  FPED::FilteredMultipedigree f2(get_mped(*f.multipedigree()));

  cout << "Filtering informative members...";

  FPED::FilterResults results = FPED::MPFilterer::
      add_family_filtered_by_members(f2, f, member_informative<typename FTYPE::member_type>());

  f2.construct();

  cout << "Done." << endl << endl;

  // Dump it to screen
  
  print_multipedigree(cout, f2);  

  assert(results.get_included_member_count() + results.get_excluded_member_count() == f.offspring_count() + 2);

  cout << "=======================================" << endl << endl;
  
  // Test3: Create a filtered_multipedigree with members who are informative
  //        and those who are structurally informative included

  FPED::FilteredMultipedigree f3(get_mped(*f.multipedigree()));

  cout << "Filtering based on family having valid sibling pairs...";

  results = FPED::MPFilterer::add_family_filtered_by_members
              (f3, f, FPED::is_family_sib_pair_inf(f,member_informative<typename FTYPE::member_type>()));

  f3.construct();

  cout << "Done." << endl << endl;

  // Dump it to screen
  
  print_multipedigree(cout, f3);  

  assert(results.get_included_member_count() + results.get_excluded_member_count() == f.offspring_count() + 2);

  cout << "=======================================" << endl << endl;
}


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
    ped_reader.add_trait_field("INFORMATIVE");

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

    // Create our RefMultiPedigree
    
    SAGE::RPED::RefMultiPedigree P;
    
    // Add a single trait "INFORMATIVE"
    
    P.info().add_binary_trait("INFORMATIVE");
    
    SAGE::RPED::RefTraitInfo& inf_trait = P.info().trait_info(0);
    
    inf_trait.set_string_affected_code("a");
    inf_trait.set_string_unaffected_code("b");
    
    // Load the data

    load_data(argv[1], argv[2], P);

    // Run Tests using RefMultiPedigree
    
    test_filter_all(P);
    
    // Create a FilteredMultipedigree and run through the tests again, using 
    // it as the source

    FPED::FilteredMultipedigree FP(P);
    
    FPED::MPFilterer::add_multipedigree(FP, P);
  
    FP.construct();
    
    test_filter_all(FP);
    
    return 0;
}

