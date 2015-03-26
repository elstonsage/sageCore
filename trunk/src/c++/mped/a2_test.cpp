#include "mped/sp.h"
#include "mped/mp.h"

using namespace SAGE::MPED;

typedef pedigree<SAGE::MPED::no_info>                       pedigree_type;
typedef multipedigree<SAGE::MPED::no_info>                  multipedigree_type;
typedef multipedigree_type::pedigree_iterator               mp_iter;
using namespace std;

#include <iomanip>

std::ostream& operator<<(std::ostream& o, SexCode c)
{
  switch(c)
  {
    case SEX_MALE    : o << "MALE";      break;
    case SEX_XMALE   : o << "XMALE";     break;
    case SEX_FEMALE  : o << "FEMALE";    break;
    case SEX_XFEMALE : o << "XFEMALE";   break;
    case SEX_MISSING : o << "MISSING";   break;
    case SEX_ARB     : o << "ARBITRARY"; break;
    default          : o << "UNKNOWN!";  break;
  }
  
  return o;
}

void    test_sex_code();

void    test0();
void    test1();
void    test2();
void    test3();
void    test4();
void    test5();
void    test6();
void    test7();
void    test8();
void    test9();

// Error tests
void    test_bad_sibship();
void    test_bad_marriage();
void    test_bad_lineage();

// Mother father test
void    test_get_moth_fath();

string  m01("m01");
string  m02("m02");
string  m03("m03");
string  m04("m04");
string  m05("m05");
string  m06("m06");
string  m07("m07");
string  m08("m08");
string  m09("m09");
string  m10("m10");
string  m11("m11");
string  m12("m12");
string  m13("m13");
string  m14("m14");
string  m15("m15");
string  m16("m16");
string  m17("m17");
string  m18("m18");
string  m19("m19");
string  m20("m20");

typedef pedigree_type::family_iterator      f_iter;
typedef pedigree_type::mate_iterator        a_iter;
typedef pedigree_type::member_iterator      m_iter;
typedef pedigree_type::offspring_iterator   o_iter;
typedef pedigree_type::parent_iterator      p_iter;
typedef pedigree_type::progeny_iterator     y_iter;
typedef pedigree_type::sibling_iterator     s_iter;
typedef pedigree_type::subpedigree_iterator u_iter;


//----------------------------------------------------------------------------
//  Function:   print()
//
//  Purpose:    Display semi-comprehensize data about pedigree.
//----------------------------------------------------------------------------
//
void    print(const string& tname, const pedigree_type& P)
{
    typedef pedigree_type::subpedigree_const_iterator   u_iter;
    typedef pedigree_type::member_const_iterator        m_iter;
    typedef pedigree_type::family_const_iterator        f_iter;
    typedef pedigree_type::offspring_const_iterator     o_iter;

    std::cout << std::endl;
    std::cout << "#######################################################";
    std::cout << std::endl << "For test: " << tname << std::endl;

    //- First print all of the members of the pedigree.
    //
    std::cout << std::endl << "  Members of: " << P.name() << std::endl;

    m_iter  mf = P.member_begin();
    m_iter  ml = P.member_end();

    for (;  mf != ml;  ++mf)
    {
        std::cout << "              " << mf->name() 
                  << " (" << mf->subpedigree()->name() << ")" << std::endl;
    }


    //- Next print a listing of each family in the pedigree.  For each
    //  family listed, also list the parents and children in that family.
    //  
    std::cout << std::endl << " Families of: " << P.name() << std::endl;

    f_iter  ff = P.family_begin();
    f_iter  fl = P.family_end();

    for (;  ff != fl;  ++ff)
    {
        std::cout << "     parents:";
        std::cout << " " << ff->name1() << " " << ff->name2() << std::endl;
        std::cout << "    children:";
    
        o_iter  of = ff->offspring_begin();
        o_iter  ol = ff->offspring_end();

        for (;  of != ol;  ++of)
        {
            std::cout << " " << of->name();
        }
        std::cout << std::endl;
    }
}


void    print2(const string& tname, const pedigree_type& P)
{
    typedef pedigree_type::subpedigree_const_iterator   u_iter;
    typedef pedigree_type::member_const_iterator        m_iter;
    typedef pedigree_type::family_const_iterator        f_iter;
    typedef pedigree_type::offspring_const_iterator     o_iter;

    std::cout << std::endl;
    std::cout << "#######################################################";
    std::cout << std::endl << "For test: " << tname << std::endl;

    u_iter  uf = P.subpedigree_begin();
    u_iter  ul = P.subpedigree_end();

    for (;  uf != ul;  ++uf)
    {
        //- First print all of the members of the pedigree.
        //
        std::cout << std::endl << "  Members of subped: " 
                  << uf->name() << std::endl;
        
        m_iter  mf = uf->member_begin();
        m_iter  ml = uf->member_end();

        for (;  mf != ml;  ++mf)
        {
            std::cout << "                     ";
            std::cout << mf->name() << std::endl;
        }

        //- Next print a listing of each family in the subpedigree.  For each
        //  family listed, also list the parents and children in that family.
        //  
        std::cout << std::endl << " Families of subped: " 
                  << uf->name() << std::endl;

        f_iter  ff = uf->family_begin();
        f_iter  fl = uf->family_end();

        for (;  ff != fl;  ++ff)
        {
            std::cout << "     parents:";
            std::cout << " " << ff->name1() << " " << ff->name2() << std::endl;
            std::cout << "    children:";
    
            o_iter  of = ff->offspring_begin();
            o_iter  ol = ff->offspring_end();

            for (;  of != ol;  ++of)
            {
                std::cout << " " << of->name();
            }
            std::cout << std::endl;
        }
    }
}

template <class Iter> void  
iprint(const char* s, Iter f, Iter l)
{
    std::cout << std::setw(16) << std::left << s << std::right << ":";

    for (;  f != l;  ++f)
    {
        std::cout << " " << f->name();
    }
    std::cout << std::endl;
}

void   test_sex_code()
{
  assert(SEX_MALE != SEX_XMALE);
  assert(SEX_MALE != SEX_FEMALE);
  assert(SEX_MALE != SEX_XFEMALE);
  assert(SEX_MALE != SEX_MISSING);
  assert(SEX_MALE != SEX_ARB);

  assert(SEX_XMALE != SEX_FEMALE);
  assert(SEX_XMALE != SEX_XFEMALE);
  assert(SEX_XMALE != SEX_MISSING);
  assert(SEX_XMALE != SEX_ARB);

  assert(SEX_FEMALE != SEX_XFEMALE);
  assert(SEX_FEMALE != SEX_MISSING);
  assert(SEX_FEMALE != SEX_ARB);

  assert(SEX_XFEMALE != SEX_MISSING);
  assert(SEX_XFEMALE != SEX_ARB);

  assert(SEX_MISSING != SEX_ARB);

  assert( is_male(SEX_MALE));
  assert( is_male(SEX_XMALE));
  assert(!is_male(SEX_FEMALE));
  assert(!is_male(SEX_XFEMALE));
  assert(!is_male(SEX_MISSING));
  assert(!is_male(SEX_ARB));

  assert(!is_female(SEX_MALE));
  assert(!is_female(SEX_XMALE));
  assert( is_female(SEX_FEMALE));
  assert( is_female(SEX_XFEMALE));
  assert(!is_female(SEX_MISSING));
  assert(!is_female(SEX_ARB));

  assert(!is_sex_unknown(SEX_MALE));
  assert(!is_sex_unknown(SEX_XMALE));
  assert(!is_sex_unknown(SEX_FEMALE));
  assert(!is_sex_unknown(SEX_XFEMALE));
  assert( is_sex_unknown(SEX_MISSING));
  assert( is_sex_unknown(SEX_ARB));

  assert(get_effective_sex(SEX_MALE)     == SEX_MALE    &&
         get_effective_sex(SEX_XMALE)    == SEX_MALE);
  assert(get_effective_sex(SEX_FEMALE)   == SEX_FEMALE  &&
         get_effective_sex(SEX_XFEMALE)  == SEX_FEMALE);
  assert(get_effective_sex(SEX_MISSING)  == SEX_MISSING &&
         get_effective_sex(SEX_ARB)      == SEX_MISSING);
}


//----------------------------------------------------------------------------
//  Function:   test0()
//
//  Scenario:   simple 2-parent lineage; iteration test
//----------------------------------------------------------------------------
//
void    test0()
{
    pedigree_type           P("P0");
    const pedigree_type&    CP = P;

    P.add_member(m01);
    P.add_member(m02);
    P.add_member(m03);
    P.add_member(m04);
    P.add_member(m05);
    P.add_member(m06);
    P.add_member(m07);
    P.add_member(m08);
    P.add_member(m09);
    P.add_member(m10);
    P.add_member(m11);
    P.add_member(m12);
    P.add_member(m13);
    P.add_member(m14);
    P.add_member(m15);
    P.add_member(m16);
    P.add_member(m17);
    P.add_member(m18);
    P.add_member(m19);
    P.add_member(m20);
    
    P.add_lineage(m11, m01, m02);
    P.add_lineage(m12, m01, m02);

    P.add_lineage(m13, m02, m03);
    P.add_lineage(m14, m02, m03);
    P.add_lineage(m15, m02, m03);

    P.add_lineage(m16, m03, m04);
    P.add_lineage(m17, m03, m04);
    P.add_lineage(m18, m03, m04);
    P.add_lineage(m19, m03, m04);
 
    P.add_lineage(m20, m03, m05);

    P.add_lineage(m03, m07, m08);
    P.add_lineage(m09, m07, m08);

    P.build();

    print("test0", P);

    SAGE::MPED::pedigree_base::subpedigree_iterator         sf = P.subpedigree_begin();
    SAGE::MPED::pedigree_base::subpedigree_const_iterator   csf = CP.subpedigree_begin();

    SAGE::MPED::family_base*        pf  = P.family_find(m03, m04);
    const SAGE::MPED::family_base*  cpf = CP.family_find(m03, m04);

    SAGE::MPED::member_base*        pm = P.member_find(m03);
    const SAGE::MPED::member_base*  cpm = CP.member_find(m03);

    std::cout << std::endl << "PEDIGREE ITERATION" << std::endl;

    iprint("members",          P  . member_begin      (),    P  . member_end      ());
    iprint("const members",    CP . member_begin      (),    CP . member_end      ());
    iprint("families",         P  . family_begin      (),    P  . family_end      ());
    iprint("const families",   CP . family_begin      (),    CP . family_end      ());
    iprint("subpeds",          P  . subpedigree_begin (),    P  . subpedigree_end ());
    iprint("const subpeds",    CP . subpedigree_begin (),    CP . subpedigree_end ());
    iprint("unconns",          P  . unconnected_begin (),    P  . unconnected_end ());
    iprint("const unconns",    CP . unconnected_begin (),    CP . unconnected_end ());
    iprint("parents 03",       P  . parent_begin      (m03), P  . parent_end      ());
    iprint("const parents 03", P  . parent_begin      (m03), P  . parent_end      ());

    std::cout << std::endl << "SUBPED ITERATION (first)" << std::endl;

    iprint("members",        sf  -> member_begin(), sf  -> member_end());
    iprint("const members",  csf -> member_begin(), csf -> member_end());
    iprint("families",       sf  -> family_begin(), sf  -> family_end());
    iprint("const families", csf -> family_begin(), csf -> family_end());

    std::cout << std::endl << "FAMILY ITERATION (m03/m04)" << std::endl;

    iprint("offspring",       pf  -> offspring_begin (), pf  -> offspring_end ());
    iprint("const offspring", cpf -> offspring_begin (), cpf -> offspring_end ());
    iprint("parents",         pf  -> parent_begin    (), pf  -> parent_end    ());
    iprint("const parents",   cpf -> parent_begin    (), cpf -> parent_end    ());

    std::cout << std::endl << "MEMBER ITERATION (m03)" << std::endl;

    iprint("mates",           pm  -> mate_begin      (),                    pm  -> mate_end      ());
    iprint("const mates",     cpm -> mate_begin      (),                    cpm -> mate_end      ());
    iprint("offspring",       pm  -> offspring_begin (*P.member_find(m04)), pm  -> offspring_end ());
    iprint("const offspring", cpm -> offspring_begin (*P.member_find(m04)), cpm -> offspring_end ());
    iprint("offspring",       pm  -> offspring_begin (m04),                 pm  -> offspring_end ());
    iprint("const offspring", cpm -> offspring_begin (m04),                 cpm -> offspring_end ());
    iprint("parents",         pm  -> parent_begin    (),                    pm  -> parent_end    ());
    iprint("const parents",   cpm -> parent_begin    (),                    cpm -> parent_end    ());
    iprint("progeny",         pm  -> progeny_begin   (),                    pm  -> progeny_end   ());
    iprint("const progeny",   cpm -> progeny_begin   (),                    cpm -> progeny_end   ());
    iprint("siblings",        pm  -> sibling_begin   (),                    pm  -> sibling_end   ());
    iprint("const siblings",  cpm -> sibling_begin   (),                    cpm -> sibling_end   ());
 
}


//----------------------------------------------------------------------------
//  Function:   test1()
//
//  Scenario:   single 2-parent lineage
//----------------------------------------------------------------------------
//
void    test1()
{   
    pedigree_type   P("P1");

    P.add_member(m01);
    P.add_member(m02);
    P.add_member(m03);
    P.add_lineage(m03, m02, m01);

    P.build();

    print("test1", P);
}


//----------------------------------------------------------------------------
//  Function:   test2()
//
//  Scenario:   multiple 2-parent lineages, 
//              multiple matched 1-parent lineages
//----------------------------------------------------------------------------
//
void    test2()
{   
    pedigree_type   P("P2");

    P.add_member(m01);
    P.add_member(m02);
    P.add_member(m03);
    P.add_member(m04);
    P.add_member(m05);
    P.add_member(m06);
    P.add_lineage(m03, m02, m01);
    P.add_lineage(m04, m01, m02);
    P.add_lineage(m05, m01);
    P.add_lineage(m05, m02);
    P.add_lineage(m06, m01);
    P.add_lineage(m06, m02);

    P.build();

    print2("test2", P);
}


//----------------------------------------------------------------------------
//  Function:   test3()
//
//  Scenario:   matched 1-parent lineages split by multiple sibships
//----------------------------------------------------------------------------
//
void    test3()
{   
    pedigree_type   P("P3");

    P.add_member(m01);
    P.add_member(m02);
    P.add_member(m03);
    P.add_member(m04);
    P.add_member(m05);
    P.add_member(m06);
    P.add_lineage(m03, m01);
    P.add_sibship(m03, m04);
    P.add_sibship(m04, m05);
    P.add_sibship(m05, m06);
    P.add_lineage(m06, m02);

    P.build();

    print("test3", P);
}


//----------------------------------------------------------------------------
//  Function:   test4()
//
//  Scenario:   single 2-parent lineage, multiple sibships
//----------------------------------------------------------------------------
//
void    test4()
{   
    pedigree_type   P("P4");

    P.add_member(m01);
    P.add_member(m02);
    P.add_member(m03);
    P.add_member(m04);
    P.add_member(m05);
    P.add_member(m06);
    P.add_member(m07);
    P.add_member(m08);
    P.add_lineage(m03, m01, m02);
    P.add_sibship(m06, m05);
    P.add_sibship(m05, m04);
    P.add_sibship(m04, m03);
    P.add_lineage(m08, m01, m07);

    P.build();

    print("test4", P);

    pedigree_type::member_pointer   m  = P.member_find(m05);
    p_iter      pf = m->parent_begin();
    p_iter      pl = m->parent_end();

    std::cout << std::endl << " parents of " << m->name() << ":";
    for (;  pf != pl;  ++pf)
    {
        std::cout << " " << pf->name();
    }
    std::cout << std::endl;

    m  = P.member_find(m01);
    pf = m->parent_begin();
    pl = m->parent_end();

    std::cout << " parents of " << m->name() << ":";
    for (;  pf != pl;  ++pf)
    {
        std::cout << " " << pf->name();
    }
    std::cout << std::endl;

    m = P.member_find(m04);

    s_iter  sf = m->sibling_begin();
    s_iter  sl = m->sibling_end();

    std::cout << "siblings of " << m->name() << ":";
    for (;  sf != sl;  ++sf)
    {
        std::cout << " " << sf->name();
    }
    std::cout << std::endl;

    m = P.member_find(m01);
    
    a_iter  af = m->mate_begin();
    a_iter  al = m->mate_end();

    std::cout << "   mates of " << m->name() << ":";
    for (;  af != al;  ++af)
    {
        std::cout << " " << af->mate().name();
    }
    std::cout << std::endl;

    m  = P.member_find(m02);
    af = m->mate_begin();
    al = m->mate_end();

    std::cout << "   mates of " << m->name() << ":";
    for (;  af != al;  ++af)
    {
        std::cout << " " << af->mate().name();
    }
    std::cout << std::endl;

    m  = P.member_find(m07);
    af = m->mate_begin();
    al = m->mate_end();

    std::cout << "   mates of " << m->name() << ":";
    for (;  af != al;  ++af)
    {
        std::cout << " " << af->mate().name();
    }
    std::cout << std::endl;

    m  = P.member_find(m03);
    af = m->mate_begin();
    al = m->mate_end();

    std::cout << "   mates of " << m->name() << ":";
    for (;  af != al;  ++af)
    {
        std::cout << " " << af->mate().name();
    }
    std::cout << std::endl;

    m = P.member_find(m01);
    y_iter  yf = m->progeny_begin();
    y_iter  yl = m->progeny_end();

    std::cout << " progeny of " << m->name() << ":";
    for (; yf != yl;  ++yf)
    {
        std::cout << " " << yf->name();
    }
    std::cout << std::endl;
}


//----------------------------------------------------------------------------
//  Function:   test5()
//
//  Scenario:   multiple 2-parent lineages
//----------------------------------------------------------------------------
//
void    test5()
{   
    pedigree_type   P("P5");

    P.add_member(m01);
    P.add_member(m02);
    P.add_member(m03);
    P.add_member(m04);
    P.add_member(m05);
    P.add_member(m06);
    P.add_member(m07);
    P.add_lineage(m03, m01, m02);
    P.add_lineage(m06, m04, m05);
    P.add_lineage(m07, m05, m04);

    P.build();

    print("test5", P);
    print2("test5", P);

    P.add_member(m08);
    P.add_lineage(m08, m03, m06);

    P.build();

    print("test5", P);
    print2("test5", P);
}


//----------------------------------------------------------------------------
//  Function:   test6()
//
//  Scenario:   total offspring iteration
//----------------------------------------------------------------------------
//
void    test6()
{   
    pedigree_type   P("P6");

    P.add_member(m01);
    P.add_member(m02);
    P.add_member(m03);
    P.add_member(m04);
    P.add_member(m05);
    P.add_member(m06);
    P.add_lineage(m03, m02, m01);
    P.add_lineage(m04, m01, m02);
    P.add_lineage(m05, m01);
    P.add_lineage(m05, m02);
    P.add_lineage(m06, m01);
    P.add_lineage(m06, m02);

    P.build();

    print("test6", P);
    print2("test6", P);

    pedigree_type::member_pointer   m = P.member_find(m01);
    y_iter  yf = m->progeny_begin();
    y_iter  yl = m->progeny_end();

    std::cout << "progeny: ";
    for (; yf != yl;  ++yf)
    {
        std::cout << " " << yf->name();
    }
    std::cout << std::endl;


    const pedigree_type&    CP = P;

    CP.subpedigree_index(0);
    CP.member_index(0);
    CP.family_index(0);

    P.subpedigree_index(0);
    P.member_index(0);
    P.family_index(0);

//    P.subpedigree_index_swap(0, 1);
    P.member_index_swap(0, 1);
//    P.family_index_swap(0, 1);
}


void test7()
{
    multipedigree_type  P;
    mp_iter it;
    string  p1("P7A");
    string  p2("P7B");

    it = P.add_member(p1, m01);
    it = P.add_member(p1, m02);
    it = P.add_member(it, m03);
    it = P.add_member(it, m04);
    it = P.add_member(it, m05);
    it = P.add_member(p1, m06);
    it = P.add_lineage(p1, m03, m02, m01);
    it = P.add_lineage(it, m04, m01, m02);
    it = P.add_lineage(p1, m05, m01);
    it = P.add_lineage(it, m05, m02);
    it = P.add_lineage(p1, m06, m01);
    it = P.add_lineage(it, m06, m02);

    P.build();

    print("test 7", *it);

    const multipedigree_type&   CP = P;

    CP.pedigree_index(0);
    P.pedigree_index(0);

    CP.member_index(0);
    P.member_index(0);

    it = P.add_member(p2, m01);
    it = P.add_member(p2, m02);
    it = P.add_member(it, m03);
    it = P.add_lineage(p2, m03, m02, m01);

    P.build();
    print("test 7b", *it);

    std::cout << P.pedigree_index(0).name() << std::endl;
    std::cout << P.pedigree_index(1).name() << std::endl;
    P.pedigree_index_swap(0, 1);
    std::cout << P.pedigree_index(0).name() << std::endl;
    std::cout << P.pedigree_index(1).name() << std::endl;

    std::cout << P.member_index(0).name() << std::endl;
    std::cout << P.member_index(1).name() << std::endl;
    P.member_index_swap(0, 1);
    std::cout << P.member_index(0).name() << std::endl;
    std::cout << P.member_index(1).name() << std::endl;
}

void    test_bad_sibship()
{   
    pedigree_type   P("P1");
    
    P.add_sibship("", "sib1missing");
    P.add_sibship("sib2missing", "");
    
    P.add_member(m01);
    P.add_member(m02);
    P.add_member(m11, SEX_MALE);
    P.add_member(m12, SEX_FEMALE);
    P.add_member(m13, SEX_MALE);
    P.add_member(m14, SEX_FEMALE);

    P.add_lineage(m01, m11, m12);
    P.add_lineage(m02, m13, m14);
    P.add_sibship(m01, m02);

    P.build();

    assert(P.error_count() == 3);
    
    pedigree_type::error_iterator err = P.error_begin();
    
    assert(err->state == error_info::bad_sibship);
    assert(err->name1 == "");
    assert(err->name2 == "sib1missing");
    assert(err->name3 == "");

    ++err;

    assert(err->state == error_info::bad_sibship);
    assert(err->name1 == "sib2missing");
    assert(err->name2 == "");
    assert(err->name3 == "");

    ++err;

    assert(err->state == error_info::bad_sibship);
    assert(err->name1 == "m01");
    assert(err->name2 == "m02");
    assert(err->name3 == "");
}

void    test_bad_marriage()
{   
    pedigree_type   P("P1");
    
    P.add_marriage("", "1missing");
    P.add_marriage("2missing", "");
    
    P.add_member(m01);
    P.add_member(m11, SEX_MALE);
    P.add_member(m12, SEX_FEMALE);

    P.add_lineage(m01, m11, m12);

    P.build();

    assert(P.error_count() == 2);
    
    pedigree_type::error_iterator err = P.error_begin();
    
    assert(err->state == error_info::bad_marriage);
    assert(err->name1 == "");
    assert(err->name2 == "1missing");
    assert(err->name3 == "");

    ++err;

    assert(err->state == error_info::bad_marriage);
    assert(err->name1 == "2missing");
    assert(err->name2 == "");
    assert(err->name3 == "");
}

void    test_bad_lineage()
{   
    pedigree_type   P("P1");
    
    P.add_lineage("", "1missing");
    P.add_lineage("2missing", "");
    P.add_lineage("", "1missing", "B");
    P.add_lineage("2missing", "", "B");
    P.add_lineage("3missing", "A", "");
    
    P.add_member(m01);
    P.add_member(m11);
    P.add_member(m12);

    P.add_lineage(m01, m01);
    P.add_lineage(m01, m01, m12);
    P.add_lineage(m11, m12, m11);

    P.build();

    assert(P.error_count() == 8);
    
    pedigree_type::error_iterator err = P.error_begin();
    
    assert(err->state == error_info::bad_lineage);
    assert(err->name1 == "");
    assert(err->name2 == "1missing");
    assert(err->name3 == "");

    ++err;

    assert(err->state == error_info::bad_lineage);
    assert(err->name1 == "2missing");
    assert(err->name2 == "");
    assert(err->name3 == "");

    ++err;

    assert(err->state == error_info::bad_lineage);
    assert(err->name1 == "");
    assert(err->name2 == "1missing");
    assert(err->name3 == "B");

    ++err;

    assert(err->state == error_info::bad_lineage);
    assert(err->name1 == "2missing");
    assert(err->name2 == "");
    assert(err->name3 == "B");

    ++err;

    assert(err->state == error_info::bad_lineage);
    assert(err->name1 == "3missing");
    assert(err->name2 == "A");
    assert(err->name3 == "");

    ++err;

    assert(err->state == error_info::bad_lineage);
    assert(err->name1 == "m01");
    assert(err->name2 == "m01");
    assert(err->name3 == "");

    ++err;

    assert(err->state == error_info::bad_lineage);
    assert(err->name1 == "m01");
    assert(err->name2 == "m01");
    assert(err->name3 == "m12");

    ++err;

    assert(err->state == error_info::bad_lineage);
    assert(err->name1 == "m11");
    assert(err->name2 == "m12");
    assert(err->name3 == "m11");
}

void    test_bad_gender()
{   
    pedigree_type   P("P1");
    
    P.add_member(m01, SEX_MALE);

    P.add_member(m01, SEX_MALE);
    P.add_member(m01, SEX_XMALE);
    P.add_member(m01, SEX_FEMALE);
    P.add_member(m01, SEX_XFEMALE);

    P.add_member(m02, SEX_XMALE);

    P.add_member(m02, SEX_MALE);
    P.add_member(m02, SEX_XMALE);
    P.add_member(m02, SEX_FEMALE);
    P.add_member(m02, SEX_XFEMALE);

    P.add_member(m03, SEX_FEMALE);

    P.add_member(m03, SEX_MALE);
    P.add_member(m03, SEX_XMALE);
    P.add_member(m03, SEX_FEMALE);
    P.add_member(m03, SEX_XFEMALE);

    P.add_member(m04, SEX_XFEMALE);

    P.add_member(m04, SEX_MALE);
    P.add_member(m04, SEX_XMALE);
    P.add_member(m04, SEX_FEMALE);
    P.add_member(m04, SEX_XFEMALE);

    P.build();
    
    assert(P.error_count() == 6);
    
    pedigree_type::error_iterator err = P.error_begin();
    
    assert(err->state == error_info::bad_gender);
    assert(err->name1 == m01);
    assert(err->name2 == "male");
    assert(err->name3 == "female");

    ++err;

    assert(err->state == error_info::bad_gender);
    assert(err->name1 == m02);
    assert(err->name2 == "male");
    assert(err->name3 == "female");

    ++err;

    assert(err->state == error_info::bad_gender);
    assert(err->name1 == m03);
    assert(err->name2 == "female");
    assert(err->name3 == "male");

    ++err;

    assert(err->state == error_info::bad_gender);
    assert(err->name1 == m03);
    assert(err->name2 == "male");
    assert(err->name3 == "female");

    ++err;

    assert(err->state == error_info::bad_gender);
    assert(err->name1 == m04);
    assert(err->name2 == "female");
    assert(err->name3 == "male");

    ++err;

    assert(err->state == error_info::bad_gender);
    assert(err->name1 == m04);
    assert(err->name2 == "male");
    assert(err->name3 == "female");

    ++err;

}

void    test_same_sex_marriage()
{   
    pedigree_type   P("P1");
    
    P.add_member(m01, SEX_MALE);
    P.add_member(m02, SEX_MALE);
    P.add_member(m03, SEX_MALE);
    
    P.add_lineage(m01, m02, m03);
    
    P.add_member(m04, SEX_FEMALE);
    P.add_member(m05, SEX_FEMALE);
    P.add_member(m06, SEX_FEMALE);
    
    P.add_lineage(m04, m05, m06);
    
    P.add_member(m07, SEX_MALE);
    P.add_member(m08, SEX_XMALE);
    P.add_member(m09, SEX_MALE);
    
    P.add_lineage(m07, m08, m09);

    P.add_member(m10, SEX_FEMALE);
    P.add_member(m11, SEX_XFEMALE);
    P.add_member(m12, SEX_FEMALE);
    
    P.add_lineage(m10, m11, m12);

    P.build();
    
    assert(P.error_count() == 4);

    pedigree_type::error_iterator err = P.error_begin();

    assert(err->state == error_info::same_sex_marriage);
    assert(err->name1 == m02);
    assert(err->name2 == m03);
    assert(err->name3 == "male");

    ++err;

    assert(err->state == error_info::same_sex_marriage);
    assert(err->name1 == m05);
    assert(err->name2 == m06);
    assert(err->name3 == "female");

    ++err;

    assert(err->state == error_info::same_sex_marriage);
    assert(err->name1 == m08);
    assert(err->name2 == m09);
    assert(err->name3 == "male");

    ++err;

    assert(err->state == error_info::same_sex_marriage);
    assert(err->name1 == m11);
    assert(err->name2 == m12);
    assert(err->name3 == "female");

    ++err;
}

void    test_arb_assignment()
{   
    pedigree_type   P("P1");
    
    P.add_member(m01, SEX_ARB);
    P.add_member(m02, SEX_ARB);
    P.add_member(m03, SEX_ARB);
    P.add_member(m04, SEX_ARB);
    P.add_member(m05, SEX_ARB);
    P.add_member(m06, SEX_ARB);
    P.add_member(m07, SEX_ARB);
    P.add_member(m08, SEX_ARB);
    
    P.add_lineage(m05, m01, m02);
    P.add_lineage(m06, m02, m03);
    P.add_lineage(m07, m03, m04);
    P.add_lineage(m08, m04, m01);

    P.add_member(m10, SEX_ARB);
    P.add_member(m11, SEX_ARB);
    P.add_member(m12, SEX_ARB);

    P.add_lineage(m10, m11, m12);
    
    P.build();
    
    for(size_t i = 0; i != P.member_count(); ++i)
      cout << P.member_index(i).name() << " -> " << P.member_index(i).get_detailed_sex() << endl;
    
    assert(P.error_count() == 0);

}

void    test_arb_bad_loop()
{   
    pedigree_type   P("P1");
    
    P.add_member(m01, SEX_ARB);
    P.add_member(m02, SEX_ARB);
    P.add_member(m03, SEX_ARB);

    P.add_member(m05, SEX_ARB);
    P.add_member(m06, SEX_ARB);
    P.add_member(m07, SEX_ARB);
    
    P.add_lineage(m05, m01, m02);
    P.add_lineage(m06, m02, m03);
    P.add_lineage(m07, m03, m01);

    P.add_member(m11, SEX_MISSING);
    P.add_member(m12, SEX_MISSING);
    P.add_member(m13, SEX_MISSING);

    P.add_member(m15, SEX_ARB);
    P.add_member(m16, SEX_ARB);
    P.add_member(m17, SEX_ARB);
    
    P.add_lineage(m15, m11, m12);
    P.add_lineage(m16, m12, m13);
    P.add_lineage(m17, m13, m11);

    P.build();
    
    for(size_t i = 0; i != P.member_count(); ++i)
      cout << P.member_index(i).name() << " -> " << P.member_index(i).get_detailed_sex() << endl;
    
    assert(P.error_count() == 5);

    pedigree_type::error_iterator err = P.error_begin();
    
    assert(err->state == error_info::no_sex_parents);
    
    ++err;

    assert(err->state == error_info::no_sex_parents);
    
    ++err;

    assert(err->state == error_info::no_sex_parents);
    
    ++err;

    assert(err->state == error_info::bad_marriage_loop);
    assert(err->name1 == m01);
    assert(err->name2 == "");
    assert(err->name3 == "");
    
    ++err;
    
    assert(err->state == error_info::bad_marriage_loop);
    assert(err->name1 == m11);
    assert(err->name2 == "");
    assert(err->name3 == "");
}

  
void    test_gender_inferred()
{
    cout << "Test Gender Inferred: " << endl;
  
    pedigree_type   P("P1");
    
    P.add_member(m01, SEX_MALE);
    P.add_member(m02, SEX_ARB);
    P.add_member(m03, SEX_FEMALE);
    P.add_member(m04, SEX_ARB);
    P.add_member(m05, SEX_FEMALE);
    P.add_member(m06, SEX_ARB);
    P.add_member(m07, SEX_MISSING);
    P.add_member(m08, SEX_MISSING);
  
    P.add_member(m10, SEX_ARB);
    P.add_member(m11, SEX_ARB);
    P.add_member(m12, SEX_ARB);
    P.add_member(m13, SEX_ARB);
    P.add_member(m14, SEX_ARB);
  
    P.add_lineage(m10, m01, m02);
    P.add_lineage(m11, m03, m04);
    P.add_lineage(m12, m05, m06);
    P.add_lineage(m13, m06, m07);
    P.add_lineage(m14, m07, m08);

    P.build();

    for(size_t i = 0; i != P.member_count(); ++i)
      cout << P.member_index(i).name() << " -> " << P.member_index(i).get_detailed_sex() << endl;
    
    assert(P.error_count() == 0);
    
    assert(P.warning_count() == 2);
}
  
void    test_no_sex_parents()
{   
    pedigree_type   P("P1");
    
    P.add_member(m01);
    P.add_member(m11);
    P.add_member(m12);

    P.add_lineage(m01, m11, m12);

    P.build();

    assert(P.error_count() == 1);
    
    pedigree_type::error_iterator err = P.error_begin();
    
    assert(err->state == error_info::no_sex_parents);
    assert(err->name1 == m11);
    assert(err->name2 == m12);
    assert(err->name3 == "");
}

void test_get_moth_fath()
{
    pedigree_base P("P1");
  
    P.add_member(m01,SEX_FEMALE);
    P.add_member(m02,SEX_MALE);
    P.add_member(m03,SEX_MISSING);
    P.add_member(m04,SEX_FEMALE);
    P.add_member(m05,SEX_MALE);
    P.add_member(m06,SEX_MISSING);
    P.add_member(m07,SEX_MISSING);
    P.add_member(m08,SEX_MISSING);

    P.add_member(m11);
    P.add_member(m12);
    P.add_member(m13);
    P.add_member(m14);

    P.add_lineage(m11, m01, m02);
    P.add_lineage(m12, m03, m04);
    P.add_lineage(m13, m05, m06);
    P.add_lineage(m14, m07, m08);

    P.build();
    
    member_base* mem = 0;
    family_base* fam = 0;
    
    mem = P.member_find(m11);
    fam = mem->family();
    
    assert(mem && fam);
    
    assert(mem->get_mother() && mem->get_mother()->name() == m01);
    assert(mem->get_father() && mem->get_father()->name() == m02);
    assert(fam->get_mother() && fam->get_mother()->name() == m01);
    assert(fam->get_father() && fam->get_father()->name() == m02);

    mem = P.member_find(m12);
    fam = mem->family();
    
    assert(mem && fam);
    
    assert(mem->get_mother() && mem->get_mother()->name() == m04);
    assert(mem->get_father() && mem->get_father()->name() == m03);
    assert(fam->get_mother() && fam->get_mother()->name() == m04);
    assert(fam->get_father() && fam->get_father()->name() == m03);
    
    mem = P.member_find(m13);
    fam = mem->family();
    
    assert(mem && fam);
    
    assert(mem->get_mother() && mem->get_mother()->name() == m06);
    assert(mem->get_father() && mem->get_father()->name() == m05);
    assert(fam->get_mother() && fam->get_mother()->name() == m06);
    assert(fam->get_father() && fam->get_father()->name() == m05);

    mem = P.member_find(m14);
    fam = mem->family();
    
    assert(mem && fam);
    
    assert(!mem->get_mother());
    assert(!mem->get_father());
    assert(!fam->get_mother());
    assert(!fam->get_father());
}
  
int main()
{
    test_sex_code();
  
    test0();
    test1();
    test2();
    test3();
    test4();
    test5();
    test6();
    test7();
    
    test_bad_sibship();
    test_bad_marriage();
    test_bad_lineage();
    test_bad_gender();
    test_same_sex_marriage();

    test_arb_assignment();
    test_arb_bad_loop();
    
    test_gender_inferred();
    
    test_no_sex_parents();
    
    test_get_moth_fath();

    return 0;
}

