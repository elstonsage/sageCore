#include <iomanip>

#ifdef _MSC_VER
    #include <app/SAGEconfig.h>
    #pragma hdrstop
#endif

#include "mlocus/phmodel.h"

using namespace SAGE::MLOCUS;

using std::cout;
using std::endl;
using std::boolalpha;
using std::copy;
using std::fixed;
using std::left;
using std::right;
using std::setw;
using std::setprecision;
using std::showpoint;

//- Allele names.
//
extern std::string  a, b, c, d;

//- Phenotype names.
//
std::string  red("red"), blue("blue"), green("green");
std::string  orange("orange"), white("white"), black("black");
std::string  gray("gray"), yellow("yellow"), purple("purple"), violet("violet");

void test(const phenotype_model& m, const char* str)
{
    ptrdiff_t   diff;
    size_t        count, i;
    bool        good;

    cout << "************************************" << endl;
    cout << "Testing phenotype model " << str << endl;

    diff  = m.phenotype_end() - m.phenotype_begin();
    count = m.phenotype_count();
    if (count == (size_t) diff)
        cout << "..........phenotype count OK" << endl;
    else
        cout << "..........phenotype count FAIL" << endl;

    phenotype_model::phenotype_iterator  pf = m.phenotype_begin();       
    phenotype_model::phenotype_iterator  pl = m.phenotype_end();     

    for (i = 0, good = true;  pf != pl;  ++pf, ++i)
    {
        if (*pf != m.get_phenotype(pf->id()))
        {
            cout << "..........phenotype index FAIL at: " << i << endl;
            good = false;
        }

        if (pf->id() != m.get_phenotype_id(pf->name()))
        {
            cout << ".........phenotype lookup FAIL at: " << i << endl;
            cout << "  pf->id: " << pf->id() << "  lookup->id:" << m.get_phenotype_id(pf->name()) << endl;
            good = false;
        }

        if (m.get_phenotype(pf->name()).id() != m.get_phenotype_id(pf->name()))
        {
            cout << ".........phenotype lookup FAIL at: " << i << endl;
            good = false;
        }
    }
    if (good)
    {
        cout << ".........phenotype lookup OK" << endl;
    }
    cout << endl;
}


void print(const phenotype_model& m, const char* str)
{
    cout << "************************************" << endl;
    cout << "Printing phenotype model " << str << endl;
    cout << "  has generated phenotypes: " << boolalpha << m.has_generated_phenotypes() << endl;
    cout << "  has external phenotypes:  " << boolalpha << m.has_external_phenotypes() << endl;

    cout << " phenotypes: ";
    phenotype_model::phenotype_iterator  pf = m.phenotype_begin();
    phenotype_model::phenotype_iterator  pl = m.phenotype_end();
    for (;  pf != pl;  ++pf)
    {
        cout << " " << pf->name();
    }
    cout << endl << endl;
}


void gg(genotype_model& gm)
{
    phenotype_model  m(gm);

    test(m, "GG - symmetric");
    print(m, "GG - symmetric");
}


void hh(genotype_model& gm)
{
    phenotype_model  m;

    m.add_phenotype(blue);
    m.add_phenotype(orange);
    m.add_phenotype(green);
    m.add_phenotype(red);
    m.add_phenotype(purple);
    m.add_phenotype(black);
    m.add_phenotype(gray);
    m.add_phenotype(yellow);
    m.add_phenotype(white);
    m.add_phenotype(violet);

    test(m, "HH - Non-Symmetric");
    print(m, "HH - Non-Symmetric");
}

void ii(genotype_model& gm)
{
    phenotype_model  m(gm);

    m.add_phenotype(blue);
    m.add_phenotype(orange);
    m.add_phenotype(green);
    m.add_phenotype(red);
    m.add_phenotype(purple);
    m.add_phenotype(black);
    m.add_phenotype(gray);
    m.add_phenotype(yellow);
    m.add_phenotype(white);
    m.add_phenotype(violet);

    test(m, "II - both sym and non-symmetric");
    print(m, "II - both sym and non-symmetric");
}

void test_x_linked_phenotypes()
{
  genotype_model gm;
  
  gm.set_model_type(X_LINKED);
  
  gm.add_allele("A", 0.2);
  gm.add_allele("B", 0.5);
  gm.add_allele("C", 0.3);
  
  phenotype_model m(gm);

  assert(m.get_phenotype_id("A/~Y") == m.get_phenotype_id("~Y/A"));
  assert(m.get_phenotype_id("A/~Y") != NPOS);

  assert(m.get_phenotype_id("A<~Y") == m.get_phenotype_id("~Y>A"));
  assert(m.get_phenotype_id("A<~Y") != NPOS);

  assert(m.get_phenotype_id("~Y<A") == NPOS);
  assert(m.get_phenotype_id("A>~Y") == NPOS);
  assert(m.get_phenotype_id("~Y/~Y") == NPOS);
  assert(m.get_phenotype_id("~Y<~Y") == NPOS);
  assert(m.get_phenotype_id("~Y>~Y") == NPOS);

  assert(m.get_phenotype_id("A/~Y") == m.get_phenotype_id("A/*missing"));
  assert(m.get_phenotype_id("A/~Y") == m.get_phenotype_id("*missing/A"));
  assert(m.get_phenotype_id("A/~Y") == m.get_phenotype_id("A/A(male)"));

  assert(m.get_phenotype_id("A<~Y") == m.get_phenotype_id("A<*missing"));
  assert(m.get_phenotype_id("A<~Y") == m.get_phenotype_id("*missing>A"));
  assert(m.get_phenotype_id("A<~Y") == m.get_phenotype_id("A<A(male)"));
  assert(m.get_phenotype_id("A<~Y") == m.get_phenotype_id("A>A(male)"));

  assert(m.get_phenotype_id("A/~Y") != m.get_phenotype_id("~Y/B"));
  assert(m.get_phenotype_id("A/~Y") != m.get_phenotype_id("B/~Y"));
  assert(m.get_phenotype_id("A<~Y") != m.get_phenotype_id("B<~Y"));
  assert(m.get_phenotype_id("~Y>A") != m.get_phenotype_id("~Y>B"));

  assert(m.get_phenotype_id("~Y/B") != NPOS);
  assert(m.get_phenotype_id("B/~Y") != NPOS);
  assert(m.get_phenotype_id("B<~Y") != NPOS);
  assert(m.get_phenotype_id("~Y>B") != NPOS);

  assert(m.get_phenotype_id("~Y<B") == NPOS);
  assert(m.get_phenotype_id("B>~Y") == NPOS);
}
void test_y_linked_phenotypes()
{
  genotype_model gm;
  
  gm.set_model_type(Y_LINKED);
  
  gm.add_allele("A", 0.2);
  gm.add_allele("B", 0.5);
  gm.add_allele("C", 0.3);
  
  phenotype_model m(gm);

  assert(m.get_phenotype_id("A/~X") == m.get_phenotype_id("~X/A"));
  assert(m.get_phenotype_id("A/~X") != NPOS);

  assert(m.get_phenotype_id("A>~X") == m.get_phenotype_id("~X<A"));
  assert(m.get_phenotype_id("A>~X") != NPOS);

  assert(m.get_phenotype_id("~X>A") == NPOS);
  assert(m.get_phenotype_id("A<~X") == NPOS);

  assert(m.get_phenotype_id("~X/~X") != NPOS);
  assert(m.get_phenotype_id("~X<~X") != NPOS);
  assert(m.get_phenotype_id("~X>~X") != NPOS);

  assert(m.get_phenotype_id("A/~X") == m.get_phenotype_id("A/*missing"));
  assert(m.get_phenotype_id("A/~X") == m.get_phenotype_id("A/A"));

  assert(m.get_phenotype_id("A>~X") == m.get_phenotype_id("A>*missing"));
  assert(m.get_phenotype_id("A>~X") == m.get_phenotype_id("A<A"));
  assert(m.get_phenotype_id("A>~X") == m.get_phenotype_id("A>A"));

  assert(m.get_phenotype_id("A/~X") != m.get_phenotype_id("~X/B"));
  assert(m.get_phenotype_id("A/~X") != m.get_phenotype_id("B/~X"));
  assert(m.get_phenotype_id("A>~X") != m.get_phenotype_id("B>~X"));
  assert(m.get_phenotype_id("~X<A") != m.get_phenotype_id("~X<B"));

  assert(m.get_phenotype_id("~X/B") != NPOS);
  assert(m.get_phenotype_id("B/~X") != NPOS);
  assert(m.get_phenotype_id("B>~X") != NPOS);
  assert(m.get_phenotype_id("~X<B") != NPOS);

  assert(m.get_phenotype_id("~X>B") == NPOS);
  assert(m.get_phenotype_id("B<~X") == NPOS);
}

