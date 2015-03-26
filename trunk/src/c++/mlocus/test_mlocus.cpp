#include <iomanip>

#include <string>

#include "LSF/LSFinit.h"

#define protected public

#include "mlocus/genotype.h"

using namespace SAGE::MLOCUS;

// Genotype special tests

void test_sex_type();

// Phenotype Model Tests (in test2.cpp):

void    gg(genotype_model& gm);
void    hh(genotype_model& gm);
void    ii(genotype_model& gm);

void test_x_linked_phenotypes();
void test_y_linked_phenotypes();

// Penetrance Model Tests (in test3.cpp):

void    aaa(genotype_model& gm);
void    bbb(genotype_model& gm);
void    ccc(genotype_model& gm);
void    ggg(genotype_model& gm);
void    hhh(genotype_model& gm);
void    iii(genotype_model& gm);
void    jjj(genotype_model& gm);
void    kkk(genotype_model& gm);
void    lll(genotype_model& gm);
void    mmm(genotype_model& gm);

void    test_sexed_penetrance();
void    test_sex_type_change();
void    test_dynamic_alleles();

// File i/o test

void test_file();

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
std::string  a("a"), b("b"), c("c"), d("d");

std::string  newseps(":{}");

void test(const genotype_model& m, const char* str)
{
    cout << "************************************" << endl;
    cout << "Testing genotype model " << str << endl;

    genotype_model_test(cout, m);
}


void print(const genotype_model& m, const char* str)
{
    cout << "************************************" << endl;
    cout << "Printing genotype model " << str << endl;

    genotype_model_print(cout, m);
    cout << endl;
}


void e()
{
    genotype_model  m;

    m.add_allele(a, 0.25);
    test(m, "E");
    print(m, "E");

    m.add_allele(b, 0.50);
    m.add_allele(d, 1.00);
    m.add_allele(c, 0.75);
    test(m, "E");
    print(m, "E");
}


void f()
{
    genotype_model  m;

    m.add_allele(a, 0.25);
    m.add_allele(b, 0.50);
    m.add_allele(c, 0.75);
    m.add_allele(d, 1.00);
    test(m, "F");
    print(m, "F");

    m.mark_for_remap(c);
    m.mark_for_remap(d);
    m.remap();

    test(m, "F");
    print(m, "F");
}


void g()
{
    genotype_model  m;

    m.add_allele(a, 0.25);
    m.add_allele(b, 0.50);
    m.add_allele(c, 0.75);
    m.add_allele(d, 1.00);
    test(m, "G");
    print(m, "G");

    genotype_model  m2(m);
    test(m2, "G");
    print(m2, "G");

    genotype_model  m3;

    m3 = m;
    test(m3, "G");
    print(m3, "G");

    m2.set_unphased_separator(':');
    test(m2, "G");
    print(m2, "G");

    m3.set_forward_separator('{');
    test(m3, "G");
    print(m3, "G");

    m.set_separators(newseps);
    test(m, "G");
    print(m, "G");
}


void h()
{
    genotype_model  m;

    m.add_allele(a, 0.25);
    m.add_allele(b, 0.50);
    m.add_allele(c, 0.75);
    m.add_allele(d, 1.00);

    gg(m);
    hh(m);
    ii(m);
    
    test_x_linked_phenotypes();
    test_y_linked_phenotypes();
}

void _i()
{
    genotype_model  m;

    m.add_allele(a, 0.25);
    m.add_allele(b, 0.50);
    m.add_allele(c, 0.75);
    m.add_allele(d, 1.00);

    aaa(m);
    bbb(m);
    ccc(m);
    ggg(m);
    hhh(m);
    iii(m);
    jjj(m);
    kkk(m);
    lll(m);
    mmm(m);
    
    test_sexed_penetrance();
    test_sex_type_change();
    test_dynamic_alleles();
}

void test_sex_type()
{
  // Create one with no type (ie, autosomal)
  genotype_model g;
  
  g.add_allele("A", 0.5);
  g.add_allele("B", 0.5);
  
  assert(g.allele_count() == 2);
  assert(g.unphased_genotype_count() == 3);
  assert(g.phased_genotype_count() == 4);
  
  assert(!g.get_sex_specific_allele().is_valid());
  
  assert(!g.is_x_linked());
  assert(!g.is_y_linked());
  assert( g.is_autosomal());
  
  test(g, "Autosomal");
  print(g, "Autosomal");

  g.set_model_type(X_LINKED);
  
  assert(g.allele_count() == 2);
  assert(g.unphased_genotype_count() == 5);
  assert(g.phased_genotype_count() == 6);
  
  assert(g.get_sex_specific_allele().is_valid());
  assert(g.get_sex_specific_allele().name() == "~Y");
  assert(g.get_sex_specific_allele().frequency() == 1.0);
  
  assert( g.get_sex_specific_allele().is_null_y_allele());
  assert(!g.get_sex_specific_allele().is_null_x_allele());
  
  assert( g.is_x_linked());
  assert(!g.is_y_linked());
  assert(!g.is_autosomal());

  test(g, "X-Linked");
  print(g, "X-Linked");

  g.set_model_type(Y_LINKED);
  
  assert(g.allele_count() == 2);
  assert(g.unphased_genotype_count() == 3);
  assert(g.phased_genotype_count() == 3);
  
  assert(g.get_sex_specific_allele().is_valid());
  assert(g.get_sex_specific_allele().name() == "~X");
  assert(g.get_sex_specific_allele().frequency() == 1.0);
  
  assert(!g.get_sex_specific_allele().is_null_y_allele());
  assert( g.get_sex_specific_allele().is_null_x_allele());

  assert(!g.is_x_linked());
  assert( g.is_y_linked());
  assert(!g.is_autosomal());

  test(g, "Y-Linked");
  print(g, "Y-Linked");

  g.set_model_type(AUTOSOMAL);
  
  assert(g.allele_count() == 2);
  assert(g.unphased_genotype_count() == 3);
  assert(g.phased_genotype_count() == 4);
  
  assert(!g.get_sex_specific_allele().is_valid());
  
  assert(!g.is_x_linked());
  assert(!g.is_y_linked());
  assert( g.is_autosomal());

  test(g, "Autosomal again");
  print(g, "Autosomal again");

}

void test_remap_sexed()
{
  // Create one with no type (ie, autosomal)
  genotype_model g;
  genotype_model remapped;
  
  g.add_allele("A", 0.5);
  g.add_allele("B", 0.5);
  g.add_allele("C", 0.5);
  
  g.mark_for_remap("A");
  g.mark_for_remap("B");
  
  g.remap(remapped);
  
  test(remapped, "Remapped Autosomal");
  print(remapped, "Remapped Autosomal");
  
  g.set_model_type(X_LINKED);

  g.remap(remapped);
  
  test(remapped, "Remapped X-Linked");
  print(remapped, "Remapped X-Linked");

  g.set_model_type(Y_LINKED);

  g.remap(remapped);
  
  test(remapped, "Remapped Y-Linked");
  print(remapped, "Remapped Y-Linked");
}

template <typename MGENO, typename FGENO>
void test_child_genotype_set
  (const MGENO& g1, const FGENO& g2,
   size_t exp_before_reduce,
   size_t exp_after_reduce)
{
    child_genotype_set cgs(g1, g2);

    assert(cgs.size() == exp_before_reduce);

    cgs.reduce();
    
    assert(cgs.size() == exp_after_reduce);

    child_genotype_set cgs2(g1, g2, true);

    assert(cgs.size() == exp_after_reduce);

    cgs.reduce();
    
    assert(cgs.size() == exp_after_reduce);
}

void test_child_genotype_set()
{
    genotype_model m;

    m.add_allele(a, 0.25);
    m.add_allele(b, 0.50);
    m.add_allele(c, 0.75);
    m.add_allele(d, 1.00);

    // Test with homozygous parents
    test_child_genotype_set(*m.phased_genotype_begin(), 
                            *m.phased_genotype_begin(),
                            4, 1);

    // Test with non-homozygous parents
    test_child_genotype_set(*(++m.phased_genotype_begin()), 
                            *(++m.phased_genotype_begin()),
                            4, 4);

    // Test with parents with different genotypes
    test_child_genotype_set(*(++m.phased_genotype_begin()), 
                            *m.phased_genotype_begin(),
                            4, 2);
                            
    // Test with homozygous parents
    test_child_genotype_set(*m.unphased_genotype_begin(), 
                            *m.unphased_genotype_begin(),
                            4, 1);

    // Test with non-homozygous parents
    test_child_genotype_set(*(++m.unphased_genotype_begin()), 
                            *(++m.unphased_genotype_begin()),
                            4, 4);

    // Test with parents with different genotypes
    test_child_genotype_set(*(++m.unphased_genotype_begin()), 
                            *m.unphased_genotype_begin(),
                            4, 2);
}
    
int
main()
{
  LSFInit();

  size_t i = (size_t) -1;

  cout << i << ' ' << (i / 2) << endl;

    e();
    f();
    g();
    h();
    test_sex_type();
    test_remap_sexed();
    
    _i(); // Runs penetrance tests
    
    test_child_genotype_set();

    test_file();

    return 0;
}
