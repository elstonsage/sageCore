#include <iomanip>

#ifdef _MSC_VER
    #include <app/SAGEconfig.h>
    #pragma hdrstop
#endif

#include "mlocus/imodel.h"

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

using namespace SAGE::MLOCUS;


//- Allele names.
//
extern std::string  a, b, c, d;

//- Phenotype names.
//
extern std::string  red, blue, green;
extern std::string  orange, white, black;
extern std::string  gray, yellow, purple, violet;

void test(const penetrance_model& m, const char* str)
{
  inheritance_model_test(cout, m, str);
}


void print(const penetrance_model& m, const char* str)
{
  inheritance_model_print(cout, m, str);
}


void aaa(genotype_model& gm)
{
    penetrance_model  m;

    test(m, "AAA - No Models");
    print(m, "AAA - No Models");

    m.add_allele(a, 0.1);    

    test(m, "AAA - Allele A");
    print(m, "AAA - Allele A");

    m.add_allele(b, 0.1);

    test(m, "AAA - Allele B");
    print(m, "AAA - Allele B");

    m.add_allele(a, 0.1);    

    test(m, "AAA - Allele A again");
    print(m, "AAA - Allele A again");

    m.add_allele(c, 0.1, true);    

    test(m, "AAA - Allele C");
    print(m, "AAA - Allele C");

    m.add_allele(d, 0.1, true, true);    

    test(m, "AAA - Allele D");
    print(m, "AAA - Allele D");

}
void bbb(genotype_model& gm)
{
    penetrance_model  m(gm), m2(gm, true);

    test(m, "BBB - complete Genotype Model, no phenotypes");
    print(m, "BBB - complete Genotype Model, no phenotypes");

    test(m2, "BBB - Both Models, no Penetrances");
    print(m2, "BBB - Both Models, no Penetrances");
}
void ccc(genotype_model& gm)
{
    penetrance_model  m(gm,true,true);

    test(m, "CCC - Both Models, Penetrances");
    print(m, "CCC - Both Models, Penetrances");
}

void ggg(genotype_model& gm)
{
    penetrance_model  m(gm, true, true);

    test(m, "GGG - codominant");
    print(m, "GGG - codominant");
   
    m.add_penetrance(1.0, "a/a", "a/b");

    test(m, "GGG - noncodominant");
    print(m, "GGG - noncodominant");

    m.add_penetrance(0.0, "a/a", "a/b");

    test(m, "GGG - after reset");
    print(m, "GGG - after reset");

    // Test strict flag

    size_t aa = m.get_phenotype_id("a/a");
    size_t ab = m.get_unphased_genotype("a/b").get_id();

    m.set_phenotype_strict(aa, false);

    test(m, "GGG - codominant aa not strict");
    print(m, "GGG - codominant aa not strict");

    m.add_unphased_penetrance(1.0, aa, ab);

    test(m, "GGG - codominant ignoring aa");
    print(m, "GGG - codominant ignoring aa");

    m.set_phenotype_strict(aa, true);

    test(m, "GGG - aa strict again");
    print(m, "GGG - aa strict again");

    m.remove_unphased_penetrance(aa, ab);

    test(m, "GGG - after reset");
    print(m, "GGG - after reset");

    m.add_unphased_penetrance(1.0, aa, ab);

    test(m, "GGG - noncodominant aa");
    print(m, "GGG - noncodominant aa");

    m.set_phenotype_strict(aa, false);

    test(m, "GGG - codominant without aa");
    print(m, "GGG - codominant without aa");

    m.remove_unphased_penetrance(aa, ab);

    test(m, "GGG - after reset");
    print(m, "GGG - after reset");

    m.set_phenotype_strict(aa, true);

    test(m, "GGG - codominant despite aa");
    print(m, "GGG - codominant despite aa");
}

void hhh(genotype_model& gm)
{
    penetrance_model  m(gm);

    test(m, "HHH - Just the Genotypes");
    print(m, "HHH - Just the Genotypes");

    std::string  g0(a + ';' + a);
    std::string  g1(a + ';' + b);
    std::string  g2(a + ';' + c);
    std::string  g3(a + ';' + d);
    std::string  g4(b + ';' + b);
    std::string  g5(b + ';' + c);
    std::string  g6(b + ';' + d);
    std::string  g7(c + ';' + c);
    std::string  g8(c + ';' + d);
    std::string  g9(d + ';' + d);

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

    test(m, "HHH - Add Phenotypes");
    print(m, "HHH - Add Phenotypes");

    m.add_penetrance(1.0, blue,   g0, Unphased, ';');
    m.add_penetrance(1.0, orange, g1, Unphased, ';');
    m.add_penetrance(1.0, green,  g2, Unphased, ';');
    m.add_penetrance(1.0, red,    g3, Unphased, ';');
    m.add_penetrance(1.0, purple, g4, Unphased, ';');
    m.add_penetrance(1.0, black,  g5, Unphased, ';');
    m.add_penetrance(1.0, gray,   g6, Unphased, ';');
    m.add_penetrance(1.0, yellow, g7, Unphased, ';');
    m.add_penetrance(1.0, white,  g8, Unphased, ';');
    m.add_penetrance(1.0, violet, g9, Unphased, ';');

    test(m, "HHH - codominant");
    print(m, "HHH - codominant");
   
    m.make_consistent();

    test(m, "HHH - codominant with phase");
    print(m, "HHH - codominant with phase");
}


void iii(genotype_model& gm)
{
    penetrance_model  m(gm);
    std::string  g0(a + ';' + a);
    std::string  g1(a + ';' + b);
    std::string  g2(a + ';' + c);
    std::string  g3(a + ';' + d);
    std::string  g4(b + ';' + b);
    std::string  g5(b + ';' + c);
    std::string  g6(b + ';' + d);
    std::string  g7(c + ';' + c);
    std::string  g8(c + ';' + d);
    std::string  g9(d + ';' + d);

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

    m.add_penetrance(1.0, red,    g0, Unphased, ';');
    m.add_penetrance(1.0, blue,   g1, Unphased, ';');
    m.add_penetrance(1.0, blue,   g2, Unphased, ';');
    m.add_penetrance(1.0, green,  g2, Unphased, ';');
    m.add_penetrance(1.0, orange, g3, Unphased, ';');
    m.add_penetrance(1.0, white,  g1, Unphased, ';');
    m.add_penetrance(1.0, white,  g4, Unphased, ';');
    m.add_penetrance(1.0, white,  g7, Unphased, ';');
    m.add_penetrance(1.0, black,  g5, Unphased, ';');
    m.add_penetrance(1.0, gray,   g6, Unphased, ';');
    m.add_penetrance(1.0, yellow, g7, Unphased, ';');
    m.add_penetrance(1.0, yellow, g1, Unphased, ';');
    m.add_penetrance(1.0, purple, g8, Unphased, ';');
    m.add_penetrance(1.0, violet, g0, Unphased, ';');
    m.add_penetrance(1.0, violet, g6, Unphased, ';');
    m.add_penetrance(1.0, violet, g9, Unphased, ';');

    test(m, "III - noncodominant");
    print(m, "III - noncodominant");

    penetrance_model  m2, m1(m);

    test(m1, "III1");
    print(m1, "III1");

    m2 = m;
    test(m2, "III2");
    print(m2, "III2");

    m2.make_consistent();

    test(m2, "III3");
    print(m2, "III3");
}

void jjj(genotype_model& gm)
{
    gm.set_missing_allele_name("*");

    penetrance_model  m(gm);

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

    m.add_penetrance(1.0, black,  "*/*");
    m.add_penetrance(1.0, blue,   "*<*");
    m.add_penetrance(1.0, gray,   "*/a");
    m.add_penetrance(1.0, white,  "b/*");
    m.add_penetrance(1.0, yellow, "b<*");
    m.add_penetrance(1.0, violet, "*<a");

    test(m, "JJJ - using missing allele");
    print(m, "JJJ - using missing allele");

    m.make_consistent();

    test(m, "JJJ - using missing allele");
    print(m, "JJJ - using missing allele");
}


void kkk(genotype_model& gm)
{
    gm.set_missing_allele_name("*");

    penetrance_model  m(gm);

    m.add_phenotype(blue);  
    m.add_phenotype(orange);
    m.add_phenotype(green); 
    m.add_phenotype(red);   
    m.add_phenotype(purple);
    m.add_phenotype(gray);  
    m.add_phenotype(yellow);
    m.add_phenotype(white); 
    m.add_phenotype(violet);
    m.add_phenotype(black); 

    m.add_penetrance(1.0, black,  "*/*");
    m.add_penetrance(1.0, blue,   "*<*");
    m.add_penetrance(1.0, gray,   "*/a");
    m.add_penetrance(1.0, violet, "*<a");
    m.add_penetrance(1.0, white,  "b/*");
    m.add_penetrance(1.0, yellow, "b<*");

    m.add_penetrance(0.0, black,  "a/b");
    m.add_penetrance(0.5, black,  "a/c");
    m.add_penetrance(0.0, blue,   "b<c");
    m.add_penetrance(0.0, gray,   "c/c");

    test(m, "KKK - using missing allele plus override");
    print(m, "KKK - using missing allele plus override");

    size_t id1 = m.get_phenotype_id(gray);
    size_t id2 = m.get_phenotype_id(orange);

    m.copy_penetrance(id1, id2);

    test(m, "KKK - copy gray to orange");
    print(m, "KKK - copy gray to orange");

    id2 = m.get_phenotype_id(white);

    m.copy_penetrance(id1, id2);

    test(m, "KKK - copy gray to white");
    print(m, "KKK - copy gray to white");

    id2 = m.get_phenotype_id(blue);

    m.copy_penetrance(id1, id2, true);

    test(m, "KKK - gray override blue");
    print(m, "KKK - gray override blue");
}

void lll(genotype_model& gm)
{
    penetrance_model  m(gm, true, true);

    test(m, "LLL - codominant");
    print(m, "LLL - codominant");
   
    m.mark_for_remap(c);
    m.mark_for_remap(d);
    cout << "REMAPPING: " << m.phenotype_count() << endl;

    m.remap();

    cout << "done." << m.phenotype_count() << endl;

    test(m, "LLL - after remap");
    print(m, "LLL - after remap");
  
}

void mmm(genotype_model& gm)
{
    penetrance_model  m(gm, true, true);

    test(m, "MMM - before copy");
    print(m, "MMM - before copy");
   
    penetrance_model n(m);

    cout << "Uniquefying" << endl;

    n.mark_for_remap(c);

    test(n, "MMM - after copy");
    print(n, "MMM - after copy");
  
}

void test_sexed_penetrance()
{
    genotype_model gm;
    
    gm.set_model_type(X_LINKED);
  
    penetrance_model  m(gm);
    
    m.add_allele(a, 0.1);    
    m.add_allele(b, 0.1);
    m.add_allele(a, 0.1);    
    m.add_allele(c, 0.1, true);    
    m.add_allele(d, 0.1, true, true);    

    test(m, "X-Linked Penetrance Test");
    print(m, "X-Linked Penetrance Test");
    
    // Expand and make sure the marker still has the right penetrances.
    m.add_allele("e", 0.1, true, true);    

    test(m, "X-Linked Penetrance Expansion Test");
    print(m, "X-Linked Penetrance Expansion Test");

    // Remap a few alleles and test again
    
    m.mark_for_remap(a);
    m.mark_for_remap(b);
    m.mark_for_remap(c);
    m.mark_for_remap("~Y"); // Note, cannot remap sex allele
    m.remap();
    
    test(m, "X-Linked Penetrance Remap Test");
    print(m, "X-Linked Penetrance Remap Test");

    gm.set_model_type(Y_LINKED);
  
    penetrance_model m2(gm);
  
    m2.add_allele(a, 0.1);    
    m2.add_allele(b, 0.1);
    m2.add_allele(a, 0.1);    
    m2.add_allele(c, 0.1, true);    
    m2.add_allele(d, 0.1, true, true); 
    
    test(m2, "Y-Linked Penetrance Test");
    print(m2, "Y-Linked Penetrance Test");
    
    // Expand and make sure the marker still has the right penetrances.
    m2.add_allele("e", 0.1, true, true);    

    test(m2, "Y-Linked Penetrance Expansion Test");
    print(m2, "Y-Linked Penetrance Expansion Test");
    // Remap a few alleles and test again
    
    m2.mark_for_remap(a);
    m2.mark_for_remap(b);
    m2.mark_for_remap(c);
    m2.mark_for_remap("~X"); // Note, cannot remap sex allele
    m2.remap();
    
    test(m2, "Y-Linked Penetrance Remap Test");
    print(m2, "Y-Linked Penetrance Remap Test");
}

void 
test_sex_type_change()
{
    genotype_model gm;
    
    penetrance_model  m(gm);
    
    m.add_allele(a, 0.1);    
    m.add_allele(b, 0.1);
    m.add_allele(a, 0.1);    
    m.add_allele(c, 0.1, true);    
    m.add_allele(d, 0.1, true, true);    

    test(m, "Before Sex Change");
    print(m,"Before Sex Change");

    m.set_model_type(X_LINKED);
    
    test(m, "X-Linked Sex Change");
    print(m,"X-Linked Sex Change");

    m.set_model_type(AUTOSOMAL);
    
    test(m, "Autosomal Sex Change");
    print(m,"Autosomal Sex Change");
}

void
test_dynamic_allele_add
    (penetrance_model& m, 
     const std::string& genotype,
     size_t expected_allele_count)
{
  m.add_genotype_dynamically(genotype);
  
  assert(m.allele_count() == expected_allele_count);
}

void test_dynamic_alleles_autosomal()
{
  penetrance_model m;
  
  m.gmodel().set_dynamic_alleles(false);
  
  // Test with non-dynamic marker.  Nothing here should add alleles
  
  test_dynamic_allele_add(m, "FOO",        0); // Note invalid genotype
  test_dynamic_allele_add(m, "A/A",        0); // One new allele
  test_dynamic_allele_add(m, "A/B",        0); // New allele second
  test_dynamic_allele_add(m, "C/B",        0); // New allele first
  test_dynamic_allele_add(m, "D/E",        0); // Two new alleles
  test_dynamic_allele_add(m, "*missing/F", 0); // One new, first missing
  test_dynamic_allele_add(m, "G/*missing", 0); // One new, second missing

  // Now we try with an actual dynamic marker

  m.gmodel().set_dynamic_alleles(true);
  
  test_dynamic_allele_add(m, "FOO",        0); // Note invalid genotype
  test_dynamic_allele_add(m, "A/A",        1); // One new allele
  test_dynamic_allele_add(m, "A/B",        2); // New allele second
  test_dynamic_allele_add(m, "C/B",        3); // New allele first
  test_dynamic_allele_add(m, "D/E",        5); // Two new alleles
  test_dynamic_allele_add(m, "*missing/F", 6); // One new, first missing
  test_dynamic_allele_add(m, "G/*missing", 7); // One new, second missing
}

void test_dynamic_alleles_x_linked()
{
  penetrance_model m;
  m.set_model_type(X_LINKED);
  
  m.gmodel().set_dynamic_alleles(false);
  
  // Test with non-dynamic marker.  Nothing here should add alleles
  
  test_dynamic_allele_add(m, "FOO",        0); // Note invalid genotype
  test_dynamic_allele_add(m, "A/A",        0); // One new allele
  test_dynamic_allele_add(m, "A/B",        0); // New allele second
  test_dynamic_allele_add(m, "C/B",        0); // New allele first
  test_dynamic_allele_add(m, "D/E",        0); // Two new alleles
  test_dynamic_allele_add(m, "*missing/F", 0); // One new, first missing
  test_dynamic_allele_add(m, "G/*missing", 0); // One new, second missing
  test_dynamic_allele_add(m, "H/~Y",       0); // One new, second ~Y
  test_dynamic_allele_add(m, "~Y/I",       0); // One new, first ~Y
  test_dynamic_allele_add(m, "~Y/~Y",      0); // Double ~Y

  // Now we try with an actual dynamic marker

  m.gmodel().set_dynamic_alleles(true);
  
  test_dynamic_allele_add(m, "FOO",        0); // Note invalid genotype
  test_dynamic_allele_add(m, "A/A",        1); // One new allele
  test_dynamic_allele_add(m, "A/B",        2); // New allele second
  test_dynamic_allele_add(m, "C/B",        3); // New allele first
  test_dynamic_allele_add(m, "D/E",        5); // Two new alleles
  test_dynamic_allele_add(m, "*missing/F", 6); // One new, first missing
  test_dynamic_allele_add(m, "G/*missing", 7); // One new, second missing
  test_dynamic_allele_add(m, "H/~Y",       8); // One new, second ~Y
  test_dynamic_allele_add(m, "~Y/I",       9); // One new, first ~Y
  test_dynamic_allele_add(m, "~Y/~Y",      9); // Double ~Y
}
void test_dynamic_alleles_y_linked()
{
  penetrance_model m;
  m.set_model_type(Y_LINKED);
  
  m.gmodel().set_dynamic_alleles(false);
  
  // Test with non-dynamic marker.  Nothing here should add alleles
  
  test_dynamic_allele_add(m, "FOO",        0); // Note invalid genotype
  test_dynamic_allele_add(m, "A/A",        0); // One new allele
  test_dynamic_allele_add(m, "A/B",        0); // New allele second
  test_dynamic_allele_add(m, "C/B",        0); // New allele first
  test_dynamic_allele_add(m, "D/E",        0); // Two new alleles
  test_dynamic_allele_add(m, "*missing/F", 0); // One new, first missing
  test_dynamic_allele_add(m, "G/*missing", 0); // One new, second missing
  test_dynamic_allele_add(m, "H/~X",       0); // One new, second ~X
  test_dynamic_allele_add(m, "~X/I",       0); // One new, first ~X
  test_dynamic_allele_add(m, "~X/~X",      0); // Double ~X

  // Now we try with an actual dynamic marker

  m.gmodel().set_dynamic_alleles(true);
  
  test_dynamic_allele_add(m, "FOO",        0); // Note invalid genotype
  test_dynamic_allele_add(m, "A/A",        1); // One new allele
  test_dynamic_allele_add(m, "A/B",        2); // New allele second
  test_dynamic_allele_add(m, "C/B",        3); // New allele first
  test_dynamic_allele_add(m, "D/E",        5); // Two new alleles
  test_dynamic_allele_add(m, "*missing/F", 6); // One new, first missing
  test_dynamic_allele_add(m, "G/*missing", 7); // One new, second missing
  test_dynamic_allele_add(m, "H/~X",       8); // One new, second ~X
  test_dynamic_allele_add(m, "~X/I",       9); // One new, first ~X
  test_dynamic_allele_add(m, "~X/~X",      9); // Double ~X
}

void test_dynamic_alleles()
{
  test_dynamic_alleles_autosomal();
  test_dynamic_alleles_x_linked();
  test_dynamic_alleles_y_linked();
}


