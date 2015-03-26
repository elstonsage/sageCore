#include "lvec/inheritance_vector.h"
#include "lvec/iv_generator.h"
#include "lvec/codom_ivgen.h"
#include "lvec/dgraph.h"
#include "lvec/fixed_bit_calculator.h"

using namespace SAGE;

void test_meiosis_map(ostream& o, const meiosis_map& mm)
{
  o << "Meiosis_Map test on Pedigree: " << mm.get_pedigree()->name() << endl;

  o << mm.get_pedigree()->member_count() << endl;
  o << mm.founder_count() << endl;
  o << mm.nonfounder_count() << endl;
  o << mm.founder_mask() << endl;
  o << mm.nonfounder_mask() << endl;

  for(int count = 0; count < (int)mm.get_pedigree()->member_count(); count++)
  {
    o << setw(3) << left << mm.member(count)->name().c_str() << ' ' 
      << right;

    o << setw(3) << (signed) mm.mother_meiosis(count) << ' '
      << setw(3) << (signed) mm.father_meiosis(count) << ' ';

    o << setw(3) << (signed) mm.mother_index(count) << ' '
      << setw(3) << (signed) mm.father_index(count) << ' ';

    o << setw(8) << (signed) mm.mother_mask(count) << ' '
      << setw(8) << (signed) mm.father_mask(count);

    o << endl;
  }

  for(int count = 0; count < (int)mm.founder_count(); ++count)
    o << (signed) mm.mask(count) << endl;
}

void test_inheritance_vector(ostream& o, inheritance_vector& iv)
{
  const meiosis_map& mm = iv.get_meiosis_map();

  cout << "IV test on pedigree: " << mm.get_pedigree()->name() << endl;

  int max = iv.num_equivalence_classes();

  if(max > 100) max = 10;

  for(inheritance_vector::storage_type s = 0; s < (size_t)max; ++s)
  {
    o << "Forwards:  ";
    iv.set_equivalence_class(s);
    for( ; !iv.isEnd(); ++iv)
    {
      o << iv.get_equivalence_class() << ": "
        << iv.get_founders() << '-' << iv.get_nonfounders() << "\t";
    }

    o << endl << "Backwards: ";

    for(--iv ; !iv.isBegin(); --iv)
    {
      o << iv.get_equivalence_class() << ": "
        << iv.get_founders() << '-' << iv.get_nonfounders() << "\t";
    }

    o << endl;
  }

  iv.set_iteration_method(inheritance_vector::BW);
 
  o << endl << endl << "Bitwise" << endl << endl; 

  max = 1 << (mm.nonfounder_count() + mm.founder_count());
  if(max > 1000)
    max = 1 << (mm.nonfounder_count() + 1);
  
  for(int x = 0; x < max && !iv.isEnd(); ++iv, ++x)
  {
      o << iv.get_equivalence_class() << ": "
        << iv.get_founders() << '-' << iv.get_nonfounders() << endl;
  }

  o << "Backwards" << endl;

  for(--iv ; !iv.isBegin(); --iv)
  {
      o << iv.get_equivalence_class() << ": "
        << iv.get_founders() << '-' << iv.get_nonfounders() << endl;
  }

  o << endl << "Setting" << endl << endl;

  iv.meiosis(0) = true;

  o << iv.get_founders() << ' ' << iv.get_nonfounders() << ' '
    << iv.get_equivalence_class() << endl;

  iv.meiosis(meiosis_map::meiosis_bits) = true;

  o << iv.get_founders() << ' ' << iv.get_nonfounders() << ' '
    << iv.get_equivalence_class() << endl;

  iv.meiosis(1) = true;

  o << iv.get_founders() << ' ' << iv.get_nonfounders() << ' '
    << iv.get_equivalence_class() << endl;

  iv.meiosis(0) = false;

  o << iv.get_founders() << ' ' << iv.get_nonfounders() << ' '
    << iv.get_equivalence_class() << endl;

  iv.meiosis(0) = false;

  o << iv.get_founders() << ' ' << iv.get_nonfounders() << ' '
    << iv.get_equivalence_class() << endl;
}

class TestIVA : public iv_acceptor
{
public:

  virtual void accept(equivalence_class e, equivalence_class dontcare, double p)
  { for(int i = 1; i <= (int)e; i*=2)
      cout << (int)(bool) (e & i);
    cout << " = " << p << endl; }
};

void test_iv_generator(ostream& o, const meiosis_map& mm, const MLOCUS::inheritance_model& pm)
{
  o << "IVGenerator test on pedigree: " << mm.get_pedigree()->name() 
    << endl;

  boost::shared_ptr<iv_acceptor> iva(new TestIVA());

  iv_generator ivg(iva);

  if(!ivg.build(mm, pm))
  {
    o << "Pedigree " << mm.get_pedigree()->name()  << " is bad!" << endl;
  }
}

void test_codom_iv_generator(ostream& o, const meiosis_map& mm, const MLOCUS::inheritance_model& pm)
{
  o << "Codom IVGenerator test on pedigree: " << mm.get_pedigree()->name() 
    << endl;

  boost::shared_ptr<iv_acceptor> iva(new TestIVA());

  codominant_iv_generator ivg(iva);

  if(!ivg.build(mm, pm))
  {
    o << "Pedigree " << mm.get_pedigree()->name()  << " is bad!" << endl;
  }
}

void out(meiosis_map::storage_type u)
{
  for(meiosis_map::storage_type i = 1; i <= u; i <<= 1)
    cout << (bool) (u & i);
  cout << '\t';
}

void test_descent_graph(ostream& o, inheritance_vector& iv)
{
  const meiosis_map& mm = iv.get_meiosis_map();

  cout << endl << "Descent Graph Iteration Testing on peidgree: "
       << mm.get_pedigree()->name()
       << endl << endl;

  descent_graph dg(iv);

  int max = iv.num_equivalence_classes();

  if(max > 100) max = 10;

  for(inheritance_vector::storage_type s = 0; s < (size_t)max; ++s)
  {
    iv.set_equivalence_class(s);
    for( ; !iv.isEnd(); ++iv)
    {
      dg.move_to(iv);
      out(iv.get_founders());
      out(iv.get_nonfounders());
      out(iv.get_equivalence_class());
      for(int j = 0; j < (int)mm.get_pedigree()->member_count(); ++j)
        cout << dg.alleles(j).first << 'x' << dg.alleles(j).second << ' ';
      cout << endl;
    }

  }
}

void test_fixed_bit_calculator(ostream& o, const FPED::Subpedigree& sped, const MLOCUS::inheritance_model& pm)
{
  o << "fixed_bit_calculator test on pedigree: " << sped.name() 
    << endl;

  fixed_bit_calculator fbc(sped, pm);

  const fixed_bit_container& fb = fbc.get_fixed_bit_container();

  o << "marker " << pm.name() << endl;

  fb.dump(o);
}

