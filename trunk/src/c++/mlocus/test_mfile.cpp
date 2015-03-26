#include "mlocus/mfile.h"

void tf1()
{
  cout << "============================" << endl;

  SAGE::MLOCUS::InheritanceModelFile i;

  SAGE::MLOCUS::inheritance_model_map m;

  i.input(m, "locus");
  i.output(m, cout);
}

void tf2()
{
  cout << "============================" << endl;

  SAGE::MLOCUS::InheritanceModelFile i;

  SAGE::MLOCUS::inheritance_model_map m;

  i.set_genotype_verbose_output(10);
  i.set_marker_verbose_output(5);

  i.input(m, "locus");
  i.output(m, cout);
}

void tf3()
{
  cout << "============================" << endl;

  SAGE::MLOCUS::InheritanceModelFile i;

  SAGE::MLOCUS::inheritance_model_map m;

  i.set_genotype_verbose_output(10);
  i.set_marker_verbose_output(5);

  i.input(m, "locus.typ");

  SAGE::MLOCUS::inheritance_model_map::index_iterator j = m.index_begin();

  for( ; j != m.index_end(); ++j)
    inheritance_model_print(cout, j->second, "INPUT TEST");
}
void tf4()
{
  cout << "============================" << endl;

  SAGE::MLOCUS::InheritanceModelFile i;

  SAGE::MLOCUS::inheritance_model_map m;

  i.set_genotype_verbose_output(10);
  i.set_marker_verbose_output(5);

  i.input(m, "inval.typ");
}
void tf5()
{
  cout << "============================" << endl;

  SAGE::MLOCUS::InheritanceModelFile i;

  SAGE::MLOCUS::inheritance_model_map m;

  i.set_genotype_verbose_output(10);
  i.set_marker_verbose_output(5);

  i.input(m, "locus.pen");
}

void test_file()
{
  tf1();
  tf2();
  tf3();
  tf5();
  tf4();
}
