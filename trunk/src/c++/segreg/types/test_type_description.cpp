#include "segreg/types/TypeDescription.h"

using namespace SAGE;
using namespace SEGREG;

int main()
{
  TypeDescription t;
  
  t.set_name("Name");
  assert(t.get_name() == "Name");
  
  t.set_description("Desc");
  assert(t.get_description() == "Desc");
  
  assert(t.get_state_count() == 0);
  
  t.reserve_state_count(10);
  
  assert(t.get_state_count() == 0);
  
  t.add_state("State1");
  t.add_state("State2");
  t.add_state("State3");
  t.add_state("State4");

  assert(t.get_state_count() == 4);
  assert(t.get_state(0).get_index() == 0 && t.get_state(0).get_name() == "State1");
  assert(t.get_state(1).get_index() == 1 && t.get_state(1).get_name() == "State2");
  assert(t.get_state(2).get_index() == 2 && t.get_state(2).get_name() == "State3");
  assert(t.get_state(3).get_index() == 3 && t.get_state(3).get_name() == "State4");

  TypeDescription t2(t);

  assert(t2.get_name()        == "Name");
  assert(t2.get_description() == "Desc");
  assert(t2.get_state_count() == 4);
  assert(t2.get_state(0).get_index() == 0 && t2.get_state(0).get_name() == "State1");
  assert(t2.get_state(1).get_index() == 1 && t2.get_state(1).get_name() == "State2");
  assert(t2.get_state(2).get_index() == 2 && t2.get_state(2).get_name() == "State3");
  assert(t2.get_state(3).get_index() == 3 && t2.get_state(3).get_name() == "State4");

  t2.add_state("State5");
  assert(t2.get_name()        == "Name");
  assert(t2.get_description() == "Desc");
  assert(t2.get_state_count() == 5);
  assert(t2.get_state(0).get_index() == 0 && t2.get_state(0).get_name() == "State1");
  assert(t2.get_state(1).get_index() == 1 && t2.get_state(1).get_name() == "State2");
  assert(t2.get_state(2).get_index() == 2 && t2.get_state(2).get_name() == "State3");
  assert(t2.get_state(3).get_index() == 3 && t2.get_state(3).get_name() == "State4");
  assert(t2.get_state(4).get_index() == 4 && t2.get_state(4).get_name() == "State5");

  assert(t.get_name()        == "Name");
  assert(t.get_description() == "Desc");
  assert(t.get_state_count() == 4);
  assert(t.get_state(0).get_index() == 0 && t.get_state(0).get_name() == "State1");
  assert(t.get_state(1).get_index() == 1 && t.get_state(1).get_name() == "State2");
  assert(t.get_state(2).get_index() == 2 && t.get_state(2).get_name() == "State3");
  assert(t.get_state(3).get_index() == 3 && t.get_state(3).get_name() == "State4");

  t = t2;

  assert(t.get_name()        == "Name");
  assert(t.get_description() == "Desc");
  assert(t.get_state_count() == 5);
  assert(t.get_state(0).get_index() == 0 && t.get_state(0).get_name() == "State1");
  assert(t.get_state(1).get_index() == 1 && t.get_state(1).get_name() == "State2");
  assert(t.get_state(2).get_index() == 2 && t.get_state(2).get_name() == "State3");
  assert(t.get_state(3).get_index() == 3 && t.get_state(3).get_name() == "State4");
  assert(t.get_state(4).get_index() == 4 && t.get_state(4).get_name() == "State5");

  OUTPUT::Section s  = convert_to_output(t);
  OUTPUT::Section s2 = convert_to_output(t2);

  std::cout << s << s2;
  
  std::cout << t << t2;

  return 0;
}

