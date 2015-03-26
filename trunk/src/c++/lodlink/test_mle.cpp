//============================================================================
// File:      test_mle.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   9/4/2   Created.  -djb
//            11/7/6  Removed synchronization check upon conversion to
//                    use of maxfunapi (no longer relevant).  -djb                                                   
//                                                                          
// Notes:     Tests the lodlink mle_sub_model class.
//    
//                                                                      
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include <fstream>
#include "lodlink/mle_sub_model.h"

using namespace std;
using namespace SAGE;
using namespace LODLINK;

void  check_set_functions(ofstream& out);
// void  check_synchronization(ofstream& out);

int main(int argc, char* argv[])
{

  ofstream out;
  out.open("testmle.out");
  
  if(! out)
  {
    cout << "Cannot open output file: testmle.out.  Exiting ..." << endl;
    exit(EXIT_FAILURE);
  }

  check_set_functions(out);
//  check_synchronization(out);

  exit(EXIT_SUCCESS);
}

void
check_set_functions(ofstream& out)
{
  out << "----------  Default instantiation  ----------\n\n";  
  mle_sub_model  inst;
  out << inst;
  
  out << "----------  set(sex-specific, no alpha)  ----------\n\n";
  inst.set(true, false);
  out << inst;
  
  out << "----------  set_strict_limits()  ----------\n\n";
  inst.set_strict_limits();
  out << inst;
  
  out << "----------  set_relaxed_limits()  ----------\n\n";
  inst.set_relaxed_limits();
  out << inst;
  
  out << "----------  set(not sex-specific, alpha)  ----------\n\n";
  inst.set(false, true);
  out << inst;
  
  out << "----------  fix_alpha()  ----------\n\n";
  inst.fix_alpha();
  out << inst;
  
  out << "----------  unfix_alpha()  ----------\n\n";
  inst.unfix_alpha();
  out << inst;
  
  out << "----------  set(sex-specific, alpha)  --------\n\n";
  inst.set(true, true);
  out << inst;
  
  out << "----------  constrain_thetas()  ----------\n\n";
  inst.constrain_thetas();
  out << inst;
  
  out << "----------  unconstrain_thetas()  ----------\n\n";
  inst.unconstrain_thetas();
  out << inst;
  
  out << "----------  set_male_theta(.3)  -----------\n\n";
  inst.set_male_theta(.3);
  out << inst;
  out << "male theta " << inst.theta(LODLINK::male) << endl;
  out << "female_theta " << inst.theta(LODLINK::female) << endl;
  
  out << "----------  set_female_theta(.3)  ----------\n\n";
  inst.set_female_theta(.3);
  out << inst;
  out << "male theta " << inst.theta(LODLINK::male) << endl;
  out << "female_theta " << inst.theta(LODLINK::female) << endl;
  
  out << "----------  reset()  ----------\n\n";
  inst.reset();
  out << inst;
  
  out << "----------  set_average_theta(.3)  ----------\n\n";
  inst.set_average_theta(.3);
  out << inst;  
  out << "male theta " << inst.theta(LODLINK::male) << endl;
  out << "female_theta " << inst.theta(LODLINK::female) << endl;
}

/* NOT APPLICABLE W. MAXFUNAPI
class synch_checker : private mle_sub_model
{
  public:
    using mle_sub_model::constrain_thetas;
  
    synch_checker(ofstream& out, vector<double>& mvs, bool ss, bool ua)
          : mle_sub_model(ss, ua), my_out(out), my_mvs(mvs)
    {}
    
    void  check()
    {
      my_out << "----------  Before synchronization  ----------\n\n";
      print_parameters();
      
      synchronize(my_mvs.begin());
      
      my_out << "----------  After synchronization  ----------\n" << endl;
      print_parameters();
    }
    
  private:
    void  print_parameters();
  
    ofstream&  my_out; 
    vector<double>&  my_mvs;
};

void
check_synchronization(ofstream& out)
{
  vector<double>  maxfun_values;
  maxfun_values.push_back(.2);
  maxfun_values.push_back(.4);
  maxfun_values.push_back(.7);
 
  synch_checker  checker(out, maxfun_values, true, true);
  checker.constrain_thetas();
  checker.check();
}



void
synch_checker::print_parameters()
{
  my_out << "Maxfun values:               " << my_mvs[0] << ", "
                                            << my_mvs[1] << ", "
                                            << my_mvs[2] << endl;
  my_out << "Sub-model parameter values:  " << my_parameters[MALE].value << ", "
                                            << my_parameters[FEMALE].value << ", "
                                            << my_parameters[ALPHA_TWO].value << endl;
  my_out << "Sub-model values:            " << male_theta() << ", "
                                            << female_theta() << ", "
                                            << alpha() << endl;
}
*/


