#ifndef MAXFUN_JOINT_TEST_H
#define MAXFUN_JOINT_TEST_H

#include "maxfunapi/Results.h"

namespace SAGE   {
namespace MAXFUN {

class JointTest
{
  public:
  
    JointTest(const Results & results1, const Results & results2)
    {
      my_H0_name  = "H0 " + results1.getFunctionName() + " " + results1.getSequenceName();
      my_H1_name  = "H1 " + results2.getFunctionName() + " " + results2.getSequenceName();

      my_H0_valid = results1.getConverged() && !results1.getWasSkipped();
      my_H1_valid = results2.getConverged() && !results2.getWasSkipped();

      my_H0_val   = results1.getFinalFunctionValue();
      my_H1_val   = results2.getFinalFunctionValue();

      my_degrees_of_freedom = abs(results1.getNumOfIndependentParams() - results2.getNumOfIndependentParams());
      my_comp_val           = 2 * (abs(my_H0_val - my_H1_val)),
      my_p_value            = chdtrc(my_degrees_of_freedom, my_comp_val);

      if(SAGE::isnan(my_p_value)) 
        my_p_value = 1.0;
    }
  
    JointTest(const JointTest & other) :
      my_H0_name            (other.my_H0_name),
      my_H1_name            (other.my_H1_name),
      my_H0_valid           (other.my_H0_valid),
      my_H1_valid           (other.my_H1_valid),
      my_H0_val             (other.my_H0_val),
      my_H1_val             (other.my_H1_val),
      my_degrees_of_freedom (other.my_degrees_of_freedom),
      my_comp_val           (other.my_comp_val),
      my_p_value            (other.my_p_value)
    {
    }
    
    JointTest & operator= (const JointTest & other)
    {
      if(this != &other)
      {
        my_H0_name            = my_H0_name;
        my_H1_name            = my_H1_name;
        my_H0_valid           = my_H0_valid;
        my_H1_valid           = my_H1_valid;
        my_H0_val             = my_H0_val;
        my_H1_val             = my_H1_val;
        my_degrees_of_freedom = my_degrees_of_freedom;
        my_comp_val           = my_comp_val;
        my_p_value            = my_p_value;
      }
      
      return *this;
    }
  
    OUTPUT::Table summarizeAsTable() const
    {
      OUTPUT::Table t("Joint test");

      OUTPUT::TableRow row0 = OUTPUT::TableRow() << my_H0_name << my_H0_val;

      if(!my_H0_valid)
        row0 << OUTPUT::String("(possible non-convergence)");

      t << row0;

      OUTPUT::TableRow row1 = OUTPUT::TableRow() << my_H1_name << my_H1_val;
  
      if(!my_H1_valid)
        row1 << OUTPUT::String("(possible non-convergence)");

      t << row1;

      t.insertBlankRow();

      t << (OUTPUT::TableRow() << OUTPUT::String("2 * |H0 - H1|")      <<      my_comp_val)
        << (OUTPUT::TableRow() << OUTPUT::String("Degrees of freedom") << (int)my_degrees_of_freedom)
        << (OUTPUT::TableRow() << OUTPUT::String("P-value")            <<      my_p_value);

      return t;
    }

  private:

    std::string my_H0_name;
    std::string my_H1_name;

    bool my_H0_valid;
    bool my_H1_valid;
    
    double my_H0_val;
    double my_H1_val;

    size_t my_degrees_of_freedom;
    double my_comp_val;
    double my_p_value;

};

} // End namespace MAXFUN
} // End namespace SAGE

#endif
