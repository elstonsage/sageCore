#include "segreg/RegPenetranceCalculator.h"

namespace SAGE {
namespace SEGREG {

//======================================================================
//
//  Penetrance functions start here.......
//
//======================================================================


double
RegPenetranceCalculator::get_penetrance
    (penetrance_info          i,
     const PenetranceContext& context) const
{
  if(!my_member_calc.is_member_valid(*i.member))
    return 1.0;

  double exp_sd     = my_member_calc.get_expected_sd       (*i.member,i.genotype);
  double exp_var    = my_member_calc.get_expected_variance (*i.member,i.genotype);
  double std_indiv  = my_member_calc.get_ascertained_standardization   (*i.member,i.genotype);

  // Calculate z

  double z = std_indiv;
  double w = 1.0;

  // Add family related adjustments, if needed
  const PenetranceContextStatistics& stats = context.my_stats;

  if(context.includes_family())
  {
    // Perform parental adjustments
    z += stats.get_parental_mean_adjustment();
    w += my_corr_adjs.get_parental_var_adjustment(stats.get_mother_informativity(),
                                                  stats.get_father_informativity());

    // If model is class D, we add sibling effects

    if(my_model_class == model_D && context.includes_siblings())
    {
      size_t prior_sib_count = get_prior_sib_count(*i.member);

      z += stats.get_sibling_mean_adjustment(prior_sib_count, context);
      w += my_corr_adjs.get_sibling_var_adjustment(stats.get_mother_informativity(),
                                                   stats.get_father_informativity(),
                                                   prior_sib_count);

    }
  }

  // Add spousal adjustments if necessary

  if(context.incorporate_spousal_adjustments(*i.member))
  {
    z += stats.get_spousal_mean_adjustment(); 
    w += my_corr_adjs.get_spousal_var_adjustment();
  }
  
  z *= exp_sd;
  w *= exp_var;

  // Calculate the penetrance based on z and w

  if(SAGE::isnan(z))
  {
/*    cout << std_indiv << endl;
    if(context.includes_family())
    {
      cout << "Fam Adj: " << context.get_parental_mean_adjustment() << endl;

      if(my_model_class == model_D && context.includes_siblings())
      {
        size_t prior_sib_count = get_prior_sib_count(*i.member);

        cout << "Sib Adj: " << context.get_sibling_mean_adjustment(prior_sib_count)
             << endl;
      }
    }

    if(context.incorporate_spousal_adjustments(*i.member))
    {
      cout << "Sps Adj: " << context.get_spousal_mean_adjustment() << endl;
    }

    cout << "EXP STD: " << exp_sd << endl;
*/
    exit(1);
  }

//  cout << z << endl;

  return calculate_continuous_penetrance(z,w, my_member_calc.get_member_class(*i.member));
}


double
RegPenetranceCalculator::get_penetrance(penetrance_info         i) const
{
  if(!my_member_calc.is_member_valid(*i.member))
    return 1.0;

  double std_indiv  = my_member_calc.get_ascertained_standardization (*i.member,i.genotype);
  double exp_sd     = my_member_calc.get_expected_sd                 (*i.member,i.genotype);
  double exp_var    = my_member_calc.get_expected_variance           (*i.member,i.genotype);
  double z          = std_indiv * exp_sd;
  double w          = exp_var;

  if(SAGE::isnan(z))
  {
/*    cout << i.member << endl;
    cout << i.member->name() << endl;
    cout << i.member->pedigree() << endl;
    cout << i.member->pedigree()->name() << endl;
    cout << "STD IND: " << std_indiv << endl;

    cout << "EXP STD: " << exp_sd << endl;

    cout << i.member->family() << endl;

    cout << i.member->family()->name() << endl;
*/
    exit(1);
  }


  return calculate_continuous_penetrance(z,w, my_member_calc.get_member_class(*i.member));
}

}
}
