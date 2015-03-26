//============================================================================
// File:      likelihood.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   12/2/2 created        -djb
//                                                                          
// Notes:     inline implementation of likelihood calculators.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

/*
inline void
check_parameters(const parameter_vector& theta, const mle_sub_model& mle)
{
  const MAXFUN::ParameterMgr* parameter_manager = mle.parameter_mgr();
  assert(parameter_manager);
  
  size_t  parameter_count = static_cast<size_t>(parameter_manager->getParamCount());
  assert(parameter_count == theta.size());
  
  for(size_t i = 0; i < parameter_count; ++i)
  {
    assert(theta[i] == (*parameter_manager)(i));
  }
}
*/

//============================================================================
// IMPLEMENTATION:  subped_calculator
//============================================================================
//
inline
subped_calculator::subped_calculator(peeler& p)
      : my_peeler(p), my_unlinked_likelihood(QNAN), 
        unlinked_likelihood_cached(false)
{
  nfe = 0;
}

inline double
subped_calculator::evaluate(parameter_vector& theta)
{
  check_parameters(theta, my_peeler.tcalc().mle());

  /*
  const MAXFUN::ParameterMgr* parameter_manager = my_peeler.tcalc().mle().parameter_mgr();
  assert(parameter_manager);
  
  size_t  parameter_count = static_cast<size_t>(parameter_manager->getParamCount());
  assert(parameter_count == theta.size());
  
  for(size_t i = 0; i < parameter_count; ++i)
  {
    assert(theta[i] == (*parameter_manager)(i));
  }
  */
  
  /*
  const vector<maxfun_parameter>&  parameters = my_peeler.tcalc().mle().parameters();

  size_t  i = 0;  
  vector<maxfun_parameter>::const_iterator  iter;
  for(iter = parameters.begin(); iter != parameters.end(); ++iter, ++i)
  {
    assert(i < theta.size() && theta[i] == iter->value);
  }
  */
  
  ++nfe;
  
  return likelihood().get_log();
}

inline int
subped_calculator::update_bounds(parameter_vector& theta)
{
  return 0;
}


//============================================================================
// IMPLEMENTATION:  ped_calculator
//============================================================================
//
inline
ped_calculator::ped_calculator(const FPED::Pedigree& ped, const mle_sub_model& mle,
                                 size_t trait, size_t marker)
      : my_ped(ped), my_mle(mle), my_trait(trait), my_marker(marker),
        my_unlinked_likelihood(QNAN), unlinked_likelihood_cached(false)
{
  nfe = 0;
}

inline double
ped_calculator::evaluate(parameter_vector& theta)
{
  check_parameters(theta, my_mle);

  /*
  const vector<maxfun_parameter>&  parameters = my_mle.parameters();
  
  size_t  i = 0;
  vector<maxfun_parameter>::const_iterator  iter;
  for(iter = parameters.begin(); iter != parameters.end(); ++iter, ++i)
  {
    assert(i < theta.size() && theta[i] == iter->value);
  }
  */

  ++nfe;
  
  return  likelihood().get_log();
}

inline int
ped_calculator::update_bounds(parameter_vector& theta)
{
  return 0;
}


//============================================================================
// IMPLEMENTATION:  group_calculator
//============================================================================
//
inline
group_calculator::group_calculator(const group& g, const FPED::FilteredMultipedigree& mped, 
                                   const mle_sub_model& mle, size_t trait, size_t marker)
      : my_mle(mle), my_trait(trait), my_marker(marker)
{
  build_group(g, mped);
  nfe = 0;
}

inline double
group_calculator::evaluate(parameter_vector& theta)
{
  check_parameters(theta, my_mle);
  
  /*
  const vector<maxfun_parameter>&  parameters = my_mle.parameters();
  
  size_t  i = 0;
  vector<maxfun_parameter>::const_iterator  iter;
  for(iter = parameters.begin(); iter != parameters.end(); ++iter, ++i)
  {
    assert(i < theta.size() && theta[i] == iter->value);
  }
  */

  ++nfe;
  
  return  likelihood().get_log();
}

inline int
group_calculator::update_bounds(parameter_vector& theta)
{
  return 0;
}


//============================================================================
// IMPLEMENTATION:  mped_calculator
//============================================================================
//
inline
mped_calculator::mped_calculator(const FPED::FilteredMultipedigree& mped, const mle_sub_model& mle,
                                 size_t trait, size_t marker)
      : my_mped(mped), my_mle(mle), my_trait(trait), my_marker(marker),
        my_unlinked_likelihood(QNAN), unlinked_likelihood_cached(false)
{
  nfe = 0;
}

inline double
mped_calculator::evaluate(parameter_vector& theta)
{
  check_parameters(theta, my_mle);
  
  /*
  const vector<maxfun_parameter>&  parameters = my_mle.parameters();
  
  // cout << endl;
  
  size_t  i = 0;
  vector<maxfun_parameter>::const_iterator  iter;
  for(iter = parameters.begin(); iter != parameters.end(); ++iter, ++i)
  {
    // cout << "theta[" << i << "] = " << theta[i] << endl;
    assert(i < theta.size() && (theta[i] == iter->value));  // || ((SAGE::isnan(theta[i]) && SAGE::isnan(iter->value))));
  }
  */

  double  like = likelihood().get_log();
  // cout << "ln likelihood = " << like << endl;

  ++nfe;
  
  return like;
}

inline int
mped_calculator::update_bounds(parameter_vector& theta)
{
  return 0;
}




      

