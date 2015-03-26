//======================================================================
//
//  File:  Model.cpp
//
//  Author:  Stephen Gross
//
//  Copyright 2002 R. C. Elston
//======================================================================

#include "ageon/Model.h"

namespace SAGE {
namespace AO {

//======================================================================
//
//  Model() CONSTRUCTOR
//
//======================================================================
Model::Model()
{
  reset();
}

//======================================================================
//
//  copy(...)
//
//======================================================================
void
Model::copy(const Model & other)
{
  lambda1_mxid             = other.lambda1_mxid;
  lambda2_mxid             = other.lambda2_mxid;

  AO_id                    = other.AO_id;
  AE_id                    = other.AE_id;
  aff_id                   = other.aff_id;
  class_id                 = other.class_id;

  my_affectedness_trait    = other.my_affectedness_trait;
  my_age_of_exam_trait     = other.my_age_of_exam_trait;
  my_age_of_onset_trait    = other.my_age_of_onset_trait;
  my_allow_averaging       = other.my_allow_averaging;
  my_pool_class            = other.my_pool_class;
  my_class_trait           = other.my_class_trait;
  my_debug                 = other.my_debug;
  my_debug_cfg             = other.my_debug_cfg;
  my_epsilon               = other.my_epsilon;
  my_parameter_mgr         = other.my_parameter_mgr;
  my_min_denominator       = other.my_min_denominator;
  my_num_of_classes        = other.my_num_of_classes;
  my_ofilename             = other.my_ofilename;
  my_pooling_cfg           = other.my_pooling_cfg;
  my_title                 = other.my_title;
  my_transf_sub_model      = other.my_transf_sub_model;
  my_truncate              = other.my_truncate;
  my_use_adjustment        = other.my_use_adjustment;
  my_class_type_map        = other.my_class_type_map;
  my_valid                 = other.my_valid;
}

//======================================================================
//
//  Model(...) COPY CONSTRUCTOR
//
//======================================================================
Model::Model(const Model & other)
{
  copy(other);
}

//======================================================================
//
//  operator=(...)
//
//======================================================================
Model &
Model::operator=(const Model & other)
{
  copy(other);

  return *this;
}

//======================================================================
//
//  reset()
//
//======================================================================
void
Model::reset()
{
  my_debug                        = false;
  my_allow_averaging              = false;
  my_pool_class                   = false;
  my_valid                        = false;
  my_use_adjustment               = AO_DEFAULT_USE_ADJUSTMENT;
  my_title                        = AO_DEFAULT_ANALYSIS_TITLE;
  my_affectedness_trait           = "";
  my_age_of_onset_trait           = "";
  my_age_of_exam_trait            = "";
  my_ofilename                    = "";
  my_class_trait                  = AO_CLASS_TYPE;
  my_num_of_classes               = AO_DEFAULT_NUM_OF_CLASSES;
  my_epsilon                      = AO_DEFAULT_EPSILON;
  my_min_denominator              = 1.0;
  my_num_of_classes               = AO_DEFAULT_NUM_OF_CLASSES;
  my_truncate                     = AO_DEFAULT_USE_TRUNCATION;

  my_parameter_mgr.reset();

  my_parameter_mgr.addGroup("Susceptibility intercepts");
  my_parameter_mgr.addGroup("Susceptibility covariates");
  my_parameter_mgr.addGroup("Mean intercept");
  my_parameter_mgr.addGroup("Mean covariates");
  my_parameter_mgr.addGroup("Variance intercept");
  my_parameter_mgr.addGroup("Variance covariates");
  my_parameter_mgr.addGroup("Transformation");

  my_parameter_mgr.addParameter("Mean intercept",
      "Mean intercept", 
      MAXFUN::Parameter::INDEPENDENT, 
      AO_DEFAULT_MEAN,
      AO_LOWER_BOUND_MEAN,
      AO_UPPER_BOUND_MEAN);

  my_parameter_mgr.getParameter("Mean intercept", "Mean intercept").setNameAbbr("Mean int.");

  my_parameter_mgr.addParameter("Variance intercept",
      "Variance intercept", 
      MAXFUN::Parameter::INDEPENDENT, 
      AO_DEFAULT_VARIANCE,
      AO_LOWER_BOUND_VARIANCE,
      AO_UPPER_BOUND_VARIANCE);

  my_parameter_mgr.getParameter("Variance intercept", "Variance intercept").setNameAbbr("Var. int.");

  lambda1_mxid = my_parameter_mgr.addParameter(
      "Transformation",
      "Lambda1", 
      AO_DEFAULT_LAMBDA1_TYPE,
      AO_DEFAULT_LAMBDA1,
      AO_LOWER_BOUND_LAMBDA1,
      AO_UPPER_BOUND_LAMBDA1);

  lambda2_mxid = my_parameter_mgr.addParameter(
      "Transformation",
      "Lambda2", 
      AO_DEFAULT_LAMBDA2_TYPE,
      AO_DEFAULT_LAMBDA2,
      AO_LOWER_BOUND_LAMBDA2,
      AO_UPPER_BOUND_LAMBDA2);

}

//======================================================================
//
//  update_after_parse()
//
//======================================================================
void
Model::update_after_parse()
{
  for( size_t i = 0; i < my_num_of_classes; i++ )
  {
    std::ostringstream class_name;
    
    if( get_class_trait() == AO_CLASS_TYPE )
      class_name << PoolingCfg::getCode(i);
    else
      class_name << "Class " << i;

    my_parameter_mgr.addParameter("Susceptibility intercepts", 
                                  class_name.str(),
                                  MAXFUN::Parameter::INDEPENDENT, 
                                  AO_DEFAULT_GEN_SUSCEPT,
                                  AO_LOWER_BOUND_GEN_SUSCEPT,
                                  AO_UPPER_BOUND_GEN_SUSCEPT);
  }
}

//======================================================================
//
//  update()
//
//======================================================================
int
Model::update(vector<double> & params, int t)
{
  // 1. If susceptibilities are equal, update the dependent ones:

  if(SusceptibilitiesEqual(t))
  {
    MAXFUN::ParameterIterator p = my_parameter_mgr.getParamBegin ("Susceptibility intercepts");
        
    double val = p->getCurrentEstimate();

    p++;
         
    for(; p != my_parameter_mgr.getParamEnd("Susceptibility intercepts"); p++)
      p->setCurrentEstimate(val);
  }

  // 2. Return success:

  return 0;
}

//======================================================================
//
//  setupSample(...)
//
//======================================================================
void
Model::setupSample(SAMPLING::PartitionedMemberDataSample & sample)
{
  if( get_pool_class() )
  {
    set<size_t> pool_vals;
    for( size_t j = 0; j < 6; ++j )
    {
      size_t new_j = my_pooling_cfg.getPool(j);
      pool_vals.insert(new_j);
      my_class_type_map[j] = new_j;
    }
  }
  else
  {
    for( size_t j = 0; j < 6; ++j )
    {
      my_class_type_map[j] = j;
    }
  }

#if 0
    cout << "my_num_of_classes = " << my_num_of_classes << endl;
    cout << "my_class_type_map" << endl;
    map<size_t, size_t>::const_iterator mi = my_class_type_map.begin();
    for( ; mi != my_class_type_map.end(); ++mi )
      cout << mi->first << " " << mi->second << endl;
#endif

  update_after_parse();

  aff_id   = sample.importField(get_affectedness_trait(), "CORE_TRAITS", "affectedness");
  AO_id    = sample.importField(get_age_of_onset_trait(), "CORE_TRAITS", "age of onset");
  AE_id    = sample.importField(get_age_of_exam_trait(),  "CORE_TRAITS", "age at exam");
  class_id = sample.importField(get_class_trait(),        "CORE_TRAITS", "classification");

  unsigned long flags = get_allow_averaging() ? SAMPLING::Field::ALLOW_AVERAGING : SAMPLING::Field::NO_FLAGS;

  sample.addGroup("Susceptibility covariates");

  for(MAXFUN::ParameterConstIterator p  = GetParameterMgr().getParamBegin ("Susceptibility covariates");
                                     p != GetParameterMgr().getParamEnd   ("Susceptibility covariates"); p++)
  {
    sample.importField(p->getName(), "Susceptibility covariates", p->getName(), flags);
  }

  sample.addGroup("Mean covariates");

  for(MAXFUN::ParameterConstIterator p  = GetParameterMgr().getParamBegin ("Mean covariates");
                                     p != GetParameterMgr().getParamEnd   ("Mean covariates"); p++)
  {
    sample.importField(p->getName(), "Mean covariates", p->getName(), flags);
  }

  sample.addGroup("Variance covariates");

  for(MAXFUN::ParameterConstIterator p  = GetParameterMgr().getParamBegin ("Variance covariates");
                                     p != GetParameterMgr().getParamEnd   ("Variance covariates"); p++)
  {
    sample.importField(p->getName(), "Variance covariates", p->getName(), flags);
  }

  // Validate individuals:

  AO::Validator validator;

  validator.setMultiPedigree(sample.getMultipedigree());

  sample.finalizeData(validator);

  sample.finalizeUserCreatedData();

  // Classify individuals:

#if 0
  cout << "sample.getIndividualCount() = " << sample.getIndividualCount() << endl;
  cout << "sample.getPartitionCount() = " << sample.getPartitionCount() << endl;
  cout << "sample.getTotalIndividualCount() = " << sample.getTotalIndividualCount() << endl;
  sample.dumpPartitions();
  sample.dumpTraitValues();
#endif

  for(size_t j = 0; j < sample.getTotalIndividualCount(); ++j)
  {
    const FPED::Member& ind = sample.getMultipedigree().member_index(j);

    size_t c_type = getClassType(sample, ind);
#if 0
    cout << "add ind " << j
         << " " << ind.pedigree()->name()
         << " " << ind.name()
         << " c_type = " << c_type << endl;
#endif
    if( get_class_trait() == AO_CLASS_TYPE && c_type > 5 )
      c_type = 6;

    sample.addIndividual(j, c_type);
  }

#if 0
  cout << "sample.getIndividualCount() = " << sample.getIndividualCount() << endl;
  cout << "sample.getPartitionCount() = " << sample.getPartitionCount() << endl;
  cout << "sample.getTotalIndividualCount() = " << sample.getTotalIndividualCount() << endl;
  sample.dumpPartitions();
  sample.dumpTraitValues();
#endif

  // Generate partitioned fields (PartitionedMemberDataSample)
  sample.generatePartitionedFields();

#if 0
  cout << "sample.getIndividualCount() = " << sample.getIndividualCount() << endl;
  cout << "sample.getPartitionCount() = " << sample.getPartitionCount() << endl;
  cout << "sample.getTotalIndividualCount() = " << sample.getTotalIndividualCount() << endl;
  sample.dumpPartitions();
  sample.dumpTraitValues();
#endif
}

//==========================================================================
//
//  getClassType(...)
//
//==========================================================================
size_t
Model::getClassType(
  const SAMPLING::PartitionedMemberDataSample & sample,
  const FPED::Member & ind)
{
  // 0. Set up local variables:

  size_t class_type = 0;

  // 1. If this is a user-defined classification, fetch it:

  if( get_class_trait() != AO_CLASS_TYPE )
  {
    size_t class_trait_num = sample.getMultipedigree().info().trait_find(get_class_trait());

    if( ind.pedigree()->info().trait_missing(ind.index(), class_trait_num) )
      class_type = my_num_of_classes;
    else
    {
      size_t class_val = (size_t)(ind.pedigree()->info().trait(ind.index(), class_trait_num));
      class_type = my_class_type_map[class_val];
    }
  }

  // 2. Otherwise, use the default classification system:

  else
  {

  // 2.1. If this individual is a founder, return 0. The individual won't be
  //      used in the analysis, anyway...
  //      NOTE TO SELF: Check to make sure that founders aren't included in class statistics counts

    if( MPED::mp_utilities::is_founder(ind) )
    {
      class_type = 6;
#if 0
      cout << "founder ";
#endif
    }

  // 2.2. If this individual is a non-founder...

    else
    {
      size_t parent1_idx       = ind.parent1()->mpindex(),
             parent2_idx       = ind.parent2()->mpindex();
      double d_aff_parent1     = sample.getAdjValue(parent1_idx, aff_id),
             d_aff_parent2     = sample.getAdjValue(parent2_idx, aff_id);

#if 0
      size_t aff_num = sample.getMultipedigree().info().trait_find(get_affectedness_trait());

      size_t parent1_idx       = ind.parent1()->index(),
             parent2_idx       = ind.parent2()->index();
      double d_aff_parent1     = ind.pedigree()->info().trait(parent1_idx, aff_num),
             d_aff_parent2     = ind.pedigree()->info().trait(parent2_idx, aff_num);
#endif
#if 0
      cout << " p1 = " << ind.parent1()->name() << " " << d_aff_parent1
           << " p2 = " << ind.parent2()->name() << " " << d_aff_parent2;
#endif

      bool   parent1_aff_known = !SAGE::isnan(d_aff_parent1),
             parent2_aff_known = !SAGE::isnan(d_aff_parent2),
             parent1_aff       = parent1_aff_known ? d_aff_parent1 : false,
             parent2_aff       = parent2_aff_known ? d_aff_parent2 : false;

           if( !parent1_aff_known && !parent2_aff_known)                     class_type = 0;
      else if((!parent1_aff_known &&  parent2_aff_known &&  parent2_aff) ||
              (!parent2_aff_known &&  parent1_aff_known &&  parent1_aff))    class_type = 1;
      else if((!parent1_aff_known &&  parent2_aff_known && !parent2_aff) ||
              (!parent2_aff_known &&  parent1_aff_known && !parent1_aff))    class_type = 2;
      else if(  parent1_aff_known &&  parent1_aff       &&
                parent2_aff_known &&  parent2_aff  )                         class_type = 3;
      else if(( parent1_aff_known &&  parent1_aff &&
                parent2_aff_known && !parent2_aff) ||
              ( parent1_aff_known && !parent1_aff &&
                parent2_aff_known &&  parent2_aff))                          class_type = 4;
      else if(  parent1_aff_known && !parent1_aff &&
                parent2_aff_known && !parent2_aff)                           class_type = 5;
    }
  }

  // X. Return reassigned class type:
#if 0
  cout << " class type = " << class_type << " ";
#endif

  if( get_class_trait() != AO_CLASS_TYPE )
    return class_type;

  return my_pooling_cfg.getPool(class_type);
}

//======================================================================
//
//  calculate_initial_est(...)
//
//======================================================================
void
Model::calculate_initial_est(const SAMPLING::PartitionedMemberDataSample & sample, int t)
{
  // 1. Calculate default mean_base and var_base:

  double exp_variance = sample.getField(AO_id).getVariance (),
         exp_mean     = sample.getField(AO_id).getMean     ();

  my_parameter_mgr.getParameter("Mean intercept",     "Mean intercept")    . setInitialEstimate(exp_mean);
  my_parameter_mgr.getParameter("Variance intercept", "Variance intercept"). setInitialEstimate(exp_variance);

  if( SusceptibilitiesFree(t) )
  {
    MAXFUN::ParameterIterator p = my_parameter_mgr.getParamBegin("Susceptibility intercepts");
    for( ; p != my_parameter_mgr.getParamEnd("Susceptibility intercepts"); p++ )
    {
      p->setInitialType(sample.getPartitionValidIndividualCount(p->getGroupIndex()) < 5 ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT);
    }
  }
  else
  {
    MAXFUN::ParameterIterator p = my_parameter_mgr.getParamBegin("Susceptibility intercepts");

    p->setInitialType(MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL);

    p++;

    for( ; p != my_parameter_mgr.getParamEnd("Susceptibility intercepts"); p++ )
      p->setInitialType(MAXFUN::Parameter::DEPENDENT);
  }

  // 2. Check those initial values for errors:

  if(SAGE::isnan(my_parameter_mgr.getParameter("Mean intercept", "Mean intercept").getInitialEstimate()))
  {
    cout << endl
         << "Error: Sample set does not yield a valid initial mean estimate. Exiting."
         << endl;

    exit(0);
  }
        
  if(SAGE::isnan(my_parameter_mgr.getParameter("Variance intercept", "Variance intercept").getInitialEstimate()))
  {
    cout << endl
         << "Error: Sample set does not yield a valid initial variance estimate. Exiting."
         << endl;
 
    exit(0);
  }
}

//======================================================================
//
//  add_trait(...)
//
//======================================================================
void
Model::add_trait(trait_type           t, 
                 string               name,
                 bool                 new_initial_est,
                 double               initial_est,
                 bool                 new_fixed,
                 bool                 fixed)
{
  // 2. Check to see if we are just modifying existing traits in the model:

  MAXFUN::Parameter* existing_param = NULL;

  string comp_name = toUpper(name);

  if(comp_name == toUpper(GetParameterMgr().getParameter("Transformation", "Lambda1").getName()))
    existing_param = &GetParameterMgr().getParameter("Transformation", "Lambda1");

  if(comp_name == toUpper(GetParameterMgr().getParameter("Transformation", "Lambda2").getName()))
    existing_param = &GetParameterMgr().getParameter("Transformation", "Lambda2");

  if(existing_param != NULL)
  {
    if(new_initial_est) existing_param->setInitialEstimate (initial_est);
    if(new_fixed)       existing_param->setInitialType     (fixed ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT);
  }

  // Have not found that the trait to be added already exists in the model. Hence,
  // we must add a new trait to one of the three vectors of traits.

  else
  {
    string group_name;

    switch(t)
    {
      case type_suscept : group_name = "Susceptibility covariates"; break;
      case type_mean    : group_name = "Mean covariates"; break;
      case type_var     : group_name = "Variance covariates"; break;
    }

    my_parameter_mgr.addParameter(  
        group_name, 
        name,
        MAXFUN::Parameter::INDEPENDENT,
        AO_DEFAULT_COVARIATE,
        AO_LOWER_BOUND_COVARIATE,
        AO_UPPER_BOUND_COVARIATE);
  
  }
}

//======================================================================
//
//  verify()
//
//======================================================================
void
Model::verify()
{
  my_valid = true;
}

void
Model::dump_model(ostream& o) const
{
  o << endl
    << "Analysis title: " << get_title() << endl
    << endl
    << "affectedness trait: " << get_affectedness_trait() << endl
    << "age_of_onset trait: " << get_age_of_onset_trait() << endl
    << "age_of_exam trait:  " << get_age_of_exam_trait() << endl
    << endl;

  if( get_allow_averaging() )
    o << "allow_averaging = true" << endl;
  else
    o << "allow_averaging = false" << endl;    

  if( get_pool_class() )
  {
    o << "pool_class = true " << endl;

    if( my_class_type_map.size() )
    {
      map<size_t, size_t>::const_iterator i = my_class_type_map.begin();
      for( ; i != my_class_type_map.end(); ++i )
      {
        if( i->first != i->second )
          o << "  class " << my_pooling_cfg.getCode(i->first)
            << " pooled with class " << my_pooling_cfg.getCode(i->second)
            << endl;
      }
    }

  }
  else
    o << "pool_class = false " << endl;

  o << endl;
}

} // End namespace AO
} // End namespace SAGE
