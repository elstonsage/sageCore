//============================================================================
//
//  File:	AnalysisOutput.cpp
//
//  Author:	Stephen Gross
//
//  Copyright 2002 R. C. Elston
//============================================================================

#include "maxfunapi/OutputFormatter.h"
#include "ageon/AnalysisOutput.h"

namespace SAGE {
namespace AO   {

//============================================================================
//
//  AnalysisOutput() CONSTRUCTOR
//
//============================================================================
AnalysisOutput::AnalysisOutput(const SAMPLING::PartitionedMemberDataSample & sample) :
        my_sample              (sample)
{ 
  my_models        .resize(num_of_analysis_types);
  my_maxfun_datas  .resize(num_of_analysis_types);
  my_Results       .resize(num_of_analysis_types);
  my_model_traits  .resize(num_of_analysis_types);
}

//============================================================================
//
//  AnalysisOutput() COPY CONSTRUCTOR
//
//============================================================================
AnalysisOutput::AnalysisOutput(const AnalysisOutput & other) :
  my_sample(other.my_sample)
{ 
  my_maxfun_datas        = other.my_maxfun_datas;
  my_Results             = other.my_Results;
  my_models              = other.my_models;
  my_model_traits        = other.my_model_traits;
}

//============================================================================
//  get_model(...)
//============================================================================
const Model & AnalysisOutput::get_model(int t) const { return my_models[t]; }

//============================================================================
//  get_models()
//============================================================================
const vector<Model> & AnalysisOutput::get_models() const { return my_models; }

//============================================================================
//  get_maxfun_data(...)
//============================================================================
const Maxfun_Data & AnalysisOutput::get_maxfun_data(int t) const { return my_maxfun_datas[t]; }

//============================================================================
//  get_Results(...)
//============================================================================
const MAXFUN::Results & AnalysisOutput::get_Results(int t) const { return *my_Results[t]; }

//============================================================================
//  get_sample()
//============================================================================
const SAMPLING::PartitionedMemberDataSample & AnalysisOutput::get_sample() const { return my_sample; }

//============================================================================
//  input(...) #1
//============================================================================
void AnalysisOutput::input(const Model & mod, int t) { my_models[t] = mod; }

//============================================================================
//  input(...) #2
//============================================================================
void AnalysisOutput::input(const Maxfun_Data & data, int t) { my_maxfun_datas[t] = data; }

//============================================================================
//  input(...) #2B
//============================================================================
void AnalysisOutput::input(const MAXFUN::Results & results, int t) { my_Results[t] = new MAXFUN::Results(results); }

//============================================================================
//  input(...) #3
//============================================================================
void AnalysisOutput::input(const MemberCovariateCalculator & mcc, int t) 
{ 
  ExtraOutput::populateModelTraits(mcc, my_model_traits[t]);
}

//============================================================================
//
//  generate_output(...)
//
//============================================================================
OUTPUT::Section
AnalysisOutput::generate_output(bool detailed)
{
  OUTPUT::Section s;
  s << generateSampleHeader()
    << generateHeader()
    << generateClassSystem()
    << generateClasses()
    << MAXFUN::OutputFormatter::convertEstimates(*my_Results[NO_TRUNCATION | SUSCEPTIBILITIES_EQUAL], detailed)
    << MAXFUN::OutputFormatter::convertEstimates(*my_Results[NO_TRUNCATION | SUSCEPTIBILITIES_FREE],  detailed);

  if(detailed)
  {
    s << MAXFUN::OutputFormatter::convertEstimates(*my_Results[USE_TRUNCATION | SUSCEPTIBILITIES_EQUAL], detailed)
      << MAXFUN::OutputFormatter::convertEstimates(*my_Results[USE_TRUNCATION | SUSCEPTIBILITIES_FREE],  detailed);
  }

  s << MAXFUN::JointTest(*my_Results[NO_TRUNCATION | SUSCEPTIBILITIES_FREE], *my_Results[NO_TRUNCATION | SUSCEPTIBILITIES_EQUAL]).summarizeAsTable();

  if(detailed)
  {
    s << MAXFUN::JointTest(*my_Results[USE_TRUNCATION | SUSCEPTIBILITIES_FREE], *my_Results[USE_TRUNCATION | SUSCEPTIBILITIES_EQUAL]).summarizeAsTable()
      << MAXFUN::OutputFormatter::convertMatrix(*my_Results[NO_TRUNCATION  | SUSCEPTIBILITIES_EQUAL])
      << MAXFUN::OutputFormatter::convertMatrix(*my_Results[NO_TRUNCATION  | SUSCEPTIBILITIES_FREE])
      << MAXFUN::OutputFormatter::convertMatrix(*my_Results[USE_TRUNCATION | SUSCEPTIBILITIES_EQUAL])
      << MAXFUN::OutputFormatter::convertMatrix(*my_Results[USE_TRUNCATION | SUSCEPTIBILITIES_FREE]);
  }

  return s;
}

//============================================================================
//
//  generateSampleHeader()
//
//============================================================================
OUTPUT::Table 
AnalysisOutput::generateSampleHeader()
{
  size_t all_num = 0;
  size_t val_num = 0;
  for( size_t i = 0; i < my_sample.getPartitionCount()-1; ++i )
  {
    all_num += my_sample.getIndividualCount(i);
    val_num += my_sample.getPartitionInvalidIndividualCount(i);
  }

  return OUTPUT::Table("SAMPLE DESCRIPTION")
         << (OUTPUT::TableRow() << "Number of pedigrees in dataset"         << my_sample.getMultipedigree().info().get_source_rped()->pedigree_count())
         << (OUTPUT::TableRow() << "Number of analyzable pedigrees"         << my_sample.getMultipedigree().pedigree_count())
         <<  OUTPUT::Table::INSERT_BLANK_ROW()
         << (OUTPUT::TableRow() << "Number of individuals in dataset"       << my_sample.getMultipedigree().info().get_source_rped()->member_count())
         << (OUTPUT::TableRow() << "Number of sibs in dataset"              << all_num)
         << (OUTPUT::TableRow() << "Number of analyzable sibs"              << all_num-val_num)
         << (OUTPUT::TableRow() << "Number of sibs with missing age values" << val_num);
}

//============================================================================
//
//  generateHeader()
//
//============================================================================
OUTPUT::Table 
AnalysisOutput::generateHeader()
{
  return OUTPUT::Table("MODEL DESCRIPTION")
         << (OUTPUT::TableRow() << "Title"              << my_models[0].get_title              ())
         << (OUTPUT::TableRow() << "Affectedness trait" << my_models[0].get_affectedness_trait ())
         << (OUTPUT::TableRow() << "Age-of-onset trait" << my_models[0].get_age_of_onset_trait ())
         << (OUTPUT::TableRow() << "Age-at-exam trait"  << my_models[0].get_age_of_exam_trait  ());

  OUTPUT::Table t("MODEL DESCRIPTION");

  t << (OUTPUT::TableRow() << "Title"              << my_models[0].get_title              ())
    << (OUTPUT::TableRow() << "Affectedness trait" << my_models[0].get_affectedness_trait ())
    << (OUTPUT::TableRow() << "Age-of-onset trait" << my_models[0].get_age_of_onset_trait ())
    << (OUTPUT::TableRow() << "Age-at-exam trait"  << my_models[0].get_age_of_exam_trait  ());

  if(  my_models[0].get_allow_averaging() )
    t << (OUTPUT::TableRow() << "Allow_averaging" << "TRUE");
}

//======================================================
//
//  generateClassSystem()
//
//======================================================
OUTPUT::Table
AnalysisOutput::generateClassSystem()
{
  OUTPUT::Table t("CLASSIFICATION SYSTEM");
  
  if( my_models[0].get_pool_class() )
  {
    t << OUTPUT::Table::INSERT_ROW_MESSAGE("Class(es) pooled.");

    map<size_t, size_t>::const_iterator mi = my_models[0].get_class_type_map().begin();
    for( ; mi != my_models[0].get_class_type_map().end(); ++mi )
    {
      if( mi->first != mi->second )
      {
        std::ostringstream mesg;

        mesg << "- class " << my_models[0].get_pooling_cfg().getCode(mi->first)
             << " pooled with " << my_models[0].get_pooling_cfg().getCode(mi->second);

        t << (OUTPUT::TableRow() << "  " << mesg.str());
      }
    }

    t << (OUTPUT::TableRow() << "\n")
      << (OUTPUT::TableRow() << "??" << "Both parents are unknown.")
      << (OUTPUT::TableRow() << "?A" << "One of the parents is unknown, the other is affected.")
      << (OUTPUT::TableRow() << "?U" << "One of the parents is unknown, the other is unaffected.")
      << (OUTPUT::TableRow() << "AA" << "Both parents are affected.")
      << (OUTPUT::TableRow() << "AU" << "One of the parents is affected, the other is unaffected.")
      << (OUTPUT::TableRow() << "UU" << "Both parents are unaffected.");

    return t;
  }
  else
  {
    return t << OUTPUT::Table::INSERT_ROW_MESSAGE("Using default classification system:")

             << (OUTPUT::TableRow() << "??" << "Both parents are unknown.")
             << (OUTPUT::TableRow() << "?A" << "One of the parents is unknown, the other is affected.")
             << (OUTPUT::TableRow() << "?U" << "One of the parents is unknown, the other is unaffected.")
             << (OUTPUT::TableRow() << "AA" << "Both parents are affected.")
             << (OUTPUT::TableRow() << "AU" << "One of the parents is affected, the other is unaffected.")
             << (OUTPUT::TableRow() << "UU" << "Both parents are unaffected.");
  }

}

//============================================================================
//
//  generateClasses()
//
//============================================================================
OUTPUT::Section
AnalysisOutput::generateClasses()
{
#if 0
  cout << my_sample.getTotalIndividualCount() << " "
       << my_sample.getValidIndividualCount() << " "
       << my_sample.getInvalidIndividualCount() << " "
       << my_sample.getIndividualCount() << " "
       << my_sample.getPartitionCount() << endl;
#endif

  OUTPUT::Section s("CLASS STATISTICS");

  for( size_t i = 0; i < my_sample.getPartitionCount(); ++i )
  {
    if( i >= my_models[0].num_of_classes() )
      continue;

#if 0
  cout << i << " "
       << my_sample.getIndividualCount(i) << " "
       << my_sample.getPartitionValidIndividualCount(i) << " "
       << my_sample.getPartitionInvalidIndividualCount(i) << endl;
#endif

    OUTPUT::Table t;

    std::ostringstream class_name;
    class_name << "CLASS ";

    // If this is the default  classification system, append the meaningful name to the table title:
    if( my_models[0].get_class_trait() == AO_CLASS_TYPE )
      class_name << PoolingCfg::getCode(i);
    else
      class_name << i;
    
    t.setTitle(class_name.str());
    
    if(    !my_sample.getIndividualCount(i)
        || !my_sample.getPartitionValidIndividualCount(i) )
    {
      t << OUTPUT::NamedString("Note", "Contains no valid individuals for analysis.");
    }
    else
    {
      //size_t total_cnt = my_sample.getIndividualCount(i);
      size_t total_used_cnt = my_sample.getPartitionValidIndividualCount(i);
      size_t total_affected_cnt = 0;

      SampleInfo AO_all;
      SampleInfo AO_affected;
      SampleInfo AO_affected_1(1);

      SampleInfo AE_all;
      SampleInfo AE_affected;
      SampleInfo AE_unaffected;
      SampleInfo AE_unaffected_1(1);

      vector<size_t>::const_iterator j = my_sample.getIndividualBegin(i);
      for( ; j != my_sample.getIndividualEnd(i); ++j )
      {
	if( !my_sample.isValid(*j) )
          continue;

        if( !SAGE::isnan(my_sample.getAdjValue(*j, "CORE_TRAITS", "age of onset")) )
          AO_all.add(my_sample.getAdjValue(*j, "CORE_TRAITS", "age of onset"));

        if( !SAGE::isnan(my_sample.getAdjValue(*j, "CORE_TRAITS", "age at exam")) )
          AE_all.add(my_sample.getAdjValue(*j, "CORE_TRAITS", "age at exam"));

        // If the individual has an affectedness trait, and it's positive, increment the affected count:
        if( my_sample.getField("CORE_TRAITS", "affectedness").isAdjValuePresent(*j) == true )
        {
          if( my_sample.getAdjValue(*j, "CORE_TRAITS", "affectedness") == 1 )
          {
            total_affected_cnt++;

            if( !SAGE::isnan(my_sample.getAdjValue(*j, "CORE_TRAITS", "age of onset")) )
            {
              AO_affected.add(my_sample.getAdjValue(*j, "CORE_TRAITS", "age of onset"));
              AO_affected_1.add(my_sample.getAdjValue(*j, "CORE_TRAITS", "age of onset"));
            }

            if( !SAGE::isnan(my_sample.getAdjValue(*j, "CORE_TRAITS", "age at exam")) )
              AE_affected.add(my_sample.getAdjValue(*j, "CORE_TRAITS", "age at exam"));
          }
          else
          {
            if( !SAGE::isnan(my_sample.getAdjValue(*j, "CORE_TRAITS", "age at exam")) )
            {
              AE_unaffected.add(my_sample.getAdjValue(*j, "CORE_TRAITS", "age at exam"));
              AE_unaffected_1.add(my_sample.getAdjValue(*j, "CORE_TRAITS", "age at exam"));
            }
          }
        }

      } // End of loop-across-individuals
      
      //t << (OUTPUT::TableRow() << "TOTAL NUMBER OF INDIVIDUALS"            << total_cnt);
      t << (OUTPUT::TableRow() << "TOTAL NUMBER OF INDIVIDUALS USED IN ANALYSIS" << total_used_cnt);
      /*
      // AO for all
      t << OUTPUT::Table::INSERT_BLANK_ROW()
        << (OUTPUT::TableRow() << "NUMBER OF INDIVIDUALS WITH AN AGE OF ONSET"   << AO_all.count())
        << (OUTPUT::TableRow() << "MEAN OF AGE OF ONSET"                         << AO_all.mean());

      OUTPUT::TableRow ao_all_variance_row;

      ao_all_variance_row << "VARIANCE OF AGE OF ONSET";

      if( !SAGE::isnan(AO_all.variance()) )
        ao_all_variance_row << AO_all.variance();
      else if( AO_all.min() == AO_all.max() )
        ao_all_variance_row << 0.0;
      else
        ao_all_variance_row << OUTPUT::UnavailableCell();
      
      t << ao_all_variance_row;
      */
      // AO for affected_1
      t << (OUTPUT::TableRow() << "NUMBER OF AFF INDIVIDUALS WITH AN AGE OF ONSET" << AO_affected_1.count())
        << (OUTPUT::TableRow() << "MEAN OF AGE OF ONSET"                           << AO_affected_1.mean());

      OUTPUT::TableRow ao_aff_1_variance_row;

      ao_aff_1_variance_row << "VARIANCE OF AGE OF ONSET";

      if( !SAGE::isnan(AO_affected_1.variance()) )
        ao_aff_1_variance_row << AO_affected_1.variance();
      else if( AO_affected_1.min() == AO_affected_1.max() )
        ao_aff_1_variance_row << 0.0;
      else
        ao_aff_1_variance_row << OUTPUT::UnavailableCell();
      
      t << ao_aff_1_variance_row;

      // affected counts
      t //<< OUTPUT::Table::INSERT_BLANK_ROW()
        << (OUTPUT::TableRow() << "NUMBER OF INDIVIDUALS AFFECTED"     << total_affected_cnt)
        //<< (OUTPUT::TableRow() << "PROPORTION OF INDIVIDUALS AFFECTED(TOTAL)" << ((double)total_affected_cnt / (double)total_cnt))
        << (OUTPUT::TableRow() << "PROPORTION OF INDIVIDUALS AFFECTED" << ((double)total_affected_cnt / (double)total_used_cnt));

      /*
      // AO for affected
      t << (OUTPUT::TableRow() << "NUMBER OF AFF INDIVIDUALS WITH AN AGE OF ONSET" << AO_affected.count())
        << (OUTPUT::TableRow() << "MEAN OF AGE OF ONSET"                           << AO_affected.mean());

      OUTPUT::TableRow ao_aff_variance_row;

      ao_aff_variance_row << "VARIANCE OF AGE OF ONSET (divided by n)";

      if( !SAGE::isnan(AO_affected.variance()) )
        ao_aff_variance_row << AO_affected.variance();
      else if( AO_affected.min() == AO_affected.max() )
        ao_aff_variance_row << 0.0;
      else
        ao_aff_variance_row << OUTPUT::UnavailableCell();
      
      t << ao_aff_variance_row;
      */
      // AE for all
      /*
      t << (OUTPUT::TableRow() << "NUMBER OF INDIVIDUALS WITH AN AGE OF EXAM"   << AE_all.count())
        << (OUTPUT::TableRow() << "MEAN OF AGE OF EXAM"                         << AE_all.mean());

      OUTPUT::TableRow ae_all_variance_row;

      ae_all_variance_row << "VARIANCE OF AGE OF EXAM";

      if( !SAGE::isnan(AE_all.variance()) )
        ae_all_variance_row << AE_all.variance();
      else if( AE_all.min() == AE_all.max() )
        ae_all_variance_row << 0.0;
      else
        ae_all_variance_row << OUTPUT::UnavailableCell();
      
      t << ae_all_variance_row;

      // AE for affected
      t << (OUTPUT::TableRow() << "NUMBER OF AFF INDIVIDUALS WITH AN AGE OF EXAM" << AE_affected.count())
        << (OUTPUT::TableRow() << "MEAN OF AGE OF EXAM"                           << AE_affected.mean());

      OUTPUT::TableRow ae_aff_variance_row;

      ae_aff_variance_row << "VARIANCE OF AGE OF EXAM";

      if( !SAGE::isnan(AE_affected.variance()) )
        ae_aff_variance_row << AE_affected.variance();
      else if( AE_affected.min() == AE_affected.max() )
        ae_aff_variance_row << 0.0;
      else
        ae_aff_variance_row << OUTPUT::UnavailableCell();
      
      t << ae_aff_variance_row;
      */
      // AE for unaffected_1
      t //<< OUTPUT::Table::INSERT_BLANK_ROW()
        << (OUTPUT::TableRow() << "MEAN OF AGE OF EXAM OF THE UNAFFECTED"           << AE_unaffected_1.mean());

      OUTPUT::TableRow ae_unaff_1_variance_row;

      ae_unaff_1_variance_row << "VARIANCE OF AGE OF EXAM OF THE UNAFFECTED";

      if( !SAGE::isnan(AE_unaffected_1.variance()) )
        ae_unaff_1_variance_row << AE_unaffected_1.variance();
      else if( AE_unaffected.min() == AE_unaffected_1.max() )
        ae_unaff_1_variance_row << 0.0;
      else
        ae_unaff_1_variance_row << OUTPUT::UnavailableCell();
      
      t << ae_unaff_1_variance_row;

      // AE for unaffected
      /*
      t << (OUTPUT::TableRow() << "NUMBER OF INDIVIDUALS WITH AN AGE OF EXAM" << AE_unaffected.count())
        << (OUTPUT::TableRow() << "MEAN OF AGE OF EXAM"                       << AE_unaffected.mean());

      OUTPUT::TableRow ae_unaff_variance_row;

      ae_unaff_variance_row << "VARIANCE OF AGE OF EXAM OF THE UNAFFECTED (divided by n)";

      if( !SAGE::isnan(AE_unaffected.variance()) )
        ae_unaff_variance_row << AE_unaffected.variance();
      else if( AE_unaffected.min() == AE_unaffected.max() )
        ae_unaff_variance_row << 0.0;
      else
        ae_unaff_variance_row << OUTPUT::UnavailableCell();
      
      t << ae_unaff_variance_row;
      */
    }
  
    s << t;
    
  }

  return s;
}

} // End namespace AO
} // End namespace SAGE
