#include "func/AdjustedTraitCreator.h"
#include <iomanip>

namespace SAGE {
namespace FUNC {

bool
AdjustedTraitCreator::createAdjustedTrait(RPED::RefMultiPedigree& rmp,
                cerrorstream errors, const FunctionParser::TraitData& par_data)
{
  AdjustedTraitCreator  c(rmp, errors);

  return  c.internal_createAdjustedTrait(par_data);
}


AdjustedTraitCreator::AdjustedTraitCreator(RPED::RefMultiPedigree& rmp, cerrorstream errors) 
  : my_rmp(rmp), my_errors(errors), my_summary_stats(1)
{
  my_reference_trait = (size_t)(-1);
  my_bin_size = 0;
  my_class_stats.clear();
  my_class_traits.clear();
}

bool
AdjustedTraitCreator::internal_createAdjustedTrait(const FunctionParser::TraitData& par_data)
{
  int err_code = parseExpressionString(par_data);

  if(err_code != 0)
  {
    my_errors << priority(error) << "Fatal errors encountered in function expression ('"
              << par_data.expr << "')" << endl;

    if(err_code == 1)
    {
      my_errors << "Expression malformed." << endl;
    }
    else if(err_code == 2)
    {
      my_errors << "At least one of the traits referred to doesn't seem to exist in your data." << endl;
    }
    else if(err_code == 3)
    {
      my_errors << "The bin size argument is misspecified." << endl;
    }

    exit(1);
  }

  // Calculate the descriptive statistics for the trait, class and non-class-based.
  //
  err_code = calculateStats(par_data);

  if(err_code != 0)
  {
    my_errors << priority(error) << "Fatal errors encountered in function expression ('"
              << par_data.expr << "')" << endl;

    if(err_code == 1)
    {
      my_errors << "Not enough valid values present to perform binning algorithm." << endl;
    }

    exit(1);
  }

  // Now actually create the trait and populate it with values:
  //
  createIndividualNewTraitValues(par_data);

  return  true;
}

void
AdjustedTraitCreator::createIndividualNewTraitValues(const FunctionParser::TraitData& par_data)
{
  // Add the trait to the RefMPedInfo and to all the RefPedInfo's:
  //
  size_t  new_trait = addTrait(par_data);

  // Loop across all individuals and populate the new trait value for each:
  //
  for(RPED::PedigreeIterator p_iter = my_rmp.pedigree_begin(); p_iter != my_rmp.pedigree_end(); ++p_iter)
  {
    for(RPED::MemberIterator m_iter = p_iter->member_begin(); m_iter != p_iter->member_end(); ++m_iter)
    {
      double  reference_value = p_iter->info().trait(m_iter->index(), my_reference_trait);

      // Skip if missing.
      //
      if(SAGE::isnan(reference_value))
      {
        p_iter->info().set_trait(m_iter->index(), new_trait, std::numeric_limits<double>::quiet_NaN());
        continue;
      }

      double  mean  = my_summary_stats.mean();
      double  stdev = my_summary_stats.standard_deviation();

      if(my_adj_type == CLASS_BASED)
      {
        // Find out which class this individual belongs to:
        //
        size_t class_trait_num = (size_t)(-1);
        for(TraitVector::iterator c = my_class_traits.begin(); c!= my_class_traits.end(); ++c)
        {
          if(p_iter->info().trait(m_iter->index(), *c) == 1)
          {
            class_trait_num = *c;
            break;
          }
        }

        // The class could not be identified, tag this individual as missing and move on:
        //
        if(class_trait_num == (size_t)(-1))
        {
          p_iter->info().set_trait(m_iter->index(), new_trait, std::numeric_limits<double>::quiet_NaN());
          continue;
        }

        mean  = my_class_stats[class_trait_num].binned_mean;
        stdev = my_class_stats[class_trait_num].binned_stdev;
      }

      // Adjust value and add it
      //
      if(my_function_type == MEAN_ADJ)
      {
        p_iter->info().set_trait(m_iter->index(), new_trait, reference_value - mean);
      }
      else if(my_function_type == VAR_ADJ)
      {
        p_iter->info().set_trait(m_iter->index(), new_trait, reference_value / stdev);
      }
      else if(my_function_type == Z_SCORE)
      {
        p_iter->info().set_trait(m_iter->index(), new_trait, (reference_value - mean) / stdev);
      }
       
    } // End member loop

  } // End pedigree loop
}

size_t 
AdjustedTraitCreator::addTrait(const FunctionParser::TraitData& par_data)
{
  size_t trait_num = 0;

  // Add trait to RefMPedInfo:
  //
  if(par_data.binary)
  {
    trait_num = my_rmp.info().add_binary_trait(par_data.trait_name, par_data.usage);

    my_rmp.info().trait_info(trait_num).set_string_affected_code(par_data.affected);
    my_rmp.info().trait_info(trait_num).set_numeric_affected_code(str2doub(par_data.affected));
    my_rmp.info().trait_info(trait_num).set_string_unaffected_code(par_data.unaffected);
    my_rmp.info().trait_info(trait_num).set_numeric_unaffected_code(str2doub(par_data.unaffected));
  }
  else // par_data is continuous
  {
    trait_num = my_rmp.info().add_continuous_trait(par_data.trait_name, par_data.usage);
  }

  my_rmp.info().trait_info(trait_num).set_string_missing_code(par_data.missing);
  my_rmp.info().trait_info(trait_num).set_numeric_missing_code(str2doub(par_data.missing));

  // Now add the trait to all the RefPedInfo's:
  //
  for(RPED::PedigreeIterator p_iter = my_rmp.pedigree_begin(); p_iter != my_rmp.pedigree_end(); ++p_iter)
    p_iter->info().resize_traits(p_iter->info().trait_count() + 1);

  // Return new trait number:
  //
  return  trait_num;
}

int
AdjustedTraitCreator::parseExpressionString(const FunctionParser::TraitData& data)
{
  if(data.verbose)
    my_errors << priority(information) << "Processing '" << data.expr << "'..." << endl;

  // Verify well-formed string
  //
  vector<string> results = UTIL::RegexUtils::matchPattern("^\\s*(\\w+)\\s*\\((.+)\\)", data.expr);

  if(results.size() != 3)
    return 1;

  // Check process name:
  //
  string process_name = results[1];
  if(process_name == "mean_adj")
  {
    my_function_type = MEAN_ADJ;
  }
  else if(process_name == "var_adj")
  {
    my_function_type = VAR_ADJ;
  }
  else if(process_name == "z_score")
  {
    my_function_type = Z_SCORE;
  }
  else
  {
    return 1;
  }

  // Extract function arguments:
  //
  vector<string> function_arguments;

  UTIL::StringUtils::splitMultiDelimitedString(results[2], " ,", function_arguments);

  // Process function arguments:
  //
  if(function_arguments.size() == 0) // No arguments? 
  {
    return 1;
  }
  else if(function_arguments.size() == 1) // One argument = no class-based adjustment
  {
    // Set the adjustment type:
    //
    my_adj_type = NON_CLASS_BASED;

    // Make sure trait exists:
    //
    if(my_rmp.info().trait_exists(function_arguments[0]))
      my_reference_trait = my_rmp.info().trait_find(function_arguments[0]);
    else
      return 2;

    return 0;
  }
  else if(function_arguments.size() == 2) // Two arguments = malformed!
  {
    return 1;
  }
  else // At least three arguments = class-based adjustment
  {
    // Set the adjustment type:
    //
    my_adj_type = CLASS_BASED;

    // Reference trait 
    //
    if(my_rmp.info().trait_exists(function_arguments[0]))
      my_reference_trait = my_rmp.info().trait_find(function_arguments[0]);
    else
      return 2;      

    // Make sure bin size is well-formed:
    //
    if(UTIL::RegexUtils::doesPatternMatch("^\\d+$", function_arguments[1]))
      my_bin_size = atoi(function_arguments[1].c_str());
    else
      return 3;

    // Class traits.
    //
    for(size_t i = 2; i < function_arguments.size(); ++i)
    {
      if(my_rmp.info().trait_exists(function_arguments[i]))
        my_class_traits.push_back(my_rmp.info().trait_find(function_arguments[i]));
      else
        return 2;
    }

    return 0;
  }
}

int
AdjustedTraitCreator::calculateStats(const FunctionParser::TraitData& data)
{
  // Populate summary SampleInfo & class-based TraitStats's with values:
  //
  populateSampleInfos(data);

  if(my_adj_type == CLASS_BASED)
    if(performBinning(data) == false)
      return  1;

  return 0;
}

void
AdjustedTraitCreator::populateSampleInfos(const FunctionParser::TraitData& data)
{
  if(data.verbose)
    my_errors << priority(information) << "Calculating non-binned means & stdevs..." << endl;

  // Loop across the class trait list and stick the numbers onto the class-based TraitStat's:
  //
  for(TraitVector::iterator c = my_class_traits.begin(); c != my_class_traits.end(); ++c)
  {
    if(data.verbose)
      my_errors << priority(information) << my_rmp.info().trait_info(*c).name() << ": ";

    for(RPED::PedigreeIterator p_iter = my_rmp.pedigree_begin(); p_iter != my_rmp.pedigree_end(); ++p_iter)
    {
      for(RPED::MemberIterator m_iter = p_iter->member_begin(); m_iter != p_iter->member_end(); ++m_iter)
      {
        double reference_value = p_iter->info().trait(m_iter->index(), my_reference_trait);

        if(! SAGE::isnan(reference_value) && (p_iter->info().trait(m_iter->index(), *c) == 1))
        {
          if(data.verbose)
            my_errors << priority(information) << reference_value << " ";

          my_class_stats[*c].sinfo.add(reference_value);
        }
      }
    }
    
    // - Class must contain at least two members.
    //
    if(my_class_stats[*c].sinfo.count() < 2)
    {
      my_errors << flush;
      my_errors << priority(fatal) << "Attempting to create class adjusted variable with fewer than two data "
                << "members in a class.  Aborting program ..." << endl;
      exit(1);    
    }

    // - Divide by n-1 when calculating standard deviation.
    //
    my_class_stats[*c].sinfo.sample_adjustment(1);
  
    my_class_stats[*c].binned_mean  = my_class_stats[*c].sinfo.mean();
    my_class_stats[*c].binned_stdev = my_class_stats[*c].sinfo.standard_deviation();

    if(data.verbose)
      my_errors << priority(information) << "mean = "   << my_class_stats[*c].binned_mean
                                     << " stdev = " << my_class_stats[*c].binned_stdev << endl;
  }
  

  // Loop across all individuals and add valid values to the summary sampleinfo:
  //
  if(data.verbose)
    my_errors << priority(information) << "Summary stats: ";

  for(RPED::PedigreeIterator p_iter = my_rmp.pedigree_begin(); p_iter != my_rmp.pedigree_end(); ++p_iter)
  {
    for(RPED::MemberIterator m_iter = p_iter->member_begin(); m_iter != p_iter->member_end(); ++m_iter)
    {
      double reference_value = p_iter->info().trait(m_iter->index(), my_reference_trait);

      if(! SAGE::isnan(reference_value))
      {
        if(data.verbose)
          my_errors << priority(information) << reference_value << " ";

        my_summary_stats.add(reference_value);
      }
    }
  }

  if(data.verbose)
    my_errors << priority(information) << "mean = "   << my_summary_stats.mean()
                                   << " stdev = " << my_summary_stats.standard_deviation() << endl;
}

bool 
AdjustedTraitCreator::performBinning(const FunctionParser::TraitData& data)
{
  if(data.verbose)
    my_errors << priority(information) << "Calculating binned means & stdevs..." << endl;

  // If bin_size == 0, then no binning is necessary:
  //
  if(my_bin_size == 0.0)
  {
    if(data.verbose)
      my_errors << priority(information) << "No minimum bin size required" << endl;

    for(TraitVector::iterator c = my_class_traits.begin(); c != my_class_traits.end(); ++c)
    {
      my_class_stats[*c].binned_mean  = my_class_stats[*c].sinfo.mean();
      my_class_stats[*c].binned_stdev = my_class_stats[*c].sinfo.standard_deviation();
    }

    return true;
  }

  // Calculate the start/endpoints of the class-based proportions.
  //
  typedef std::map<size_t, std::pair<double, double> > ClassEndpointsMap;

  ClassEndpointsMap endpoints;

  double current_x = 0.0;

  if(data.verbose)
    my_errors << priority(information) << "Bin endpoints: ";

  for(TraitVector::iterator c = my_class_traits.begin(); c != my_class_traits.end(); ++c)
  {
    double proportion = (SAGE::isnan(my_class_stats[*c].sinfo.standard_deviation()) || SAGE::isnan(my_class_stats[*c].sinfo.mean())) ?
                        0.0 : (double)my_class_stats[*c].sinfo.count() / (double)my_bin_size;

    endpoints[*c].first  = current_x;
    endpoints[*c].second = current_x + proportion;

    if(data.verbose)
      my_errors << priority(information) << my_rmp.info().trait_info(*c).name() << ": "
                << endpoints[*c].first << " to " << endpoints[*c].second << " ";

    current_x += proportion;
  }

  if(data.verbose)
    my_errors << endl;

  double final_x = current_x;

  // Make sure there is at least 100% of a bin is present in the data altogether:
  //
  if(data.verbose)
    my_errors << priority(information) << "Total bin proportion? " << final_x << endl;

  if(final_x < 1.0)
    return false;

  // Loop through class traits, fleshing out the classes having insufficient data.
  //
  for(TraitVector::iterator c = my_class_traits.begin(); c != my_class_traits.end(); ++c)
  {
    // Skip this class if it has enough people:
    //
    if(my_class_stats[*c].sinfo.count() >= my_bin_size)
      continue;

    if(data.verbose)
      my_errors << priority(information) << "Binning " << my_rmp.info().trait_info(*c).name() 
                << " with " << my_class_stats[*c].sinfo.count() << " of " << my_bin_size 
                << " individuals" << endl;

    // Calculate the centerpoint and left/right bounds for this class's proportion:
    //
    double  center_point = (endpoints[*c].first + endpoints[*c].second) / 2;
    double  left_bound   = center_point - 0.5;
    double  right_bound  = center_point + 0.5;

    if(left_bound < 0)
    {
      right_bound += -left_bound;
      left_bound   =  0.0;
    }
    else if(right_bound > final_x)
    {
      left_bound  -= right_bound - final_x;
      right_bound  = final_x;
    }

    // Now loop through all classes, calculate their respective proportions,
    // and add the proportion of the statistic in question to the binned statistic:
    //
    my_class_stats[*c].binned_mean  = 0.0;
    my_class_stats[*c].binned_stdev = 0.0;
    for(TraitVector::iterator t = my_class_traits.begin(); t != my_class_traits.end(); ++t)
    {
      double  x1 = endpoints[*t].first;
      double  x2 = endpoints[*t].second;
      double  proportion = 0.0;

      enum EndPointStatus { LEFT_OF_REGION, WITHIN_REGION, RIGHT_OF_REGION };

      EndPointStatus x1_status = x1 < left_bound ? LEFT_OF_REGION : (x1 <= right_bound ? WITHIN_REGION : RIGHT_OF_REGION),
                     x2_status = x2 < left_bound ? LEFT_OF_REGION : (x2 <= right_bound ? WITHIN_REGION : RIGHT_OF_REGION);

      // Is it ENTIRELY to the left of the region?
      //
      if(x2_status == LEFT_OF_REGION)
        proportion = 0.0;

      // Is it PARTIALLY to the left of the region?
      //
      else if(x1_status == LEFT_OF_REGION && x2_status == WITHIN_REGION)
        proportion = x2 - left_bound;

      // Is it WHOLLY ENCLOSED within the region?
      //
      else if(x1_status == WITHIN_REGION && x2_status == WITHIN_REGION)
        proportion = x2 - x1;

      // Is it PARTIALLY to the right of the region?
      //
      else if(x1_status == WITHIN_REGION && x2_status == RIGHT_OF_REGION)
        proportion = right_bound - x1;

      // Is it ENTIRELY to the right of the region?
      //
      else if(x1_status == RIGHT_OF_REGION)
        proportion = 0.0;

      // Add the proportion of the two statistics:
      //
      if(data.verbose)
      {
        // Store precision so we don't muss the stream
        size_t p = my_errors.precision();
        
        my_errors << priority(information) 
                  << std::setprecision(1) << std::fixed << (proportion * 100.0) << "% of " 
                  << my_rmp.info().trait_info(*t).name() << ", ";

        my_errors << setprecision(p);
      }

      my_class_stats[*c].binned_mean  += proportion * my_class_stats[*t].sinfo.mean();
      my_class_stats[*c].binned_stdev += proportion * my_class_stats[*t].sinfo.standard_deviation();
    }

    if(data.verbose)
      my_errors << priority(information)
                << "MEAN changed from "   << my_class_stats[*c].sinfo.mean() << " to " << my_class_stats[*c].binned_mean 
                << " STDEV changed from " << my_class_stats[*c].sinfo.standard_deviation() << " to " << my_class_stats[*c].binned_stdev 
                << endl;

  } // End of class trait loop

  return  true;
}


} // End namespace FUNC
} // End namespace SAGE
