#include "error/internal_error.h"
#include <iostream>
#include "sampling/Field.h"

namespace SAGE {
namespace SAMPLING {

//=====================================================
//
//  CONSTRUCTOR #1
//
//=====================================================
Field::Field()
{
  my_orig_name   = "";
  my_field_name  = "";
  my_group_name  = "";
  my_field_id    = FieldID();
  my_group_id    = GroupID();
  my_input_flags = 0;
  
  // Clear sample info and flags
  my_sample_info        . clear();
  my_adj_value_presents . clear();

  my_valid           = false;
  my_finalized       = false;
  my_user_defined    = false;
}

Field::Field(size_t size)
{
  my_orig_name  = "";
  my_field_name = "";
  my_group_name = "";
  my_field_id   = FieldID();
  my_group_id   = GroupID();
  my_input_flags = 0;

  my_adj_value_presents.resize(size, std::numeric_limits<double>::quiet_NaN());
  
  initializeForImport(size);
}

void 
Field::initializeForImport(size_t size)
{
  // Clear value vectors using the swap trick to minimize memory use
  std::vector<double> temp(size, std::numeric_limits<double>::quiet_NaN());
  
  my_orig_values     . swap( temp );

  temp.resize(size, std::numeric_limits<double>::quiet_NaN());

  my_adj_values      . swap( temp );
  
  // Clear sample info and flags
  my_sample_info     . clear();

  my_valid           = false;
  my_finalized       = false;
  my_user_defined    = false;
}

//=====================================================
//
//  COPY CONSTRUCTOR
//
//=====================================================
Field::Field(const Field & other)
{
  copy(other);
}

//=====================================================
//
//  operator=
//
//=====================================================
Field& 
Field::operator=(const Field & other)
{
  if(this != &other)
    copy(other);

  return *this;
}

//=====================================================
//
//  copy
//
//=====================================================
void
Field::copy(const Field & other)
{
  my_adj_values         = other.my_adj_values;
  my_adj_value_presents = other.my_adj_value_presents;
  my_field_id           = other.my_field_id;
  my_field_name         = other.my_field_name;
  my_finalized          = other.my_finalized;
  my_group_id           = other.my_group_id;
  my_group_name         = other.my_group_name;
  my_input_flags        = other.my_input_flags;
  my_orig_name          = other.my_orig_name;
  my_orig_values        = other.my_orig_values;
  my_replaced_values    = other.my_replaced_values;
  my_sample_info        = other.my_sample_info;
  my_user_defined       = other.my_user_defined;
  my_valid              = other.my_valid;  
}

//=====================================================
//
//  setOrigValue
//
//=====================================================
int
Field::setOrigValue(size_t i, double val)
{
  // Check finalization status:
  if(my_finalized)
    return 1;

  // Change value:
  my_orig_values[i] = val;

  // Return success:
  return 0;
}

//=====================================================
//
//  setAdjValue
//
//=====================================================
int
Field::setAdjValue(size_t i, double val)
{
  // Change value:
  my_adj_values[i] = val;

  // Return success:
  return 0;
}

//=====================================================
//
//  finalizeData
//
//=====================================================
int
Field::finalizeData()
{
  // Update statistics (mean, variance, etc.) for original values:
  updateStats();

  // Copy original values to replaced values vector:
  my_replaced_values = my_orig_values;

  // Replace the missing values if requested:
  if(getInputFlags() & Field::ALLOW_AVERAGING)
    replaceMissing();

  // Copy replaced values to adjusted values:
  my_adj_values = my_replaced_values;

  // Adjust the adjusted values vector as necessary:
  if(getInputFlags() & Field::MEAN_ADJUST)
    meanAdjustValues();

  if(getInputFlags() & Field::STDEV_ADJUST)
    stdevAdjustValues();

  // Calculate the presence attribute for all individuals:
  for(size_t i = 0; i < my_adj_value_presents.size(); ++i)
    my_adj_value_presents[i] = !SAGE::isnan(getAdjValue(i));

  // Set the finalized status:
  my_finalized = true;

  // Return success:
  return 0;
}

//=====================================================
//
//  updateStats()
//
//=====================================================
int
Field::updateStats()
{
  // 0. Check my_finalized:

	if(my_finalized)
	  return 1;

  // 1. Check for valid population size:

        if(my_orig_values.size() == 0)
        {
          my_valid = false;
          return 0;
        }

  // 1. Set validity to false:

	my_valid = false;

  // 2. Create sampleinfo object:

	for(std::vector<double>::const_iterator val = my_orig_values.begin(); val != my_orig_values.end(); val++)
	  if(!SAGE::isnan(*val) && finite(*val))
	    my_sample_info.add(*val);

  // 3. Check population_size and sample_size before calculating descriptive statistics:

	if(my_sample_info.count() == 0)
	  my_valid = false;
	else
	  my_valid = true;

  // 4. Return success:

	return 0;
}

//=====================================================
//
//  replaceMissing
//
//=====================================================
int 
Field::replaceMissing()
{
  if(my_finalized)
    return 1;

  for(std::vector<double>::iterator val_itr = my_replaced_values.begin(); val_itr != my_replaced_values.end(); ++val_itr)
    if(SAGE::isnan(*val_itr) || !finite(*val_itr))
      *val_itr = getMean();

  return 0;
}
     
//=====================================================
//
//  meanAdjustValues
//
//=====================================================
int 
Field::meanAdjustValues()
{
  if(my_finalized)
    return 1;

  for(std::vector<double>::iterator val_itr = my_adj_values.begin(); val_itr != my_adj_values.end(); ++val_itr)
    if(!SAGE::isnan(*val_itr) && finite(*val_itr))
      *val_itr -= getMean();

  return 0;
}

//=====================================================
//
//  stdevAdjustValues
//
//=====================================================
int 
Field::stdevAdjustValues()
{
  if(my_finalized)
    return 1;

  for(std::vector<double>::iterator val_itr = my_adj_values.begin(); val_itr != my_adj_values.end(); ++val_itr)
    if(!SAGE::isnan(*val_itr) && finite(*val_itr))
      *val_itr /= getStdev();

  return 0;
}

} // End namespace SAMPLING
} // End namespace SAGE
