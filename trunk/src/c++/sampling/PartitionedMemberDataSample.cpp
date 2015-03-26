#include "sampling/PartitionedMemberDataSample.h"

namespace SAGE {
namespace SAMPLING {

//================================================================
//  CONSTRUCTOR
//================================================================
PartitionedMemberDataSample::PartitionedMemberDataSample(const FPED::Multipedigree& rmp, 
                                                                        cerrorstream& err) 
    : MemberDataSample(rmp, err), IndividualPartition()
{}

//================================================================
//  COPY CONSTRUCTOR
//================================================================
PartitionedMemberDataSample::PartitionedMemberDataSample(const PartitionedMemberDataSample& other) 
    : MemberDataSample(other), IndividualPartition(other)
{
  my_valid_counts       = other.my_valid_counts;
  my_invalid_counts     = other.my_invalid_counts;
  my_partitioned_fields = other.my_partitioned_fields;
}

//================================================================
//  Destructor
//================================================================
PartitionedMemberDataSample::~PartitionedMemberDataSample()
{}

//================================================================
//  reset()
//================================================================
void
PartitionedMemberDataSample::reset()
{
  IndividualPartition::clear();
  MemberDataSample::reset();
  
  my_valid_counts.clear();
  my_invalid_counts.clear();
  my_partitioned_fields.clear();
}
    
//================================================================
//  generatePartitionedFields()
//================================================================
bool 
PartitionedMemberDataSample::generatePartitionedFields()
{
  // Clear the partitioned fields
  my_partitioned_fields.clear();

  // Create partitioned field entries

  for(FieldConstIterator field_itr = getFieldBegin(); field_itr != getFieldEnd(); ++field_itr)
  {
    // Create an entry for this field in the my_partitioned_fields map and keep
    // a reference to the field vector
    pair<string, string> composite_name = std::make_pair(field_itr->getGroupName(), field_itr->getFieldName());

    vector<Field>& field_vector = my_partitioned_fields[composite_name];
    
    field_vector.resize(getPartitionCount());

    // Iterate through each partition and create a sub-field with the members
    // of it with this field's data
    for(size_t partition_id = 0; partition_id < getPartitionCount(); ++partition_id)
    {
      // Create the sub-Field by copying the complete one
      field_vector[partition_id] = *field_itr;
      
      // Re-initialize the sub-Field (clear the data) to get ready for the new values
      field_vector[partition_id].initializeForImport(getIndividualCount(partition_id));

      // Import the new values from the complete Field
      size_t sub_index = 0;

      for(vector<size_t>::const_iterator i  = getIndividualBegin (partition_id);
                                         i != getIndividualEnd   (partition_id); ++i, ++sub_index)
      {
        field_vector[partition_id].setOrigValue (sub_index, field_itr->getOrigValue (*i));
        field_vector[partition_id].setAdjValue  (sub_index, field_itr->getAdjValue  (*i));
      }

      // Update our statistics based upon the new imported data
      field_vector[partition_id].updateStats();

    } // End of partition loop

  } // End of field loop

  // Update partition counts:
  updatePartitionCounts();

  // Return success
  return true;
}
  
//==============================================================
//  getPartitionedField()
//==============================================================
const Field& 
PartitionedMemberDataSample::getPartitionedField(const string& group_name, 
                                                 const string& field_name, size_t partition_idx) const
{
  return my_partitioned_fields.find(std::make_pair(group_name, field_name))->second[partition_idx];
}

//=================================================================
//  getPartitionValidIndividualCount()
//=================================================================
int PartitionedMemberDataSample::getPartitionValidIndividualCount(size_t partition_idx) const
{ 
  return  my_valid_counts[partition_idx]; 
}

//=================================================================
//  getPartitionInvalidIndividualCount()
//=================================================================

int PartitionedMemberDataSample::getPartitionInvalidIndividualCount(size_t partition_idx) const
{ 
  return  my_invalid_counts[partition_idx]; 
}

//================================================================
//  updatePartitionCounts()
//================================================================
void
PartitionedMemberDataSample::updatePartitionCounts()
{
  my_valid_counts.resize(getPartitionCount(), 0);
  my_invalid_counts.resize(getPartitionCount(), 0);

  for(size_t i = 0; i < getTotalIndividualCount(); ++i)
  {
    if(isValid(i)) 
    {
      my_valid_counts[getPartition(i)]++;
    }
    else           
    {
      my_invalid_counts[getPartition(i)]++;
    }
  }
}

} // End namespace SAMPLING
} // End namespace SAGE
