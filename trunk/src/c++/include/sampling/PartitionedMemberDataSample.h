#ifndef PARTITIONED_MEMBER_DATA_SAMPLE_H
#define PARTITIONED_MEMBER_DATA_SAMPLE_H

#include "fped/fped.h"
#include "sampling/MemberDataSample.h"
#include "sampling/IndividualPartition.h"

namespace SAGE {
namespace SAMPLING {

/** \brief Combines the functionality of the MemberDataSample with the IndividualPartition
  *
  * \par Introduction
  *
  * Let's say that in addition to tracking trait data, you want to track classification-based
  * trait data. In that case, you can use the PartitionedMemberDataSample.
  *
  * \par Getting started
  *
  * To begin with, you use this class in the same way that you would use the
  * MemberDataSample:
  *
  * \code
  * sample.importField(p->getName(), "Susceptibility covariates", p->getName(), flags);
  * sample.finalizeData();
  * sample.finalizeUserCreatedData();
  * \endcode
  *
  * Having finished importing/finalizing your PartitionedMemberDataSample, we move on to
  * the partitioning stage. First, you must classify all the individuals in the set:
  *
  * Note that in the following example, classifyIndividual() is assumed to be a function that
  * you have already written. It should return an unsigned integer, indicating the classification
  * code of the individual.
  *
  * \code
  * for(int j = 0; j < sample.getTotalIndividualCount(); ++j)
  *   sample.addIndividual(j, classifyIndividual(sample.getIndividual(j)));
  * \endcode
  *
  * Now that you have classified all the individuals, you can generate the partitioned fields.
  * Basically, the ParititionedMemberDataSample will take the existing fields it stores, and for
  * each field it will generate k partitioned fields (where k is the number of classifications).
  * Please note that the adjustment applied to each partitioned field (mean, standardization) will
  * be based on the mean and standard deviation of the \b entire field, not those of the partitioned
  * group:
  *
  * \code
  * sample.generatePartitionedFields();
  * \endcode
  *
  * Now, in addition to accessing Field data via the getField() function, you can access partitioned
  * Field data via the getPartitionedField() function.
  */
class PartitionedMemberDataSample : public MemberDataSample, public IndividualPartition
{
public:
    PartitionedMemberDataSample(const FPED::Multipedigree & mp, cerrorstream & err = sage_cerr);
    PartitionedMemberDataSample(const PartitionedMemberDataSample & other);
    ~PartitionedMemberDataSample();
    
    void  reset();
  
  //@}

  /// @name Partitioned field generation
  //@{

    ///
    /// Generates the information for the partitioned fields.
    bool generatePartitionedFields();
  
  //@}
  
  /// @name Individual counts
  //@{
  
    ///
    /// Returns the number of valid individuals in the indicated partition.
    /// \param partition_idx The index number of the partition in question 
    int getPartitionValidIndividualCount (size_t partition_idx) const;
    
    ///
    /// Returns the number of invalid individuals in the indicated partition.
    /// \param partition_idx The index number of the partition in question 
    int getPartitionInvalidIndividualCount (size_t partition_idx) const;
  
  //@}

  /// @name Partitioned field extraction
  //@{
  
    ///
    /// Returns the partitioned Field.
    /// \param group_name The name of the group to which the Field belongs
    /// \param field_name The name of the field
    /// \param partition_idx The index number of the partition requested
    const Field & getPartitionedField(const string & group_name, const string & field_name, size_t partition_idx) const;

  //@}

private:

    PartitionedMemberDataSample& operator= (const PartitionedMemberDataSample &) { return *this; }

    ///
    /// \internal
    /// Updates the my_valid_counts and my_invalid_counts vectors.
    void updatePartitionCounts();

    vector<size_t> my_valid_counts;
    vector<size_t> my_invalid_counts;
    
    std::map<pair<string, string>, vector<Field> > my_partitioned_fields;
};

} // End namespace SAMPLING
} // End namespace SAGE

#endif
