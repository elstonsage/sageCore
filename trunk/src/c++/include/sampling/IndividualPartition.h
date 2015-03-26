#ifndef INDIVIDUALPARTITION_H
#define INDIVIDUALPARTITION_H

#include <vector>
#include <sys/types.h>
#include <iostream>

namespace SAGE {
namespace SAMPLING {

/** \brief Stores information on which individuals belong to which partitions (disjoint classes)
  *
  * \par Introduction
  *
  * The IndividualPartition is a mechanism for tracking the classification number of a set
  * of individuals. The basic idea is that you have a number of partitions (n > 0), and a set
  * of individuals. Each individual belongs to one and only one partition.
  */
class IndividualPartition
{
public:
    IndividualPartition();
    IndividualPartition(const IndividualPartition& other);
    ~IndividualPartition();
    IndividualPartition& operator=(const IndividualPartition& other);

  //@}

  /// @name Data manipulation
  //@{

      ///
      /// To speed up processing, you can set the number of individuals in advance.
      /// \param individual_count The total number of individuals
      void setIndividualCount(size_t individual_count);

      ///
      /// To speed up processing, you can set the number of partitions in advance.
      /// \param partition_count The number of partitions
      void setPartitionCount(size_t partition_count);

      ///
      /// Adds an entry for the numbered individual, belonging to the indicated partition
      /// \param individual_idx The index number of the individual
      /// \param partition_idx The index number of the partition to which the individual belongs
      /// \retval true Individual successfully added
      /// \retval false Individual \b not successfully added
      bool addIndividual(size_t individual_idx, size_t partition_idx);
      
      void  clear();

  //@}
  
  /// @name Data extraction
  //@{
  
      ///
      /// Returns the total number of individuals.
      size_t getIndividualCount() const;

      ///
      /// Returns the number of individuals in the indicated partition.
      /// \param partition_idx The index number of the partition in question
      size_t getIndividualCount(size_t partition_idx) const;

      ///
      /// Returns the number of partitions.
      size_t getPartitionCount() const;

      ///
      /// Returns \c true if the indicated partition has no members, \c false otherwise.
      /// \param partition_idx The index number of the partition in question
      bool isEmpty(size_t partition_idx) const;
      
      ///
      /// Returns the partition index number of the indicated individual.
      /// \param individual_idx The index number of the individual in question
      size_t getPartition(size_t individual_idx) const;
      
  //@}
    
  /// @name Individual traversal
  //@{
  
      ///
      /// Returns a const begin iterator for the numbered partition.
      /// \param partition_idx The index number of the partition in question
      std::vector<size_t>::const_iterator getIndividualBegin(size_t partition_idx) const;

      ///
      /// Returns a const end iterator for the numbered partition.
      /// \param partition_idx The index number of the partition in question
      std::vector<size_t>::const_iterator getIndividualEnd(size_t partition_idx) const;

  //@}

  /// @name Debugging
  //@{
  
    ///
    /// \internal
    /// Dumps a list of the members of each group.
    void dumpPartitions() const;
    
  //@}

private:
    void copy(const IndividualPartition& other);

    // Data members
    std::vector<size_t>  my_individuals;               // Index = individual_idx, value = partition_idx
    std::vector<std::vector<size_t> >  my_partitions;      // Top-level index = partition_idx, 
                                                           // bottom-level index = member of partition, 
                                                           // value = individual_idx
  
};

//=====================
//  INLINE FUNCTIONS
//=====================

inline
IndividualPartition::IndividualPartition()
{
  my_individuals.clear();
  my_partitions.clear();
}
    
inline
IndividualPartition::IndividualPartition(const IndividualPartition & other)
{
  copy(other);
}
    
inline
IndividualPartition::~IndividualPartition()
{
//cout << "IndividualPartition::~IndividualPartition()" << endl;
}

inline IndividualPartition & 
IndividualPartition::operator= (const IndividualPartition & other)
{
  copy(other);
  
  return *this;
}

inline void 
IndividualPartition::setIndividualCount(size_t individual_count)
{
  my_individuals.resize(individual_count);
}

inline void 
IndividualPartition::setPartitionCount(size_t partition_count)
{
  my_partitions.resize(partition_count);
}

inline bool 
IndividualPartition::addIndividual(size_t individual_idx, size_t partition_idx)
{
  // Resize the individual vector if necessary:

  if(individual_idx >= my_individuals.size())
    my_individuals.resize(individual_idx + 1);
    
  // Resize the partitions vector if necessary:

  if(partition_idx >= my_partitions.size())
    my_partitions.resize(partition_idx + 1);

  // Stick data in:
  
  my_individuals[individual_idx] = partition_idx;
  my_partitions[partition_idx].push_back(individual_idx);

  return true;
}

inline void
IndividualPartition::clear()
{
  my_individuals.clear();
  my_partitions.clear();
}

inline size_t 
IndividualPartition::getIndividualCount() const
{
  return my_individuals.size();
}

inline size_t 
IndividualPartition::getIndividualCount(size_t partition_idx) const
{
  return my_partitions[partition_idx].size();
}

inline size_t
IndividualPartition::getPartitionCount() const
{
  return my_partitions.size();
}

inline bool 
IndividualPartition::isEmpty(size_t partition_idx) const
{
  return my_partitions[partition_idx].size() != 0;
}

inline size_t 
IndividualPartition::getPartition(size_t individual_idx) const
{
  return my_individuals[individual_idx];
}
      
inline void 
IndividualPartition::copy (const IndividualPartition & other)
{
  my_individuals = other.my_individuals;
  my_partitions  = other.my_partitions;
}

inline std::vector<size_t>::const_iterator 
IndividualPartition::getIndividualBegin(size_t partition_idx) const
{
  return my_partitions[partition_idx].begin();
}

inline std::vector<size_t>::const_iterator 
IndividualPartition::getIndividualEnd(size_t partition_idx) const
{
  return my_partitions[partition_idx].end();
}

inline void 
IndividualPartition::dumpPartitions() const
{
  for(size_t i = 0; i < my_partitions.size(); ++i)
  {
    cout << "Partition #" << i << ":" << endl;

    for(vector<size_t>::const_iterator j = getIndividualBegin(i); j != getIndividualEnd(i); ++j)
      cout << *j << " ";

    cout << endl << endl;
  }
}
    
} // End namespace SAMPLING
} // End namespace SAGE

#endif
