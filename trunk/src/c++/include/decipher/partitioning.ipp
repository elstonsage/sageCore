//============================================================================
// File:      partitioning.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   4/18/7 created        -djb
//                                                                          
// Notes:     Inline implementation of partitioner class.
//                                                                          
// Copyright (c) 2007 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


//============================================================================
// IMPLEMENTATION:  partitioner
//============================================================================
//
inline
partitioner::partitioner(const instructions& instr, 
                         const FilteredMultipedigree& filtered_mped, cerrorstream& errors)
      : my_instructions(instr), my_filtered_mped(filtered_mped), my_errors(errors)
{}

inline bool
partitioner::designated_by_user(member_iterator m_iter) const
{
  bool designated = false;

  if( my_instructions.family_rep != (size_t)(-1) )
  {
    if(    my_instructions.rep_field_type == instructions::CONTINUOUS_FIELD
        || my_instructions.rep_field_type == instructions::BINARY_FIELD )
    {
      double trait_value = m_iter->info().trait(my_instructions.family_rep);

      if( !SAGE::isnan(trait_value) )
        designated = (trait_value == my_instructions.family_rep_value.dbl);
    }
    else if( my_instructions.rep_field_type == instructions::STRING_FIELD )
    {
      string string_value = m_iter->info().get_string(my_instructions.family_rep);

      if( !string_value.empty() )
        designated = (string_value == my_instructions.family_rep_value.str);
    }
  }

  return  designated;
}

inline size_t
partitioner::get_genotyped_count(member_iterator m_iter, const locus_group& loci) const
{
  size_t genotype_cnt = 0;
  
  for( size_t m = 0; m < loci.size(); ++m )
  {
    size_t mi = loci[m].first;
    const MLOCUS::inheritance_model& minfo = *(loci[m].second);

    if( !m_iter->info().phenotype_missing(mi, minfo) )
      ++genotype_cnt;
  }

  return genotype_cnt;
}

inline bool  
partitioner::member_missing_partition_value(const member_iterator& m_iter, 
                                         const partition_data& partition) const
{
  assert(partition.valid());
  
  bool  missing = false;
  
  if(partition.type == instructions::CONTINUOUS_FIELD ||
     partition.type == instructions::BINARY_FIELD       )
  {
    missing = SAGE::isnan(m_iter->info().trait(partition.field_index));
  }
  else if(partition.type == instructions::STRING_FIELD)
  {
    missing = m_iter->info().get_string(partition.field_index).empty();
  }  
  
  return  missing;
}

inline pair< bool, pair<string, string> >
partitioner::member_belongs(member_iterator m_iter, const partition_vector& partitions) const
{
  assert(partitions.size() == 2);

  pair<bool, string> partition0 = member_belongs_to_partition(m_iter, partitions[0]);
  pair<bool, string> partition1 = member_belongs_to_partition(m_iter, partitions[1]);

  bool belong = partition0.first && partition1.first;

  pair<string, string> sub_pop_keys = make_pair(partition0.second, partition1.second);

  return make_pair(belong, sub_pop_keys);
}

inline pair<bool, string>
partitioner::member_belongs_to_partition(member_iterator m_iter,
                                      const partition_data& partition) const
{
  if(    partition.type == instructions::CONTINUOUS_FIELD
      || partition.type == instructions::BINARY_FIELD )
  {
    double trait_value = m_iter->info().trait(partition.field_index);

    return check_double_partition(trait_value, partition.sub_pops);
  }
  else if( partition.type == instructions::STRING_FIELD )
  {
    string string_value = m_iter->info().get_string(partition.field_index);

    return check_string_partition(string_value, partition.sub_pops);
  }

  return make_pair(true, NULL_STRING);
}

inline pair<bool, string>
partitioner::check_double_partition(double trait_value, const std::map<string, instructions::value>& sub_pops) const
{
  std::map<string, instructions::value>::const_iterator pi = sub_pops.begin();
  for( ; pi != sub_pops.end(); ++pi )
  {
    if( trait_value == pi->second.dbl )
      return make_pair(true, pi->first);
  }

  return make_pair(false, NULL_STRING);
}

inline pair<bool, string>
partitioner::check_string_partition(string string_value, const std::map<string, instructions::value>& sub_pops) const
{
  std::map<string, instructions::value>::const_iterator pi = sub_pops.begin();
  for( ; pi != sub_pops.end(); ++pi )
  {
    if( string_value == pi->second.str )
      return make_pair(true, pi->first);
  }

  return make_pair(false, NULL_STRING);
}





