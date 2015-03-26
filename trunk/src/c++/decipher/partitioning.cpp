//============================================================================
// File:      partitioning.cpp                      
//                                                                          
// Author:    Dan Baechle                                     
//                                                                          
// History:   4/18/7 - created.                         djb
//                                                                          
// Notes:     Implementation of analysis class.   
//               
//                                                                          
// Copyright (c) 2007 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "decipher/partitioning.h"

namespace SAGE
{

namespace DECIPHER
{

//============================================================================
// IMPLEMENTATION:  partitioner
//============================================================================
//
void
partitioner::get_genotype_cnt_sorted_members(std::map< size_t, vector< member_iterator >, greater<size_t> >& g_sorted_members,
                                          subpedigree_iterator sped, const locus_group& loci) const
{
  member_iterator iter = sped->member_begin();
  for( ; iter != sped->member_end(); ++iter )
  {
    size_t g_cnt = get_genotyped_count(iter, loci);

    if( designated_by_user(iter) )
    {
      if( g_cnt > 0 )
        g_sorted_members[g_cnt].push_back(iter);
      else
      {
        my_errors << priority(warning)
                  << "Individual " << iter->name() << " in pedigree " << iter->pedigree()->name()
                  << " is designated as family representative, but is not genotyped in the "
                  << "haplotyping region.  Ignoring designation ..." << endl;
      }
    }
    else if( my_instructions.family_rep == (size_t)(-1) && g_cnt > 0 )
    {
      g_sorted_members[g_cnt].push_back(iter);
    }
  } // end of member

  return;
}

pair<member_iterator, pair<string, string> >
partitioner::get_best_family_rep(subpedigree_iterator sped, const locus_group& loci,
                              const partition_vector& partitions) const
{
  std::map<size_t, vector<member_iterator>, greater<size_t> > g_sorted_members;

  instructions::value val;
  get_genotype_cnt_sorted_members(g_sorted_members, sped, loci);

#if 0
  cout << "g_sorted_members :" << endl;
  std::map<size_t, vector<member_iterator>, greater<size_t> >::const_iterator gi = g_sorted_members.begin();
  for( ; gi != g_sorted_members.end(); ++gi )
  {
    cout << "genotyped count = " << gi->first << endl;

    vector<member_iterator>::const_iterator m_iter = gi->second.begin();

    for( ;  m_iter != gi->second.end(); ++ m_iter )
    {
      cout << "\t" << (*m_iter)->pedigree()->name() << ":" << (*m_iter)->name();
    }
    cout << endl;
  }
  cout << endl;
#endif

  if( g_sorted_members.begin() != g_sorted_members.end() )
  {
    std::map<size_t, vector<member_iterator>, greater<size_t> >::const_iterator g_iter = g_sorted_members.begin();

    for( ; g_iter != g_sorted_members.end(); ++g_iter )
    {
      vector<member_iterator>::const_iterator m_iter = g_iter->second.begin();

      for( ;  m_iter != g_iter->second.end(); ++ m_iter )
      {
        pair< bool, pair<string, string> > m_belong = member_belongs(*m_iter, partitions);

        if( m_belong.first )
          return make_pair(*m_iter, m_belong.second);

      }
    }
  }

  return make_pair(sped->member_end(), make_pair(NULL_STRING, NULL_STRING));
}

// - Divide data into subpopulations per user instructions.  Covers unrelated,
//   family and pool cases.
//
void
partitioner::partition_population(sub_pop_directory& sub_pop_members, 
                               const locus_group& loci, const partition_vector& partitions) const
{
  init_sub_pop_directory(sub_pop_members, partitions);
  
  switch(my_instructions.analysis_unit)
  {  
    case instructions::FAMILY_REP:
    {
      pedigree_iterator ped = my_filtered_mped.pedigree_begin();
      for( ; ped != my_filtered_mped.pedigree_end(); ++ped )
      {
        subpedigree_iterator sped = ped->subpedigree_begin();
        for( ; sped != ped->subpedigree_end(); ++sped )
        {
          pair<member_iterator, pair<string, string> > best_rep 
                    = get_best_family_rep(sped, loci, partitions);

          if( best_rep.first != sped->member_end() )
          {
            string key0 = best_rep.second.first;
            string key1 = best_rep.second.second;

            sub_pop_members[key0][key1].push_back(&(*(best_rep.first)));
          }
          else
          {
            my_errors << priority(information)
                      << "No valid family representative found for constituent pedigree "
                      << "in pedigree " << sped->pedigree()->name() << " containing member " 
                      << sped->member_begin()->name()
                      << ".  Skipping this constituent pedigree ..." << endl;
          }
        } // end of subpedigree
        
        find_singletons(sub_pop_members, ped, loci, partitions);      
       
      } // end of pedigree

  #if 0
    cout << "Related case:" << endl;
    dump_sub_pops_map(sub_pop_members);
  #endif
  
      break;
    }

    case instructions::EACH_INDIVIDUAL:
    {
      pedigree_iterator ped = my_filtered_mped.pedigree_begin();
      for( ; ped != my_filtered_mped.pedigree_end(); ++ped )
      {
        member_iterator m_iter = ped->member_begin();
        for( ; m_iter != ped->member_end(); ++m_iter )
        {
          size_t g_cnt = get_genotyped_count(m_iter, loci);
          
          if( g_cnt < 1 )
          {
            /*  This error message could be an annoyance as often data has many untyped people
            my_errors << priority(information)
                      << "Individual " << m_iter->name() << " in pedigree " << m_iter->pedigree()->name()
                      << " is not genotyped in the haplotyping region.  Ignoring this individual ..." << endl;
            */

            continue;
          }

          pair< bool, pair<string, string> > m_belong = member_belongs(m_iter, partitions);
          
          if( m_belong.first )
          {
            string key0 = m_belong.second.first;
            string key1 = m_belong.second.second;

            sub_pop_members[key0][key1].push_back(&(*m_iter));
          }
        } // end of member
      }
  #if 0
    cout << "Unrelated case:" << endl;
    dump_sub_pops_map(sub_pop_members);
  #endif
  
      break;
    }
    
    case instructions::POOL:
    {
      pedigree_iterator ped = my_filtered_mped.pedigree_begin();
      for( ; ped != my_filtered_mped.pedigree_end(); ++ped )
      {
        member_iterator m_iter = ped->member_begin();
        for( ; m_iter != ped->member_end(); ++m_iter )
        {
          pair< bool, pair<string, string> > m_belong = member_belongs(m_iter, partitions);
          
          if( m_belong.first )
          {
            string key0 = m_belong.second.first;
            string key1 = m_belong.second.second;

            sub_pop_members[key0][key1].push_back(&(*m_iter));
          }
        } // end of member
      }
  #if 0
    cout << "Pool case:" << endl;
    dump_sub_pops_map(sub_pop_members);
  #endif
  
      break;
    }    
    
    default:
      assert(false);
  }
}

template<typename M> void
partitioner::init_sub_pop_directory(std::map<string, std::map<string, M> >& sub_pop_members,
                                 const partition_vector& partitions) const
{
  bool  outer_partition_valid = partitions[0].valid();
  bool  inner_partition_valid = partitions[1].valid();
  
  M  mems;
  
  if(! inner_partition_valid)
  {
    assert(! outer_partition_valid);
    
    sub_pop_members[""][""] = mems;    // No subpopulations specified.
  }
  else if(! outer_partition_valid)
  {
    std::map<string, instructions::value>::const_iterator s1 = partitions[1].sub_pops.begin();
    for( ; s1 != partitions[1].sub_pops.end(); ++ s1 )
    {
      string key1 = s1->first;

      sub_pop_members[""][key1] = mems;
    }  
  }
  else
  {
    std::map<string, instructions::value>::const_iterator s0 = partitions[0].sub_pops.begin();
    for( ; s0 != partitions[0].sub_pops.end(); ++ s0 )  
    {
      string key0 = s0->first;

      std::map<string, instructions::value>::const_iterator s1 = partitions[1].sub_pops.begin();
      for( ; s1 != partitions[1].sub_pops.end(); ++ s1 )
      {
        string key1 = s1->first;

        sub_pop_members[key0][key1] = mems;
      }
    }
  }
}

// - Add all singletons in the given pedigree to the directory.
//
void  
partitioner::find_singletons(sub_pop_directory& sub_pop_members, 
                          pedigree_iterator& ped, const locus_group& loci,
                          const partition_vector& partitions) const
{
  member_iterator  m_iter     = ped->unconnected_begin();
  member_iterator  m_end_iter = ped->unconnected_end();
  for(; m_iter != m_end_iter; ++m_iter)
  {
    size_t g_cnt = get_genotyped_count(m_iter, loci);

    if( g_cnt < 1 )
    {
      /*  This error message could be an annoyance as often data has many untyped people
      my_errors << priority(information)
                << "Individual " << m_iter->name() << " in pedigree " << m_iter->pedigree()->name()
                << " is not genotyped in the haplotyping region.  Ignoring this individual ..." << endl;
      */

      continue;
    }

    pair< bool, pair<string, string> > m_belong = member_belongs(m_iter, partitions);

    if( m_belong.first )
    {
      string key0 = m_belong.second.first;
      string key1 = m_belong.second.second;

      sub_pop_members[key0][key1].push_back(&(*m_iter));
    }        
  }  
}


// - Add all singletons in the given pedigree to the directory.
//
void  
partitioner::find_singletons(founders_sub_pop_directory& sub_pop_members, 
                          pedigree_iterator& ped, const locus_group& loci,
                          const partition_vector& partitions) const
{
  member_iterator  m_iter     = ped->unconnected_begin();
  member_iterator  m_end_iter = ped->unconnected_end();
  for(; m_iter != m_end_iter; ++m_iter)
  {
    size_t g_cnt = get_genotyped_count(m_iter, loci);

    if( g_cnt < 1 )
    {
      /*  This error message could be an annoyance as often data has many untyped people
      my_errors << priority(information)
                << "Individual " << m_iter->name() << " in pedigree " << m_iter->pedigree()->name()
                << " is not genotyped in the haplotyping region.  Ignoring this individual ..." << endl;
      */

      continue;
    }

    pair< bool, pair<string, string> > m_belong = member_belongs(m_iter, partitions);

    if( m_belong.first )
    {
      string key0 = m_belong.second.first;
      string key1 = m_belong.second.second;

      sub_pop_members[key0][key1].second.push_back(&(*m_iter));
    }        
  }  
}

// - Do all constituent pedigree founders belong to the same subpopulation?
//   If so, return true and the subpopulation names.  If a founder has missing
//   or contradictary partition information, constituent pedigree is inconsistent 
//   and false with null strings for names is returned.
//
pair<bool, pair<string, string> >
partitioner::partition_data_consistent(subpedigree_iterator& sped, 
                                    const partition_vector& partitions) const
{
  bool    consistent = false;
  string  outer_name = "";
  string  inner_name = "";

  size_t  partition_count = valid_partition_count(partitions);
  
  if(partition_count == 0)
  {
    consistent = true;
  }
  else
  {
    member_iterator  m_iter     = sped->member_begin();
    member_iterator  m_end_iter = sped->member_end();
    for(; m_iter != m_end_iter; m_iter++)
    {
      if(m_iter->parent1() != 0)     // Consider founders only.
      {
        continue;
      }
      
      if(partition_count == 1)
      {
        if(member_missing_partition_value(m_iter, partitions[1]))
        {
          continue;
        }
      
        pair<bool, string>  member_sub_pop_info1 = member_belongs_to_partition(m_iter, partitions[1]);
        if(member_sub_pop_info1.first)
        {
          if(inner_name.empty())
          {
            inner_name = member_sub_pop_info1.second;
            consistent = true;
          }
          else if(member_sub_pop_info1.second != inner_name)
          {
            consistent = false;
            outer_name = "";
            inner_name = "";
            break;
          }
        }
        else
        {
          consistent = false;
          outer_name = "";
          inner_name = "";
          break;          
        }      
      }
      else if(partition_count == 2)
      {
        bool  missing0 = member_missing_partition_value(m_iter, partitions[0]);
        bool  missing1 = member_missing_partition_value(m_iter, partitions[1]);
        
        if(missing0 && missing1)
        {
          continue;
        }

        if(! (missing0 || missing1))
        {
          pair<bool, string>  member_sub_pop_info0 = member_belongs_to_partition(m_iter, partitions[0]);                
          pair<bool, string>  member_sub_pop_info1 = member_belongs_to_partition(m_iter, partitions[1]);        
        
          if(member_sub_pop_info0.first && member_sub_pop_info1.first)
          {
            if(outer_name.empty())
            {
              outer_name = member_sub_pop_info0.second;
              if(! inner_name.empty())
              {
                consistent = true;
              }
            }
            else if(member_sub_pop_info0.second != outer_name)
            {
              consistent = false;
              outer_name = "";
              inner_name = "";
              break;
            }                        
          
            if(inner_name.empty())
            {
              inner_name = member_sub_pop_info1.second;
              if(! outer_name.empty())
              {
                consistent = true;
              }
            }
            else if(member_sub_pop_info1.second != inner_name)
            {
              consistent = false;
              outer_name = "";
              inner_name = "";
              break;
            }
          }
          else
          {
            consistent = false;
            outer_name = "";
            inner_name = "";
            break;
          }      
        }
        else if(missing0)
        {
          pair<bool, string>  member_sub_pop_info1 = member_belongs_to_partition(m_iter, partitions[1]);        
        
          if(member_sub_pop_info1.first)
          {
            if(inner_name.empty())
            {
              inner_name = member_sub_pop_info1.second;
              if(! outer_name.empty())
              {
                consistent = true;
              }
            }
            else if(member_sub_pop_info1.second != inner_name)
            {
              consistent = false;
              outer_name = "";
              inner_name = "";
              break;
            }
          }
          else
          {
            consistent = false;
            outer_name = "";
            inner_name = "";
            break;
          }        
        }
        else
        {
          pair<bool, string>  member_sub_pop_info0 = member_belongs_to_partition(m_iter, partitions[0]);        
        
          if(member_sub_pop_info0.first)
          {     
            if(outer_name.empty())
            {
              outer_name = member_sub_pop_info0.second;
              if(! inner_name.empty())
              {
                consistent = true;
              }
            }
            else if(member_sub_pop_info0.second != inner_name)
            {
              consistent = false;
              outer_name = "";
              inner_name = "";
              break;
            }
          }
          else
          {
            consistent = false;
            outer_name = "";
            inner_name = "";
            break;
          }                
        }
      }
      else
      {
        assert(false);
      }
    }
  }
  
  return  make_pair(consistent, make_pair(outer_name, inner_name));
}

// - Return a founder with at least some genotype information.  If none are
//   found, return a null pointer.
//
member  
partitioner::get_informative_founder(subpedigree_iterator& sped, const locus_group& loci) const
{
  member_iterator  m_iter     = sped->member_begin();
  member_iterator  m_end_iter = sped->member_end();
  for(; m_iter != m_end_iter; ++m_iter)
  {
    if(MPED::mp_utilities::is_founder(m_iter) && get_genotyped_count(m_iter, loci) > 0)     
    {
      return  &(*m_iter);
    }
  }
  
  return  0;
}

// - Divide data into subpopulations per user instructions for founder pool
//   case.
//
void      
partitioner::partition_population(founders_sub_pop_directory& sub_pop_members, 
                               const locus_group& loci, const partition_vector& partitions) const
{
  init_sub_pop_directory(sub_pop_members, partitions);
  
  pedigree_iterator ped = my_filtered_mped.pedigree_begin();
  for( ; ped != my_filtered_mped.pedigree_end(); ++ped )
  {
    subpedigree_iterator sped = ped->subpedigree_begin();
    for( ; sped != ped->subpedigree_end(); ++sped )
    {
      pair<bool, pair<string, string> >  sub_pop_info = partition_data_consistent(sped, partitions);
      if(sub_pop_info.first)
      {
        string  key0 = sub_pop_info.second.first;
        string  key1 = sub_pop_info.second.second;
        member  m = get_informative_founder(sped, loci);
        if(m)
        {
          sub_pop_members[key0][key1].first.push_back(m);
          continue;
        }
        else
        {
          my_errors << priority(information)
                    << "Constituent pedigree "
                    << "in pedigree " << sped->pedigree()->name() << " containing member " 
                    << sped->member_begin()->name() << " contains no genotyped founders.  "
                    << "Looking for a family representative ..." << endl;        
        }
      }
      else
      {
        my_errors << priority(information)
                  << "Partition information inconsistent for constituent pedigree "
                  << "in pedigree " << sped->pedigree()->name() << " containing member " 
                  << sped->member_begin()->name()
                  << ".  Looking for a family representative ..." << endl;
      }
        
      pair<member_iterator, pair<string, string> > best_rep 
                = get_best_family_rep(sped, loci, partitions);

      if( best_rep.first != sped->member_end() )
      {
        string key0 = best_rep.second.first;
        string key1 = best_rep.second.second;

        sub_pop_members[key0][key1].second.push_back(&(*(best_rep.first)));
      }
      else
      {
        my_errors << priority(information)
                  << "No valid family representative found for constituent pedigree "
                  << "in pedigree " << sped->pedigree()->name() << " containing member " 
                  << sped->member_begin()->name()
                  << ".  Skipping this constituent pedigree ..." << endl;
      }
    } 

    find_singletons(sub_pop_members, ped, loci, partitions);
  }
  
#if 0
  cout << "Founder case :" << endl;
  dump_sub_pops_map(sub_pop_members);
#endif
  
}

void
partitioner::dump_sub_pops_map(const sub_pop_directory& sub_pops) const
{
  sub_pop_directory::const_iterator si1 = sub_pops.begin();
  for( ; si1 != sub_pops.end(); ++si1 )
  {
    std::map<string, vector<member> >::const_iterator si2 = si1->second.begin();
    for( ; si2 != si1->second.end(); ++si2 )
    {
      cout << "Sub pop : '" << si1->first << "', '" << si2->first << "'" << endl;

      vector<member>::const_iterator mi = si2->second.begin();
      for( ; mi != si2->second.end(); ++mi )
        cout << "\t" << (*mi)->pedigree()->name() << ":" << (*mi)->name();
      cout << endl;
    }
  }
}

// - For founders.
//
void
partitioner::dump_sub_pops_map(const founders_sub_pop_directory& sub_pops) const
{
  founders_sub_pop_directory::const_iterator si1 = sub_pops.begin();
  for( ; si1 != sub_pops.end(); ++si1 )
  {
    std::map<string, pair<vector<member>, vector<member> > >::const_iterator si2 = si1->second.begin();
    for( ; si2 != si1->second.end(); ++si2 )
    {
      cout << "\nSub pop : '" << si1->first << "', '" << si2->first << "'" << endl;

      cout << "Use founders:" << endl;
      vector<member>::const_iterator mi = si2->second.first.begin();
      for( ; mi != si2->second.first.end(); ++mi )
        cout << "\t" << (*mi)->pedigree()->name() << ":" << (*mi)->name();
      cout << "\n" << endl;
      
      cout << "Reps and singletons:" << endl;
      mi = si2->second.second.begin();
      for( ; mi != si2->second.second.end(); ++mi )
        cout << "\t" << (*mi)->pedigree()->name() << ":" << (*mi)->name();
      cout << endl;
    }
  }
}

}
}

