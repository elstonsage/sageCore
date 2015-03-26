//****************************************************************************
//* File:      pairsetdata.cpp                                               *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0                                                   *
//*                                                                          *
//* Notes:     This source file stores the present inter or intra classes'   *
//*            pairset in multipedigrees and also keeps correlation and      *
//*            standard errors for each pairset.                             *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "fcor/setweight.h"
#include "fcor/pairsetcountinfo.h"
#include "fcor/pairsetdata.h"

using namespace std;

namespace SAGE {
namespace FCOR {

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     PairSetData                                                  ~
// ~                                                                         ~
// ~ Purpose:   Store present inter or intra class pair sets in pedigrees.   ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// ---------------------------------------------------------------------------
// Out-of-Line Implementation of PairSetData
// ---------------------------------------------------------------------------

PairSetData::PairSetData()
{
  my_parser     = NULL;
  my_main_built = false;
}

void
PairSetData::build_pairset(const FcorParser& fp)
{
  RelationMatrix  relmatrix(fp.get_multi_pedigree());

  build_pairset(relmatrix, fp);
}

void
PairSetData::build_pairset(const RelationMatrix& relmatrix,
                           const FcorParser&     fp)
{
  my_parser = &fp;

  sub_pairset_map sp_map;
  sub_count_map   sc_map;

  const vector<name_index_type>& traits = my_parser->get_trait_list();

  // build pair count map, per each pedigree.
  //
  for( size_t k = 0; k < my_parser->get_pedigree_count(); ++k )
  {
    const RPED::RefPedigree& sp = my_parser->get_multi_pedigree()->pedigree_index(k);

    size_t rank = sp.member_count();

    const RPED::RefPedInfo& trait_info = sp.info();

    vector<bool> informative(rank, false);

    // member with a valid trait value become informative.
    //
    for(size_t i = 0; i < rank; ++i)
      for(size_t t = 0; t < traits.size(); ++t)
        if( finite( trait_info.trait(i, traits[t].second) ) )
        {
          informative[i] = true;
          break;
        }

    // make pairs only with informative members.
    //
    for( size_t i = 0; i < rank ; ++i )
    {
      if(!informative[i])
        continue;

      for( size_t j = 0; j <= i; ++j )
      {
        if(!informative[j])
          continue;

        const_pedigree_member_pair m_pair( &sp.member_index(i), &sp.member_index(j) );

        CompleteRelationType c = relmatrix.get_complete_relationship(m_pair);

        c.normalize(m_pair, gender_less());

        if( sp_map.find(c) == sp_map.end() )
        {
          sp_map[c].resize(my_parser->get_pedigree_count());
          sc_map[c] = 0;
        }

        sp_map[c][k].push_back(WeightedMemberPair(m_pair, WEIGHT_COUNT));
        sc_map[c] = sc_map[c] + 1;

        if( is_intraclass(c) && !is_selfclass(c) )
        {
          const_pedigree_member_pair rp( m_pair.second, m_pair.first );
          sp_map[c][k].push_back(WeightedMemberPair(rp, WEIGHT_COUNT));
        }
      }
    }
  }

  // save pairs.
  sub_pairset_const_iterator sp_begin = sp_map.begin();
  sub_pairset_const_iterator sp_end   = sp_map.end();

  PairSetCountInfo count_info;

  for( ; sp_begin != sp_end; ++sp_begin )
  {
    //if( sc_map[sp_begin->first] < 3 )
    //  continue;

    sc_map[sp_begin->first] = my_subtype_pairset.size();

    my_subtype_pairset.push_back(sp_begin->second);
    my_subtypes.push_back(sp_begin->first);

    pairset_info p_info;
    p_info.name  = long_relationship_name(sp_begin->first.relationship);
    p_info.gname = long_relationship_gender_name(sp_begin->first);
    p_info.rcode = gender_lists(sp_begin->first);
    p_info.distance1 = sp_begin->first.relationship.distance1();
    p_info.distance2 = sp_begin->first.relationship.distance2();

    if( is_selfclass(sp_begin->first) )
      p_info.type = SELF;
    else if( is_unrelated(sp_begin->first) )
      p_info.type = UNRELATED;
    else if( is_mate(sp_begin->first) )
      p_info.type = MATE;
    else if( is_intraclass(sp_begin->first) )
      p_info.type = INTRA;
    else
      p_info.type = OTHER;

    count_info.build_member_info(sp_begin->second, (p_info.type == INTRA));

    p_info.total_pair_count = count_info.get_total_pair_count();
    p_info.distinctive_pair_count = count_info.get_distinctive_pair_count();

    my_subtype_info.push_back(p_info);

#if 0
  cout << p_info.name << " ";

  if( p_info.type == SELF )
    cout << "self" << endl;
  else if( p_info.type == UNRELATED )
    cout << "unrelated" << endl;
  else if( p_info.type == MATE )
    cout << "mate" << endl;  
  else if( p_info.type == INTRA )
    cout << "intra" << endl;
  else
    cout << "inter" << endl;
#endif
  }

  set_pair_weight(my_subtype_pairset);
  //view_pairset(my_subtype_pairset, my_subtype_info, cout);

  build_main_pairset();
}

void
PairSetData::build_main_pairset()
{
  if( !my_subtype_pairset.size() )
    return;

  for( size_t s = 0; s < my_subtypes.size(); ++s )
    my_pairset_group.insert(make_pair(my_subtypes[s].relationship, s));

  //if( my_parser->get_analysis_options().pairset != SUBTYPES )
  //{
    main_to_sub_type_const_iterator     rel_group = my_pairset_group.type_begin();
    main_to_sub_type_const_iterator rel_group_end = my_pairset_group.type_end();

    PairSetCountInfo count_info;

    for( ; rel_group != rel_group_end; ++rel_group )
    {
      pairset_by_pedigree_type main_set;
      main_set.resize(my_parser->get_pedigree_count());

      // per each subtype.
      //
      for( size_t s = 0; s < rel_group->second.size(); ++s )
      {
        size_t k = rel_group->second[s];

        // per each pedigree.
        //
        for( size_t p = 0; p < my_subtype_pairset[k].size(); ++p )
        {
          const pairset_type& sub_set = my_subtype_pairset[k][p];

          // get a memeber.
          //
          for( size_t i = 0; i < sub_set.size(); ++i )
          {
            main_set[p].push_back(sub_set[i]);
          }
        }

        // if maintype is intraclass while subtype is not, need to include
        // reverse set.
        if(     is_intraclass(rel_group->first)
            && !is_selfclass(rel_group->first)
            && !is_intraclass(my_subtypes[k]) )
        {
          for( size_t p = 0; p < my_subtype_pairset[k].size(); ++p )
          {
            const pairset_type& sub_set = my_subtype_pairset[k][p];

            // get a memeber.
            //
            for( size_t i = 0; i < sub_set.size(); ++i )
            {
              WeightedMemberPair reverse_pair(WEIGHT_COUNT);

              reverse_pair.member_pair.first  = sub_set[i].member_pair.second;
              reverse_pair.member_pair.second = sub_set[i].member_pair.first;
              reverse_pair.pair_weight        = sub_set[i].pair_weight;

              main_set[p].push_back(reverse_pair);
            }
          }
        }
      }
      my_maintype_pairset.push_back(main_set);

      pairset_info p_info;
      p_info.name  = p_info.gname = long_relationship_name(rel_group->first);
      p_info.rcode = relationship_number(rel_group->first);
      p_info.distance1 = rel_group->first.distance1();
      p_info.distance2 = rel_group->first.distance2();

      if( is_intraclass(rel_group->first) )
        p_info.type = INTRA;
      else if( is_selfclass(rel_group->first) )
        p_info.type = SELF;
      else if( is_unrelated(rel_group->first) )
        p_info.type = UNRELATED;
      else if( is_mate(rel_group->first) )
        p_info.type = MATE;
      else
        p_info.type = OTHER;

      count_info.build_member_info(main_set, (p_info.type == INTRA));

      p_info.total_pair_count = count_info.get_total_pair_count();
      p_info.distinctive_pair_count = count_info.get_distinctive_pair_count();

      my_maintype_info.push_back(p_info);
    }

    set_pair_weight(my_maintype_pairset);
    //view_pairset(my_maintype_pairset, my_maintype_info, cout);

    my_main_built = true;
  //}
}

void
PairSetData::set_pair_weight(pairset_vector& p_vector)
{
  if( my_parser == NULL )
    return;

  SetWeight pair_weight;

  for( size_t s = 0; s < p_vector.size(); ++s )
    for( size_t w = 0; w < WEIGHT_COUNT; ++w )
      pair_weight.set_weight(p_vector[s], (weight_type)w);
}

void
PairSetData::view_pairset(const pairset_vector&      p_vector,
                          const pairset_info_vector& p_info,
                          ostream&                   out) const
{
  for( size_t r = 0; r < p_vector.size(); ++r )
  {
    out << p_info[r].name << endl;

    out << "  pedigree_count [n=" << setw(3) << p_vector[r].size() << "]: " << endl;

    for( size_t i = 0; i < p_vector[r].size(); ++i )
    {
      pairset_type::const_iterator        p = p_vector[r][i].begin();
      pairset_type::const_iterator pair_end = p_vector[r][i].end();

      for( size_t k = 0; p != pair_end; ++p, ++k )
      {
        if( k && k % 4 == 0 )
          out << endl;
        out << setprecision(3);
        out.setf(ios_base::fixed,ios_base::floatfield);
        out << '(' << setw(2) << p->member_pair.first->pedigree()->name()
            << ':' << setw(2) << p->member_pair.first->name()  << '-' << setw(2) << p->member_pair.first->index()
            << ',' << setw(2) << p->member_pair.second->name() << '-' << setw(2) << p->member_pair.second->index();
        for( size_t w = 0; w < p->pair_weight.size(); ++w )
          out << '|' << setw(5) << p->pair_weight[w];
        out << ") ";
      }
      out << endl;
    }
    out << "  total pair count : " << p_info[r].total_pair_count << endl
        << "  distintive pair count : "  << p_info[r].distinctive_pair_count << endl << endl;
  }
}

} // end of namespace FCOR
} // end of namespace SAGE
