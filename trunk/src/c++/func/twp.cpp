//============================================================================
// File:      twp.cpp (trimming-winsorization process)
//
// Author:    Kai He
//
//
// History:   12/2003 Initial version
//
// Notes:
//
// Copyright (c) 2002 R.C. Elston
// All Rights Reserved
//============================================================================

#include "func/twp.h"
#include <iostream>

namespace SAGE {
namespace FUNC {

//----------------------------------------------------------------------------
//
// default constructor
//
//----------------------------------------------------------------------------

trim_winsor_process::trim_winsor_process
()
{ }

//----------------------------------------------------------------------------
//
// constructor with one argument
//
//----------------------------------------------------------------------------

trim_winsor_process::trim_winsor_process
(cerrorstream& err):my_errors(err)
{
  my_people.clear();
  my_sort_people.clear();
  my_trait = "";
  my_gamma = 0;
  my_missing_members = 0;
}

//----------------------------------------------------------------------------
//
// copy constructor
//
//----------------------------------------------------------------------------

trim_winsor_process::trim_winsor_process
(const trim_winsor_process& rhs)
{
  my_errors = rhs.my_errors;
  my_people = rhs.my_people;
  my_sort_people = rhs.my_sort_people;
  my_trait  = rhs.my_trait;
  my_gamma  = rhs.my_gamma;
  my_missing_members= rhs.my_missing_members;
}

//----------------------------------------------------------------------------
//
// assignment constructor
//
//----------------------------------------------------------------------------

trim_winsor_process& trim_winsor_process::operator=
(const trim_winsor_process& rhs)
{
  if(this != &rhs)
  {
     my_errors = rhs.my_errors;
     my_people = rhs.my_people;
     my_sort_people = rhs.my_sort_people;
     my_trait  = rhs.my_trait;
     my_gamma  = rhs.my_gamma;
     my_missing_members= rhs.my_missing_members;
  }
  return *this;
}

//----------------------------------------------------------------------------
//
// destructor
//
//----------------------------------------------------------------------------

trim_winsor_process::~trim_winsor_process()
{}

//----------------------------------------------------------------------------
//
// check_parenthese(...)
//
//----------------------------------------------------------------------------

size_t 
trim_winsor_process::check_parenthese(string seq)
{
  size_t c1,c2; c1=c2=0;
  for(size_t i=0; i< seq.size(); ++i)
  {
      if(seq[i]=='(')         c1++;
      if(seq[i]==')')         c2++;
  }
  if(c1==c2)
     return 1;
  else if(c1>c2)
     return 2;
  else if(c1<c2)
     return 3;
  else
     return 0;
}

//----------------------------------------------------------------------------
//
// remove parenthese(...)
//
//----------------------------------------------------------------------------

string 
trim_winsor_process::deparenthese(string str)
{
  string str2;
  size_t first = str.find('(',0);
  size_t last  = str.find(')',0);
  for(size_t i=first+1;i<last;++i)
      str2+=str[i];
  return str2;
}

//----------------------------------------------------------------------------
//
// parsing the expression in function for twp
//
//----------------------------------------------------------------------------

void 
trim_winsor_process::parse_twp_expression(RPED::RefMultiPedigree& mp, string expr)
{
  if(expr[0]=='+'||expr[0]=='-'||expr[0]=='*'||expr[0]=='/')
  {
     cout<<"The first character is invalid! please add it."<<endl;
     exit(1);
  }

  if(check_parenthese(expr)!=1)
  {
     cout<<endl
         <<"One or more parenthese(s) is(are) missing in the expression of "
         <<"the function block(s), check and try again!"<<endl;
     exit(1);
  }

  // get my_process

  string process_name;
  size_t first = expr.find('(',0);
  for(size_t i=0; i<first;++i)
      process_name+=expr[i];

  my_process = process_name;

  // remove parentheses

  string str = deparenthese(expr);
  parameter_list plst;
  plst.clear();

  if(str.size() > 0)
  {
    char* ch  = &str[0];
    char* pch = strtok (ch,",");

    while (pch != NULL)
    {
      plst.push_back(pch);
      pch = strtok (NULL, " ,");
    }
  }

  my_trait = *plst.begin();

  size_t count = 0;  
  for(size_t ti = 0; ti < mp.info().trait_count(); ++ti)
  {
     string name = mp.info().trait_info(ti).name();

     if(my_trait == name)
  
        ++count;
  }
  if(count < 1)
        write_error_message(invalid_trait, my_trait);

  plst.pop_front();

  const string &st = *plst.begin();

  my_gamma = str2doub(st);

  if( my_gamma < 0 || my_gamma > 1)

     write_error_message(invalid_gamma, "gamma");
  
}

//----------------------------------------------------------------------------
//
// adding new trait number for added trait(s) in pedigree data
//
//----------------------------------------------------------------------------

size_t 
trim_winsor_process::add_trait_number(RPED::RefMPedInfo& mi, parser_data par_data)
{
  size_t trait_no;
  const string&  name = par_data.trait_name;

  if(par_data.binary)
  {
    trait_no = mi.add_binary_trait(name, par_data.usage);
    mi.trait_info(trait_no).set_string_missing_code(par_data.missing);
    mi.trait_info(trait_no).set_numeric_missing_code(str2doub(par_data.missing));
    mi.trait_info(trait_no).set_string_affected_code(par_data.affected);
    mi.trait_info(trait_no).set_numeric_affected_code(str2doub(par_data.affected));
    mi.trait_info(trait_no).set_string_unaffected_code(par_data.unaffected);
    mi.trait_info(trait_no).set_numeric_unaffected_code(str2doub(par_data.unaffected));
  }
  else
  {
    trait_no = mi.add_continuous_trait(name, par_data.usage);
    mi.trait_info(trait_no).set_string_missing_code(par_data.missing);
    mi.trait_info(trait_no).set_numeric_missing_code(str2doub(par_data.missing));
  }

  return trait_no;
}

//----------------------------------------------------------------------------
//
// create twp data as new trait in pedigree data
//
//----------------------------------------------------------------------------

void
trim_winsor_process::create_twp(RPED::RefMultiPedigree& mp, parser_data pd)
{
  ///// parsing associated allele with given marker

  parse_twp_expression(mp, pd.expr);

  ///// build data for sorting and positioning

  build_data(mp,pd);

  ///// process trimming

  if(my_process == "trim")
  {
     sort_positioning_data();
     trim_data();
     create_twp_trait(mp,pd);
  }

  ///// process winsorization

  if(my_process == "winsor")
  {
     sort_positioning_data();
     winsor_data();
     create_twp_trait(mp,pd);
  }

}

//----------------------------------------------------------------------------
//
// building data structure for twp
//
//----------------------------------------------------------------------------

void
trim_winsor_process::build_data(RPED::RefMultiPedigree& mp, parser_data pd)
{
  my_people.clear();
  RPED::RefMPedInfo& mpinfo = mp.info();
  size_t trait_number = mpinfo.trait_find(my_trait);
  pedigree_iterator ped;
  for(ped = mp.pedigree_begin(); ped!= mp.pedigree_end(); ++ped)
  {
    member_iterator    mem;
    RPED::RefPedInfo& pinfo = ped->info();
    for(mem = ped->member_begin(); mem != ped->member_end(); ++mem)
    {
      double value = pinfo.trait(mem->index(),trait_number);

      if(SAGE::isnan(value))
      {
         my_missing_members++;
         my_errors << priority(error) 
	           << "An individual "
                   << mem->name()
                   << " in "
                   << ped->name()
                   <<" has a missing value."<<endl;
      }
      else
      {
        Person ps;
        ps.mem = mem;
        ps.ped = ped;
        ps.trait = value;
        ps.position = 0;
        my_people.push_back(ps);
      }
    }
  }
}

//----------------------------------------------------------------------------
//
// sorting and positioning member data
//
//----------------------------------------------------------------------------

void
trim_winsor_process::sort_positioning_data()
{
  my_pos_map.clear();
  
  vector<Person>::iterator it  = my_people.begin();
  for(; it!=my_people.end();++it)
     my_trait_map.insert(pair(it->trait,*it));

  trait_map_iterator tmi = my_trait_map.begin();
  for(size_t i = 1; tmi != my_trait_map.end(); ++tmi, ++i)
  {
     tmi->second.position = i;
     my_pos_map.insert(make_pair(i,tmi->second));
  }
}

//----------------------------------------------------------------------------
//
// trimming data
//
//----------------------------------------------------------------------------

void
trim_winsor_process::trim_data()
{
  size_t g = my_trait_map.size() * (size_t)my_gamma;

  size_t p = my_trait_map.size() - g;

  pos_map_iterator tmi = my_pos_map.begin();

  for(; tmi != my_pos_map.end(); ++tmi)
  {
     if(tmi->second.position < g)

        my_pos_map.erase(tmi);

     if(tmi->second.position > p+1)

        my_pos_map.erase(tmi);
  } 
        
}

//----------------------------------------------------------------------------
//
// winsorizing data
//
//----------------------------------------------------------------------------

void
trim_winsor_process::winsor_data()
{
  size_t g = my_trait_map.size() * (size_t)my_gamma;

  size_t p = my_trait_map.size() - g;

  trait_map_iterator tmi = my_trait_map.begin();

  double value_g = 0.0;

  double value_p = 0.0;

  for(; tmi != my_trait_map.end(); ++tmi)
  {
     if(tmi->second.position == g)
        value_g = tmi->second.trait;
     if(tmi->second.position == p+1)
        value_p = tmi->second.trait;
  }

  pos_map_iterator pmi = my_pos_map.begin();
  for(; pmi != my_pos_map.end(); ++pmi)
  {
     if(pmi->first < g)
        pmi->second.trait = value_g;
     if(pmi->first > p+1)
        pmi->second.trait = value_p;
  }
}

//----------------------------------------------------------------------------
//
// create twp data as new trait in pedigree data
//
//----------------------------------------------------------------------------

void
trim_winsor_process::create_twp_trait(RPED::RefMultiPedigree& mp, parser_data pd)
{
  RPED::RefMPedInfo&  mi = mp.info();
  size_t        trait_no;
  trait_no = add_trait_number(mi, pd);

  RPED::RefMultiPedigree::pedigree_iterator   p_iter;
  for(p_iter = mp.pedigree_begin(); p_iter != mp.pedigree_end(); ++p_iter)
  {
    RPED::RefPedInfo&   ped_info = p_iter->info();
    ped_info.resize_traits(ped_info.trait_count() + 1);

    RPED::RefPedigree::member_iterator  m_iter;
    for(m_iter=p_iter->member_begin();m_iter!=p_iter->member_end();++m_iter)
    {
      size_t m_index = m_iter->index();

      pos_map_iterator it = my_pos_map.begin();
      for(; it!=my_pos_map.end();++it)
      {
          if(it->second.ped == p_iter && it->second.mem->index() == m_index)
             ped_info.set_trait(m_index,
                                trait_no,
                                doub2str(it->second.trait),
                                mi.trait_info(trait_no));

      }

    }
  }
}

//----------------------------------------------------------------------------
//
// dump and test member data
//
//----------------------------------------------------------------------------

void 
trim_winsor_process::dump_members()
{
  cout<<endl;
  cout<<"The non-sorted and non-positioning data for each member"<<endl;
  cout<<"my_process:"<<get_process()<<endl;
  cout<<"my_trait:"  <<get_trait()<<endl;
  cout<<"my_gamma:"  <<get_gamma()<<endl;
  cout<<"my_people:" <<get_members()<<endl;
  cout<<"my_missing:"<<get_missing_members()<<endl;

  cout<<setfill(' ')<<setw(10)
      <<"Pid"
      <<setfill(' ')<<setw(10)
      <<"id"
      <<setfill(' ')<<setw(10)
      <<"trait"
      <<setfill(' ')<<setw(10)
      <<"position"<<endl;  

  vector<Person>::iterator it = my_people.begin();
  for(; it!=my_people.end();++it) 
  {
    cout<<setfill(' ')<<setw(10)
        <<it->ped->name()
        <<setfill(' ')<<setw(10)  
        <<it->mem->name()
        <<setfill(' ')<<setw(10)  
        <<it->trait
	<<setfill(' ')<<setw(10)
        <<it->position<<endl;
  }
  cout<<endl;

}

//----------------------------------------------------------------------------
//
// dump and test sorted member data
//
//----------------------------------------------------------------------------

void 
trim_winsor_process::dump_sort_members()
{
  cout<<endl;
  cout<<"The sorted and positioned data for each member"<<endl;
  cout<<"my_process:"<<get_process()<<endl;
  cout<<"my_trait:"  <<get_trait()<<endl;
  cout<<"my_gamma:"  <<get_gamma()<<endl;
  cout<<"my_people:" <<my_trait_map.size()<<endl;
  cout<<"my_missing:"<<2 * (my_people.size()*my_gamma - 1)<<endl;

  cout<<setfill(' ')<<setw(10)
      <<"Pid"
      <<setfill(' ')<<setw(5)
      <<"id"
      <<setfill(' ')<<setw(10)
      <<"ps.trait"
      <<setfill(' ')<<setw(10)
      <<"position"<<endl;  

  trait_map_iterator it = my_trait_map.begin();
  for(; it!=my_trait_map.end();++it) 
  {
    cout<<setfill(' ')<<setw(10)
        <<it->second.ped->name()
        <<setfill(' ')<<setw(5)  
        <<it->second.mem->name()
        <<setfill(' ')<<setw(10)
        <<it->second.trait     
	<<setfill(' ')<<setw(10)
        <<it->second.position<<endl;
  }
  cout<<endl;
}

//----------------------------------------------------------------------------
//
// dump and test sorted and positioned member data
//
//----------------------------------------------------------------------------

void 
trim_winsor_process::dump_pos_members()
{
  cout<<endl;
  cout<<"The positioned data for each member"<<endl;
  cout<<"my_process:"<<get_process()<<endl;
  cout<<"my_people:" <<my_pos_map.size()<<endl;

  cout<<setfill(' ')<<setw(10)
      <<"Pid"
      <<setfill(' ')<<setw(5)
      <<"id"
      <<setfill(' ')<<setw(10)
      <<"ps.trait"
      <<setfill(' ')<<setw(10)
      <<"position"<<endl;  

  pos_map_iterator it = my_pos_map.begin();
  for(; it!=my_pos_map.end();++it) 
  {
    cout<<setfill(' ')<<setw(10)
        <<it->second.ped->name()
        <<setfill(' ')<<setw(5)  
        <<it->second.mem->name()
        <<setfill(' ')<<setw(10)
        <<it->second.trait     
	<<setfill(' ')<<setw(10)
        <<it->second.position<<endl;
  }
  cout<<endl;
}

//----------------------------------------------------------------------------
//
// writing error message
//
//----------------------------------------------------------------------------

void 
trim_winsor_process::write_error_message(error_message msg, const string token)
{
  switch(msg)
  {
     case invalid_trait :
          my_errors<<priority(error)<<token<<" is invalid because it doesn't"
                                           <<" exist in the pedigree."<<endl;
	  break;

     case duplicate_trait :
          my_errors<<priority(error)<<token<<" is existing trait in the pedigree."
                                           <<" Duplicate."<<endl;
          break;
     case invalid_gamma :
          my_errors<<priority(error)<<token<<"The second parameter, gamma, should be"
                                    <<" decimal number. 0 <= gamma <= 1"<<endl;
          break;
     default :;
  }
}

} // End namespace FUNC
} // End namespace SAGE
