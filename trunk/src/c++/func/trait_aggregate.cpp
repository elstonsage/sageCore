//===========================================================================
// File:      trait_aggregate.cpp
//                                                                          
// Author:    Kai He
//                                                                          
// History:   11/01/02
//                                                                          
// Notes:     Implementation of the following classes:
//            trait_aggregate
//                                                                          
// Copyright (c) 2000 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "func/trait_aggregate.h"
#include <iomanip>

using namespace std;
using namespace SAGE;

size_t add_trait_no(RPED::RefMPedInfo& mi, const FunctionParser::TraitData& par_data)
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

//============================================================================
//
//----------------------------------------------------------------------------
size_t check_parenthese(std::string seq)
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
//============================================================================
//
//----------------------------------------------------------------------------
string deparenthese(string str)
{   
  string str2;              
  size_t first = str.find('(',0); 
  size_t last  = str.find(')',0);
  for(size_t i=first+1;i<last;++i)
      str2+=str[i]; 
  return str2;    
}

//============================================================================
//
//----------------------------------------------------------------------------
void
trait_aggregate::process_aggregation(RPED::RefMultiPedigree& mp, string st1, string st2)
{
  RPED::RefMPedInfo&  mi = mp.info();
  const std::string& cst = st1;
  const std::string& vst = st2;
  my_sinfo.clear();

  RPED::RefMultiPedigree::pedigree_iterator   p_iter;
  for(p_iter = mp.pedigree_begin(); p_iter != mp.pedigree_end(); ++p_iter)
  {
    size_t m_index;
    size_t trait_no1,    trait_no2;
    double trait_value1, trait_value2, trait_value;

    RPED::RefPedInfo&   ped_info = p_iter->info();
    RPED::RefPedigree::member_iterator  m_iter, m_iter2;
    for(m_iter=p_iter->member_begin(); m_iter!=p_iter->member_end();++m_iter)
    {
      my_v_record.clear();
      my_record.mem_value.clear();
      my_traits.clear();

      m_index      = m_iter->index();
      trait_no1    = mi.trait_find(cst);
      trait_no2    = mi.trait_find(vst);
      trait_value1 = ped_info.trait(m_index, trait_no1);
      trait_value2 = ped_info.trait(m_index, trait_no2);

      if( (trait_no1 > trait_no2 && trait_value1 != 0) ||
          (trait_no1 < trait_no2 && trait_value2 != 0)  )
      {
         
         trait_value = trait_value1 * trait_value2;
         my_trait_no   = trait_no1;
         my_trait_name = mi.trait_info(trait_no1).name();
         my_record.trait_name=my_trait_name;
         my_record.trait_no  =my_trait_no;
         my_record.mem_value.insert(pair<string,double>(m_iter->name(),trait_value));
         my_v_record.push_back(my_record);
         my_traits.push_back(trait_value);
         my_sinfo.add(trait_value);
      }
      //print_out();
    }
  } 

}

//============================================================================
//
//----------------------------------------------------------------------------
string
trait_aggregate::parse_expression(RPED::RefMultiPedigree& mp, string express)
{
  string st;
  if(express.find("aggregate")!= string::npos)
     st = parse_aggregation(mp,express);

  if(express.find("merge")!=string::npos)
     st = parse_merge(express);

  return st;

}

//============================================================================
//
//----------------------------------------------------------------------------
void
trait_aggregate::create_adjusted_trait
(RPED::RefMultiPedigree& mp, const FunctionParser::TraitData& par_data)
{
  if(par_data.expr.find("aggregate")!= string::npos)
  {
     create_class_trait(mp, par_data);
     //print_out();
     dump();
  }
  if(par_data.expr.find("merge")!=string::npos)
     merge_class_trait(mp,par_data);

}
//============================================================================
//
//----------------------------------------------------------------------------
void
trait_aggregate::create_class_trait
(RPED::RefMultiPedigree& mp, const FunctionParser::TraitData& par_data)
{
  const string&  name = par_data.trait_name;
  const FunctionParser::constant_vector& constants = par_data.constants;
  unsigned int  time_limit = par_data.time_limit;

  const string& expr = parse_expression(mp, par_data.expr);

  //write_test_msg(expr);

  //dump();

  function my_func(my_errors);


  if(! (name.size() && expr.size()))
  {
    my_func.write_error_msg(SAGE::function::missing_info_msg, "", "");
    return;
  }

  if(my_func.exists(mp, name))
  {
    return;
  }
  
  function::python calculator(my_errors);
  if(! my_func.set_constants(calculator, constants, time_limit))
  {
    return;
  }
  
  if(! my_func.find_variables_used(mp, expr, calculator))
  {
    my_func.write_error_msg(SAGE::function::bad_syntax_msg, "", expr);
    return;
  }
  
  if(! calculator.compile(expr))
  {
    my_func.write_error_msg(SAGE::function::bad_syntax_msg, "", expr);
    return;
  }
  
  RPED::RefMPedInfo&  mi = mp.info();

  size_t   trait_no;

  trait_no = add_trait_no(mi,par_data);

  /*
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
  */

  SampleInfo sf;

  bool  runtime_error_written = false;
  RPED::RefMultiPedigree::pedigree_iterator   p_iter;
  for(p_iter = mp.pedigree_begin(); p_iter != mp.pedigree_end(); ++p_iter)
  {
    RPED::RefPedInfo&   ped_info = p_iter->info();
    ped_info.resize_traits(ped_info.trait_count() + 1);
    size_t ind_count = 0;
    
    RPED::RefPedigree::member_iterator  m_iter;
    for(m_iter=p_iter->member_begin();m_iter!=p_iter->member_end();++m_iter)
    {
      size_t    m_index = m_iter->index();
      double    trait_value;

      if(my_func.set_inputs(calculator, mi, ped_info, m_index))  
      {
        if(expr.find("class")!=std::string::npos)
        {
           if(is_in_class(mp,ped_info,expr,m_index)==true )
           {
        	trait_value = calculator.run(expr, time_limit);
        	if(SAGE::isnan(trait_value) && (! runtime_error_written))
           		runtime_error_written = true;
       		sf.add(trait_value);
       		ind_count++;
       		ped_info.set_trait(m_index,trait_no,
	        	           doub2str(trait_value),
		        	   mi.trait_info(trait_no));
	   }
           else
           {
                ped_info.set_trait(m_index,trait_no,
				   doub2str(numeric_limits<double>::quiet_NaN()),
                                   mi.trait_info(trait_no));
           }
	}
        else
	{
                trait_value = calculator.run(expr, time_limit);
                if(SAGE::isnan(trait_value) && (! runtime_error_written))
                        runtime_error_written = true;
                sf.add(trait_value);
                ind_count++;
                ped_info.set_trait(m_index,trait_no,
                                   doub2str(trait_value),
                                   mi.trait_info(trait_no));
	}
      } 
    }   

    if(ind_count == 0)
    {
          my_errors << priority(error)
                    << "ZERO member in the class "
                    << "'" << mi.trait_info(trait_no).name() //expr
                    << " or " << expr
                    << "', classify was failed. " << endl;    
          runtime_error_written = true;
    }
    if(ind_count == 1)
    {
          my_errors << priority(error)
                    << "Only one member in the class "
                    << "'" << mi.trait_info(trait_no).name() //expr
                    << " or " << expr
                    << "', trait adjustment in this class was failed. " << endl;
          runtime_error_written = true;
    }
    if(ind_count > 1)
    {
      my_errors<<priority(error)
      <<setfill(' ')<<setw(9)<<mi.trait_info(trait_no).name()
      <<setfill(' ')<<setw(10)<<setprecision(5)<<sf.max()          
      <<setfill(' ')<<setw(10)<<sf.min()
      <<setfill(' ')<<setw(05)<<sf.count()
      <<setfill(' ')<<setw(12)<<sf.mean()
      <<setfill(' ')<<setw(12)<<sf.variance()
      <<setfill(' ')<<setw(12)<<sf.standard_deviation()<<endl;
      runtime_error_written = true;
    }
  }  
}

//============================================================================
//
//----------------------------------------------------------------------------
void   
trait_aggregate::merge_class_trait
(RPED::RefMultiPedigree& mp, const FunctionParser::TraitData& par_data)
{
  //const string&  name = par_data.trait_name;

  string str = par_data.expr;

  char* pch;

  list_string lst;

  string str2 = deparenthese(str);

  char* ch = &str2[0];

  pch = strtok (ch,",");  

  while (pch != NULL)
  {
    lst.push_back(pch);
    pch = strtok (NULL, " ,.");
  }

  RPED::RefMPedInfo&  mi = mp.info();

  size_t   trait_no;

  trait_no = add_trait_no(mi, par_data);

  /*
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
    trait_no = mi.add_continuous_trait(name,par_data.usage);
    mi.trait_info(trait_no).set_string_missing_code(par_data.missing);
    mi.trait_info(trait_no).set_numeric_missing_code(str2doub(par_data.missing));
  }
  */

  SampleInfo sf;
  bool  runtime_error_written = false;
  RPED::RefMultiPedigree::pedigree_iterator   p_iter;
  for(p_iter = mp.pedigree_begin(); p_iter != mp.pedigree_end(); ++p_iter)
  {
    RPED::RefPedInfo&   ped_info = p_iter->info();
    ped_info.resize_traits(ped_info.trait_count() + 1);
    size_t ind_count = 0;

    RPED::RefPedigree::member_iterator  m_iter;
    for(m_iter=p_iter->member_begin();m_iter!=p_iter->member_end();++m_iter)
    {
      size_t m_index = m_iter->index();
      size_t trait_no2; 
      double trait_value;
      for(list_string_iterator it=lst.begin(); it!=lst.end();++it)
      {
         const string& trait_name = *it;
         trait_no2   = mi.trait_find(trait_name);
         if(SAGE::isnan(trait_value) && (! runtime_error_written))
            runtime_error_written = true;
         if(!ped_info.trait_missing(m_index, trait_no2) )
         {
            trait_value = ped_info.trait(m_index, trait_no2);
            //cout<<"trait_value="<<trait_value<<endl;
            ind_count++;
            sf.add(trait_value);
            ped_info.set_trait(m_index,trait_no,
                               doub2str(trait_value),
                               mi.trait_info(trait_no));
         }
      }
    }
    //dump(sf,mi.trait_info(trait_no).name());
  }

}

//============================================================================
//
//----------------------------------------------------------------------------
string
trait_aggregate::parse_merge(string express)
{
  string seq = express;

  if(seq[0]=='+'||seq[0]=='-'||seq[0]=='*'||seq[0]=='/')
  {
     cout<<"The first character is invalid! please add it."<<endl;
     exit(1);
  }

  if(check_parenthese(seq)!=1)
  {
     cout<<"The parenthese (is)are missing or duplicating."<<endl;
     exit(1);
  }

  list sq1,sq2,sq3;
  for(size_t i=0; i<express.size(); ++i)
      sq1.push_back(express[i]);

  list_iterator f1, f2;
  string str1, str2, str3, str4;

  f1 = find(sq1.begin(),sq1.end(),'(');
  f2 = find(sq1.begin(),sq1.end(),')');

  for(list_iterator i = f1; i != f2; i++)
  {
      if(*i!='(')
         str1+=*i;
  }
  for(size_t i=0;i<str1.size();++i)
  {
      if(str1[i]==',')
         str1[i]='+';
  }

  return str1;
}

//============================================================================
//
//----------------------------------------------------------------------------
string
trait_aggregate::parse_aggregation(RPED::RefMultiPedigree& mp, string express)
{
  string seq     = express;
  string subseq  = "aggregate";
  string subseq2 = "isMissing";

  if(seq[0]=='+'||seq[0]=='-'||seq[0]=='*'||seq[0]=='/')
  {
     cout<<"The first character is invalid! please add it."<<endl;                  
     exit(1);
  }

  if(check_parenthese(seq)!=1)
  {
     cout<<endl
         <<"One or more parenthese(s) is(are) missing in the expression of "
         <<"the function block(s), check and try again!"<<endl;
     exit(1);
  }

  list sq1,sq2,sq3;
  for(int i=0; i<seq.size(); ++i)
      sq1.push_back(seq[i]);
  for(int j = 0; j < subseq.size(); ++j)
      sq2.push_back(subseq[j]);

  list_iterator f_end, f1, f2;
  string str1,str2,str3,str4;
  f_end = find_end(sq1.begin(),sq1.end(),sq2.begin(),sq2.end() );
  for(list_iterator i = f_end; i != sq1.end(); i++)
      sq3.push_back(*i), str1+=*i;
  f1 = find(sq3.begin(),sq3.end(),'(');
  f2 = find(sq3.begin(),sq3.end(),')');
  for(list_iterator i = sq1.begin(); i != f_end; i++)
      str2+=*i;
  for(list_iterator i=sq3.begin(); i!=f1; ++i)
      str3+=*i;
  for(list_iterator i = f1; i != sq3.end(); ++i)
      str4+=*i;

  string q1,q2,q3;
  for(int k=1; k<str4.size()-1; ++k) 
      q1+=str4[k];
  int pos = q1.find(',',0);
  for(int k=0; k<pos; ++k) 
      q2+=q1[k];    
  for(int k=pos+1; k<q1.size(); ++k)
      q3+=q1[k];

  process_aggregation(mp,q2,q3);

  
  double dr;

  if(str3=="aggregate_mean")      dr = my_sinfo.mean();

  else if(str3=="aggregate_var")  dr = my_sinfo.variance();

  else  my_errors<<priority(error)<<"Invalid function in the expression '"
   	         <<str3<<"'"
                 <<endl;

  str2 += doub2str(dr,0,-1,0,'\0');

  /*
  if(check_parenthese(str2)==2)
     str2 += ')';
  if(check_parenthese(str2)==3)
     str2.insert(0,'(');
  */

  return str2;

}

//============================================================================
//
//----------------------------------------------------------------------------
bool
trait_aggregate::is_in_class(RPED::RefMultiPedigree& mp, RPED::RefPedInfo& ped_info, 
			     string seq, size_t m_index)
{
  size_t i, start, end;
  string str,str2,str3;
  
  if(seq.find('(')!=string::npos)
  {
	start = seq.find('(',0);
	end   = seq.find(')',0);
  }
  for(i=start+1; i!=end;++i)
	  str+=seq[i];
  
  for(i=0;i<seq.size();++i)
  {
      if(isalnum(str[i])||str[i]=='_' ||str[i]=='-')
	     str2+=str[i];
      else   break;
  }

  /*
  for(i=0;i<6;++i)
  {
	if(stra[i]==str2)
	  str3.assign(str2);
  }
  */

  size_t trait_no    = mp.info().trait_find(str2);
  double trait_value = ped_info.trait(m_index, trait_no);
  if(trait_value == 1)
     return true;
  else
     return false;

}

//============================================================================
//
//----------------------------------------------------------------------------
void 
trait_aggregate::dump()
{
  cout<<endl
      <<"SampleInfo data dump: "<<endl
      <<"-----------------------------------------------------"<<endl;
  for(std::vector<double>::iterator it = my_traits.begin();
  				    it!= my_traits.end(); ++it)    
  cout<<*it<<" ";
  cout<<endl
      <<setfill(' ')<<setw(10)<<"Max" 
      <<setfill(' ')<<setw(10)<<"Min" 
      <<setfill(' ')<<setw(10)<<"Count" 
      <<setfill(' ')<<setw(10)<<"Mean"
      //<<setfill(' ')<<setw(10)<<"Sum" 
      <<setfill(' ')<<setw(10)<<"Var" 
      <<setfill(' ')<<setw(10)<<"S.D"<<endl;        
      //<<setfill(' ')<<setw(10)<<"S.E"    
      //<<setfill(' ')<<setw(10)<<"V_Coeff"<<endl; 
  cout<<my_sinfo.max()  
      <<setfill(' ')<<setw(10)<<my_sinfo.min()  
      <<setfill(' ')<<setw(10)<<my_sinfo.count()                
      <<setfill(' ')<<setw(10)<<my_sinfo.mean()  
      <<setfill(' ')<<setw(10)<<my_sinfo.variance()                
      <<setfill(' ')<<setw(10)<<my_sinfo.standard_deviation()<<endl<<endl;                

}
//============================================================================
//
//----------------------------------------------------------------------------
void   
trait_aggregate::dump(SampleInfo sf, string name)
{
    my_errors<<priority(error)
    <<setfill(' ')<<setw(9)<<name  //mi.trait_info(trait_no).name()
    <<setfill(' ')<<setw(10)<<setprecision(5)<<sf.max()
    <<setfill(' ')<<setw(10)<<sf.min()    
    <<setfill(' ')<<setw(05)<<sf.count()  
    <<setfill(' ')<<setw(12)<<sf.mean()   
    <<setfill(' ')<<setw(12)<<sf.variance()
    <<setfill(' ')<<setw(12)<<sf.standard_deviation()<<endl;


}
//============================================================================
//
//----------------------------------------------------------------------------
void 
trait_aggregate::print_out()
{
   for(v_record_iterator it =my_v_record.begin();
		         it!=my_v_record.end();++it)
   {
       cout<<"trait_name: "<<it->trait_name;
       cout<<"  trait_no: "<<it->trait_no<<endl;

       for(miter it2 =it->mem_value.begin();
                 it2!=it->mem_value.end();++it2)
       {
           cout<<"Name: "<<it2->first;
           cout<<"  value: "<<it2->second;
           cout<<endl;
       }
       cout<<endl;
   }

}
//============================================================================
//
//----------------------------------------------------------------------------
void
trait_aggregate::write_error_msg(const string& message)
{
  my_errors << message << endl;
}
void          
trait_aggregate::write_test_msg(const string&  msg)
{
  my_errors << priority(error) << msg << endl;
}
//============================================================================
//
//----------------------------------------------------------------------------
