//===================================================================
//
//  File:       mp_utilities.ipp
//
//  Author:     Kai He
//
//===================================================================
namespace SAGE {
namespace MPED {

inline int 
mp_utilities::nuclear_family_count(member_const_pointer mcp)
{
/*
  uint num = 0;
  if( mcp->parent1()==NULL&&mcp->mate_count()==1 ) num=1;
  if( mcp->parent1()==NULL&&mcp->mate_count()>1  ) num=mcp->mate_count();
  if( mcp->parent1()!=NULL&&mcp->mate_count()>0  ) num=mcp->mate_count()+1;
  if( mcp->parent1()!=NULL&&mcp->mate_count()==0 ) num=1;
*/
  uint num = mcp->mate_count() + (int)(mcp->parent1()!=NULL);
  
  return num;
}
inline int
mp_utilities::nuclear_family_count(member_pointer mp)
{
  return nuclear_family_count( (member_const_pointer)mp );
}
inline int
mp_utilities::nuclear_family_count(member_const_iterator mci)
{
  return nuclear_family_count(&*mci);
}
inline int
mp_utilities::nuclear_family_count(member_iterator mi)
{
  return nuclear_family_count(&*mi);
}
//===================================================================
inline int 
mp_utilities::connector_count(family_const_pointer fcp)
{
  int num=0;

  if ( is_connector(fcp->parent1())==true )
      ++num;
  if ( is_connector(fcp->parent2())==true )  
      ++num;
/*
  for(member_const_iterator mi =(fcp->subpedigree())->member_begin(); 
                            mi!=(fcp->subpedigree())->member_end(); ++mi)
  {
     if( fcp->parent1()==mi->parent1()&&fcp->parent2()==mi->parent2() )
     {   if(is_connector(mi)==true)	++num;    }
  }
*/

  for(offspring_const_iterator oi = fcp->offspring_begin();
                               oi!= fcp->offspring_end(); ++oi)
  {
     if( is_connector(&*oi)==true ) ++num;
  }

  return num;
}
inline int
mp_utilities::connector_count(family_pointer fp)
{
  return connector_count( (family_const_pointer)fp );
}
inline int
mp_utilities::connector_count(family_const_iterator fci)          
{
  return connector_count( &*fci );
}
inline int
mp_utilities::connector_count(family_iterator fi)          
{
  return connector_count( &*fi );
}
//===================================================================    

inline bool
mp_utilities::is_founder(member_const_pointer mcp)
{   
    
    if(mcp->parent1()==NULL)// || mcp->parent2()==0)
        return true;
    else
        return false;
}
inline bool     
mp_utilities::is_founder(member_pointer mp)           
{
  return is_founder( (member_const_pointer)mp );
}
inline bool     
mp_utilities::is_founder(member_const_iterator mci)
{
  return is_founder(&*mci);
}
inline bool     
mp_utilities::is_founder(member_iterator mi)           
{
  return is_founder(&*mi);
}
inline bool     
mp_utilities::is_founder(const member_base& mi)           
{
  return is_founder(&mi);
}

//===================================================================

inline bool
mp_utilities::is_connector(member_const_pointer mcp)
{
  int num = nuclear_family_count(mcp);

  if(num>1)  return true;
  else	     return false;

}
inline bool
mp_utilities::is_connector(member_pointer mp)
{
    return is_connector( (member_const_pointer)mp );
}
inline bool
mp_utilities::is_connector(member_const_iterator mci)
{
    return is_connector(&*mci);
}
inline bool
mp_utilities::is_connector(member_iterator mi)
{      
    return is_connector(&*mi);
}

//===================================================================
inline bool 
mp_utilities::in_family(member_const_pointer mcp, family_const_pointer fcp)
{
/*
  if(mcp->name()==fcp->name1()     || 
     mcp->name()==fcp->name2()     ||
     (mcp->parent1()==fcp->parent1()&&mcp->parent2()==fcp->parent2()) )
     return true;
  else
     return false;
*/
  member_const_pointer p1=fcp->parent1();

  member_const_pointer p2=fcp->parent2();

  return mcp==p1 || mcp==p2 || (mcp->parent1()==p1 && mcp->parent2()==p2 );

}
inline bool
mp_utilities::in_family(member_pointer mp, family_pointer fp)
{
  return in_family( (member_const_pointer)mp, (family_const_pointer)fp );
}
inline bool
mp_utilities::in_family(member_const_iterator mci, family_const_iterator fci)
{
  return in_family(&*mci, &*fci);
}
inline bool
mp_utilities::in_family(member_iterator mi, family_iterator fi)
{
  return in_family(&*mi, &*fi);
}
//===================================================================
inline bool
mp_utilities::is_terminal_family(family_const_pointer fcp)
{
/*
  if(connector_count(fcp) <= 1)
     return true;
  else
     return false;     
*/

  return connector_count(fcp)<=1;
}
inline bool
mp_utilities::is_terminal_family(family_pointer fp)
{
  return is_terminal_family((family_const_pointer)fp);
}
inline bool
mp_utilities::is_terminal_family(family_const_iterator fci)
{
  return is_terminal_family(&*fci);
}
inline bool
mp_utilities::is_terminal_family(family_iterator fi)
{
  return is_terminal_family(&*fi);
}
//===================================================================

} // End namespace MPED
} // End namespace SAGE
