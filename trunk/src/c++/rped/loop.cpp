/////////////////////////////////////////////////////////////
// LoopChecker:  function object used to perform pedigree  //
//               loops (marriage loops, non-marriage loops)//
//               checking.                                 //
//     (implementation code )                              //
// Author : Qing Sun                                       //
//                                                         //
// History: version 1.0, August 8, 1996                    //
//          version 1.2  modified to use new pedigree class//
//                       , Oct,29,1998                     //
//                                                         //
//          bug fix in check_nmloop()                      //
//                       June 11, 2002  -djb               //
//                                                         //
// Copyright (c) 1996 R.C. Elston                          //
//                                                         //
/////////////////////////////////////////////////////////////

#include <limits.h>

#include "rped/loop.h"
#include "error/errormanip.h"

namespace SAGE {
namespace RPED {

void LoopChecker::reset()
{
  my_marriage_loops = 0;
  my_non_marriage_loops = 0;
}

//If check_pedigree fails, return false, otherwise return true.
bool LoopChecker::check_pedigree(const RefPedigree* ped, bool verbose)
{
  reset();

  if(!ped)  
  {
   errors << priority(error) << "Pedigree does not exist." << endl;
   return false;
  }

  //for each subpedigree 
  RefPedigree::subpedigree_const_iterator sp;
  for(sp=ped->subpedigree_begin(); sp != ped->subpedigree_end(); ++sp)
  { 
    if(sp->member_count() < 4)  continue;

    if(!check_mloop(&*sp, verbose) || !check_nmloop(&*sp))
      return false;
  }

  if( marriage_loops() && non_marriage_loops() )    
    my_non_marriage_loops -= my_marriage_loops;    

  return true;
}

//************  checking marriage loop  ****************
bool LoopChecker::check_mloop(const RefSubpedigree* subp, bool verbose)
{
  if(!subp)
    return false;

  if(subp->member_count() < 4)
    return true;

  //family_set[i]=  0 : unvisited
  //                1 : visited once
  //                2 : visited twice
  //                3 : loop cut at this marriage

  vector<uint>  family_set( subp->family_count(), 0);
  set< int,less<int> > loop_members;
  
  for(size_t i=0; i< subp->family_count(); ++i)
  {
    RefPedigree::family_const_pointer family = &subp->family_index(i);

    if(family->parent1()->mate_count() == 1 && 
       family->parent2()->mate_count() == 1 )
      continue;
  
    //not a loop
    if(!get_spouse_set(i,family_set,subp) )
    {
      for(size_t j=0; j<family_set.size(); ++j)
        if(family_set[j] != 3)  family_set[j] = 0;
          continue;
    }
    else   //loop exist
    {
      ++my_marriage_loops;

      //get mmebers in the loop into loop_member
      for(size_t l=0; l<family_set.size(); ++l)
      {
        if(family_set[l]==2)
        { 
          loop_members.insert( subp->family_index(l).parent1()->index());
          loop_members.insert( subp->family_index(l).parent2()->index());
        } 
    
       if(family_set[l] != 3 )      family_set[l] = 0;
      }
    }

   if(verbose)
   {
     errors << "A marriage loop exists in pedigree '" << subp->name()
            << "' among individuals:";

     set<int,less<int> >::const_iterator p;
     for(p = loop_members.begin(); p != loop_members.end(); ++p)
       errors << "  '" << subp->member_index(*p).name() << "'";

     loop_members.clear();
    }
  }

  return true;
}

bool LoopChecker::get_spouse_set( int pivot_fam, vector<uint>& fam_set,
                                  const RefSubpedigree* subp) const
{
  RefPedigree::member_const_pointer  p1, p2;

  p1 = subp->family_index(pivot_fam).parent1();
  p2 = subp->family_index(pivot_fam).parent2();

  deque<RefPedigree::member_const_pointer >  que;

  if     ( p1->mate_count() > 1 ) que.push_back(p1);
  else if( p2->mate_count() > 1)  que.push_back(p2);
  else                            return false;

  set<uint,less<uint> >     spouse_set ;

  while(!que.empty())
  {
    p1 = que.front();
    spouse_set.insert(p1->index());
    que.pop_front();

    size_t loop_flag = 0;
    RefPedigree::mate_const_iterator         miter;
    for(miter = p1->mate_begin(); miter != p1->mate_end(); ++miter)
    {
      p2 = &(miter->mate());

      if(p2->mate_count() <= 1)  continue;

      RefPedigree::family_const_pointer  family;
      family = subp->pedigree()->family_find(p1->name(),p2->name());
     
      if(!family) // recover with error
        continue;

      size_t fam_index = family->subindex();   

      if(fam_set[fam_index] == 3) 
        continue;   //has been cut from a previous loop

      if(spouse_set.find(p2->index()) == spouse_set.end()) 
      {
        spouse_set.insert(p2->index());
        que.push_back(p2); 
      }

      if(fam_set[fam_index]++ == 1)
        loop_flag++;

      if(loop_flag >= 2)
      { 
        fam_set[pivot_fam] = 3; 
        return true; 
      }      
    }
  }

  return false;
}

// - This function is mis-named.  It actually calculates a 
//   lower_bound for ALL loops, both marriage and non-marriage.  
//   -djb  6/11/02
//
bool LoopChecker::check_nmloop(const RefSubpedigree* subp)
{
  int Np ; //# of parents with at least one offspring
  int Npp; //# of pairs of parents with at least one offspring
  int Npc; //# of individuals who have both parents and offspring
  int fnc; //# of marriages without child.

  Np=Npp=Npc=fnc=0;
  Npp = subp->family_count();
  RefPedigree::member_const_iterator i;
  RefPedigree::family_const_iterator f;

  //count # of marriages without child.  
  for(f = subp->family_begin(); f != subp->family_end(); ++f)
   if( f->offspring_count() == 0)  fnc++;

  Npp -= fnc;

  for(i = subp->member_begin(); i != subp->member_end(); ++i)
    if(i->offspring_count() > 0) 
    {
      ++Np;
      if( i->parent1() && i->parent2() )
        ++Npc;
    }

  // - Replaced '=' w. '+='.  -djb  6-11-02
  //
  my_non_marriage_loops += ( Npp - (Np-Npc-1) );

  return true;
}

} // End namespace RPED
} // End namespace SAGE
