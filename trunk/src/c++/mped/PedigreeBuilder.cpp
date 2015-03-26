#include "mped/spbase.h"
 
namespace SAGE {
namespace MPED {

//----------------------------------------------------------------------------
//  Function:   add_lineage()
//
//  Purpose:    This member function creates a lineage meta-relationship
//              given the name of the child and the name of one parent.
//              It performs validation of all names.
//----------------------------------------------------------------------------
//
bool
PedigreeBuilder::add_lineage(const string& child, const string& parent)
{
    PedigreeBuilder::name_pair   names(parent);

    if (child.size() > 0  &&  parent.size() > 0  &&  child != parent)
    {
        add_lineage(child, name_pair(parent));
        return true;
    }
    else
    {
        my_errors.push_back(error_info(error_info::bad_lineage, child, parent));
        return false;
    }
}
//----------------------------------------------------------------------------
//  Function:   add_lineage()
//
//  Purpose:    This member function creates a lineage meta-relationship
//              given the name of the child and the names of both parents.
//              It performs validation of all names.
//----------------------------------------------------------------------------
//
bool
PedigreeBuilder::add_lineage
(const string& child, const string& parent1, const string& parent2)
{
    if (child.size() > 0  &&  parent1.size() > 0  &&  parent2.size() > 0  &&
        child != parent1  &&  child != parent2  &&  parent1 != parent2)
    {
        add_lineage(child, name_pair(parent1, parent2));
        return true;
    }
    else
    {
        my_errors.push_back(error_info(error_info::bad_lineage, child, parent1, parent2));
        return false;
    }
}

//----------------------------------------------------------------------------
//  Function:   add_marriage()
//
//  Purpose:    This member function adds a marriage meta-relationship, given
//              the names of both mates.
//----------------------------------------------------------------------------
//
bool
PedigreeBuilder::add_marriage(const string& mate1, const string& mate2)
{
    if (mate1.size() > 0  &&  mate2.size() > 0)
    {
        my_marriages.push_back( PedigreeBuilder::name_pair(mate1, mate2) );
        return true;
    }
    else
    {
        error_info  err(error_info::bad_marriage, mate1, mate2);
        my_errors.push_back(err);
        return false;
    }
}

//----------------------------------------------------------------------------
//  Function:   add_sibship()
//
//  Purpose:    This member function adds a sibling meta-relationship, given
//              the names of both siblings.
//----------------------------------------------------------------------------
//
bool
PedigreeBuilder::add_sibship(const string& sib1, const string& sib2)
{
    if (sib1.size() > 0  &&  sib2.size() > 0)
    {
        my_sibships.push_back( name_pair(sib1, sib2) );
        return true;
    }
    else
    {
        my_errors.push_back(error_info(error_info::bad_sibship, sib1, sib2));
        return false;
    }
}

//----------------------------------------------------------------------------
//  Function:   add_lineage()
//
//  Purpose:    This member function adds a lineage meta-relationship
//              given the valid name of the child and a valid pair of
//              parental names.
//
//  Notes:      While this function looks complicated, it is actually the
//              result of coding the 6 entries of the following table:
//
//               \ src      (1) <p,0>                 (2) <p,q>
//            dst \-----------------------------------------------------
//                 |                           |
//       (a) <0,0> |       dst <-- src         |      dst <-- src
//                 |                           |
//                 |----------------------------------------------------
//                 |                           |
//       (b) <j,0> |    if (j != p)            |   if (j == p  ||  j == q)
//                 |        dst <-- <j,p>      |       dst <-- <p,q>
//                 |    else                   |   else
//                 |        **do nothing**     |       **error**
//                 |                           | 
//                 |----------------------------------------------------
//                 |                           |
//       (c) <j,k> |  if (j != p  &&  k != p)  |  if (j != p  &&  k != q)
//                 |        **error**          |      **error**
//                 |                           |
//
//
//              Cases (a1) and (a2):  In both of these cases, there is no
//              existing pair of parent names in the lineage map.  Therefore
//              the src pair is just copied to the destination pair.
//
//              Case (b1):  In this case, the existing pair lists one parent,
//              and the new pair lists one parent.  If the parents are 
//              different, then the existing pair can be overwritten with
//              a new pair that lists both parents.
//
//              Case (b2):  In this case, the existing pair lists one parent,
//              and the new pair lists both parents.  If either parent in the
//              new pair is equal to the single existing parent, then the
//              existing pair can be overwritten by the new pair.  Otherwise,
//              an error has occurred.
//
//              Case (c1):  In this case, the existing pair lists both
//              parents.  If the single parent listed in the new pair is
//              not equal to the two existing parents, an error has occurred.
//
//              Case (c2):  In this case, the existing pair lists both
//              parents.  If the two parents listed in the new pair are
//              not equal to the two existing parents, an error has occurred.
//----------------------------------------------------------------------------
//
void
PedigreeBuilder::add_lineage(const string& child, const name_pair& src)
{
    //- Start by obtaining a reference to the "parents" portion of the
    //  "child/parents" association.
    //
    name_pair&  dst = my_lineages[child];

    //- Next, precalculate the sizes of the four relevant strings for
    //  convenience later.
    //
    uint    dsize1 = dst.name1.size();
    uint    dsize2 = dst.name2.size();
//    uint    ssize1 = src.name1.size(); <- *Not* actually used later, so let's not do it
    uint    ssize2 = src.name2.size(); 

    //- This portion of the "if" statement represents cases (a1) and (a2).
    //
    if (dsize1 == 0  &&  dsize2 == 0)
    {
        dst = src;
    }

    //- This portion of the "if" statement represents row (b).
    //
    else if (dsize1 != 0  &&  dsize2 == 0)
    {
        if (ssize2 == 0)
        {
            if (dst.name1 != src.name1)
            {
                dst = PedigreeBuilder::name_pair(dst.name1, src.name1);
            }
        }
        else
        {
            if (dst.name1 == src.name1  ||  dst.name1 == src.name2)
            {
                dst = src;
            }
            else
            {
                //- Error condition
            }
        }
    }

    //- This portion of the "if" statement represents row (c).
    //
    else if (dsize1 != 0  &&  dsize2 != 0)
    {
        if (ssize2 == 0)
        {
            if (dst.name1 != src.name1  &&  dst.name1 != src.name2)
            {
                //- Error condition
            }
        }
        else
        {
            if (dst.name1 != src.name1  ||  dst.name2 != src.name2)
            {
                //- Error condition
            }
        }
    }

    //- If we get to this point, something very bad has happened.
    //
    else if (dsize1 == 0  &&  dsize2 != 0)
    {
        //- Fatal error condition (data corruption)
    }
}

void
PedigreeBuilder::build_pedigree(pedigree_base& ped)
{
    process_marriages(ped);
    process_lineages(ped);
    process_sibships(ped);
    
    build_subpedigrees(ped);
    build_indices(ped);
    
    infer_sexes(ped);

    check_parental_missing_sexes(ped);

    check_parental_sex_consistency(ped);
    
    assign_arbitrary_sexes(ped);
    
    test_marriage_loop_consistency(ped);
    
    set_family_mother_father(ped);

    cleanup();
}


//----------------------------------------------------------------------------
//  Function:   process_marriage()
//
//  Purpose:    This member function processes marriage meta-relationship
//              information.
//----------------------------------------------------------------------------
//
void
PedigreeBuilder::process_marriages(pedigree_base& ped)
{
    marriage_list::iterator     mf = my_marriages.begin();
    marriage_list::iterator     ml = my_marriages.end();

    //- Iterate over each unused marriage in the list.
    //
    for (;  mf != ml;  ++mf)
    {
        if (mf->state != name_pair::unused) continue;

        //- Determine whether parents exist.
        //
        member_id   par1 = ped.lookup_member(mf->name1);
        member_id   par2 = ped.lookup_member(mf->name2);

        //- If the described family already exists, we can mark this
        //  marriage meta-info as used.
        //
        if ( ped.lookup_family(par1, par2) )
        {
            mf->state = name_pair::used;
        }

        //- Otherwise, we will try to add a new family.  In order to 
        //  create a new family, both parents must exist.  If a new
        //  family is created, then we can mark this meta-info as used.
        //
        else if (par1  &&  par2)
        {
            ped.add_family(par1, par2, 0);
            mf->state = name_pair::used;
        }
    }
}


//----------------------------------------------------------------------------
//  Function:   process_sibships()
//
//  Purpose:    This member function siblings marriage meta-relationship
//              information.
//----------------------------------------------------------------------------
//
void
PedigreeBuilder::process_sibships(pedigree_base& ped)
{
    uint    count;

    //- Because sibships are not inserted in any particular order, multiple
    //  passes over the list of sibhsip meta-info may be necessary in order
    //  process every relevant meta-relationship;
    //
    do
    {
        sibship_list::iterator      sf = my_sibships.begin();
        sibship_list::iterator      sl = my_sibships.end();

        //- Iterate over each unmarked sibship in the list.
        //
        for (count = 0;  sf != sl;  ++sf)
        {
            if (sf->state != name_pair::unused) continue;

            //- First get the member IDs of the siblings.
            //
            member_id   sib1 = ped.lookup_member(sf->name1);
            member_id   sib2 = ped.lookup_member(sf->name2);

            //- In order to link them, both siblings must exist.
            //
            if (sib1  &&  sib2)
            {
                //- Next, get the ID of each sib's family.
                //
                family_id   fam1 = sib1->family();
                family_id   fam2 = sib2->family();

                //- If sib "A" is linked to a family and sib "B" is
                //  unlinked, then sib "B" is linked to the family of 
                //  sib "A" (and vice versa).  In both cases, the 
                //  relationship meta-info is marked as used.
                //
                if (fam1  &&  !fam2)
                {
                    fam1->add_offspring(sib2);
                    sf->state = name_pair::used;
                    ++count;
                }
                else if (!fam1  &&  fam2)
                {
                    fam2->add_offspring(sib1);
                    sf->state = name_pair::used;
                    ++count;
                }

                //- In the case that both sibs are already linked to the same 
                //  family, there's nothing to do except mark the relationship 
                //  meta-info as used.
                //
                else if (fam1  &&  fam2  &&  fam1 == fam2) 
                {
                    sf->state = name_pair::used;
                    ++count;
                }

                //- Otherwise, something untoward has happened and should
                //  be reported.
                else
                {
                    my_errors.push_back(error_info(error_info::bad_sibship, sib1->name(), sib2->name()));
                }
            }
        }
    }
    while (count > 0);
}


//----------------------------------------------------------------------------
//  Function:   process_lineages()
//
//  Purpose:    This member function processes lineage meta-relationship
//              information.
//----------------------------------------------------------------------------
//
void
PedigreeBuilder::process_lineages(pedigree_base& ped)
{
    lineage_map::iterator       lf = my_lineages.begin();
    lineage_map::iterator       ll = my_lineages.end();

    //- Iterate over every unused lineage in the map.
    //
    for (;  lf != ll;  ++lf)
    {
        if (lf->second.state != name_pair::unused) continue;

        //- Use references to the data members of each association in
        //  the map.  While not strictly necessary, this will aid in
        //  debugging.
        //
        const string&   child   = lf->first;
        const string&   parent1 = lf->second.name1;
        const string&   parent2 = lf->second.name2;
        name_pair&      parents = lf->second;

        //- Perform initial lookup for relevant members.
        //
        member_id   kid1 = ped.lookup_member(child);
        member_id   par1 = ped.lookup_member(parent1);
        member_id   par2 = ped.lookup_member(parent2);

        //- At the very least, a child and one parent must be specified
        //  (and exist).
        //
        if (kid1  &&  par1)
        {
            family_id   fam1;

            //- If the second parent is given in this map association,
            //  then the lineage is considered fully specified.
            //
            if (par2)
            {
                fam1 = ped.lookup_family(parent1, parent2);

                if (fam1)
                {
                    fam1->add_offspring(kid1);
                }
                else
                {
                    fam1 = ped.add_family(par1, par2, kid1);
                }

                //- In either case, we have used the lineage meta-info, and
                //  so it can be marked as used.
                //
                parents.state = name_pair::used;
            }

            //- Otherwise, the lineage is considered partially specified.
            //  We must construct a sib-chain and look for the other
            //  parent.  If we find the other parent, then we can add the
            //  sibs to the family.
            //
            else
            {
                sib_chain   chain;

                //- If a proper sib-chain can be constructed, the ID of
                //  the second (unknown) parent will be returned.
                //
                par2 = build_sib_chain(ped, kid1, par1, chain);

                if (par2)
                {
                    //- If a family for this pair of parents does not yet
                    //  exist, then it has to be created.
                    //
                    fam1 = ped.lookup_family(par1, par2);
                    
                    if (!fam1)
                    {
                        fam1 = ped.add_family(par1, par2, kid1);
                    }

                    //- Now we iterate over the list of siblings, linking
                    //  each one to the family.
                    //
                    sib_chain::iterator     cf = chain.begin();
                    sib_chain::iterator     cl = chain.end();

                    for (;  cf != cl;  ++cf)
                    {   
                        fam1->add_offspring(*cf);
                    }

                    //- Since a sib chain was found and linked, the current 
                    //  lineage meta-info can be marked as used.
                    //
                    parents.state = name_pair::used;
                }
            }
        }
    }
}

//----------------------------------------------------------------------------
//  Function:   build_sib_chain()
//
//  Purpose:    This member function builds a chain of siblings, given
//              a sibling id and a parent id.
//----------------------------------------------------------------------------
//
member_id
PedigreeBuilder::build_sib_chain
  (pedigree_base& ped,
   member_id kid1, 
   member_id par1, 
   sib_chain& chain)
{
    sib_chain       stack;

    name_set        parents;
    member_id       par2(0);
    
    //- Start by pushing the first known sibling onto the stack.  At the end
    //  of this loop, we should have a sib-chain for 'kid1', and a list of 
    //  iterators to all meta-info entries used to create the sib-chain.
    //
    for (stack.push_back(kid1);  stack.size() > 0;)
    {
        member_id   curr_sib, new_sib;

        //- Remove the sib from the top of the stack.  Since only members
        //  that exist are placed onto the stack, we can add this sib to the
        //  sib chain.
        //
        curr_sib = stack.back();
        stack.pop_back();
        chain.push_back(curr_sib);

        sibship_list::iterator  sf = my_sibships.begin();
        sibship_list::iterator  sl = my_sibships.end();

        //- Iterate over all unused entries in the sibship meta-info list.
        //
        for (;  sf != sl;  ++sf)
        {
            if (sf->state != name_pair::unused) continue;

            //- If the first name in the pair matches the current sib's name,
            //  then we lookup the second name in the pair (and vice versa).
            //  If neither matches, then the lookup result is set to null.
            //
            if (sf->name1 == curr_sib->name())
            {
                new_sib = ped.lookup_member(sf->name2);
            }
            else if (sf->name2 == curr_sib->name())
            {
                new_sib = ped.lookup_member(sf->name1);
            }
            else
            {
                new_sib = 0;
            }

            //- A non-null lookup result means that an existing sibling was
            //  found.  We can therefore push the new sib onto the stack, 
            //  and mark the meta-info as possibly used.
            //
            if (new_sib)
            {
                stack.push_back(new_sib);
                sf->state = name_pair::possibly_used;
            }
        }
    }

    //- Next, we have to iterate over the sib-chain to determine the set of
    //  candidate parents.
    //
    sib_chain::iterator     cf = chain.begin();
    sib_chain::iterator     cl = chain.end();

    //- We know that the chain has at least one entry, and we know one parent
    //  for that entry.  We also know that the known sibling has a lineage
    //  with the known parent, and that the known sibling is the first entry
    //  in the chain.  Therefore, we can add the known parent's name to the
    //   set, and advance the initial sib-chain iterator by one.
    //
    for (parents.insert(par1->name()), ++cf;  cf != cl;  ++cf)
    {
        //- The first step is to look in the lineage map to see if a lineage
        //  entry exists for the current sibling.
        //
        lineage_map::iterator   entry = my_lineages.find((*cf)->name());

        //- If such an entry exists, then the non-empty parent names
        //  are added to the set of parents.
        //
        if (entry != my_lineages.end())
        {
            const string&   parent1 = entry->second.name1;
            const string&   parent2 = entry->second.name2;

            if (parent1.size() > 0)
            {
                parents.insert(parent1);
            }
            if (parent2.size() > 0)
            {
                parents.insert(parent2);
            }

            //- Since we have used a lineage meta-relationship, we have
            //  to mark it as possibly used.
            //
            entry->second.state = name_pair::possibly_used;
        }
    }

    //- Next, check the set of possible parents.  If the number of parent
    //  names is not equal to 2, then the sib-chain cannot be valid.  In
    //  addition, both parents must exist for the chain to be valid.
    //
    if (parents.size() == 2)
    {
        name_set::iterator  ip1 = parents.begin();
        name_set::iterator  ip2 = ++ip1;

        if ((*ip1)  ==  par1->name())
        {
            par2 = ped.lookup_member(*ip2);
        }
        else
        {
            par2 = ped.lookup_member(*ip1);
        }
    }

    //- In any event, the sibship and lineage meta-info entries must be 
    //  iterated through and any entries marked as "possibly used" must 
    //  be marked with the value given by "new_state".  The value of 
    //  "new_state" depends on the outcome of the search for a second 
    //  valid parent.
    //
    name_pair::UseStatus    new_state;

    if (par2)
    {
        new_state = name_pair::used;
    }
    else
    {
        new_state = name_pair::unused;
        chain.clear();
    }

    sibship_list::iterator  sf = my_sibships.begin();
    sibship_list::iterator  sl = my_sibships.end();

    for (;  sf != sl;  ++sf)
    {
        if (sf->state == name_pair::possibly_used)
        {
            sf->state = new_state;
        }
    }

    lineage_map::iterator   mf = my_lineages.begin();
    lineage_map::iterator   ml = my_lineages.end();

    for (;  mf != ml;  ++mf)
    {
        if (mf->second.state == name_pair::possibly_used)
        {
            mf->second.state = new_state;
        }
    }

    //- If a valid sib-chain could be constructed, then the ID of the second
    //  parent is returned.  Otherwise, a null ID is returned, indicating
    //  failure.
    //
    return par2;
}

//----------------------------------------------------------------------------
//  Function:   build_subpedigrees()
//
//  Purpose:    This member function scans the member list looking for 
//              unconnected members, and attempting to connect them to,
//              and/or build, subpedigrees.
//----------------------------------------------------------------------------
//
void
PedigreeBuilder::build_subpedigrees(pedigree_base& ped)
{
    family_id       f;
    member_id       p1, p2, s1, m1;
    subped_id       sp, old_sp;
    subped_set      spset;

    fset_iterator   ff = ped.my_families.begin();
    fset_iterator   fl = ped.my_families.end();

    //- Start by looping over all nuclear families.
    //
    for (;  ff != fl;  ++ff)
    {
        f = *ff;

        //- We're only interested in families that are still unconnected
        //  to any subpedigree.
        //
        if (f->subpedigree() == ped.my_unconnecteds)
        {
            spset.clear();
            p1 = f->parent1();
            p2 = f->parent2();
            s1 = f->offspring();
            
            //- If parents are connected, add their subpedigree ID to 
            //  the temporary set of subpedigree IDs.
            //
            if (p1->subpedigree() != ped.my_unconnecteds)
            {
                spset.insert(p1->subpedigree());
            }

            if (p2->subpedigree() != ped.my_unconnecteds)
            {
                spset.insert(p2->subpedigree());
            }

            //- If any child is connected, add its subpedigree ID to 
            //  the temporary set of subpedigree IDs.
            //
            while (s1)
            {
                if (s1->subpedigree() != ped.my_unconnecteds)
                {
                    spset.insert(s1->subpedigree());
                }
                s1 = s1->siblings();
            }

            //- The next step depends on how many different subpedigrees
            //  are encountered.  If no-one in the family is connected,
            //  then we allocate a new subpedigree and connect everyone
            //  to it.
            //
            if (spset.size() == 0)
            {
                mark_family(f, ped.add_subped());
            }

            //- If exactly one person is connected, we connect everyone
            //  in the family to this subpedigree.
            //
            else if (spset.size() == 1)
            {
                mark_family(f, *spset.begin());
            }

            //- If more than one person is connected, we connect everyone
            //  in the family to the first subpedigree, then we reset 
            //  everyone else in the other subpedegrees to belong to the
            //  first one.
            //
            else if (spset.size() >= 2)
            {
                sset_iterator   sf = spset.begin();
                sset_iterator   sl = spset.end();
                
                sp = *sf;
                mark_family(f, sp);

                for (++sf;  sf != sl;  ++sf)
                {
                    old_sp = *sf;
                    mark_all(ped, old_sp, sp);
                    ped.my_subpeds.erase(old_sp);
                    ped.deallocate_subped(old_sp);
                }
            }
        }
    }

    //- Next, we look for any members that were added to an offspring
    //  chain *after* their nuclear family was linked to a subpedigree.
    //
    mset_iterator   mf = ped.my_members.begin();
    mset_iterator   ml = ped.my_members.end();

    for (;  mf != ml;  ++mf)
    {
        m1 = *mf;

        if (m1->family())
        {
            sp = m1->family()->subpedigree();

            if (m1->subpedigree() != sp)
            {
                m1->set_subped(sp);
            }
        }
    }
}

//----------------------------------------------------------------------------
//  Function:   mark_family()
//
//  Purpose:    This member sets the subpedigree ID of the family "F" and
//              its parents and offspring to be "SP".
//----------------------------------------------------------------------------
//
void
PedigreeBuilder::mark_family(family_id f, subped_id sp)
{
    if (f)
    {
        f->set_subped(sp);
        f->parent1()->set_subped(sp);
        f->parent2()->set_subped(sp);

        member_id   s = f->offspring();

        while (s)
        {
            s->set_subped(sp);
            s = s->siblings();
        }
    }
}


//----------------------------------------------------------------------------
//  Function:   mark_all()
//
//  Purpose:    This member function sets the subpedigree ID of all members
//              and families that have the old ID to be that of the new ID.
//----------------------------------------------------------------------------
//
void
PedigreeBuilder::mark_all(pedigree_base& ped, subped_id old_id, subped_id new_id)
{
    //- Iterate through the set of all families, updating any with the
    //  old subpedigree ID to the new subpedigree ID.
    //
    fset_iterator   ff = ped.my_families.begin();
    fset_iterator   fl = ped.my_families.end();

    for (;  ff != fl;  ++ff)
    {
        if ((*ff)->subpedigree() == old_id)
        {
            (*ff)->set_subped(new_id);
        }
    }

    //- Iterate through the set of all members, updating any with the
    //  old subpedigree ID to the new subpedigree ID.
    //
    mset_iterator   mf = ped.my_members.begin();
    mset_iterator   ml = ped.my_members.end();

    for (;  mf != ml;  ++mf)
    {
        if ((*mf)->subpedigree() == old_id)
        {
            (*mf)->set_subped(new_id);
        }
    }
}

//----------------------------------------------------------------------------
//  Function:   build_indices()
//
//  Purpose:    This member function builds the vector of pointers used
//              for indexing into the pedigree.
//----------------------------------------------------------------------------
//
void
PedigreeBuilder::build_indices(pedigree_base& ped)
{
    uint    i;

    //- Populate the subpedigree index array.  Note that the member
    //  and family index arrays of each subpedigree are cleared in
    //  this loop.
    //
    sset_iterator   sf = ped.my_subpeds.begin();
    sset_iterator   sl = ped.my_subpeds.end();

    i = ped.my_subped_index.size();
    ped.my_subped_index.resize(ped.my_subpeds.size());

    for (;  sf != sl;  ++sf)
    {
        if( (*sf)->index() != (uint) -1 )
          continue;

        ped.my_subped_index[i] = *sf;
        (*sf)->set_index(i++);
//        (*sf)->clear_index_arrays();
    }

    //- Populate the family index array.  Note that each family is
    //  added to the relevant subpedigree as the set of families is
    //  being traversed.
    //
    fset_iterator   ff = ped.my_families.begin();
    fset_iterator   fl = ped.my_families.end();

    i = ped.my_family_index.size();
    ped.my_family_index.resize(ped.my_families.size());

    for (;  ff != fl;  ++ff)
    {
        if( (*ff)->index() != (uint) -1 )
          continue;

        ped.my_family_index[i] = *ff;
        (*ff)->set_index(i++);
        (*ff)->subpedigree()->add_family(*ff);
    }

    //- Populate the member index array.  Note that each member is
    //  added to the relevant subpedigree as the set of members is
    //  being traversed.
    //
    mset_iterator   mf = ped.my_members.begin();
    mset_iterator   ml = ped.my_members.end();

    i = ped.my_member_index.size();
    ped.my_member_index.resize(ped.my_members.size());

    for (;  mf != ml;  ++mf)
    {
        if( (*mf)->index() != (uint) -1 )
          continue;

        ped.my_member_index[i] = *mf;
        (*mf)->set_index(i++);
        (*mf)->subpedigree()->add_member(*mf);
    }

    //- Now that the pedigree index arrays have been built,
    //  the subpedigrees must each build their own index arrays.
    //

    for (sf = ped.my_subpeds.begin(), sl = ped.my_subpeds.end();  sf != sl;  ++sf)
    {
        (*sf)->build_indices();
    }

    //- The unconnecteds subpedigree must now remove all members and
    //  and families that no longer belong to it.
    //

    ped.my_unconnecteds->cleanup();

    //- Finally, set the overall counts.
    //
    ped.my_fcount = ped.my_families.size();
    ped.my_mcount = ped.my_members.size();
    ped.my_scount = ped.my_subpeds.size();
    ped.my_ucount = ped.my_unconnecteds->member_count();
}


/// Finds matings where only one of the mating pair is sexed and infers the other.
/// This algorithm repeats this process until there are no more matings exist
/// with this condition.
void
PedigreeBuilder::infer_sexes(pedigree_base& ped)
{
  // Collect all parental pairs with at least one parent who has unknown sex.
  typedef std::list<std::pair<member_id, member_id> > ParentPairContainer;

  ParentPairContainer parent_pairs;
  
  for(family_set::iterator fam = ped.my_families.begin(); fam != ped.my_families.end(); ++fam)
  {
    member_id p1 = (*fam)->parent1();
    member_id p2 = (*fam)->parent2();
    
    if(p1->is_sex_unknown() || p2->is_sex_unknown())
    {
      parent_pairs.push_back(std::make_pair(p1,p2));
    }
  }
  
  // Keep track of changes.  If there's no change for a loop, we're done
  bool changed = true;
    
  // While ther might be inferences to make.
  while(changed && parent_pairs.size())
  {
    changed = false;
    
    // Go through the pairs, looking for parents who are sexed for inference
    for(ParentPairContainer::iterator pp = parent_pairs.begin(); pp != parent_pairs.end(); )
    {
      if(pp->first->is_sex_unknown() && pp->second->is_sex_unknown())
      {
        // Do nothing but increment the iterator, both parents unsexed
        ++pp;
      }
      else if(pp->first->is_sex_unknown())
      {
        // Second must be either male or female
        if(pp->second->is_female())
        {
          if(pp->first->my_gender == SEX_MISSING)
            my_warnings.push_back(error_info(error_info::gender_inferred, pp->first->name(), "male"));
          pp->first->set_sex(SEX_XMALE);
        }
        else
        {
          if(pp->first->my_gender == SEX_MISSING)
            my_warnings.push_back(error_info(error_info::gender_inferred, pp->first->name(), "female"));
          pp->first->set_sex(SEX_XFEMALE);
        }
        
        // Remove this pair, done
        changed = true;
        parent_pairs.erase(pp++);
      }
      else if(pp->second->is_sex_unknown())
      {
        // First must be either male or female
        if(pp->first->is_female())
        {
          if(pp->second->my_gender == SEX_MISSING)
            my_warnings.push_back(error_info(error_info::gender_inferred, pp->second->name(), "male"));
          pp->second->set_sex(SEX_XMALE);
        }
        else
        {
          if(pp->second->my_gender == SEX_MISSING)
            my_warnings.push_back(error_info(error_info::gender_inferred, pp->second->name(), "female"));
          pp->second->set_sex(SEX_XFEMALE);
        }

        // Remove this pair, done
        changed = true;
        parent_pairs.erase(pp++);
      }
      else // Both must be sexed
      {
        // Remove this pair, done.  Note that this is not 'changed', as it
        // didn't make inferences.
        //
        // Note we don't care about parents being the same sex at this stage.
        parent_pairs.erase(pp++);
      }
    }
  }
}

/// Looks for marriages where both members are the same sex and reports them to the
/// error handler.
void
PedigreeBuilder::check_parental_sex_consistency(pedigree_base& ped)
{
  // Look for parents with the same sex to report errors
  for(family_base_iterator fam = ped.family_begin(); fam != ped.family_end(); ++fam)
  {
    member_id p1 = fam->parent1();
    member_id p2 = fam->parent2();
    
    if(p1->is_male()   && p2->is_male())
    {
        my_errors.push_back(error_info(error_info::same_sex_marriage, p1->name(), p2->name(), "male"));
    }
    else if (p1->is_female() && p2->is_female())
    {
        my_errors.push_back(error_info(error_info::same_sex_marriage, p1->name(), p2->name(), "female"));
    }
  }
}

/// Finds matings where both parents are unsexed and gives errors to the 
/// error handler.
void
PedigreeBuilder::check_parental_missing_sexes(pedigree_base& ped)
{
  // Look for parent pairs which are both still missing to report them as such
  for(family_base_iterator fam = ped.family_begin(); fam != ped.family_end(); ++fam)
  {
    member_id p1 = fam->parent1();
    member_id p2 = fam->parent2();

    if(p1->my_gender == SEX_MISSING && p2->my_gender == SEX_MISSING)
    {
      my_errors.push_back(error_info(error_info::no_sex_parents, p1->name(), p2->name()));
    }
  }
}

/// Attempts to assign internally consistent sexes to individuals with the SEX_ARB
/// sex code.  These sexes are propagated to individuals connected in the marriage
/// chain.  If no consistent sexing can be performed, sexes are set to SEX_MISSING.
/// Errors will be generated later by the check_marriage_loop_consistency() function.
///
/// SEX_ARB individuals can only be mated to other SEX_ARB individuals or SEX_MISSING
/// individuals when this operation is performed, since sex inference should already 
/// be done.  If a SEX_MISSING is found, all
/// individuals in the marriage chain are assigned SEX_MISSING.  
void
PedigreeBuilder::assign_arbitrary_sexes(pedigree_base& ped)
{
  // Create a vector of pairs of bools and sexes.  This will keep track of
  // who we visited and what their sex would be to maintain consistency.
  std::vector<std::pair<bool, SexCode> > visit_vector(ped.member_count(), std::make_pair(false, SEX_MISSING));
  
  // Look for individuals with arbitrary sex
  for(member_base_iterator i = ped.member_begin(); i != ped.member_end(); ++i)
  {
    if(i->my_gender == SEX_ARB)
    {
      // Attempt to make them female, but assigning thier sex (and recursively
      // their mates and so on).  We just choose female here arbitrarily
      if(assign_arbitrary_sex(&*i, 0, SEX_XFEMALE, visit_vector))
      {
        // This assignment was successful, so go through the visited vector and
        // assign everyone in it to the sexes they have listed
        for(size_t vv_index = 0; vv_index != visit_vector.size(); ++vv_index)
        {
          if(visit_vector[vv_index].first)
            ped.member_index(vv_index).set_sex(visit_vector[vv_index].second);
        }
      }
      else
      {
        // This assignment was not successful, so go through the visited vector 
        // and assign everyone in it to missing.
        for(size_t vv_index = 0; vv_index != visit_vector.size(); ++vv_index)
        {
          if(visit_vector[vv_index].first)
            ped.member_index(vv_index).set_sex(SEX_MISSING);
        }
      }
      
      // Clear the visit vector for the next chain
      visit_vector.assign(ped.member_count(), std::make_pair(false, SEX_MISSING));
    }
  }
}

/// Attempts to assign \c sex to \c ind and propagate this to any mates of the individual.
///
/// \return \c true if the assignment works, \c false if there are problems.
bool
PedigreeBuilder::assign_arbitrary_sex(member_id ind, member_id mate, SexCode sex,
                                    std::vector<std::pair<bool, SexCode> >& visit_vector)
{
  // If there's a missing sex individual in this marriage chain, we have to label
  // everyone missing.
  if(ind->my_gender == SEX_MISSING) return false;

  // If the ind has been previously visited, we must check their previously
  // arbitrarily assigned sex versus the current one.
  if(visit_vector[ind->index()].first)
  {
    // If the sex is not the same, there is a conflict, so we can't consistently
    // assign the arbitrary sexes.  Assign everyone missing.  Error generated when
    // marriage loop checking performed.
    if(visit_vector[ind->index()].second != sex)
    {
      return false;
    }
  
    // Sexes match, and we've previously seen this one, so we can terminate
    // successfully.
    
    return true;
  }

  // At this point, we have established that the ind has not been previously 
  // visited and they have arbitrary sex.
  
  // Set them visited and their sex to sex
  visit_vector[ind->index()].first  = true;
  visit_vector[ind->index()].second = sex;  
  
  SexCode mate_sex = (sex == SEX_XMALE) ? SEX_XFEMALE : SEX_XMALE;
  
  // Iterate through all the mates of ind (except mate) and see if they can
  // be assigned to mate_sex consistently.  If they can't, we'll have to 
  // return an error.
  for(mate_base_iterator m2 = ind->mate_begin(); m2 != ind->mate_end(); ++m2)
  {
    if(&m2->mate() == mate) continue;
    
    if(!assign_arbitrary_sex(&m2->mate(), ind, mate_sex, visit_vector))
      return false;
  }
  
  // We're ok, and our mates (and the chain off of them) are as well, so we're
  // successful
  return true;
}

/// Visits marriage pairs of SEX_MISSING members, looking for marriage loops 
/// that might indicate sex inconsistencies.  A marriage loop with an odd individual
/// count will cause these errors.  Any errors detected are given to the
/// error handler.
void
PedigreeBuilder::test_marriage_loop_consistency(pedigree_base& ped)
{
  // Create a vector of pairs of bools and sexes.  This will keep track of
  // who we visited and what their sex would be to maintain consistency.
  std::vector<std::pair<bool, SexCode> > visit_vector(ped.member_count(), std::make_pair(false, SEX_MISSING));

  std::vector<bool>                      checked_status(ped.member_count(), false);
  
  // Look for individuals with arbitrary sex
  for(member_base_iterator i = ped.member_begin(); i != ped.member_end(); ++i)
  {
    // If we're unknown and haven't been checked
    if(i->is_sex_unknown() && !checked_status[i->index()])
    {
      // Attempt to label them female, but assigning thier sex (and recursively
      // their mates and so on).  We just choose female here arbitrarily
      test_marriage_loop_consistency(&*i, 0, SEX_XFEMALE, visit_vector);
      
      for(size_t vv_index = 0; vv_index != visit_vector.size(); ++vv_index)
      {
        if(visit_vector[vv_index].first)
          checked_status[vv_index] = true;
      }
      
      // Clear the visit vector for the next chain
      visit_vector.assign(ped.member_count(), std::make_pair(false, SEX_MISSING));
    }
  }
}

/// Tests an individual for consistency given a particular sex.  Note that
/// this sex is only a placeholder, and is not ever given to the individual.
/// This assignment is then propagated to mates of the individual.
///
/// \return \c true if the assignment works, \c false if there are problems.
bool
PedigreeBuilder::test_marriage_loop_consistency
    (member_id ind,
     member_id mate,
     SexCode sex,
     std::vector<std::pair<bool, SexCode> >& visit_vector)
{
  // If the ind has been previously visited, we must check their previously
  // arbitrarily assigned sex versus the current one.
  if(visit_vector[ind->index()].first)
  {
    // If the sex is not the same, there is a conflict, so we can't consistently
    // assign the arbitrary sexes.  Generate an error
    if(visit_vector[ind->index()].second != sex)
    {
      my_errors.push_back(error_info(error_info::bad_marriage_loop, ind->name()));
      return false;
    }
  
    // Sexes match, and we've previously seen this one, so we can terminate
    // successfully.
    
    return true;
  }

  // At this point, we have established that the ind has not been previously 
  // visited
  
  // Set them visited and their sex to sex
  visit_vector[ind->index()].first  = true;
  visit_vector[ind->index()].second = sex;  
  
  SexCode mate_sex = (sex == SEX_XMALE) ? SEX_XFEMALE : SEX_XMALE;
  
  // Iterate through all the mates of ind (except mate) and see if they can
  // be assigned to mate_sex consistently.  If they can't, we'll have to 
  // return an error.
  for(mate_base_iterator m2 = ind->mate_begin(); m2 != ind->mate_end(); ++m2)
  {
    if(&m2->mate() == mate) continue;
    
    if(!test_marriage_loop_consistency(&m2->mate(), ind, mate_sex, visit_vector))
      return false;
  }
  
  // We're ok, and our mates (and the chain off of them) are as well, so we're
  // successful
  return true;
}

void
PedigreeBuilder::set_family_mother_father(pedigree_base& ped)
{
  for(family_base_iterator fam = ped.family_begin(); fam != ped.family_end(); ++fam)
  {
    if(fam->parent1()->is_male() && fam->parent2()->is_female())
    {
      fam->my_mother = fam->parent2();
      fam->my_father = fam->parent1();
    }
    else if(fam->parent1()->is_female() && fam->parent2()->is_male())
    {
      fam->my_mother = fam->parent1();
      fam->my_father = fam->parent2();
    }
  }
}

//----------------------------------------------------------------------------
//  Function:   flush_build_info()
//
//  Purpose:    This member function deletes all meta-relationship information
//              and clears all STL containers used to hold and/or represent 
//              this information.
//----------------------------------------------------------------------------
//
void
PedigreeBuilder::flush_build_info()
{
    my_marriages.clear();
    my_lineages.clear();
    my_sibships.clear();
    my_errors.clear();
    my_warnings.clear();
}

//----------------------------------------------------------------------------
//  Function:   cleanup()
//
//  Purpose:    This very simple member function removes all meta-info 
//              that has been marked as 'used'.
//----------------------------------------------------------------------------
//
void
PedigreeBuilder::cleanup()
{
    PedigreeBuilder::marriage_list::iterator     mf = my_marriages.begin();
    PedigreeBuilder::marriage_list::iterator     ml = my_marriages.end();

    while (mf != ml)
    {
        if (mf->state == PedigreeBuilder::name_pair::used)
        {
            my_marriages.erase(mf++);
        }
        else
        {
            ++mf;
        }
    }

    PedigreeBuilder::lineage_map::iterator       lf = my_lineages.begin();
    PedigreeBuilder::lineage_map::iterator       ll = my_lineages.end();

    while (lf != ll)
    {
        if (lf->second.state == PedigreeBuilder::name_pair::used)
        {
            my_lineages.erase(lf++);
        }
        else
        {
            ++lf;
        }
    }

    PedigreeBuilder::sibship_list::iterator      sf = my_sibships.begin();
    PedigreeBuilder::sibship_list::iterator      sl = my_sibships.end();

    while (sf != sl)
    {
        if (sf->state == PedigreeBuilder::name_pair::used)
        {
            my_sibships.erase(sf++);
        }
        else
        {
            ++sf;
        }
    }
}

}
}

