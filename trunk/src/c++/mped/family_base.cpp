#include "mped/spbase.h"

namespace SAGE {
namespace MPED {

//============================================================================
//  IMPLEMENTATION: family_base
//============================================================================
//
family_base::family_base(subped_id SP, member_id P1, member_id P2, member_id C)
  : my_index((uint) -1), 
    my_subindex((uint) -1), 
    my_subped(SP),
    my_parent1(P1), 
    my_parent2(P2), 
    my_offspring(0), 
    my_offspring_count(0)
{
    if (P2->name() < P1->name())
    {
        std::swap(my_parent1, my_parent2);
    }

    if (C)
    {
        add_offspring(C);
    }
}


string
family_base::name() const
{
    string  fname;

    fname = my_parent1->name() + ":" + my_parent2->name();

    return fname;
}

//----------------------------------------------------------------------------
//  Function:   in_offspring()
//
//  Purpose:    This member function looks through the singly-linked list
//              of offspring to see if the argument is among them.
//----------------------------------------------------------------------------
//
bool
family_base::in_offspring(member_id C) const
{
    member_id T = my_offspring;  
    
    while (T  &&  T != C)
    {
        T = T->siblings();
    }

    return T == C;
}


//----------------------------------------------------------------------------
//  Function:   add_offspring()
//
//  Purpose:    This member function adds a new child to its singly-linked
//              list of offspring.  It makes sure all the appropriate actions
//              are taken to update the child and its parents.
//----------------------------------------------------------------------------
//
void
family_base::add_offspring(member_id C)
{
    //- We can link the child into the offspring list if it exists and
    //  is not currently a member of that list.
    //
    if (C  &&  !in_offspring(C))
    {
        //- The child must be linked to its family_base (this object).
        //
        C->set_origin(this);

        //- The child must be inserted at the head of the singly-linked list
        //  of offspring.
        //
        C->set_sibling(my_offspring);
        my_offspring = C;

        //- The number of children of each of the parents and also of this
        //  family_base must be incremented.
        //
        my_parent1->add_child();
        my_parent2->add_child();
        my_offspring_count++;
    }
}


//----------------------------------------------------------------------------
//  Function:   reset_links()
//  
//  Purpose:    This member function re-links a family to its surrounding
//              members.
//----------------------------------------------------------------------------
//
void
family_base::reset_links(subped_id SP, member_id P1, member_id P2, member_id C)
{
    my_subped = SP;

    if (P1->name() < P2->name())
    {
        my_parent1 = P1;
        my_parent2 = P2;
    }
    else
    {
        my_parent1 = P2;
        my_parent2 = P1;
    }

    P1->add_mate();
    P2->add_mate();

    my_offspring = 0;

    if (C)
    {
        add_offspring(C);
    }
}

}
}
