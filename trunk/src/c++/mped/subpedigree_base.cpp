#include "mped/spbase.h"

namespace SAGE {
namespace MPED {

//============================================================================
//  IMPLEMENTATION: subpedigree_base
//============================================================================
//
subpedigree_base::subpedigree_base()
  : my_index((uint) -1), my_name(), my_pedigree(0)
{}


subpedigree_base::subpedigree_base
(const string& name, pedigree_id pped)
  : my_index((uint) -1), my_name(name), my_pedigree(pped)
{}


subpedigree_base::subpedigree_base(const string& name, const no_info&)
  : my_index((uint) -1), my_name(name), my_pedigree(0)
{}


//----------
//
void
subpedigree_base::family_index_swap(uint i, uint j)
{
    std::swap(my_family_index[i], my_family_index[j]);
    my_family_index[i]->set_subindex(i);
    my_family_index[j]->set_subindex(j);
}


void
subpedigree_base::member_index_swap(uint i, uint j)
{
    std::swap(my_member_index[i], my_member_index[j]);
    my_member_index[i]->set_subindex(i);
    my_member_index[j]->set_subindex(j);
}


void
subpedigree_base::add_family(family_id f)
{
    if ( my_families.find(f) != my_families.end() )
        return;

    uint    size = my_family_index.size();
    uint    capy = my_family_index.capacity();

    if (size == 0)
    {
        my_family_index.reserve(8);
    }
    else if (size >= (capy - 1))
    {
        my_family_index.reserve(2 * capy);
    }

    my_families.insert(f);
    my_family_index.push_back(f);
}


void
subpedigree_base::add_member(member_id m)
{
    if( my_members.find(m) != my_members.end() )
      return;

    uint    size = my_member_index.size();
    uint    capy = my_member_index.capacity();

    if (size == 0)
    {
        my_member_index.reserve(16);
    }
    else if (size >= (capy - 1))
    {
        my_member_index.reserve(2 * capy);
    }

    my_members.insert(m);
    my_member_index.push_back(m);
}


void
subpedigree_base::build_indices()
{
    uint    i;

    for (i = 0;  i < my_member_index.size();  ++i)
    {
        my_member_index[i]->set_subindex(i);
    }

    for (i = 0;  i < my_family_index.size();  ++i)
    {
        my_family_index[i]->set_subindex(i);
    }
}


void
subpedigree_base::cleanup()
{
    uint        i;
    subped_id   this_id = this;
    fidx_vector new_families;
    midx_vector new_members;

    for (i = 0;  i < my_family_index.size();  ++i)
    {
        if (my_family_index[i]->subpedigree() == this_id)
        {
            new_families.push_back(my_family_index[i]);
        }
    }
    my_family_index.swap(new_families);

    for (i = 0;  i < my_member_index.size();  ++i)
    {
        if (my_member_index[i]->subpedigree() == this_id)
        {
            new_members.push_back(my_member_index[i]);
        }
    }
    my_member_index.swap(new_members);
}

}
}
