//============================================================================
//  File:       mpbase.h
//
//  Author:     Bob Steagall
//
//  History:    Version 0.90
//
//  Copyright (c) 1998 R.C. Elston
//  All Rights Reserved
//============================================================================
//
namespace SAGE {
namespace MPED {

inline uint 
multipedigree_base::pedigree_count() const
{
    return my_ped_index.size();
}

inline uint 
multipedigree_base::member_count() const
{
    return my_mem_mpindex.size();
}

//----------
//
inline const pedigree_base&
multipedigree_base::pedigree_index(uint i) const
{
    return *(my_ped_index[i]);
}

inline multipedigree_base::pedigree_const_iterator
multipedigree_base::pedigree_begin() const
{
    return pedigree_const_iterator(my_ped_index.begin());
}

inline multipedigree_base::pedigree_const_iterator
multipedigree_base::pedigree_end() const
{
    return pedigree_const_iterator(my_ped_index.end());
}

inline multipedigree_base::pedigree_const_iterator
multipedigree_base::pedigree_last() const
{
    return pedigree_const_iterator(my_last_iter);
}

inline const member_base&
multipedigree_base::member_index(uint i) const
{
    return *(my_mem_mpindex[i]);
}

//----------
//
inline pedigree_base&
multipedigree_base::pedigree_index(uint i)
{
    return *(my_ped_index[i]);
}

inline multipedigree_base::pedigree_iterator
multipedigree_base::pedigree_begin()
{
    return pedigree_iterator(my_ped_index.begin());
}

inline multipedigree_base::pedigree_iterator
multipedigree_base::pedigree_end()
{
    return pedigree_iterator(my_ped_index.end());
}

inline multipedigree_base::pedigree_iterator
multipedigree_base::pedigree_last()
{
    return my_last_iter;
}

inline member_base&
multipedigree_base::member_index(uint i)
{
    return *(my_mem_mpindex[i]);
}

//----------
//
inline const string&
multipedigree_base::last_name() const
{
    return my_last_name;
}

} // End namespace MPED
} // End namespace SAGE
