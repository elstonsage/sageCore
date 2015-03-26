namespace SAGE {
namespace MPED {

//============================================================================
//  IMPLEMENTATION: subpedigree<GI,FI,SI,PI,MI>
//============================================================================
//
template <class GI, class FI, class SI, class PI, class MI> inline
subpedigree<GI,FI,SI,PI,MI>::subpedigree
(const string& name, pedigree_base* pped)
  : subpedigree_base(name, pped)
{}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline uint
subpedigree<GI,FI,SI,PI,MI>::index() const
{
    return base::index();
}

template <class GI, class FI, class SI, class PI, class MI> inline uint
subpedigree<GI,FI,SI,PI,MI>::family_count() const
{
    return base::family_count();
}

template <class GI, class FI, class SI, class PI, class MI> inline uint
subpedigree<GI,FI,SI,PI,MI>::member_count() const
{
    return base::member_count();
}

template <class GI, class FI, class SI, class PI, class MI> inline const SI&
subpedigree<GI,FI,SI,PI,MI>::info() const
{
    return my_info;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
const string&
subpedigree<GI,FI,SI,PI,MI>::name() const
{
    return base::name();
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename subpedigree<GI,FI,SI,PI,MI>::pedigree_const_pointer
subpedigree<GI,FI,SI,PI,MI>::pedigree() const
{
    return reinterpret_cast<pedigree_const_pointer>(base::pedigree());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename subpedigree<GI,FI,SI,PI,MI>::multipedigree_const_pointer
subpedigree<GI,FI,SI,PI,MI>::multipedigree() const
{
    return reinterpret_cast<multipedigree_const_pointer>(base::multipedigree());
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline
const typename subpedigree<GI,FI,SI,PI,MI>::family_type&
subpedigree<GI,FI,SI,PI,MI>::family_index(uint i) const
{
    return reinterpret_cast<const family_type&>(base::family_index(i));
}

template <class GI, class FI, class SI, class PI, class MI> inline
const typename subpedigree<GI,FI,SI,PI,MI>::member_type&
subpedigree<GI,FI,SI,PI,MI>::member_index(uint i) const
{
    return reinterpret_cast<const member_type&>(base::member_index(i));
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename subpedigree<GI,FI,SI,PI,MI>::family_const_iterator
subpedigree<GI,FI,SI,PI,MI>::family_begin() const
{
    return family_const_iterator(base::family_begin());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename subpedigree<GI,FI,SI,PI,MI>::family_const_iterator
subpedigree<GI,FI,SI,PI,MI>::family_end() const
{
    return family_const_iterator(base::family_end());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename subpedigree<GI,FI,SI,PI,MI>::member_const_iterator
subpedigree<GI,FI,SI,PI,MI>::member_begin() const
{
    return member_const_iterator(base::member_begin());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename subpedigree<GI,FI,SI,PI,MI>::member_const_iterator
subpedigree<GI,FI,SI,PI,MI>::member_end() const
{
    return member_const_iterator(base::member_end());
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline SI&
subpedigree<GI,FI,SI,PI,MI>::info()
{
    return my_info;
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename subpedigree<GI,FI,SI,PI,MI>::pedigree_pointer
subpedigree<GI,FI,SI,PI,MI>::pedigree()
{
    return reinterpret_cast<pedigree_pointer>(base::pedigree());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename subpedigree<GI,FI,SI,PI,MI>::multipedigree_pointer
subpedigree<GI,FI,SI,PI,MI>::multipedigree()
{
    return reinterpret_cast<multipedigree_pointer>(base::multipedigree());
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename subpedigree<GI,FI,SI,PI,MI>::family_type&
subpedigree<GI,FI,SI,PI,MI>::family_index(uint i)
{
    return reinterpret_cast<family_type&>(base::family_index(i));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename subpedigree<GI,FI,SI,PI,MI>::member_type&
subpedigree<GI,FI,SI,PI,MI>::member_index(uint i)
{
    return reinterpret_cast<member_type&>(base::member_index(i));
}

template <class GI, class FI, class SI, class PI, class MI> inline void
subpedigree<GI,FI,SI,PI,MI>::family_index_swap(uint i, uint j)
{
    base::family_index_swap(i, j);
}

template <class GI, class FI, class SI, class PI, class MI> inline void
subpedigree<GI,FI,SI,PI,MI>::member_index_swap(uint i, uint j)
{
    base::member_index_swap(i, j);
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename subpedigree<GI,FI,SI,PI,MI>::family_iterator
subpedigree<GI,FI,SI,PI,MI>::family_begin()
{
    return family_iterator(base::family_begin());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename subpedigree<GI,FI,SI,PI,MI>::family_iterator
subpedigree<GI,FI,SI,PI,MI>::family_end()
{
    return family_iterator(base::family_end());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename subpedigree<GI,FI,SI,PI,MI>::member_iterator
subpedigree<GI,FI,SI,PI,MI>::member_begin()
{
    return member_iterator(base::member_begin());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename subpedigree<GI,FI,SI,PI,MI>::member_iterator
subpedigree<GI,FI,SI,PI,MI>::member_end()
{
    return member_iterator(base::member_end());
}

}
}
