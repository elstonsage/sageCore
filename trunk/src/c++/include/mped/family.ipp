namespace SAGE {
namespace MPED {

//============================================================================
//  IMPLEMENTATION: family<GI,FI,SI,PI,MI>
//============================================================================
//
template <class GI, class FI, class SI, class PI, class MI> inline uint
family<GI,FI,SI,PI,MI>::index() const
{
    return base::index();
}

template <class GI, class FI, class SI, class PI, class MI> inline uint
family<GI,FI,SI,PI,MI>::subindex() const
{
    return base::subindex();
}

template <class GI, class FI, class SI, class PI, class MI> inline uint
family<GI,FI,SI,PI,MI>::offspring_count() const
{
    return base::offspring_count();
}

template <class GI, class FI, class SI, class PI, class MI> inline const FI&
family<GI,FI,SI,PI,MI>::info() const
{
    return my_info;
}

template <class GI, class FI, class SI, class PI, class MI> inline string
family<GI,FI,SI,PI,MI>::name() const
{
    return base::name();
}

template <class GI, class FI, class SI, class PI, class MI> inline 
const string&
family<GI,FI,SI,PI,MI>::name1() const
{
    return base::name1();
}

template <class GI, class FI, class SI, class PI, class MI> inline 
const string&
family<GI,FI,SI,PI,MI>::name2() const
{
    return base::name2();
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename family<GI,FI,SI,PI,MI>::subpedigree_const_pointer
family<GI,FI,SI,PI,MI>::subpedigree() const
{
    return reinterpret_cast<subpedigree_const_pointer>(base::subpedigree());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename family<GI,FI,SI,PI,MI>::pedigree_const_pointer
family<GI,FI,SI,PI,MI>::pedigree() const
{
    return reinterpret_cast<pedigree_const_pointer>(base::pedigree());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename family<GI,FI,SI,PI,MI>::multipedigree_const_pointer
family<GI,FI,SI,PI,MI>::multipedigree() const
{
    return reinterpret_cast<multipedigree_const_pointer>(base::multipedigree());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename family<GI,FI,SI,PI,MI>::member_const_pointer
family<GI,FI,SI,PI,MI>::parent1() const
{
    return reinterpret_cast<member_const_pointer>(base::parent1());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename family<GI,FI,SI,PI,MI>::member_const_pointer
family<GI,FI,SI,PI,MI>::parent2() const
{
    return reinterpret_cast<member_const_pointer>(base::parent2());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename family<GI,FI,SI,PI,MI>::member_const_pointer
family<GI,FI,SI,PI,MI>::get_mother() const
{
    return reinterpret_cast<member_const_pointer>(base::get_mother());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename family<GI,FI,SI,PI,MI>::member_const_pointer
family<GI,FI,SI,PI,MI>::get_father() const
{
    return reinterpret_cast<member_const_pointer>(base::get_father());
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename family<GI,FI,SI,PI,MI>::offspring_const_iterator
family<GI,FI,SI,PI,MI>::offspring_begin() const
{
    return offspring_const_iterator(base::offspring_begin());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename family<GI,FI,SI,PI,MI>::offspring_const_iterator
family<GI,FI,SI,PI,MI>::offspring_end() const
{
    return offspring_const_iterator(base::offspring_end());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename family<GI,FI,SI,PI,MI>::parent_const_iterator
family<GI,FI,SI,PI,MI>::parent_begin() const
{
    return parent_const_iterator(base::parent_begin());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename family<GI,FI,SI,PI,MI>::parent_const_iterator
family<GI,FI,SI,PI,MI>::parent_end() const
{
    return parent_const_iterator(base::parent_end());
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline FI&
family<GI,FI,SI,PI,MI>::info()
{
    return my_info;
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename family<GI,FI,SI,PI,MI>::subpedigree_pointer
family<GI,FI,SI,PI,MI>::subpedigree()
{
    return reinterpret_cast<subpedigree_pointer>(base::subpedigree());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename family<GI,FI,SI,PI,MI>::pedigree_pointer
family<GI,FI,SI,PI,MI>::pedigree()
{
    return reinterpret_cast<pedigree_pointer>(base::pedigree());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename family<GI,FI,SI,PI,MI>::multipedigree_pointer
family<GI,FI,SI,PI,MI>::multipedigree()
{
    return reinterpret_cast<multipedigree_pointer>(base::multipedigree());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename family<GI,FI,SI,PI,MI>::member_pointer
family<GI,FI,SI,PI,MI>::parent1()
{
    return reinterpret_cast<member_pointer>(base::parent1());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename family<GI,FI,SI,PI,MI>::member_pointer
family<GI,FI,SI,PI,MI>::parent2() 
{
    return reinterpret_cast<member_pointer>(base::parent2());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename family<GI,FI,SI,PI,MI>::member_pointer
family<GI,FI,SI,PI,MI>::get_mother()
{
    return reinterpret_cast<member_pointer>(base::get_mother());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename family<GI,FI,SI,PI,MI>::member_pointer
family<GI,FI,SI,PI,MI>::get_father() 
{
    return reinterpret_cast<member_pointer>(base::get_father());
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename family<GI,FI,SI,PI,MI>::offspring_iterator
family<GI,FI,SI,PI,MI>::offspring_begin()
{
    return offspring_iterator(base::offspring_begin());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename family<GI,FI,SI,PI,MI>::offspring_iterator
family<GI,FI,SI,PI,MI>::offspring_end()
{
    return offspring_iterator(base::offspring_end());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename family<GI,FI,SI,PI,MI>::parent_iterator
family<GI,FI,SI,PI,MI>::parent_begin()
{
    return parent_iterator(base::parent_begin());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename family<GI,FI,SI,PI,MI>::parent_iterator
family<GI,FI,SI,PI,MI>::parent_end()
{
    return parent_iterator(base::parent_end());
}

}
}

