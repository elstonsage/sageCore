
namespace SAGE {
namespace MPED {

template <class GI, class FI, class SI, class PI, class MI> inline
member<GI,FI,SI,PI,MI>::member(const string& name, SexCode gender, const geninfo_type& g)
  : member_base(name, gender), my_info(g)
{}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline SexCode
member<GI,FI,SI,PI,MI>::get_effective_sex() const
{
    return base::get_effective_sex();
}

template <class GI, class FI, class SI, class PI, class MI> inline SexCode
member<GI,FI,SI,PI,MI>::get_detailed_sex() const
{
    return base::get_detailed_sex();
}

template <class GI, class FI, class SI, class PI, class MI> inline bool
member<GI,FI,SI,PI,MI>::is_male() const
{
    return base::is_male();
}

template <class GI, class FI, class SI, class PI, class MI> inline bool
member<GI,FI,SI,PI,MI>::is_female() const
{
    return base::is_female();
}

template <class GI, class FI, class SI, class PI, class MI> inline uint
member<GI,FI,SI,PI,MI>::mpindex() const
{
    return base::mpindex();
}

template <class GI, class FI, class SI, class PI, class MI> inline uint
member<GI,FI,SI,PI,MI>::index() const
{
    return base::index();
}

template <class GI, class FI, class SI, class PI, class MI> inline uint
member<GI,FI,SI,PI,MI>::subindex() const
{
    return base::subindex();
}

template <class GI, class FI, class SI, class PI, class MI> inline uint
member<GI,FI,SI,PI,MI>::mate_count() const
{
    return base::mate_count();
}

template <class GI, class FI, class SI, class PI, class MI> inline uint
member<GI,FI,SI,PI,MI>::offspring_count() const
{
    return base::offspring_count();
}

template <class GI, class FI, class SI, class PI, class MI> inline uint
member<GI,FI,SI,PI,MI>::sibling_count() const
{
    return base::sibling_count();
}

template <class GI, class FI, class SI, class PI, class MI> inline const GI&
member<GI,FI,SI,PI,MI>::info() const
{
    return my_info;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
const string&
member<GI,FI,SI,PI,MI>::name() const
{
    return base::name();
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::subpedigree_const_pointer
member<GI,FI,SI,PI,MI>::subpedigree() const
{
    return reinterpret_cast<subpedigree_const_pointer>(base::subpedigree());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::pedigree_const_pointer
member<GI,FI,SI,PI,MI>::pedigree() const
{
    return reinterpret_cast<pedigree_const_pointer>(base::pedigree());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::multipedigree_const_pointer
member<GI,FI,SI,PI,MI>::multipedigree() const
{
    return reinterpret_cast<multipedigree_const_pointer>(base::multipedigree());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::family_const_pointer
member<GI,FI,SI,PI,MI>::family() const
{
    return reinterpret_cast<family_const_pointer>(base::family());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::member_const_pointer
member<GI,FI,SI,PI,MI>::parent1() const
{
    return reinterpret_cast<member_const_pointer>(base::parent1());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::member_const_pointer
member<GI,FI,SI,PI,MI>::parent2() const
{
    return reinterpret_cast<member_const_pointer>(base::parent2());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::member_const_pointer
member<GI,FI,SI,PI,MI>::get_mother() const
{
    return reinterpret_cast<member_const_pointer>(base::get_mother());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::member_const_pointer
member<GI,FI,SI,PI,MI>::get_father() const
{
    return reinterpret_cast<member_const_pointer>(base::get_father());
}


//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::mate_const_iterator
member<GI,FI,SI,PI,MI>::mate_begin() const
{
    return mate_const_iterator(base::mate_begin());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::mate_const_iterator
member<GI,FI,SI,PI,MI>::mate_end() const
{
    return mate_const_iterator(base::mate_end());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::offspring_const_iterator
member<GI,FI,SI,PI,MI>::offspring_begin(const mate_const_iterator& m) const
{
    return offspring_const_iterator(base::offspring_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::offspring_const_iterator
member<GI,FI,SI,PI,MI>::offspring_begin(const member_const_iterator& m) const
{
    return offspring_const_iterator(base::offspring_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::offspring_const_iterator
member<GI,FI,SI,PI,MI>::offspring_begin(const member& m) const
{
    return offspring_const_iterator(base::offspring_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::offspring_const_iterator
member<GI,FI,SI,PI,MI>::offspring_begin(const string& m) const
{
    return offspring_const_iterator(base::offspring_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::offspring_const_iterator
member<GI,FI,SI,PI,MI>::offspring_end() const
{
    return offspring_const_iterator(base::offspring_end());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::parent_const_iterator
member<GI,FI,SI,PI,MI>::parent_begin() const
{
    return parent_const_iterator(base::parent_begin());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::parent_const_iterator
member<GI,FI,SI,PI,MI>::parent_end() const
{
    return parent_const_iterator(base::parent_end());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::progeny_const_iterator
member<GI,FI,SI,PI,MI>::progeny_begin() const
{
    return progeny_const_iterator(base::progeny_begin());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::progeny_const_iterator
member<GI,FI,SI,PI,MI>::progeny_end() const
{
    return progeny_const_iterator(base::progeny_end());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::sibling_const_iterator
member<GI,FI,SI,PI,MI>::sibling_begin() const
{
    return sibling_const_iterator(base::sibling_begin());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::sibling_const_iterator
member<GI,FI,SI,PI,MI>::sibling_end() const
{
    return sibling_const_iterator(base::sibling_end());
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline GI&
member<GI,FI,SI,PI,MI>::info()
{
    return my_info;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::subpedigree_pointer
member<GI,FI,SI,PI,MI>::subpedigree()
{
    return reinterpret_cast<subpedigree_pointer>(base::subpedigree());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::pedigree_pointer
member<GI,FI,SI,PI,MI>::pedigree()
{
    return reinterpret_cast<pedigree_pointer>(base::pedigree());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::multipedigree_pointer
member<GI,FI,SI,PI,MI>::multipedigree()
{
    return reinterpret_cast<multipedigree_pointer>(base::multipedigree());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::family_pointer
member<GI,FI,SI,PI,MI>::family()
{
    return reinterpret_cast<family_pointer>(base::family());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::member_pointer
member<GI,FI,SI,PI,MI>::parent1()
{
    return reinterpret_cast<member_pointer>(base::parent1());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::member_pointer
member<GI,FI,SI,PI,MI>::parent2()
{
    return reinterpret_cast<member_pointer>(base::parent2());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::member_pointer
member<GI,FI,SI,PI,MI>::get_mother()
{
    return reinterpret_cast<member_pointer>(base::get_mother());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::member_pointer
member<GI,FI,SI,PI,MI>::get_father()
{
    return reinterpret_cast<member_pointer>(base::get_father());
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::mate_iterator
member<GI,FI,SI,PI,MI>::mate_begin()
{
    return mate_iterator(base::mate_begin());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::mate_iterator
member<GI,FI,SI,PI,MI>::mate_end()
{
    return mate_iterator(base::mate_end());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::offspring_iterator
member<GI,FI,SI,PI,MI>::offspring_begin(const mate_const_iterator& m)
{
    return offspring_iterator(base::offspring_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::offspring_iterator
member<GI,FI,SI,PI,MI>::offspring_begin(const member_const_iterator& m)
{
    return offspring_iterator(base::offspring_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::offspring_iterator
member<GI,FI,SI,PI,MI>::offspring_begin(const member& m)
{
    return offspring_iterator(base::offspring_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::offspring_iterator
member<GI,FI,SI,PI,MI>::offspring_begin(const string& m)
{
    return offspring_iterator(base::offspring_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::offspring_iterator
member<GI,FI,SI,PI,MI>::offspring_end()
{
    return offspring_iterator(base::offspring_end());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::parent_iterator
member<GI,FI,SI,PI,MI>::parent_begin()
{
    return parent_iterator(base::parent_begin());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::parent_iterator
member<GI,FI,SI,PI,MI>::parent_end()
{
    return parent_iterator(base::parent_end());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::progeny_iterator
member<GI,FI,SI,PI,MI>::progeny_begin()
{
    return progeny_iterator(base::progeny_begin());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::progeny_iterator
member<GI,FI,SI,PI,MI>::progeny_end()
{
    return progeny_iterator(base::progeny_end());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member<GI,FI,SI,PI,MI>::sibling_iterator
member<GI,FI,SI,PI,MI>::sibling_begin()
{
    return sibling_iterator(base::sibling_begin());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename member<GI,FI,SI,PI,MI>::sibling_iterator
member<GI,FI,SI,PI,MI>::sibling_end()
{
    return sibling_iterator(base::sibling_end());
}


//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline void
member<GI,FI,SI,PI,MI>::set_sex(SexCode s)
{
    base::set_sex(s);
}

}
}

