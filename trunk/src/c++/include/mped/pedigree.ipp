namespace SAGE {
namespace MPED {

//============================================================================
//  IMPLEMENTATION: pedigree<GI,FI,SI,PI,MI>
//============================================================================
//
template <class GI, class FI, class SI, class PI, class MI>
pedigree<GI,FI,SI,PI,MI>::pedigree(const string& name)
  : pedigree_base(name), my_build_callback(0), my_freeze_callback(0)
{}

template <class GI, class FI, class SI, class PI, class MI>
pedigree<GI,FI,SI,PI,MI>::pedigree
(const string& name, multipedigree_type& mp)
  : pedigree_base(name, reinterpret_cast<multipedigree_base*>(&mp)), 
    my_build_callback(0), my_freeze_callback(0)
{}

template <class GI, class FI, class SI, class PI, class MI>
pedigree<GI,FI,SI,PI,MI>::~pedigree()
{
    clear_all();
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline uint
pedigree<GI,FI,SI,PI,MI>::index() const
{
    return base::index();
}

template <class GI, class FI, class SI, class PI, class MI> inline uint
pedigree<GI,FI,SI,PI,MI>::error_count() const
{
    return base::error_count();
}

template <class GI, class FI, class SI, class PI, class MI> inline uint
pedigree<GI,FI,SI,PI,MI>::family_count() const
{
    return base::family_count();
}

template <class GI, class FI, class SI, class PI, class MI> inline uint
pedigree<GI,FI,SI,PI,MI>::member_count() const
{
    return base::member_count();
}

template <class GI, class FI, class SI, class PI, class MI> inline uint
pedigree<GI,FI,SI,PI,MI>::subpedigree_count() const
{
    return base::subpedigree_count();
}

template <class GI, class FI, class SI, class PI, class MI> inline uint
pedigree<GI,FI,SI,PI,MI>::unconnected_count() const
{
    return base::unconnected_count();
}

template <class GI, class FI, class SI, class PI, class MI> inline const PI&
pedigree<GI,FI,SI,PI,MI>::info() const
{
    return my_info;
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::multipedigree_const_pointer
pedigree<GI,FI,SI,PI,MI>::multipedigree() const
{
    return reinterpret_cast<multipedigree_const_pointer>(base::multipedigree());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
const string&
pedigree<GI,FI,SI,PI,MI>::name() const
{
    return base::name();
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline
const typename pedigree<GI,FI,SI,PI,MI>::subpedigree_type&
pedigree<GI,FI,SI,PI,MI>::subpedigree_index(uint i) const
{
    return reinterpret_cast<const subpedigree_type&>(base::subpedigree_index(i));
}

template <class GI, class FI, class SI, class PI, class MI> inline
const typename pedigree<GI,FI,SI,PI,MI>::member_type&
pedigree<GI,FI,SI,PI,MI>::member_index(uint i) const
{
    return reinterpret_cast<const member_type&>(base::member_index(i));
}

template <class GI, class FI, class SI, class PI, class MI> inline
const typename pedigree<GI,FI,SI,PI,MI>::family_type&
pedigree<GI,FI,SI,PI,MI>::family_index(uint i) const
{
    return reinterpret_cast<const family_type&>(base::family_index(i));
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::family_const_pointer
pedigree<GI,FI,SI,PI,MI>::family_find(const string& p1, const string& p2) const
{
    return reinterpret_cast<family_const_pointer>(base::family_find(p1, p2));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::family_const_pointer
pedigree<GI,FI,SI,PI,MI>::family_find
(const member_type& p1, const member_type& p2) const
{
    return reinterpret_cast<family_const_pointer>(base::family_find(p1, p2));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::family_const_pointer
pedigree<GI,FI,SI,PI,MI>::family_find
(const member_const_iterator& p1, const member_const_iterator& p2) const
{
    return reinterpret_cast<family_const_pointer>(base::family_find(p1, p2));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::member_const_pointer
pedigree<GI,FI,SI,PI,MI>::member_find(const string& m) const
{
    return reinterpret_cast<member_const_pointer>(base::member_find(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::subpedigree_const_pointer
pedigree<GI,FI,SI,PI,MI>::subpedigree_find(const string& m) const
{
    return reinterpret_cast<subpedigree_const_pointer>(base::subpedigree_find(m));
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::error_iterator
pedigree<GI,FI,SI,PI,MI>::error_begin() const
{
    return base::error_begin();
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::error_iterator
pedigree<GI,FI,SI,PI,MI>::error_end() const
{
    return base::error_end();
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::family_const_iterator
pedigree<GI,FI,SI,PI,MI>::family_begin() const
{
    return family_const_iterator(base::family_begin());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::family_const_iterator
pedigree<GI,FI,SI,PI,MI>::family_begin(const subpedigree_type& S) const
{
    return family_const_iterator(base::family_begin(S));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::family_const_iterator
pedigree<GI,FI,SI,PI,MI>::family_begin(const subpedigree_const_iterator& S) const
{
    return family_const_iterator(base::family_begin(S));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::family_const_iterator
pedigree<GI,FI,SI,PI,MI>::family_end() const
{
    return family_const_iterator(base::family_end());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::family_const_iterator
pedigree<GI,FI,SI,PI,MI>::family_end(const subpedigree_type& S) const
{
    return family_const_iterator(base::family_end(S));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::family_const_iterator
pedigree<GI,FI,SI,PI,MI>::family_end(const subpedigree_const_iterator& S) const
{
    return family_const_iterator(base::family_end(S));
}

//-
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::mate_const_iterator
pedigree<GI,FI,SI,PI,MI>::mate_begin(const string& m) const
{
    return mate_const_iterator(base::mate_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::mate_const_iterator
pedigree<GI,FI,SI,PI,MI>::mate_begin(const member_type& m) const
{
    return mate_const_iterator(base::mate_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::mate_const_iterator
pedigree<GI,FI,SI,PI,MI>::mate_begin(const member_const_iterator& m) const
{
    return mate_const_iterator(base::mate_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::mate_const_iterator
pedigree<GI,FI,SI,PI,MI>::mate_end(const string& m) const
{
    return mate_const_iterator(base::mate_end(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::mate_const_iterator
pedigree<GI,FI,SI,PI,MI>::mate_end(const member_type& m) const
{
    return mate_const_iterator(base::mate_end(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::mate_const_iterator
pedigree<GI,FI,SI,PI,MI>::mate_end(const member_const_iterator& m) const
{
    return mate_const_iterator(base::mate_end(m));
}

//-
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::member_const_iterator
pedigree<GI,FI,SI,PI,MI>::member_begin() const
{
    return member_const_iterator(base::member_begin());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::member_const_iterator
pedigree<GI,FI,SI,PI,MI>::member_begin(const subpedigree_type& S) const
{
    return member_const_iterator(base::member_begin(S));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::member_const_iterator
pedigree<GI,FI,SI,PI,MI>::member_begin(const subpedigree_const_iterator& S) const
{
    return member_const_iterator(base::member_begin(S));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::member_const_iterator
pedigree<GI,FI,SI,PI,MI>::member_end() const
{
    return member_const_iterator(base::member_end());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::member_const_iterator
pedigree<GI,FI,SI,PI,MI>::member_end(const subpedigree_type& S) const
{
    return member_const_iterator(base::member_end(S));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::member_const_iterator
pedigree<GI,FI,SI,PI,MI>::member_end(const subpedigree_const_iterator& S) const
{
    return member_const_iterator(base::member_end(S));
}

//-
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::offspring_const_iterator
pedigree<GI,FI,SI,PI,MI>::offspring_begin(const family_type& f) const
{
    return offspring_const_iterator(base::offspring_begin(f));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::offspring_const_iterator
pedigree<GI,FI,SI,PI,MI>::offspring_begin(const family_const_iterator& f) const
{
    return offspring_const_iterator(base::offspring_begin(f));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::offspring_const_iterator
pedigree<GI,FI,SI,PI,MI>::offspring_begin
(const member_type& p1, const member_type& p2) const
{
    return offspring_const_iterator(base::offspring_begin(p1, p2));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::offspring_const_iterator
pedigree<GI,FI,SI,PI,MI>::offspring_begin
(const member_const_iterator& p1, const member_const_iterator& p2) const
{
    return offspring_const_iterator(base::offspring_begin(p1, p2));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::offspring_const_iterator
pedigree<GI,FI,SI,PI,MI>::offspring_begin(const string& p1, const string& p2) const
{
    return offspring_const_iterator(base::offspring_begin(p1, p2));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::offspring_const_iterator
pedigree<GI,FI,SI,PI,MI>::offspring_end() const
{
    return offspring_const_iterator(base::offspring_end());
}

//-
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::parent_const_iterator
pedigree<GI,FI,SI,PI,MI>::parent_begin(const string& m) const
{
    return parent_const_iterator(base::parent_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::parent_const_iterator
pedigree<GI,FI,SI,PI,MI>::parent_begin(const member_type& m) const
{
    return parent_const_iterator(base::parent_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::parent_const_iterator
pedigree<GI,FI,SI,PI,MI>::parent_begin(const member_const_iterator& m) const
{
    return parent_const_iterator(base::parent_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::parent_const_iterator
pedigree<GI,FI,SI,PI,MI>::parent_begin(const family_type& f) const
{
    return parent_const_iterator(base::parent_begin(f));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::parent_const_iterator
pedigree<GI,FI,SI,PI,MI>::parent_begin(const family_const_iterator& f) const
{
    return parent_const_iterator(base::parent_begin(f));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::parent_const_iterator
pedigree<GI,FI,SI,PI,MI>::parent_end() const
{
    return parent_const_iterator(base::parent_end());
}

//-
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::progeny_const_iterator
pedigree<GI,FI,SI,PI,MI>::progeny_begin(const string& m) const
{
    return progeny_const_iterator(base::progeny_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::progeny_const_iterator
pedigree<GI,FI,SI,PI,MI>::progeny_begin(const member_type& m) const
{
    return progeny_const_iterator(base::progeny_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::progeny_const_iterator
pedigree<GI,FI,SI,PI,MI>::progeny_begin(const member_const_iterator& m) const
{
    return progeny_const_iterator(base::progeny_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::progeny_const_iterator
pedigree<GI,FI,SI,PI,MI>::progeny_end(const string& m) const
{
    return progeny_const_iterator(base::progeny_end(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::progeny_const_iterator
pedigree<GI,FI,SI,PI,MI>::progeny_end(const member_type& m) const
{
    return progeny_const_iterator(base::progeny_end(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::progeny_const_iterator
pedigree<GI,FI,SI,PI,MI>::progeny_end(const member_const_iterator& m) const
{
    return progeny_const_iterator(base::progeny_end(m));
}

//-
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::sibling_const_iterator
pedigree<GI,FI,SI,PI,MI>::sibling_begin(const string& m) const
{
    return sibling_const_iterator(base::sibling_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::sibling_const_iterator
pedigree<GI,FI,SI,PI,MI>::sibling_begin(const member_type& m) const
{
    return sibling_const_iterator(base::sibling_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::sibling_const_iterator
pedigree<GI,FI,SI,PI,MI>::sibling_begin(const member_const_iterator& m) const
{
    return sibling_const_iterator(base::sibling_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::sibling_const_iterator
pedigree<GI,FI,SI,PI,MI>::sibling_end() const
{
    return sibling_const_iterator(base::sibling_end());
}

//-
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::subpedigree_const_iterator
pedigree<GI,FI,SI,PI,MI>::subpedigree_begin() const
{
    return subpedigree_const_iterator(base::subpedigree_begin());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::subpedigree_const_iterator
pedigree<GI,FI,SI,PI,MI>::subpedigree_end() const
{
    return subpedigree_const_iterator(base::subpedigree_end());
}

//-
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::member_const_iterator
pedigree<GI,FI,SI,PI,MI>::unconnected_begin() const
{
    return member_const_iterator(base::unconnected_begin());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::member_const_iterator
pedigree<GI,FI,SI,PI,MI>::unconnected_end() const
{
    return member_const_iterator(base::unconnected_end());
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline PI&
pedigree<GI,FI,SI,PI,MI>::info()
{
    return my_info;
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::multipedigree_pointer
pedigree<GI,FI,SI,PI,MI>::multipedigree()
{
    return reinterpret_cast<multipedigree_pointer>(base::multipedigree());
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
subpedigree<GI,FI,SI,PI,MI>&
pedigree<GI,FI,SI,PI,MI>::subpedigree_index(uint i)
{
    return reinterpret_cast<subpedigree_type&>(base::subpedigree_index(i));
}

template <class GI, class FI, class SI, class PI, class MI> inline 
member<GI,FI,SI,PI,MI>&
pedigree<GI,FI,SI,PI,MI>::member_index(uint i)
{
    return reinterpret_cast<member_type&>(base::member_index(i));
}

template <class GI, class FI, class SI, class PI, class MI> inline 
family<GI,FI,SI,PI,MI>&
pedigree<GI,FI,SI,PI,MI>::family_index(uint i)
{
    return reinterpret_cast<family_type&>(base::family_index(i));
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline void
pedigree<GI,FI,SI,PI,MI>::subpedigree_index_swap(uint i, uint j)
{
    base::subpedigree_index_swap(i, j);
}

template <class GI, class FI, class SI, class PI, class MI> inline void
pedigree<GI,FI,SI,PI,MI>::member_index_swap(uint i, uint j)
{
    base::member_index_swap(i, j);
}

template <class GI, class FI, class SI, class PI, class MI> inline void
pedigree<GI,FI,SI,PI,MI>::family_index_swap(uint i, uint j)
{
    base::family_index_swap(i, j);
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::family_pointer
pedigree<GI,FI,SI,PI,MI>::family_find(const string& p1, const string& p2)
{
    return reinterpret_cast<family_pointer>(base::family_find(p1, p2));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::family_pointer
pedigree<GI,FI,SI,PI,MI>::family_find
(const member_type& p1, const member_type& p2)
{
    return reinterpret_cast<family_pointer>(base::family_find(p1, p2));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::family_pointer
pedigree<GI,FI,SI,PI,MI>::family_find
(const member_const_iterator& p1, const member_const_iterator& p2)
{
    return reinterpret_cast<family_pointer>(base::family_find(p1, p2));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::member_pointer
pedigree<GI,FI,SI,PI,MI>::member_find(const string& m)
{
    return reinterpret_cast<member_pointer>(base::member_find(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::subpedigree_pointer
pedigree<GI,FI,SI,PI,MI>::subpedigree_find(const string& m)
{
    return reinterpret_cast<subpedigree_pointer>(base::subpedigree_find(m));
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::family_iterator
pedigree<GI,FI,SI,PI,MI>::family_begin()
{
    return family_iterator(base::family_begin());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::family_iterator
pedigree<GI,FI,SI,PI,MI>::family_begin(const subpedigree_type& S)
{
    return family_iterator(base::family_begin(S));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::family_iterator
pedigree<GI,FI,SI,PI,MI>::family_begin(const subpedigree_const_iterator& S)
{
    return family_iterator(base::family_begin(S));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::family_iterator
pedigree<GI,FI,SI,PI,MI>::family_end()
{
    return family_iterator(base::family_end());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::family_iterator
pedigree<GI,FI,SI,PI,MI>::family_end(const subpedigree_type& S)
{
    return family_iterator(base::family_end(S));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::family_iterator
pedigree<GI,FI,SI,PI,MI>::family_end(const subpedigree_const_iterator& S)
{
    return family_iterator(base::family_end(S));
}

//-
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::mate_iterator
pedigree<GI,FI,SI,PI,MI>::mate_begin(const string& m)
{
    return mate_iterator(base::mate_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::mate_iterator
pedigree<GI,FI,SI,PI,MI>::mate_begin(const member_type& m)
{
    return mate_iterator(base::mate_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::mate_iterator
pedigree<GI,FI,SI,PI,MI>::mate_begin(const member_const_iterator& m)
{
    return mate_iterator(base::mate_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::mate_iterator
pedigree<GI,FI,SI,PI,MI>::mate_end(const string& m)
{
    return mate_iterator(base::mate_end(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::mate_iterator
pedigree<GI,FI,SI,PI,MI>::mate_end(const member_type& m)
{
    return mate_iterator(base::mate_end(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::mate_iterator
pedigree<GI,FI,SI,PI,MI>::mate_end(const member_const_iterator& m)
{
    return mate_iterator(base::mate_end(m));
}

//-
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::member_iterator
pedigree<GI,FI,SI,PI,MI>::member_begin()
{
    return member_iterator(base::member_begin());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::member_iterator
pedigree<GI,FI,SI,PI,MI>::member_begin(const subpedigree_type& S)
{
    return member_iterator(base::member_begin(S));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::member_iterator
pedigree<GI,FI,SI,PI,MI>::member_begin(const subpedigree_const_iterator& S)
{
    return member_iterator(base::member_begin(S));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::member_iterator
pedigree<GI,FI,SI,PI,MI>::member_end()
{
    return member_iterator(base::member_end());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::member_iterator
pedigree<GI,FI,SI,PI,MI>::member_end(const subpedigree_type& S)
{
    return member_iterator(base::member_end(S));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::member_iterator
pedigree<GI,FI,SI,PI,MI>::member_end(const subpedigree_const_iterator& S)
{
    return member_iterator(base::member_end(S));
}

//-
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::offspring_iterator
pedigree<GI,FI,SI,PI,MI>::offspring_begin(const family_type& f)
{
    return offspring_iterator(base::offspring_begin(f));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::offspring_iterator
pedigree<GI,FI,SI,PI,MI>::offspring_begin(const family_const_iterator& f)
{
    return offspring_iterator(base::offspring_begin(f));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::offspring_iterator
pedigree<GI,FI,SI,PI,MI>::offspring_begin
(const member_type& p1, const member_type& p2)
{
    return offspring_iterator(base::offspring_begin(p1, p2));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::offspring_iterator
pedigree<GI,FI,SI,PI,MI>::offspring_begin
(const member_const_iterator& p1, const member_const_iterator& p2)
{
    return offspring_iterator(base::offspring_begin(p1, p2));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::offspring_iterator
pedigree<GI,FI,SI,PI,MI>::offspring_begin(const string& p1, const string& p2)
{
    return offspring_iterator(base::offspring_begin(p1, p2));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::offspring_iterator
pedigree<GI,FI,SI,PI,MI>::offspring_end()
{
    return offspring_iterator(base::offspring_end());
}

//-
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::parent_iterator
pedigree<GI,FI,SI,PI,MI>::parent_begin(const string& m)
{
    return parent_iterator(base::parent_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::parent_iterator
pedigree<GI,FI,SI,PI,MI>::parent_begin(const member_type& m)
{
    return parent_iterator(base::parent_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::parent_iterator
pedigree<GI,FI,SI,PI,MI>::parent_begin(const member_const_iterator& m)
{
    return parent_iterator(base::parent_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::parent_iterator
pedigree<GI,FI,SI,PI,MI>::parent_begin(const family_type& f)
{
    return parent_iterator(base::parent_begin(f));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::parent_iterator
pedigree<GI,FI,SI,PI,MI>::parent_begin(const family_const_iterator& f)
{
    return parent_iterator(base::parent_begin(f));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::parent_iterator
pedigree<GI,FI,SI,PI,MI>::parent_end()
{
    return parent_iterator(base::parent_end());
}

//-
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::progeny_iterator
pedigree<GI,FI,SI,PI,MI>::progeny_begin(const string& m)
{
    return progeny_iterator(base::progeny_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::progeny_iterator
pedigree<GI,FI,SI,PI,MI>::progeny_begin(const member_type& m)
{
    return progeny_iterator(base::progeny_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::progeny_iterator
pedigree<GI,FI,SI,PI,MI>::progeny_begin(const member_const_iterator& m)
{
    return progeny_iterator(base::progeny_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::progeny_iterator
pedigree<GI,FI,SI,PI,MI>::progeny_end(const string& m)
{
    return progeny_iterator(base::progeny_end(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::progeny_iterator
pedigree<GI,FI,SI,PI,MI>::progeny_end(const member_type& m)
{
    return progeny_iterator(base::progeny_end(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::progeny_iterator
pedigree<GI,FI,SI,PI,MI>::progeny_end(const member_const_iterator& m)
{
    return progeny_iterator(base::progeny_end(m));
}

//-
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::sibling_iterator
pedigree<GI,FI,SI,PI,MI>::sibling_begin(const string& m)
{
    return sibling_iterator(base::sibling_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::sibling_iterator
pedigree<GI,FI,SI,PI,MI>::sibling_begin(const member_type& m)
{
    return sibling_iterator(base::sibling_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::sibling_iterator
pedigree<GI,FI,SI,PI,MI>::sibling_begin(const member_const_iterator& m)
{
    return sibling_iterator(base::sibling_begin(m));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::sibling_iterator
pedigree<GI,FI,SI,PI,MI>::sibling_end()
{
    return sibling_iterator(base::sibling_end());
}

//-
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::subpedigree_iterator
pedigree<GI,FI,SI,PI,MI>::subpedigree_begin()
{
    return subpedigree_iterator(base::subpedigree_begin());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::subpedigree_iterator
pedigree<GI,FI,SI,PI,MI>::subpedigree_end()
{
    return subpedigree_iterator(base::subpedigree_end());
}

//-
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::member_iterator
pedigree<GI,FI,SI,PI,MI>::unconnected_begin()
{
    return member_iterator(base::unconnected_begin());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::member_iterator
pedigree<GI,FI,SI,PI,MI>::unconnected_end()
{
    return member_iterator(base::unconnected_end());
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline void
pedigree<GI,FI,SI,PI,MI>::add_lineage
(const string& child, const string& parent)
{
    base::add_lineage(child, parent);
}

template <class GI, class FI, class SI, class PI, class MI> inline void
pedigree<GI,FI,SI,PI,MI>::add_lineage
(const string& child, const string& parent1, const string& parent2)
{
    base::add_lineage(child, parent1, parent2);
}

template <class GI, class FI, class SI, class PI, class MI> inline void
pedigree<GI,FI,SI,PI,MI>::add_marriage
(const string& spouse1, const string& spouse2)
{
    base::add_marriage(spouse1, spouse2);
}

template <class GI, class FI, class SI, class PI, class MI> inline void
pedigree<GI,FI,SI,PI,MI>::add_member(const string& n, SexCode x)
{
    base::add_member(n, x, 0);
}

template <class GI, class FI, class SI, class PI, class MI> inline void
pedigree<GI,FI,SI,PI,MI>::add_member
(const string& n, SexCode x, const geninfo_type& i)
{
    base::add_member(n, x, &i);
}

template <class GI, class FI, class SI, class PI, class MI> inline void
pedigree<GI,FI,SI,PI,MI>::add_sibship(const string& sib1, const string& sib2)
{
    base::add_sibship(sib1, sib2);
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline void
pedigree<GI,FI,SI,PI,MI>::build()
{
    base::build();
    if (my_build_callback != 0)
        (*my_build_callback)(*this);
}

template <class GI, class FI, class SI, class PI, class MI> inline void
pedigree<GI,FI,SI,PI,MI>::clear()
{
    base::clear();
}

template <class GI, class FI, class SI, class PI, class MI> inline void
pedigree<GI,FI,SI,PI,MI>::flush_build_info()
{
    base::flush_build_info();
}

template <class GI, class FI, class SI, class PI, class MI> inline void
pedigree<GI,FI,SI,PI,MI>::freeze()
{
    base::freeze();
    if (my_freeze_callback != 0)
        (*my_freeze_callback)(*this);
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::callback_type
pedigree<GI,FI,SI,PI,MI>::set_build_callback(callback_type f)
{
    callback_type   tmp = my_build_callback;
    my_build_callback = f;
    return tmp;
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename pedigree<GI,FI,SI,PI,MI>::callback_type
pedigree<GI,FI,SI,PI,MI>::set_freeze_callback(callback_type f)
{
    callback_type   tmp = my_freeze_callback;
    my_freeze_callback = f;
    return tmp;
}

//----------------------------------------------------------------------------
//
template <class GI, class FI, class SI, class PI, class MI> subped_id
pedigree<GI,FI,SI,PI,MI>::allocate_subped
(const string& name, pedigree_base* pped)
{
    return new subpedigree_type(name, pped);
}

template <class GI, class FI, class SI, class PI, class MI> member_id
pedigree<GI,FI,SI,PI,MI>::allocate_member
(const string& name, SexCode x, const void* vp)
{
    const geninfo_type*     pinfo = static_cast<const geninfo_type*>(vp);

    if (pinfo)
    {
        return new member_type(name, x, *pinfo);
    }
    else
    {
        return new member_type(name, x, geninfo_type());
    }
}

template <class GI, class FI, class SI, class PI, class MI> family_id
pedigree<GI,FI,SI,PI,MI>::allocate_family()
{
    return new family_type();
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> void
pedigree<GI,FI,SI,PI,MI>::deallocate_subped(subped_id S)
{
    delete (static_cast<subpedigree_type*>(S));
}

template <class GI, class FI, class SI, class PI, class MI> void
pedigree<GI,FI,SI,PI,MI>::deallocate_member(member_id M)
{
    delete (static_cast<member_type*>(M));
}

template <class GI, class FI, class SI, class PI, class MI> void
pedigree<GI,FI,SI,PI,MI>::deallocate_family(family_id F)
{
    delete (static_cast<family_type*>(F));
}

}
}

