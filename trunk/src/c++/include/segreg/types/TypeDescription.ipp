#ifndef TYPES_TYPE_DESCRIPTION_H
#include "segreg/types/TypeDescription.h"
#endif

namespace SAGE   {
namespace SEGREG {

// =====================================
//
// TypeDescription::State
//
// =====================================
  
inline
TypeDescription::State::State()
  : my_index((size_t) -1),
    my_name ()
{ }

inline
TypeDescription::State::State(const State& s)
  : my_index(s.my_index),
    my_name (s.my_name )
{ }

inline
TypeDescription::State& 
  TypeDescription::State::operator=(const State& s)
{
  if(this != &s)
  {
    my_index = s.my_index;
    my_name  = s.my_name;
  }
  return *this;
}

inline
size_t
  TypeDescription::State::get_index () const
{
  return my_index;
}

inline
std::string
  TypeDescription::State::get_name  () const
{
  return my_name;
}

inline
TypeDescription::State::State(size_t index, std::string name)
  : my_index(index),
    my_name (name )
{ }

// =====================================
//
// TypeDescription
//
// =====================================
  
inline
TypeDescription::TypeDescription()
  : my_name        (),
    my_description (),
    my_states      ()
{ }

inline
TypeDescription::TypeDescription(const TypeDescription& t)
  : my_name        (t.my_name),
    my_description (t.my_description),
    my_states      (t.my_states)
{ }

inline
TypeDescription& TypeDescription::operator=(const TypeDescription& t)
{
  if(this != &t)
  {
    my_name        = t.my_name;
    my_description = t.my_description;
    my_states      = t.my_states;
  }
  
  return *this;
}

inline
void TypeDescription::set_name(std::string s)
{
  my_name = s;
}

inline
std::string TypeDescription::get_name() const
{
  return my_name;
}

inline
void TypeDescription::set_description(std::string s)
{
  my_description = s;
}

inline
std::string TypeDescription::get_description() const
{
  return my_description;
}

inline
void TypeDescription::reserve_state_count(size_t s)
{
  my_states.reserve(s);
}

inline
size_t TypeDescription::add_state(std::string s)
{
  my_states.push_back(State(my_states.size(), s));
  
  return my_states.size() - 1;
}

inline
size_t TypeDescription::get_state_count() const
{
  return my_states.size();
}

inline
const TypeDescription::State& TypeDescription::get_state(size_t s) const
{
  return my_states[s];
}

inline
TypeDescription::StateIterator TypeDescription::begin() const
{
  return my_states.begin();
}

inline
TypeDescription::StateIterator TypeDescription::end()   const
{
  return my_states.end();
}

inline std::ostream& operator<<(std::ostream& o, const TypeDescription& t)
{
  o << convert_to_output(t);
  
  return o;
}

} // End Namespace SEGREG
} // End Namespace SAGE

