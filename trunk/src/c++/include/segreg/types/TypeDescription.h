#ifndef TYPES_TYPE_DESCRIPTION_H
#define TYPES_TYPE_DESCRIPTION_H

#include "output/Output.h"
#include <ostream>
#include <vector>

namespace SAGE   {
namespace SEGREG {

class TypeDescription
{
  public:
  
    class State
    {
        friend class TypeDescription;
      
      public:
        inline State();
        inline State(const State&);
        
        inline State& operator=(const State&);
        
        inline size_t      get_index () const;
        inline std::string get_name  () const;
    
      private:
      
        inline State(size_t index, std::string name);
      
        size_t      my_index;
        std::string my_name;
    };

  private:

    typedef std::vector<State>          StateVector;
    
  public:
    typedef StateVector::const_iterator StateIterator;
    
    inline TypeDescription();
    inline TypeDescription(const TypeDescription&);
    
    inline TypeDescription& operator=(const TypeDescription&);
    
    inline void        set_name       (std::string);
    inline std::string get_name       (           ) const;

    inline void        set_description(std::string);
    inline std::string get_description(           ) const;
    
  /// \name State management
  //@{
    inline void reserve_state_count(size_t);
    
    inline size_t add_state(std::string);
    
    inline size_t get_state_count() const;
    
    inline const State& get_state(size_t) const;
    
    inline StateIterator begin() const;
    inline StateIterator end()   const;
  //@}

  private:

    std::string my_name;
    std::string my_description;
    
    StateVector my_states;
};

OUTPUT::Section convert_to_output(const TypeDescription& t);

inline std::ostream& operator<<(std::ostream&, const TypeDescription&);

} // End Namespace SEGREG
} // End Namespace SAGE

#include "segreg/types/TypeDescription.ipp"

#endif

