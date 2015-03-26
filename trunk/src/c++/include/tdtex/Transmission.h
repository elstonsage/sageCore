#ifndef TDTEX_TRANSMISSION_H
#define TDTEX_TRANSMISSION_H

#include <vector>
#include <iostream>
#include <algorithm>
#include "mlocus/genotype.h"
#include "rped/rped.h"

namespace SAGE  {
namespace TDTEX {

/// \struct Transmission
struct Transmission
{
  /// @name Constructor & operators
  //@{
  
    ///
    /// Constructor.
    /// \param transmitted The transmitted id
    /// \param not_transmitted The non-transmitted id
    Transmission(
            RPED::MemberConstPointer   _source,
      const MLOCUS::allele           & _transmitted, 
      const MLOCUS::allele           & _not_transmitted) : 
      
      source          (_source),
      transmitted     (_transmitted), 
      not_transmitted (_not_transmitted) { }

    ///
    /// Copy constructor.
    Transmission(const Transmission & other) : source(other.source), transmitted(other.transmitted), not_transmitted(other.not_transmitted) { }

    ///
    /// Operator=
    Transmission& operator=(const Transmission & other) 
    {
      if(this != &other)
      {
        source          = other.source;
        transmitted     = other.transmitted;
        not_transmitted = other.not_transmitted;
      }
      
      return *this;
    }

    ///
    /// Operator==
    bool operator==(const Transmission & other) const
    {
      if(!source || !other.source)
        return (transmitted == other.transmitted) && (not_transmitted == other.not_transmitted);
      else
        return (source == other.source) && ((transmitted == other.transmitted) && (not_transmitted == other.not_transmitted));
    }

    ///
    /// Operator!=
    bool operator!=(const Transmission & other) const {   return !(*this == other); }
    
  //@}

  RPED::member_const_pointer source;           // The 'source' of the transmission (contributing parent)
  MLOCUS::allele             transmitted;      // The identifier of the thing that WAS transmitted.
  MLOCUS::allele             not_transmitted;  // The identified of the thing that was NOT transmitted.
};

typedef std::vector<Transmission> TransmissionVector;

class TransmissionList
{
public:

  /// @name Constructors / operators
  //@{
  
    ///
    /// Constructor.
    TransmissionList();
  
    ///
    /// Copy constructor.
    TransmissionList(const TransmissionList& other);

    ///
    /// Operator=.
    TransmissionList& operator =(const TransmissionList& other);
    
  //@}
  
  /// @name Accessing the underlying list
  //@{

    ///
    /// Const access...
    const TransmissionVector & get_list() const { return my_list; }
  
    ///
    /// Non-const access...
    TransmissionVector & get_list() { return my_list; }
    
  //@}
  
  /// @name Error status
  //@{

    ///
    /// Indicates whether or not there's an error present.
    bool get_error() const { return my_error; }

    ///
    /// Sets the error indicator to e.
    void set_error(bool e) { my_error = e; }
    
  //@}

  /// @name Debugging
  //@{
  
    void dump() const;
  
  //@}

protected:
  bool my_error;

private:

  TransmissionVector my_list;
};


TransmissionList transmission_intersection(const TransmissionList& t1, const TransmissionList& t2);

std::ostream& operator<<(std::ostream& o, const TransmissionList& tl);

inline
TransmissionList::TransmissionList()
{
  set_error(false);
}

inline TransmissionList::TransmissionList(const TransmissionList& other) : my_error(other.my_error), my_list(other.my_list) { }

inline
TransmissionList& TransmissionList::operator=(const TransmissionList& other)
{
  if(this != &other)
  {
    my_error = other.my_error;
    my_list  = other.my_list;
  }
  
  return *this;
}


} // End namespace TDTEX
} // End namespace SAGE

#endif
