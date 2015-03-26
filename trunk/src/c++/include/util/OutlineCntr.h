#ifndef UTIL_OUTLINE_CNTR_H
#define UTIL_OUTLINE_CNTR_H

#include <string>
#include <sstream>

namespace SAGE {
namespace UTIL {

/// \brief Describes an 'outline' index number (ie: 5.4.0.1)
///
/// An 'outline' index number consists of a sequence of between
/// one and twenty (20) sequentially indexed integers (ie: 5,
/// 1.1.1.1.1.1.0.1.1.1.1, 63.1009.4.0.1).
class OutlineCntr
{
public:

  /// @name Constructors
  //@{

    ///
    /// Constructor.
    OutlineCntr() { reset();  }
    
    ///
    /// Copy constructor.
    OutlineCntr(const OutlineCntr & other) { my_depth = other.my_depth; for(int i = 0; i < 20; ++i) my_counters[i] = other.my_counters[i]; }

    ///
    /// Assignment operator.
    OutlineCntr& operator= (const OutlineCntr & other) { my_depth = other.my_depth; for(int i = 0; i < 20; ++i) my_counters[i] = other.my_counters[i]; return *this; }

  //@}
  
  /// @name Status
  //@{
  
    ///
    /// Returns the current depth level, where 0 is the minimum
    /// depth (ie: 0 for 4, 1 for 4.1, 2 for 4.1.45, 3 for 4.1.45.0).
    int getDepth() const { return my_depth; }

    ///
    /// Returns the current counter (ie: 5 for 3.5,
    /// 0 for 1.2.3.0).
    int getCounter() const { return my_counters[my_depth]; }
    
    ///
    /// Casts this object as a string (Equivalent to toString() ).
    operator std::string () const { return toString(); }

    ///
    /// Renders as string (ie: "4.1.0")
    /// \param start_depth By default, the string rendering starts at depth level
    /// 0 and progresses through the counters. You can, however, override this by
    /// specifying your own minimum depth level. For instance, if the value
    /// of the OutlineCntr is 4.1.2, and you invoke toString(1), you'll get
    /// "1.2" instead of "4.1.2". Please note that if start_depth is \b greater
    /// than the current depth level, this function returns an empty string.
    std::string toString(int start_depth = 0) const
    {
      if(start_depth > my_depth)
      {
        return "";
      }
      else
      {
        std::ostringstream s;
      
        for(int i = start_depth; i <= my_depth; ++i)
          s << (i > start_depth ? "." : "") << my_counters[i];
      
        return s.str();
      }
    }
  
  //@}
  
  
  /// @name Increment / decrement operations
  //@{

    ///
    /// Increases the depth by one unit (ie: goes from
    /// 5.4 to 5.4.0, or from 1 to 1.0).
    void increaseDepth()
    {
      if(my_depth < 19)
      {
        my_depth++;
        my_counters[my_depth] = 1;
      }
    }
  
    ///
    /// Decreases the depth by one unit (ie: goes from
    /// 3.8.7.4 to 3.8.7, or from 1.1.1.0.9.9 to 1.1.1.0.9).
    void decreaseDepth()
    {
      if(my_depth > 0)
      {
        my_depth--;
      }
    }

    ///
    /// Post-increments the counter at the current
    /// depth level (ie: goes from 3.4 to 3.5, or from
    /// 65.12.13 to 65.12.14).
    OutlineCntr operator++(int)
    {
      OutlineCntr next(*this);
      this->operator++();
      return next;
    }
    
    ///
    /// Pre-increment the counter at the current
    /// depth level (ie: goes from 3.4 to 3.5, or from
    /// 65.12.13 to 65.12.14).
    OutlineCntr operator++()
    {
      ++my_counters[my_depth];
      return *this;
    }
    
  //@}

private:

  void reset()
  {
    for(int i = 0; i < 20; ++i)
      my_counters[i] = 1;
      
    my_depth = 0;
  }

  int my_counters[20];
  int my_depth;
};

} // End namespace UTIL
} // End namespace SAGE


#endif
