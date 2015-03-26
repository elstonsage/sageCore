#ifndef AO_POOLING_H
#define AO_POOLING_H

#include <string>
#include <list>
#include <vector>
#include <cmath>

namespace SAGE {
namespace AO   {

class PoolingCfg
{
  public:

  /// @name Constructors / operators
  //@{
  
    PoolingCfg() { for(size_t i = 0; i < 6; ++i) my_pools[i] = i; }  
   
    PoolingCfg(const PoolingCfg & other) { for(size_t i = 0; i < 6; ++i) my_pools[i] = other.my_pools[i]; }
   
    PoolingCfg& operator=(const PoolingCfg & other)
    {
      if(this != &other)
        for(size_t i = 0; i < 6; ++i)
          my_pools[i] = other.my_pools[i];
     
      return *this;
    }
   
  //@}   
  
  /// @name Helper functions
  //@{
  
    ///
    /// Returns the classification for the given text code.
    /// ?? = 
    /// ?A =
    /// ?U =
    /// AA =
    /// AU =
    /// UU =
    /// [invalid] = -1
    ///
    /// Please note that it's case INsensitive, and order insensitive (AU = UA).
    static size_t getClass(const std::string & code)
    {
      // Check the length:
      if(code.length() != 2)
        return (size_t)-1;
        
      std::string adj_code = code;
      
      // Fix the case:
      if(adj_code[0] == 'a') adj_code[0] = 'A';
      if(adj_code[0] == 'u') adj_code[0] = 'U';
      if(adj_code[1] == 'a') adj_code[1] = 'A';
      if(adj_code[1] == 'u') adj_code[1] = 'U';
      
      // Make sure the "?" comes first:
      if(adj_code == "A?") adj_code = "?A";
      if(adj_code == "U?") adj_code = "?U";

      // Make sure "A" comes first:
      if(adj_code == "UA") adj_code = "AU";

      // Now read and return the correct code:
           if(code == "??") { return 0;          }
      else if(code == "?A") { return 1;          }
      else if(code == "?U") { return 2;          }
      else if(code == "AA") { return 3;          }
      else if(code == "AU") { return 4;          }
      else if(code == "UU") { return 5;          }
      else                  { return (size_t)-1; }
    }
  
    ///
    /// Returns the string code for the given classification.
    /// 0  = ??
    /// 1  = ?A
    /// 2  = ?U
    /// 3  = AA
    /// 4  = AU
    /// 5  = UU
    /// >5 = "" (Invalid classification)
    ///
    static std::string getCode(size_t classification)
    {
      switch(classification)
      {
        case 0: return "??";
        case 1: return "?A";
        case 2: return "?U";
        case 3: return "AA";
        case 4: return "AU";
        case 5: return "UU";
        default: return "";
      }
    }
  
  //@}

  /// @name Setting / getting pooling
  //@{
  
    ///
    /// Indicates that individuals from class \c from_class should be reassigned to \to_class .
    /// Returns \c true if successful; \c false otherewise. A \c false return indicates that
    /// the reassignment is illogical (for instance, reassigning AA to UU doesn't make sense!).
    bool setPool(size_t from_class, size_t to_class)
    {
      // Convert them, so it's easier to write:
      std::string from_code = getCode(from_class),
                  to_code   = getCode(to_class);
                 
      // Check that the conversion makes sense:
      if((from_code == "??") ||                                          // ?? can become anything
         (from_code == "?A" && (to_code == "AA" || to_code == "AU")) ||  // ?A can become AA or AU
         (from_code == "?U" && (to_code == "UU" || to_code == "AU")))    // ?U can become UU or AU
      {
        // Set the reassignment:
        my_pools[from_class] = to_class;

        // Return success:        
        return true;
      }

      // Ack! Invalid reassignment:
      else
      {
        return false;
      }
    }
  
    ///
    /// Given the class \c from_class to which an individual belongs, indicates to which
    /// class that individual should be re-assigned.  
    size_t getPool(size_t from_class) const
    {
      return my_pools[from_class];
    }
  
  //@}

  private:
  
    // At index i, where i is the classification of an individual, my_pools[i]
    // indicates the classification to which the individual SHOULD be (re)assigned.
    size_t my_pools[6];
};

} // End namespace AO
} // End namespace SAGE

#endif
