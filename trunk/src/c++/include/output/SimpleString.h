#ifndef OUTPUT_SIMPLE_STRING_H
#define OUTPUT_SIMPLE_STRING_H

#include "output/ViewPrettyPrint.h"

namespace SAGE   {
namespace OUTPUT {

class SimpleString
{
  private:
  
    struct convertDouble
    {
      static void go(double d, std::string & t)
      {
        std::ostringstream s;
        s << OUTPUT::Double(d);
        t = s.str();
      }
    };
    
    struct convertInt
    {
      static void go(int i, std::string & t)
      {
        std::ostringstream s;
        s << OUTPUT::Int(i);
        t = s.str();
      }
    };
    
    struct convertString
    {
      static void go(const std::string & s, std::string & t)
      {
        t = s;
      }
    };
    
    struct cantConvert
    {
      template<typename T>
      static void go(const T & x, std::string & t)
      {
        BOOST_STATIC_ASSERT((sizeof(T) == 0));
      }
    };
    
    

  public:

    SimpleString() : my_text("") { }
    
    template<typename T>
    SimpleString(const T & t)
    {
      using namespace boost::mpl;
      using namespace boost;

      typedef is_same        <T, double>       is_double;
      typedef is_integral    <T>               treat_as_int;
      typedef is_convertible <T, std::string>  treat_as_string;

      if_<is_double, convertDouble,             typename
          if_<treat_as_int, convertInt,         typename
              if_<treat_as_string, convertString, cantConvert>::type 
          >::type 
      >::type::go(t, my_text);
    }
    
    SimpleString(const SimpleString & other) : my_text(other.my_text) { }
    
    SimpleString& operator=(const SimpleString & other)
    {
      if(this != &other)
      {
        my_text = other.my_text;
      }
      
      return *this;
    }
    
    const std::string & toString() const { return my_text; }

  private:
  
    std::string my_text;
};

} // End namespace OUTPUT
} // End namespace SAGE


#endif
