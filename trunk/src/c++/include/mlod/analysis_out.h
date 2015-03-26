#ifndef ANALYSIS_OUT_H
#define ANALYSIS_OUT_H

#include <iostream>
#include <list>

namespace SAGE {
namespace MLOD {


class MultiStream
{
  public:
    MultiStream()                                                { }
    MultiStream(const MultiStream& m) : my_streams(m.my_streams) { }
    
    ~MultiStream() { }
    
    MultiStream& operator=(const MultiStream& m) { if(this != &m) { my_streams = m.my_streams; } return *this; }
    
    void add_ostream(std::ostream& o) { my_streams.push_back(&o); }
    
    template<class T>
    MultiStream& operator<<(const T& t)
    {
      for(std::list<std::ostream*>::iterator i = my_streams.begin();
          i != my_streams.end(); ++i)
      {
        (**i) << t;
      }
      
      return *this;
    }

    MultiStream& operator+= (std::ostream& o)
    {
      add_ostream(o);
      
      return *this;
    }
    MultiStream operator+ (std::ostream& o)
    {
      MultiStream tmp(*this);
      
      return tmp+= o;
    }
    
  private:
  
    std::list<std::ostream*> my_streams;
    
};

inline 
MultiStream operator+ (std::ostream& o1, std::ostream& o2)
{
  MultiStream s;
  
  return s + o1 + o2;
}

}
}

#endif
