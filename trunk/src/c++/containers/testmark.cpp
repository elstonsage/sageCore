#include <iostream>
#include <iomanip>
#include "marked_vector.h"

using namespace std;
using namespace SAGE;

template <class T>
void print(const marked_vector<T> &v)
{
  for(size_t i = 0; i < v.size(); ++i)
    cout << v[i] << " ";
  cout << endl;

  typename marked_vector<T>::mark_const_iterator m = v.mark_begin();

  for(size_t i = 0; i < v.size() && m != v.mark_end(); ++i)
  {
    if( i == *m )
    {
      cout << "^ ";
      ++m;
    }
    else
      cout << "  ";
  }

  cout << "size = " << v.mark_size() << endl;
  for(size_t i = 0; i < v.mark_size(); ++i)
    cout << *( v.mark_begin() + i ) << " ";
  cout << endl;
}

int main()
{
  marked_vector<char> v;

  cout << "Creating intial vector." << endl;
  v.push_back('a'); v.mark_back();
  v.push_back('b'); v.mark_back();
  v.push_back('c'); v.mark_back();
  v.push_back('d'); v.mark_back();
  v.push_back('e'); v.mark_back();
  v.push_back('f'); v.mark_back();
  print(v);

  cout << "Resizing it larger." << endl;
  v.resize(v.size()+2,'z');
  print(v);

  cout << "Resizing it smaller." << endl;
  v.resize(v.size()/2);
  print(v);

  cout << "Clearing remaining elements." << endl;
  v.clear();
  print(v);

  v.push_back('a'); v.mark_back();
  v.push_back('b');
  v.push_back('c'); v.mark_back();
  v.push_back('d');
  v.push_back('e'); v.mark_back();
  v.push_back('f');
  print(v);

  cout << "Clearing remaining elements." << endl;
  v.clear();
  print(v);

  v.push_back('a'); 
  v.push_back('b'); v.mark_back();
  v.push_back('c'); 
  v.push_back('d'); v.mark_back();
  v.push_back('e'); 
  v.push_back('f'); v.mark_back();
  print(v);

  cout << "Setting all marks." << endl;
  v.mark(v.begin()+0);
  v.mark(v.begin()+1);
  v.mark(v.begin()+2);
  v.mark(v.begin()+3);
  v.mark(v.begin()+4);
  v.mark(v.begin()+5);
  print(v);

  cout << "Unsetting marks back." << endl;
  v.unmark(v.begin()+1);
  v.unmark(v.begin()+3);
  v.unmark(v.begin()+5);
  print(v);

  cout << "Flipping all marks." << endl;
  v.flip_mark(v.begin()+0);
  v.flip_mark(v.begin()+1);
  v.flip_mark(v.begin()+2);
  v.flip_mark(v.begin()+3);
  v.flip_mark(v.begin()+4);
  v.flip_mark(v.begin()+5);
  print(v);

  cout << "Inserting g before c" << endl;
  v.insert( v.begin()+2, 'g' );
  print(v);

  cout << "Erasing g" << endl;
  v.erase( v.begin()+2 );
  print(v);

  cout << "Erasing b" << endl;
  v.erase( v.begin()+1 );
  print(v);

  vector<char> c;
  c.push_back('1'); c.push_back('2'); c.push_back('3');

  cout << "Inserting 1,2,3 before a." << endl;
  v.insert(v.begin(), c.begin(), c.end());
  print(v);

  cout << "Erasing 1,2,3." << endl;
  v.erase(v.begin(), v.begin()+3);
  print(v);

  cout << "Inserting b after a." << endl;
  v.insert(v.begin()+1, 'b');
  print(v);

  cout << "Erasing c--d" << endl;
  v.erase( v.begin()+2, v.end()-2 );
  print(v);

  cout << "Erasing e--f" << endl;
  v.erase( v.end()-2, v.end() );
  print(v);
}

