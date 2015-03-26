#include <iostream>
#include "numerics/svd.h"
#include "numerics/linalg.h"

using std::cout;
using std::endl;

int main()
{
  Matrix2D<double> a(4,4),aa,aaa;
  a(0,0) = 1.0;
  a(1,1) = 2.0;
  a(2,2) = 4.0;
  a(1,2) = 0.5;
  a(2,1) = 0.5;
  a(3,2) = 2.0;
  a(3,0) = 1.0;
  a(3,3) = 0.25;

  Matrix2D<double> u,v,w;

  singular_value_decomposition(a,u,w,v);

  cout << "A=";
  print_matrix(a,cout);
  cout << endl;

  cout << "U=";
  print_matrix(u,cout);
  cout << endl;

  cout << "W=";
  print_matrix(w,cout);
  cout << endl;

  cout << "V=";
  print_matrix(v,cout);
  cout << endl;
 
  Matrix2D<double> temp, temp2, temp3;

  Multiply(u,w, temp);
  Multiply(temp,v,temp3);

  cout << "A'=";
  print_matrix(temp3,cout);
  cout << endl;

  temp = inverse(u,w,v);

  cout << "A^-1=";
  print_matrix(temp,cout);
  cout << endl;

  Multiply(a,temp,temp2);

  cout << "I=";
  print_matrix(temp2,cout);
  cout << endl;

  cout << "det = " << determinant(u,w,v) << endl;

  a.resize(6,6);
  int x = 0;
  for(int i = 0; i < 4; ++i)
    for(int j = 0; j < i; ++j)
    {
      int y = 0;
      for(int k = 0; k < 4; ++k)
        for(int l = 0; l < k; ++l)
        {
          if( (i == k && j == l) || (i == l && j == k) )
            a(x,y) = 1.29480;
          else if( i == k || j == l || i == l || j == k )
            a(x,y) = -0.18299;
          else
            a(x,y) = -0.18047;
          ++y;
        }
      ++x;
    }

  SAGE::cholesky(a, aa);
  Transpose(aa);
  singular_value_decomposition(a,u,w,v);

  cout << "A=";
  print_matrix(a,cout);
  cout << endl;

  cout << "U=";
  print_matrix(u,cout);
  cout << endl;

  cout << "W=";
  print_matrix(w,cout);
  cout << endl;

  cout << "V=";
  print_matrix(v,cout);
  cout << endl;

  cout << "det = " << determinant(u,w,v) << endl;

  Multiply(u,w,a);
  Multiply(a,v,aaa);
  aaa -= aa;

  cout << "Z=";
  print_matrix(aaa,cout);
  cout << endl;
  
}
