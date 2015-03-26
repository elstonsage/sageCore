#include <functional>
#include <algorithm>
#include <vector>
#include <iostream>
#include <iomanip>
#include <string>

#include "numerics/fmatrix.h"

using namespace std;
using namespace SAGE;

int main1()
{
  FortranMatrix<double> m(2,2), mm;
  m(0,0) = 4.0;
  m(0,1) = 1.5;
  m(1,0) = 4.5;
  m(1,1) = 5.0;

  print_matrix(  m, cout, "m"  );
  print_matrix( -m, cout, "-m" );
  print_matrix(transpose(m), cout, "m'");
  print_matrix(m+m, cout, "m+m");
  print_matrix(m-m, cout, "m-m");
  print_matrix(m*m, cout, "m*m");

  mm = transpose(m);
  mm += m + m - m;
  mm -= m;
  mm.transpose();
  mm *= m;
  mm -= m*m;
  print_matrix(  mm, cout, "Z"  );

  m.resize(4,2);

  print_matrix(  m, cout, "m"  );
  print_matrix( -m, cout, "-m" );
  print_matrix(transpose(m), cout, "m'");
  print_matrix(m+m, cout, "m+m");
  print_matrix(m-m, cout, "m-m");
  print_matrix(m*transpose(m), cout, "m*m'");
  print_matrix(transpose(m)*m, cout, "m'*m");

  mm = m;
  mm += m + m - m;
  mm -= m;
  print_matrix(  mm, cout, "Z"  );
  mm *= transpose(m);
  print_matrix(  mm, cout, "Z"  );
  mm -= m*transpose(m);
  print_matrix(  mm, cout, "Z"  );
  return 0;
}

int main()
{
  FortranMatrix<double> m(2,2), mm(2,2), mmm;
  m(0,0) = 5.0;
  m(0,1) = 4.5;
  m(1,0) = 4.5;
  m(1,1) = 6.0;

  cout << "LU inverse(m): " << endl;
  LUinverse(m,mm);
  print_matrix( m,    cout, "   m"         );
  print_matrix( mm,   cout, "   m^1"       );
  print_matrix( mm*m, cout, "   m^1*m = I" );
  print_matrix( m*mm, cout, "   m*m^1 = I" );

  cout << "SVD inverse(m): " << endl;
  SVDinverse(m,mm);
  print_matrix(  m, cout, "m"  );
  print_matrix( mm, cout, "m^1"  );
  print_matrix( mm*m, cout, "m^1*m = I"  );
  print_matrix( m*mm, cout, "m*m^1 = I"  );

  cout << "Symmetric inverse(m): " << endl;
  SYMinverse(m,mm);
  print_matrix(  m, cout, "m"  );
  print_matrix( mm, cout, "m^1"  );
  print_matrix( mm*m, cout, "m^1*m = I"  );
  print_matrix( m*mm, cout, "m*m^1 = I"  );

  cout << "Cholesky facorization of m: " << endl;
  Cholesky(m, mm, false);
  print_matrix(  m, cout, "m"  );
  print_matrix( mm, cout, "chol(m)"  );
  print_matrix(mm*transpose(mm)-m, cout, "Z");

  cout << "Inverse of triangular matrix mm: " << endl;
  mmm = TRIinverse(mm,false);
  print_matrix(mm*mmm, cout, "I");

  cout << "Matrix square-root of m: " << endl;
  mm = Msqrt(m);
  print_matrix(  m, cout, "m"  );
  print_matrix( mm, cout, "Msqrt(m)"  );
  print_matrix(mm*mm-m, cout, "Z");

  cout << "Inverse Matrix square-root of m: " << endl;
  mm = Minverse_sqrt(m);
  print_matrix(  m, cout, "m"  );
  print_matrix( mm, cout, "Minverse_sqrt(m)"  );
  print_matrix(LUinverse(mm*mm)-m, cout, "Z");

  SVD svd(m);
  cout << "SVD Matrix square-root of m: " << endl;
  svd.sqrt(mm);
  print_matrix(mm*mm-m, cout, "Z");

  cout << "Inverse Matrix square-root of m: " << endl;
  svd.inverse_sqrt(mm);
  print_matrix(LUinverse(mm*mm)-m, cout, "Z");

  cout << "Transpose test for row/column vector optimizations" << endl;
  m.resize_fill(5,1);
  m(0,0) = 0;
  m(1,0) = 1;
  m(2,0) = 2;
  m(3,0) = 3;
  m(4,0) = 4;

  print_matrix(m, cout, "m");
  print_matrix(m.transpose(), cout, "m'");
  print_matrix(m.transpose(), cout, "m''");
  print_matrix(m.transpose(), cout, "m'''");
  print_matrix(m.transpose(), cout, "m''''");

  FortranMatrix<double> U(3,3,0), UU;
  U(0,0) = 1; U(0,1) = 2; U(0,2) = 3;
              U(1,1) = 4; U(1,2) = 5;
                          U(2,2) = 6;

  FortranMatrix<double> B(3,1), X;
  B(0,0) = 0; B(1,0) = -8; B(2,0) = 6;

  print_matrix(U, cout, "U");
  print_matrix(B, cout, "B");

  UU=U;
  X=B;
  back_substitute(UU,X);

  print_matrix(U*X-B, cout, "Z");

  B(0,0) = 6; B(1,0) = -8; B(2,0) = 1;
  U.transpose();
  UU=U;
  X=B;
  print_matrix(U, cout, "L");
  print_matrix(B, cout, "B");
  forward_substitute(UU,X);

  print_matrix(U*X-B, cout, "Z");

}
