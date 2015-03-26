#include <iostream>
#include <vector>
#include "numerics/matfunctions.h"
#include "numerics/linalg.h"

using namespace SAGE;

using std::cout;
using std::endl;

int main()
{
  Matrix2D<double> A(1,1);
  Matrix2D<double> Q,R;
  
  A(0,0) = 2.0;

  cout << "A=";
  print_matrix(A,cout);
  
  qr(A, Q, R);

  cout << endl << "Q=";
  print_matrix(Q,cout);

  cout << endl << "R=";
  print_matrix(R,cout);

  A.resize(3,3);
  
  A(0,0) = 1.0;
  A(0,1) = 0.5;
  A(0,2) = 0.5;

  A(1,0) = 0.5;
  A(1,1) = 1.0;
  A(1,2) = 0.5;

  A(2,0) = 0.5;
  A(2,1) = 0.5;
  A(2,2) = 1.0;

  cout << endl << endl;
  cout << "A=";
  print_matrix(A,cout);

  Matrix2D<double> G, AA;
  cholesky(A,G);
  cout << "G=";
  print_matrix(G,cout);
  XTX(G, AA);
  AA -= A;
  cout << "AA=";
  print_matrix(AA,cout);


  msqrt(A, AA);
  Multiply(AA,AA,G);
  G -= A;
  cout << "(A^.5)^2-A=";
  print_matrix(G,cout);
  
  
  qr(A, Q, R);

  Multiply(Q,R,AA);
  AA -= A;
  cout << "Q*R-A=";
  print_matrix(AA,cout);

  A.resize(4,4);
  
  A(0,0) = 1.0;
  A(0,1) = 0.5;
  A(0,2) = 0.5;
  A(0,3) = 0.5;

  A(1,0) = 0.5;
  A(1,1) = 1.0;
  A(1,2) = 0.5;
  A(1,3) = 3.5;

  A(2,0) = 0.5;
  A(2,1) = 0.5;
  A(2,2) = 1.0;
  A(2,3) = 2.0;

  A(3,0) = 2.5;
  A(3,1) = 1.5;
  A(3,2) = 3.0;
  A(3,3) = 1.0;

  cout << endl << endl;
  cout << "A=";
  print_matrix(A,cout);
  
  qr(A, Q, R);

  Multiply(Q,R,AA);

  Multiply(Q,R,AA);
  AA -= A;
  cout << "Q*R-A=";
  print_matrix(AA,cout);

}

