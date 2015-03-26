//============================================================================
// File:      test_matrix.cpp
//                                                                          
// Author:    Yeunjoo Song
//                                                                          
// History:   Initial version                                        Oct. 2007
//                                                                          
// Notes:     Tests the various matrix inversion methods.
//                                                                      
// Copyright (c) 2007 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "relpal/definitions.h"
#include "numerics/trimatrix.h"
#include "numerics/clapack.h"

using namespace std;
using namespace SAGE;
using namespace RELPAL;

void get_trimatrix_kron1(const trimatrix& b, trimatrix& kb)
{
  kb.resize(b.size()*b.size());

  for( size_t r1 = 0, r = 0; r1 < b.size(); ++r1 )
  {
    for( size_t r2 = 0; r2 < b.size(); ++r2, ++r )
    {
      for( size_t c1 = 0, c = 0; c1 < b.size(); ++c1 )
      {
        for( size_t c2 = 0; c2 < b.size(); ++c2, ++c )
          kb(r,c) = 0.5 * b(r1,c1) * b(r2,c2);  
      }
    }
  }

  return;
}

void get_trimatrix_kron2(const trimatrix& b, trimatrix& kb)
{
  kb.resize(b.size()*b.size());

  for( size_t r1 = 0, r = 0; r1 < b.size(); ++r1 )
  {
    for( size_t c1 = r1, c = 0; c1 < b.size(); ++c1 )
    {
      for( size_t r2 = 0; r2 < b.size(); ++r2, ++r )
      {
        for( size_t c2 = 0; c2 < b.size(); ++c2, ++c )
          kb(r,c) = 0.5 * b(r1,c1) * b(r2,c2);  
      }
    }
  }

  return;
}

void test_matrix_inversion(ofstream& out)
{
  matrix A;

  A.resize_nofill(5, 5);

  A(0,0) = 1.72424;  A(0,1) = 0.941884;  A(0,2) = 0.326347;  A(0,3) = 0.24562;  A(0,4) = 0.695199;
  A(1,0) = 0.941884; A(1,1) = 0.547951;  A(1,2) = 0.24562;   A(1,3) = 0.161439; A(1,4) = 0.452352;
  A(2,0) = 0.326347; A(2,1) = 0.24562;   A(2,2) = 1.72424;   A(2,3) = 0.941884; A(2,4) = 0.695199;
  A(3,0) = 0.24562;  A(3,1) = 0.161439;  A(3,2) = 0.941884;  A(3,3) = 0.547951; A(3,4) = 0.452352;
  A(4,0) = 0.695199; A(4,1) = 0.452352;  A(4,2) = 0.695199;  A(4,3) = 0.452352; A(4,4) = 1.72424;

/*
  A.resize_nofill(10, 10);

  A(0,0) = 1.72424;  A(0,1) = 0.941884;  A(0,2) = 0.326347;  A(0,3) = 0.24562;  A(0,4) = 0.695199;
  A(0,5) = 0.452352; A(0,6) = 0.326347;  A(0,7) = 0.24562;   A(0,8) = 0.326347; A(0,9) = 0.24562;

  A(1,0) = 0.941884; A(1,1) = 0.547951;  A(1,2) = 0.24562;   A(1,3) = 0.161439; A(1,4) = 0.452352;
  A(1,5) = 0.290614; A(1,6) = 0.24562;   A(1,7) = 0.161439;  A(1,8) = 0.24562;  A(1,9) = 0.161439;

  A(2,0) = 0.326347; A(2,1) = 0.24562;   A(2,2) = 1.72424;   A(2,3) = 0.941884; A(2,4) = 0.695199;
  A(2,5) = 0.452352; A(2,6) = 0.326347;  A(2,7) = 0.24562;   A(2,8) = 0.326347; A(2,9) = 0.24562;

  A(3,0) = 0.24562;  A(3,1) = 0.161439;  A(3,2) = 0.941884;  A(3,3) = 0.547951; A(3,4) = 0.452352;
  A(3,5) = 0.290614; A(3,6) = 0.24562;   A(3,7) = 0.161439;  A(3,8) = 0.24562;  A(3,9) = 0.161439;

  A(4,0) = 0.695199; A(4,1) = 0.452352;  A(4,2) = 0.695199;  A(4,3) = 0.452352; A(4,4) = 1.72424;
  A(4,5) = 0.941884; A(4,6) = 0.695199;  A(4,7) = 0.452352;  A(4,8) = 0.695199; A(4,9) = 0.452352;

  A(5,0) = 0.452352; A(5,1) = 0.290614;  A(5,2) = 0.452352;  A(5,3) = 0.290614; A(5,4) = 0.941884;
  A(5,5) = 0.547951; A(5,6) = 0.452352;  A(5,7) = 0.290614;  A(5,8) = 0.452352; A(5,9) = 0.290614;

  A(6,0) = 0.326347; A(6,1) = 0.24562;   A(6,2) = 0.326347;  A(6,3) = 0.24562;  A(6,4) = 0.695199;
  A(6,5) = 0.452352; A(6,6) = 1.72424;   A(6,7) = 0.941884;  A(6,8) = 0.326347; A(6,9) = 0.24562;

  A(7,0) = 0.24562;  A(7,1) = 0.161439;  A(7,2) = 0.24562;   A(7,3) = 0.161439; A(7,4) = 0.452352;
  A(7,5) = 0.290614; A(7,6) = 0.941884;  A(7,7) = 0.547951;  A(7,8) = 0.24562;  A(7,9) = 0.161439;

  A(8,0) = 0.326347; A(8,1) = 0.24562;   A(8,2) = 0.326347;  A(8,3) = 0.24562;  A(8,4) = 0.695199;
  A(8,5) = 0.452352; A(8,6) = 0.326347;  A(8,7) = 0.24562;   A(8,8) = 1.72424;  A(8,9) = 0.941884;

  A(9,0) = 0.24562;  A(9,1) = 0.161439;  A(9,2) = 0.24562;   A(9,3) = 0.161439; A(9,4) = 0.452352;
  A(9,5) = 0.290614; A(9,6) = 0.24562;  A(9,7) = 0.161439;  A(9,8) = 0.941884;  A(9,9) = 0.547951;
*/

  out << "== kron first, SVD inverse" << endl;
  matrix AkA = kron(A, A);

  matrix AkA2 = 2.0 * AkA;

  SAGE::SVD svd;

  matrix AkA_i = svd.inverse_of(AkA2);

  print_matrix(A, out, "A");
  print_matrix(AkA, out, "AkA");
  print_matrix(AkA2, out, "2AxA");
  print_matrix(AkA_i, out, "(2AxA)^-1");

  out << "== SVD inverse first, kron" << endl;

  matrix A_i = svd.inverse_of(A);

  matrix A_ikA_i = kron(A_i, A_i);  

  matrix A_ikA_i2 = 0.5 * A_ikA_i;

  print_matrix(A_i, out, "A^-1");
  print_matrix(A_ikA_i, out, "A^-1kA^-1");
  print_matrix(A_ikA_i2, out, "0.5(A^-1kA^-1)");

  out << "== symmetric inverse (matrix)" << endl;

  matrix A1 = A;

  vector<double> work;
  vector<int> w;

  //char jobz = 'V';
  char uplo = 'U';
  int n     = A1.rows();
  int lda   = A1.lda();
  int lwork = n * n;;
  int info;

  w.resize(n);
  work.resize(lwork);

  F77_CALL(dsytrf)(&uplo, &n, A1.raw_storage(), &lda, &w[0], &work[0], &lwork, &info);
  F77_CALL(dsytri)(&uplo, &n, A1.raw_storage(), &lda, &w[0], &work[0], &info);

  print_matrix(A1, out, "A^-1");

  out << "== symmetric inverse (triangle matrix)" << endl;

  trimatrix tri;
  tri.resize(A.rows());

  for( size_t i = 0; i < A1.rows(); ++i )
    for( size_t j = i; j < A1.cols(); ++j )
      tri(i, j) = A(i, j);

  print_trimatrix(tri, out, "Tri A");

  work.resize(n);

  F77_CALL(dsptrf)(&uplo, &n, tri.raw_storage(), &w[0], &info);
  F77_CALL(dsptri)(&uplo, &n, tri.raw_storage(), &w[0], &work[0], &info);

  print_trimatrix(tri, out, "Tri A^-1");

  trimatrix Ak;
  get_trimatrix_kron1(tri, Ak);
  print_trimatrix(Ak, out, "0.5(Tri A^-1)kron(Tri A^-1)");

  return;
}

void test_trimatrix_kron(ofstream& out)
{
  trimatrix b, wi, wi2;

  b.resize(3);

  b(0,0) = 1.0;
  b(0,1) = 2.0;
  b(0,2) = 3.0;
  b(1,1) = 4.0;
  b(1,2) = 5.0;
  b(2,2) = 6.0;

  get_trimatrix_kron1(b, wi);
  get_trimatrix_kron2(b, wi);

  print_trimatrix(b, out, "b");
  print_trimatrix(wi, out, "wi");
  print_trimatrix(wi2, out, "wi2");

  return;
}

int main(int argc, char* argv[])
{

  ofstream out;
  out.open("test_matrix.out");
  
  if(! out)
  {
    cout << "Cannot open output file: test_matirx.out.  Exiting ..." << endl;
    exit(EXIT_FAILURE);
  }

  test_matrix_inversion(out);

  exit(EXIT_SUCCESS);
}

