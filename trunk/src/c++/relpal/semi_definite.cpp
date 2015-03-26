//==========================================================================
//  File:       semi_definite.cpp
//
//  Author:     Yeunjoo Song
//
//  History:    Version 1.0  Initial implementation.             yjs Sep. 07
//
//  Notes:
//
//  Copyright (c) 2007 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "relpal/semi_definite.h"

namespace SAGE   {
namespace RELPAL {

//==========================================================================
//
//  Implements the base semi_definite class.
//
//==========================================================================

semi_definite::semi_definite(const matrix& U, const matrix& si, size_t t)
             : my_vector_U(U), my_sigma_i(si), my_trait_count(t), my_invalid(false)
{
  nfe = 0;

  // Initial estimate using U
  //
  my_current_estimate.resize_nofill(my_trait_count, my_trait_count);
  for( size_t i = 0, t = 0; i < my_trait_count; ++i )
  {
    for( size_t j = i; j < my_trait_count; ++j, ++t )
    {
      if( i == j  && my_vector_U(t,0) < 0.0 )
        my_current_estimate(i,j) = (-my_vector_U(t,0))/10.0;
      else
        my_current_estimate(i,j) = my_current_estimate(j,i) = my_vector_U(t,0);
    }
  }

#if 0
  print_matrix(my_vector_U, cout, "my_vector_U");
  print_matrix(my_current_estimate, cout, "(my_current_estimate(b)");
#endif

  compute_eigenvalue();
}

void
semi_definite::compute_eigenvalue()
{
#if 0
  cout << endl << "compute_eigenvalue()..." << endl;
  //print_matrix(my_current_estimate, cout, "b");
#endif

  matrix b = my_current_estimate;
  vector<double> w;

  int info = Eigen(b, w); //calculate_eigenvalue(b, w);

  my_current_e_value  = w;
  my_current_e_vector = b;
  my_current_info     = info;

  if( info != 0 )
    my_invalid = true;

  if( my_current_e_value[my_current_e_value.size()-1] < 0.0 )
    my_invalid = true;

  return;
}

void
semi_definite::initialize()
{
#if 0
  cout << endl << "initialize()..." << endl;
#endif

  matrix diag_eigen;
  diag_eigen.resize_fill(my_current_e_value.size(), my_current_e_value.size(), 0.0);

  for( size_t i = 0; i < my_current_e_value.size(); ++i )
  {
    if( my_current_e_value[i] < 0.0 )
      diag_eigen(i,i) = my_current_e_value[my_current_e_value.size()-1] / 10.;
    else
      diag_eigen(i,i) = my_current_e_value[i];
  }

  matrix bd;

  multiply(my_current_e_vector, diag_eigen, bd);
  XZT(bd, my_current_e_vector, my_initial_b);

#if 0
  print_matrix(my_current_e_vector, cout, "my_current_e_vector");
  print_matrix(diag_eigen, cout, "diag_eigen");
  print_matrix(bd, cout, "bd");
  print_matrix(my_initial_b, cout, "my_initial_b");
#endif

  return;
}

const matrix&
semi_definite::get_estimates() const
{
  return my_current_estimate;
}

bool
semi_definite::is_invalid() const
{
  return my_invalid;
}

double
semi_definite::evaluate(vector<double>& theta)
{
#if 0
  cout << nfe << " evaluate().." << endl;
  for( size_t i = 0; i < my_vector_U.rows(); ++i )
    cout << "  theta[" << i << "] = " << theta[i] << endl;
  cout << endl;
#endif

  ++nfe;

  return evaluate_derived(theta);
}

int
semi_definite::update_bounds(vector<double>& theta)
{
#if 0
  cout << nfe << " update_bounds().." << endl;
  for( size_t i = 0; i < my_vector_U.rows(); ++i )
    cout << "  theta[" << i << "] = " << theta[i] << endl;
  cout << endl;
#endif

  return update_bounds_derived(theta);
}

//==========================================================================
//
//  Implements the Cholesky decomposition method.
//
//==========================================================================

cholesky_decom::cholesky_decom(const matrix& U, const matrix& si, size_t t)
              : semi_definite(U, si, t)
{}


void
cholesky_decom::initialize()
{
  semi_definite::initialize();

  Cholesky(my_initial_b, my_current_estimate);

#if 0
  print_matrix(my_current_estimate, cout, "my_current_estimate(w)");

  matrix b;
  XZT(my_current_estimate, my_current_estimate, b);
  print_matrix(b, cout, "b");
#endif

  return;
}

double
cholesky_decom::evaluate_derived(vector<double>& theta)
{
  my_current_estimate.resize_fill(my_trait_count, my_trait_count, 0.0);

  for( size_t r = 0, i = 0; r < my_trait_count; ++r )
  {
    for( size_t c = r; c < my_trait_count; ++c, ++i )
    {
      my_current_estimate(r,c) = theta[i];
    }
  }

  matrix b_mat;
  XZT(my_current_estimate, my_current_estimate, b_mat);

  matrix b_vec;
  b_vec.resize_nofill(my_vector_U.rows(), 1);
  for( size_t r = 0, i = 0; r < my_trait_count; ++r )
  {
    for( size_t c = r; c < my_trait_count; ++c, ++i )
    {
      b_vec(i,0) = b_mat(r,c);
    }
  }

  matrix Ub = my_vector_U - b_vec;
  matrix UbTsigma_i, C;

  XTZ(Ub, my_sigma_i, UbTsigma_i);
  multiply(UbTsigma_i, Ub, C);

  double return_value = C(0,0);

#if 0
  print_matrix(my_current_estimate, cout, "my_current_estimate");
  print_matrix(b_mat, cout, "b_mat");
  print_matrix(my_vector_U, cout, "my_vector_U");
  print_matrix(b_vec, cout, "b_vec");
  print_matrix(Ub, cout, "U-b");
  print_matrix(my_sigma_i, cout, "my_sigma_i");
  print_matrix(UbTsigma_i, cout, "(U-b)'sigma^");
  print_matrix(C, cout, "(U-b)'sigma^(U-b)");
  cout << "evaluate return value = " << -return_value << endl;
#endif

  return -return_value;
}

int
cholesky_decom::update_bounds_derived(vector<double>& theta)
{
  return 0;
}

//==========================================================================
//
//  Implements the cutting plane algorithm.
//
//==========================================================================

cutting_plane::cutting_plane(const matrix& U, const matrix& si, size_t t)
             : semi_definite(U, si, t)
{}

bool
cutting_plane::is_positive_semidefinite() const
{
  if( my_current_info == 0 && my_current_e_value[0] >= 0. )
    return true;

  return false;
}

void
cutting_plane::reset_b()
{
#if 0
  cout << endl << "reset_b()... with intial_b" << endl;
#endif

  my_current_estimate = my_initial_b;

#if 0
  print_matrix(my_current_estimate, cout, "my_current_estimate");
#endif

  return;
}

void
cutting_plane::add_constraints()
{
#if 0
  cout << "adding constraints..., current v size = " << my_v.size() << endl;
  cout << "cur_e_val :" << endl;
  for( size_t i = 0; i < my_current_e_value.size(); ++i )
    cout << my_current_e_value[i] << " ";
  cout << endl;
  print_matrix(my_current_e_vector, cout, "cur_e_vec");
#endif

  for( size_t e = 0; e < my_current_e_value.size(); ++e )
  {
    //cout << "my_current_e_value[" << e << "] = " << my_current_e_value[e] << endl;
    if( my_current_e_value[e] < 0.0 )
    {
      matrix v;
      v.resize_nofill(my_current_e_vector.rows(), 1);
      //print_matrix(v, cout, "v");

      for( size_t i = 0; i < my_current_e_vector.rows(); ++i )
      {
        //cout << "v(" << i << "," << e << ")=" << my_current_e_vector(i,e) << " ";
        v(i,0) = my_current_e_vector(i,e);
      }
      //cout << endl;
      //print_matrix(v, cout, "v");

      my_v.push_back(v);
    }
  }

#if 0
  cout << "new v size = " << my_v.size() << endl;
#endif

  return;
}

double
cutting_plane::evaluate_derived(vector<double>& theta)
{
  matrix b;
  b.resize_nofill(my_vector_U.rows(), 1);

  for( size_t i = 0; i < b.rows(); ++i )
  {
     b(i,0) = theta[i];
  }

  my_current_estimate.resize_nofill(my_trait_count, my_trait_count);

  for( size_t r = 0, i = 0; r < my_trait_count; ++r )
  {
    for( size_t c = r; c < my_trait_count; ++c, ++i )
    {
      my_current_estimate(r,c) = my_current_estimate(c,r) = theta[i];
    }
  }

  matrix Ub = my_vector_U - b;
  matrix UbTsigma_i, C;

  XTZ(Ub, my_sigma_i, UbTsigma_i);
  multiply(UbTsigma_i, Ub, C);

  double return_value = C(0,0);

#if 0
  print_matrix(my_current_estimate, cout, "my_current_estimate");
  print_matrix(my_vector_U, cout, "my_vector_U");
  print_matrix(b, cout, "b");
  print_matrix(my_sigma_i, cout, "my_sigma_i");
  print_matrix(Ub, cout, "U-b");
  print_matrix(UbTsigma_i, cout, "(U-b)'sigma^");
  print_matrix(C, cout, "(U-b)'sigma^(U-b)");
  cout << "evaluate return value = " << -return_value << endl;
#endif

  return -return_value;
}

int
cutting_plane::update_bounds_derived(vector<double>& theta)
{
  matrix b;
  b.resize_nofill(my_trait_count, my_trait_count);

  for( size_t r = 0, i = 0; r < my_trait_count; ++r )
  {
    for( size_t c = r; c < my_trait_count; ++c, ++i )
    {
      b(r,c) = b(c,r) = theta[i];
    }
  }

  for( size_t i = 0; i < my_v.size(); ++i )
  {
    matrix v = my_v[i];

    matrix vTb, vTbv;

    XTZ(v, b, vTb);
    multiply(vTb, v, vTbv);

#if 0
    print_matrix(v, cout, "v");
    print_matrix(vTb, cout, "v'B");
    print_matrix(vTbv, cout, "v'Bv");
#endif

    double constraint_value = vTbv(0,0);

    if( constraint_value < 0. )
      return 1;
  }

  return 0;
}

//==========================================================================
/*
int
calculate_eigenvalue(matrix& b, vector<double>& w)
{
#if 0
  cout << endl << "compute_eigenvalue()..." << endl;
  //print_matrix(my_current_estimate, cout, "b");
#endif

  vector<double> work;

  char jobz = 'V';
  char uplo = 'U';
  int n     = b.rows();
  int lda   = b.lda();
  int lwork = 8*(3*n-1);
  int info;

  w.resize(n);
  work.resize(lwork);

  F77_CALL(dsyev)(&jobz, &uplo, &n, b.raw_storage(), &lda, &w[0], &work[0], &lwork, &info);

#if 0
  cout << "info = " << info << endl;

  if( info == 0 )
  {
    cout << "eigenvalue :" << endl;
    for( int i = 0; i < n; ++i )
      cout << w[i] << " ";
    cout << endl;
    //print_matrix(b, cout, "eigenvector");

    for( int i = 0; i < n; ++i )
      if( fabs(w[i]) < 1.0e-10 )
        w[i] = 0.0;

    cout << "adjusted eigenvalue :" << endl;
    for( int i = 0; i < n; ++i )
      cout << w[i] << " ";
    cout << endl;
  }
#endif

  return info;
}
*/

} // end of namespace RELPAL
} // end of naemspace SAGE
