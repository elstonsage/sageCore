#include "relpal/two_level_calculator.h"

namespace SAGE   {
namespace RELPAL {

relpal_score::relpal_score(size_t m)
            : calculator_base(m)
{
  reset(m);
}

relpal_score::~relpal_score()
{}

void 
relpal_score::reset(size_t m)
{
  parameters = m;

  observation_count = 0;
  cluster_count     = 0;
  svd_return_code   = 1;

  U.clear();
  sigma.clear();
  DWD.clear();

  U.resize_fill(m, 1, 0.0);
  sigma.resize_fill(m, m, 0.0);
  DWD.resize_fill(m, m, 0.0);
  BPB.resize_fill(m, m, 0.0);

  Us.resize(0);
  DWDs.resize(0);
}

void 
relpal_score::add_block(const matrix& D, const trimatrix& W, const matrix& S, trimatrix& Wi)
{
  if( !U || !sigma )
    return;

  assert( Wi.size() && D.rows() == Wi.size() && D.rows() == S.rows() );
  assert( D && D.cols() == parameters );
  assert( S && S.rows() > 0 );
/*
  trimatrix Wi;
  int info = get_tri_inverse(W, Wi);

  if( info != 0 )
  {
    U.setstate(matrix::badbit);
    return;
  }
*/
  observation_count += S.rows();
  cluster_count     += 1;

  compute_U_sigma(D, Wi, S);

  return;
}

void 
relpal_score::add_block_kron(const matrix& D, const trimatrix& b, const matrix& S)
{
  if( !U || !sigma )
    return;

  assert( D && D.cols() == parameters );
  assert( S && S.rows() > 0 );

  trimatrix Wi;
  get_tri_kron(b, Wi);

  observation_count += S.rows();
  cluster_count     += 1;

  compute_U_sigma(D, Wi, S);

  return;
}

void 
relpal_score::add_block_IBD(const matrix& B, const trimatrix& varIBD)
{
  if( !U || !BPB )
    return;

  assert( B && B.rows() == parameters );

  matrix BP, BPBk;

  XtriZ(B, varIBD, BP);

  XZT(BP, B, BPBk);

  BPB += BPBk;
#if 0
  print_matrix(BP, cout, "Bvar(IBD)");
  print_matrix(BPBk, cout, "Bvar(IBD)BT");
  print_matrix(BPB, cout, "BPB");
#endif

  return;
}

bool
relpal_score::compute()
{
  U_star.clear();
  sigma_sw.clear();
  sigma_sw_i.clear();
  sigma_sw_sqrt.clear();
  sigma_na.clear();
  sigma_na_i.clear();
  sigma_na_sqrt.clear();
  sigma_al.clear();
  sigma_al_i.clear();
  sigma_al_sqrt.clear();
  sigma_ib.clear();
  sigma_ib_i.clear();
  sigma_ib_sqrt.clear();
  T_sw.clear();
  T_na.clear();
  T_al.clear();
  T_ib.clear();
  df_na = df_sw = df_al = df_ib = 0;

  // ' : transpose
  // ^ : inverse
  //
  // U*                = (D' W^ D)^ U = (D' W^ D)^ (D' W^ S)
  //
  // sigma_naive       = (D' W^ D)^
  //
  // sigma_sandwich    = (D' W^ D)^ sigma (D' W^ D)^
  //                   = (D' W^ D)^ [(D' W^ S)(S' W^ D)] (D' W^ D)^
  //                   = (D' W^ D)^ [UU'] (D' W^ D)^
  //
  // sigma_alternative = (D' W^ D)^ [(U - E(U))(U - E(U))] (D' W^ D)^
  //              E(U) = (D' W^ D) U*
  //
  // sigma_IBD         = B var(phi) B'
#if 0
  cout << "cluster_count     = " << cluster_count << endl;
  cout << "observation_count = " << observation_count << endl;
  cout << "Us size = " << Us.size() << ", DWDs size = " << DWDs.size() << endl;
  cout << "var option = " << is_naive_var() << is_sandwich_var() << is_alternative_var() << endl;
#endif

  // Naive estimate of Var(U*)
  //
  svd.inverse_of(DWD, sigma_na);

  if( !sigma_na )
  {
    reset();
    T_na.setstate(matrix::badbit);
    return false;
  }

  multiply(sigma_na, U, U_star);

  df_na = df_sw = df_al = df_ib = U_star.rows();

#if 0
  print_matrix(U, cout, "U");
  //print_matrix(sigma, cout, "sigma");

  //print_matrix(DWD, cout, "DWD");
  //print_matrix(sigma_na, cout, "DWD^ = sigma_naive");

  print_matrix(U_star, cout, "U*");
  cout << "df = " << df_na << endl;
#endif

  if( is_naive_var() )
  {
    sigma_na_i = DWD;
    svd.inverse_sqrt(sigma_na_sqrt);

    // T_na = U*'sigma_na^U*
    //
    matrix Utsigma_na_i;

    XTZ(U_star, sigma_na_i, Utsigma_na_i);
    multiply(Utsigma_na_i, U_star, T_na);

#if 0
  print_matrix(sigma_na, cout, "sigma_na");
  print_matrix(sigma_na_i, cout, "sigma_na^-1");
  print_matrix(sigma_na_sqrt, cout, "sigma_na^0.5");
  print_matrix(T_na, cout, "T_na");
#endif
  }

  if( is_sandwich_var() )
  {
    // Robust sandwich estimate of Var(U*)
    //
    matrix DWDi_sigma;

    multiply(sigma_na, sigma, DWDi_sigma);
    multiply(DWDi_sigma, sigma_na, sigma_sw);

    int return_code = svd.compute(sigma_sw);

    if( return_code == 0 ) // Not invertable.
    {
      return false;
    }

    if( svd.positive_diagonal.size() < sigma_sw.cols() )
    {
      svd.general_inverse(sigma_sw_i);
      svd.general_sqrt(sigma_sw_sqrt);
      df_sw = svd.positive_diagonal.size();
    }
    else
    {
      svd.inverse(sigma_sw_i);
      svd.sqrt(sigma_sw_sqrt);
    }

    // T_sw = U*'sigma_sw^U*
    //
    matrix Utsigma_sw_i;

    XTZ(U_star, sigma_sw_i, Utsigma_sw_i);
    multiply(Utsigma_sw_i, U_star, T_sw);

#if 0
  print_matrix(sigma_sw, cout, "sigma_sandwich");
  print_matrix(sigma_sw_i, cout, "sigma_sw^-1");
  print_matrix(sigma_sw_sqrt, cout, "sigma_sw^0.5");
  print_matrix(T_sw, cout, "T_sw");
  cout << "df_sw = " << df_sw << endl;
#endif
  }

  if( is_alternative_var() )
  {
    // Alternative estimate of Var(U*)
    //
    matrix UEUUEUT;
    UEUUEUT.resize_fill(U.rows(), U.rows(), 0.0);
    for( size_t k = 0; k < Us.size(); ++k )
    {
      matrix EU_k, UEU_k, UEUUEUT_k;
      multiply(DWDs[k], U_star, EU_k);
      UEU_k = Us[k] - EU_k;
      XXT(UEU_k, UEUUEUT_k);
      UEUUEUT += UEUUEUT_k;
#if 0
  print_matrix(EU_k, cout, "E(U)_k");
  print_matrix(UEU_k, cout, "U-E(U)_k");
  print_matrix(UEUUEUT_k, cout, "(U-E(U))(U-E(U))'_k");
#endif
    }

#if 0
  print_matrix(UEUUEUT, cout, "(U-E(U))(U-E(U))'");
#endif

    matrix DWDi_sigma_al;

    multiply(sigma_na, UEUUEUT, DWDi_sigma_al);
    multiply(DWDi_sigma_al, sigma_na, sigma_al);

    int return_code = svd.compute(sigma_al);

    if( return_code == 0 ) // Not invertable.
    {
      return false;
    }

    if( svd.positive_diagonal.size() < sigma_al.cols() )
    {
      svd.general_inverse(sigma_al_i);
      svd.general_sqrt(sigma_al_sqrt);
      df_al = svd.positive_diagonal.size();
    }
    else
    {
      svd.inverse(sigma_al_i);
      svd.sqrt(sigma_al_sqrt);
    }

    // T_al = U*'sigma_al^U*
    //
    matrix Utsigma_al_i;

    XTZ(U_star, sigma_al_i, Utsigma_al_i);
    multiply(Utsigma_al_i, U_star, T_al);

#if 0
  print_matrix(sigma_al, cout, "sigma_alternative");
  print_matrix(sigma_al_i, cout, "sigma_al^-1");
  print_matrix(sigma_al_sqrt, cout, "sigma_al^0.5");
  print_matrix(T_al, cout, "T_al");
  cout << "df_sw = " << df_sw << endl;
#endif
  }

  if( is_IBD_var() )
  {
    // IBD variance conditional on trait as an estimate of Var(U*)
    //
    matrix DWDi_sigma;

    multiply(sigma_na, BPB, DWDi_sigma);
    multiply(DWDi_sigma, sigma_na, sigma_ib);

    int return_code = svd.compute(sigma_ib);

    if( return_code == 0 ) // Not invertable.
    {
      return false;
    }

    if( svd.positive_diagonal.size() < sigma_ib.cols() )
    {
      svd.general_inverse(sigma_ib_i);
      svd.general_sqrt(sigma_ib_sqrt);
      df_ib = svd.positive_diagonal.size();
    }
    else
    {
      svd.inverse(sigma_ib_i);
      svd.sqrt(sigma_ib_sqrt);
    }

    // T_ib = U*'sigma_ib^U*
    //
    matrix Utsigma_ib_i;

    XTZ(U_star, sigma_ib_i, Utsigma_ib_i);
    multiply(Utsigma_ib_i, U_star, T_ib);

#if 0
  print_matrix(BPB, cout, "BPB");
  print_matrix(sigma_ib, cout, "sigma_IBD");
  print_matrix(sigma_ib_i, cout, "sigma_ib^-1");
  print_matrix(sigma_ib_sqrt, cout, "sigma_ib^0.5");
  print_matrix(T_ib, cout, "T_ib");
  cout << "df_ib = " << df_ib << endl;
#endif
  }

  return true;
}

void
relpal_score::compute_U_sigma(const matrix& D, const trimatrix& Wi, const matrix& S)
{
  assert( Wi.size() && Wi.size() == D.rows()  && D.rows() == S.rows() );

  // ' : transpose
  // ^ : inverse
  //
  // DW     = D' W^
  // DWS    = D' W^ S = U
  // SWD    = S' W^ D = U'
  // DWSSWD = (D' W^ S)(S' W^ D) = sigma = UU'
  matrix DW, DWS, DWSSWD;

  // Compute U & robust var(U) = sigma
  XTtriZ(D, Wi, DW);
  multiply(DW, S, DWS);

  U += DWS;

  //transpose(DWS, SWD);
  //multiply(DWS, SWD, DWSSWD);
  XXT(DWS, DWSSWD);

  sigma += DWSSWD;

  // Compute D'W^D
  matrix DWDk;

  multiply(DW, D, DWDk);

  DWD += DWDk;

  if( is_alternative_var() )
  {
    Us.push_back(DWS);
    DWDs.push_back(DWDk);
  }

#if 0
  print_trimatrix_first10(Wi, cout, "Wi");
  print_matrix(DWS, cout, "DWS");
  //print_matrix(DWSSWD, cout, "DWSSWD");
  //print_matrix(DWDk, cout, "DWDk");
#endif

  return;
}

} // end of namespace RELPAL
} // end of namespace SAGE
