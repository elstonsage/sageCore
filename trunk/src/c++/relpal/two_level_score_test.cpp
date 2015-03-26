#include "maxfunapi/maxfunapi.h"
#include "relpal/two_level_test.h"

namespace SAGE   {
namespace RELPAL {

two_level_score_test::two_level_score_test(cerrorstream& err)
                    : two_level_base(err), my_valid_null(false), my_reliable_score(false)
{
  my_null_IBD_covariances.resize(0);
  my_IBD_covariances.resize(0);

  my_correction_na = 0.0;
  my_correction_sw = 0.0;
  my_correction_al = 0.0;
  my_correction_ib = 0.0;

  my_emp_pvalue_na = QNAN;
  my_emp_pvalue_sw = QNAN;
  my_emp_pvalue_al = QNAN;
  my_emp_pvalue_ib = QNAN;

  my_rep_count_na = 0;
  my_rep_count_sw = 0;
  my_rep_count_al = 0;
  my_rep_count_ib = 0;
}

void
two_level_score_test::compute_null_weight(const regression_model& null_model,
                                          regression_result& null_result)
{
  if( my_debug_out )
    debug_out() << endl << "# two_level_score_test::compute_null_weight()..." << endl;

  my_valid_null = false;
  my_null_weights.resize(0);
  my_residuals.resize(0);

  two_level_regression null_reg(errors);

  null_reg.set_pairs(my_pairs);
  null_reg.set_model(null_model);         
  null_reg.set_data(my_data);

  null_reg.set_data_out(my_data_out);
  null_reg.set_debug_out(my_debug_out);

  my_valid_null = null_reg.get_null_weight(my_null_weights, my_residuals);

  if( my_valid_null )
  {
    null_result.ind_beta     = null_reg.get_gls_1().beta;
    null_result.ind_variance = null_reg.get_gls_1().Variance;

    null_result.ped_beta     = null_reg.get_gls_2().beta;
    null_result.ped_variance = null_reg.get_gls_2().Variance;
  }

  return;
}

bool
two_level_score_test::do_two_level_score_test(score_test_result& score_result)
{
  // Set variance options.
  //
  my_score_2.set_naive_var(my_model.get_analysis_options().naive_variance);
  my_score_2.set_sandwich_var(my_model.get_analysis_options().sandwich_variance);
  my_score_2.set_alternative_var(my_model.get_analysis_options().alternative_variance);
  my_score_2.set_IBD_var(my_model.get_analysis_options().IBD_variance);

  compute_covariate_means();

  if( my_debug_out )
  {
    debug_out() << endl << "# do_two_level_score_test..." << endl;
    //my_model.dump_model(*my_pairs, cout);
  }

  if( my_score_2.is_IBD_var() )
  {
    compute_covariance_matrix();
    compute_covariance_matrix_given_m();
  }

  // Compute score statistic
  //
  if( !do_pedigree_level() )
    return false;

  score_result.U = get_score_2().U_star;

  if(    (my_score_2.df_sw < get_score_2().U_star.rows())
      || (my_score_2.df_al < get_score_2().U_star.rows())
      || (my_score_2.df_ib < get_score_2().U_star.rows()) )
    my_reliable_score = false;
  else
    my_reliable_score = true;

  if( my_score_2.is_naive_var() )
  {
    // Compute correction
    //
    my_correction_na = compute_correction(my_score_2.U_star, my_score_2.sigma_na_i);

    if( my_debug_out )
      debug_out() << "my_correction_naive = " << my_correction_na << endl;

    // Compute P-value
    //
    my_emp_pvalue_na = compute_empirical_p_value(my_rep_count_na,
                                                 my_score_2.df_na,
                                                 my_score_2.T_na(0,0),
                                                 my_correction_na,
                                                 my_score_2.sigma_na_i,
                                                 my_score_2.sigma_na_sqrt);

    if( my_debug_out )
      debug_out() << "my_emp_pvalue_naive = " << my_emp_pvalue_na << endl;

    score_result.sigma.push_back(my_score_2.sigma_na);
    score_result.T.push_back(my_score_2.T_na(0,0));
    score_result.correction.push_back(my_correction_na);
    score_result.emp_pvalue.push_back(my_emp_pvalue_na);
    score_result.rep_count.push_back(my_rep_count_na);
  }

  if( my_score_2.is_sandwich_var() )
  {
    my_correction_sw = compute_correction(my_score_2.U_star, my_score_2.sigma_sw_i);

    if( my_debug_out )
      debug_out() << "my_correction_sandwich = " << my_correction_sw << endl;

    my_emp_pvalue_sw = compute_empirical_p_value(my_rep_count_sw,
                                                 my_score_2.df_sw,
                                                 my_score_2.T_sw(0,0),
                                                 my_correction_sw,
                                                 my_score_2.sigma_sw_i,
                                                 my_score_2.sigma_sw_sqrt);

    if( my_debug_out )
      debug_out() << "my_emp_pvalue_sandwich = " << my_emp_pvalue_sw << endl;

    score_result.sigma.push_back(my_score_2.sigma_sw);
    score_result.T.push_back(my_score_2.T_sw(0,0));
    score_result.correction.push_back(my_correction_sw);
    score_result.emp_pvalue.push_back(my_emp_pvalue_sw);
    score_result.rep_count.push_back(my_rep_count_sw);
  }

  if( my_score_2.is_alternative_var() )
  {
    my_correction_al = compute_correction(my_score_2.U_star, my_score_2.sigma_al_i);

    if( my_debug_out )
      debug_out() << "my_correction_alternative = " << my_correction_al << endl;

    // Compute P-value
    //
    my_emp_pvalue_al = compute_empirical_p_value(my_rep_count_al,
                                                 my_score_2.df_al,
                                                 my_score_2.T_al(0,0),
                                                 my_correction_al,
                                                 my_score_2.sigma_al_i,
                                                 my_score_2.sigma_al_sqrt);

    if( my_debug_out )
      debug_out() << "my_emp_pvalue_alternative = " << my_emp_pvalue_al << endl;

    score_result.sigma.push_back(my_score_2.sigma_al);
    score_result.T.push_back(my_score_2.T_al(0,0));
    score_result.correction.push_back(my_correction_al);
    score_result.emp_pvalue.push_back(my_emp_pvalue_al);
    score_result.rep_count.push_back(my_rep_count_al);
  }

  if( my_score_2.is_IBD_var() )
  {
    my_correction_ib = compute_correction(my_score_2.U_star, my_score_2.sigma_ib_i);

    if( my_debug_out )
      debug_out() << "my_correction_IBD = " << my_correction_ib << endl;

    // Compute P-value
    //
    my_emp_pvalue_ib = compute_empirical_p_value(my_rep_count_ib,
                                                 my_score_2.df_ib,
                                                 my_score_2.T_ib(0,0),
                                                 my_correction_ib,
                                                 my_score_2.sigma_ib_i,
                                                 my_score_2.sigma_ib_sqrt);

    if( my_debug_out )
      debug_out() << "my_emp_pvalue_IBD = " << my_emp_pvalue_ib << endl;

    score_result.sigma.push_back(my_score_2.sigma_ib);
    score_result.T.push_back(my_score_2.T_ib(0,0));
    score_result.correction.push_back(my_correction_ib);
    score_result.emp_pvalue.push_back(my_emp_pvalue_ib);
    score_result.rep_count.push_back(my_rep_count_ib);
  }

  return true;
}

bool
two_level_score_test::do_pedigree_level()
{
  size_t t_count = my_model.get_trait_count();
  size_t p_count = my_model.get_ped_parameter_count();
  size_t g_count = ((t_count+1)*t_count)/2;

  my_score_2.reset(g_count);

  for( size_t sp = 0; sp < my_data.size(); ++sp )
  {
    const vector< mem_pointer >&    subped_members = my_data[sp].members;
    const TriangleMatrix< size_t >& subped_pairs   = my_data[sp].member_to_pair;

    // Construct second level matrices(S, D) for a sped.
    //
    size_t n        = subped_members.size();
    size_t nn       = my_residuals[sp].rows();
    size_t row_size = nn*nn;

    matrix S;
    S.resize_nofill(row_size, 1);

    matrix D;
    D.resize_nofill(row_size, g_count);

    matrix phi;
    phi.resize_nofill(n*n, 1);

    for( size_t i = 0, k = 0; i < nn; ++i )
    {
      for( size_t j = 0; j < nn; ++j, ++k )
      {
        double resid_val1 = my_residuals[sp](i, 0);
        double resid_val2 = my_residuals[sp](j, 0);

        double r = resid_val1 * resid_val2;
        S(k, 0) = r - my_null_weights[sp].first(i, j);
      }
    }

    for( size_t i = 0, k = 0; i < n; ++i )
    {
      mem_pointer mem_i = subped_members[i];

      for( size_t t1 = 0; t1 < t_count; ++t1 )
      {
        for( size_t j = 0; j < n; ++j )
        {
          mem_pointer mem_j = subped_members[j];

          size_t pair_index = subped_pairs(i,j);

          for( size_t t2 = 0; t2 < t_count; ++t2, ++k )
          {
            //cout << "k = " << k << ", (" << i << ", " << t1 << ", " << j << ", " << t2 << ")" << endl;
            for( size_t p = 0, g = 0; p < p_count; ++p )
            {
              if( my_model.get_ped_parameter(p).test_variable )
              {
                double a_par_val = get_ped_parameter_value(t1, t2, p, pair_index, mem_i, mem_j);

                D(k, g) = a_par_val;
                ++g;

                if( t1 == t2 )
                  phi(j * n  + i, 0) = a_par_val;
              }
            }
          }
        }
      }
    }

    if( my_debug_out && sp < 10 )
    {
      debug_out() << "subped " << sp << ", name = "
                  << subped_members[0]->pedigree()->name()
                  << subped_members[0]->subpedigree()->name() << endl;
      debug_out() << " number of individual in the analysis = " << n << endl;
      debug_out() << "S size : " << row_size << " by 1" << endl;
      print_matrix_first10(S, debug_out(), "S");
      print_matrix_first10(D, debug_out(), "D");
      print_trimatrix_first10(my_null_weights[sp].first,  debug_out(), "V");
      print_trimatrix_first10(my_null_weights[sp].second, debug_out(), "V^-1");
    }

    my_score_2.add_block_kron(D, my_null_weights[sp].second, S);

    if( my_score_2.is_IBD_var() )
    {
      matrix B;
      B.resize_fill(g_count, n*n, 0.0);

      get_B_matrix(n, my_null_weights[sp].second, S, B);
#if 0
      print_matrix(B, cout, "B");
      print_matrix(phi, cout, "IBD");

      matrix U;
      multiply(B, phi, U);
      print_matrix(U, cout, "B*IBD");
      print_trimatrix(my_IBD_covariances[sp], cout, "my_IBD_covariances[sp]");
#endif
      my_score_2.add_block_IBD(B, my_IBD_covariances[sp]);
    }
  }

  return my_score_2.compute();
}

double
two_level_score_test::compute_correction(const matrix& U_star, const matrix& sigma_star_i)
{
  if( my_debug_out )
  {
    debug_out() << endl << "# compute correction with" << endl;
    print_matrix(U_star, debug_out(), "U*");
    print_matrix(sigma_star_i, debug_out(), "sigma*^");
  }

  size_t t_count = my_model.get_trait_count();

  //double correction_value = -1.0;
  double correction_value = QNAN;

#if 1
  cholesky_decom correction_fn(U_star, sigma_star_i, t_count);

  if( correction_fn.is_invalid() )
    return correction_value;

  correction_fn.initialize();

  correction_value = maximize_correction_1(correction_fn);

#else
  cutting_plane correction_fn(U_star, sigma_star_i, t_count);

  if( correction_fn.is_positive_semidefinite() )
    correction_value = 0.0;
  else
  {
    correction_fn.initialize();

    size_t max_iteration = 10;

    for( size_t i = 0; i < max_iteration && !correction_fn.is_positive_semidefinite(); ++i )
    {
      correction_fn.reset_b();
      correction_fn.add_constraints();
      correction_value = maximize_correction_1(correction_fn, 0.0);
      correction_fn.compute_eigenvalue();
    }
  }
#endif

#if 0
  const matrix& w = correction_fn.get_estimates();

  cout << "w matrix after maximization:" << endl;
  print_matrix(w, cout, "w");
#endif

  return correction_value;
}

double
two_level_score_test::maximize_correction_1(semi_definite& correction_fn, double lb)
{
  size_t t_count = my_model.get_trait_count();
  size_t g_count = ((t_count+1)*t_count)/2;

  Maxfun maxfun(correction_fn);

  maxfun.nt() = g_count;

  const matrix& b = correction_fn.get_estimates();

  for( size_t i = 0, t = 0; i < t_count; ++i )
  {
    for( size_t j = i; j < t_count; ++j, ++t )
    {
      double lower_bound = NE_INF;
      if( i == j )
      {
        lower_bound = lb;
      }

      maxfun.thin(t)  = b(i,j);
      maxfun.thl(t)   = lower_bound;
      maxfun.thu(t)   = INF;
      maxfun.istin(t) = 1;
      maxfun.stpin(t) = 0.2;
    }
  }

  int    max_mthds[] = { 2,     5,     2,     6,     1,     3 };
  double max_epsc1[] = { 1e-3,  1e-3,  1e-4,  1e-12, 1e-3,  1e-3 };
  double max_epsc2[] = { 1e-12, 1e-12, 1e-12, 1e-12, 1e-15, 1e-15 };
  int    max_maxit[] = { 1,     20,    50,    20,    20,    20 };

  double correction_value = 0.0;

  for( size_t m = 0; m < 4; ++m )
  {
#if 0
  cout << "method " << m << endl;
#endif
    maxfun.method()  = max_mthds[m];
    maxfun.epsc1()   = max_epsc1[m];
    maxfun.epsc2()   = max_epsc2[m];
    maxfun.maxit()   = maxfun.nt() * max_maxit[m];  

    maxfun.run();

    correction_value = maxfun.value();

#if 0
    cout << "maximized function value = " << correction_value << endl;
    for( int i = 0; i < maxfun.nt(); ++i )
    {
      cout << "param " << i << " : " << maxfun.param(i) << endl;
    }
    cout << "nfe : " << maxfun.evaluations()
         << ", lfl : " << maxfun.last_error() << endl << endl;
#endif

    if( maxfun.evaluations() > 0 && maxfun.last_error() < 4 )
      break;

    for( int i = 0; i < maxfun.nt(); ++i )
      maxfun.thin(i) = maxfun.param(i);
  }

  return correction_value;
}

double
two_level_score_test::maximize_correction_2(semi_definite& correction_fn, double lb)
{
  size_t t_count = my_model.get_trait_count();

  MAXFUN::ParameterMgr pm;

  const matrix& b = correction_fn.get_estimates();

  for( size_t i = 0; i < t_count; ++i )
  {
    for( size_t j = i; j < t_count; ++j )
    {
      string name = "covariane";
      double lower_bound = NE_INF;
      if( i == j )
      {
        name = "variance";
        lower_bound = lb;
      }
      pm.addParameter("b", name, MAXFUN::Parameter::INDEPENDENT, b(i,j), lower_bound, INF, 0.2);
    }
  }

  //MAXFUN::DebugCfg    db;
  MAXFUN::DebugCfg    db(MAXFUN::DebugCfg::COMPLETE);
  MAXFUN::SequenceCfg sc(MAXFUN::SequenceCfg::DEFAULT_MAXIMIZATION, "A");
/*
  MAXFUN::SequenceCfg sc(MAXFUN::SequenceCfg::USER_DEFINED, "A");

  sc.addRunCfg(MAXFUN::RunCfg::DIRECT_WITHOUT, 1);
  sc.getLatestRunCfg().epsilon1 = 1e-3;
  sc.getLatestRunCfg().epsilon2 = 1e-15;
  sc.getLatestRunCfg().var_cov = MAXFUN::RunCfg::FINAL;
*/

  sc.setKeepBestRun(true);
  MAXFUN::Results maxfun_info = MAXFUN::Maximizer::Maximize(sc, pm, correction_fn, db);

  if( maxfun_info.getConverged() )
  {
    double correction_value = maxfun_info.getFinalFunctionValue();
#if 0
    cout << "final function value = " << correction_value << endl;
    for( size_t i = 0; i < ((t_count+1)*t_count)/2; ++i )
    {
      cout << "param " << i << " : "
           << maxfun_info.getParameterMgr().getParameter(i).getFinalEstimate() << endl;
    }
    cout << "nfe : " << maxfun_info.getIterations()
         << ", lfl : " << maxfun_info.getExitFlag() << endl << endl;
#endif
    return correction_value;
  }

  return 0.0;
}

double
two_level_score_test::compute_empirical_p_value(size_t& rep_count, size_t test_df,
                                                double T, double correction,
                                                const matrix& sigma_star_i,
                                                const matrix& sigma_star_sqrt)
{
  double T_org = T + correction;

  // If original T value <= 0.0, no need for sampling since p-value = 1.0.
  if( T_org <= 0.0 )
  {
    T_org = 0.0;

    rep_count = 0;

    return 1.0;
  }

  int df = my_score_2.U_star.rows();

  if( my_debug_out )
  {
    double pu        = 1.0 - chi_square_cdf(test_df, T_org);
    double threshold = my_model.get_pvalue_options().threshold;

    debug_out() << endl << "# compute empirical p-value with chi-square p-value = " << pu
                << ", threshold = " << threshold << endl;
  }

  size_t seed = my_model.get_pvalue_options().seed;

  if( seed == 0 )
    my_random_generator.reseed((long) time(NULL) + getpid());
  else
    my_random_generator.reseed((long)seed);

  double alpha     = my_model.get_pvalue_options().confidence;
  double width     = my_model.get_pvalue_options().width;
  double precision = pow(inv_normal_cdf((1.+alpha)/2.)/width, 2);

  if( my_debug_out )
  {
    debug_out() << "seed = " << seed
                << ", alpha = " << alpha
                << ", width = " << width
                << ", precision = " << precision << endl;
  }

  size_t min_replicates = my_model.get_pvalue_options().min_replicates;
  size_t max_replicates = my_model.get_pvalue_options().max_replicates;

  SampleInfo di_info;

  for( size_t i = 0; i < max_replicates;  ++i )
  {
    matrix C_vec;
    C_vec.resize_nofill(df, my_score_2.U_star.cols());

    double C_sq_sum = 0.0;

    for( int j = 0; j < df; ++j )
    {
      double C = get_random_normal();

      C_vec(j, 0) = C;
      C_sq_sum += (C * C);
    }

    matrix sim_U_star;
    multiply(sigma_star_sqrt, C_vec/sqrt(C_sq_sum), sim_U_star);

    if( my_debug_out )
    {
      debug_out() << endl << "replicate" << i << " : C_sq_sum = " << C_sq_sum << endl;
      print_matrix(C_vec, debug_out(), "C");

      //print_matrix(sigma_sqrt, debug_out(), "sigma^0.5");
      //print_matrix(sim_U_star, debug_out(), "sim_U*");

      matrix temp, sim_T;
      XTZ(sim_U_star, sigma_star_i, temp);
      multiply(temp, sim_U_star, sim_T);

      print_matrix(sim_T, debug_out(), "sim_T");
    }

    double sim_correction = compute_correction(sim_U_star, sigma_star_i);

    if( my_debug_out )
      debug_out() << "sim_correction = " << sim_correction << endl;

    if( isnan(sim_correction) )
      continue;

    //double Ti = sim_T(0, 0) + sim_correction;
    // no need to compute sim_T since always 1.0
    double Ti = 1.0 + sim_correction;

    if( Ti > 0.0 )
    {
      double T_ratio = T_org/Ti;
      double cdfi    = chi_square_cdf(test_df, T_ratio);

      if( !isnan(cdfi) )
        di_info += (1.0 - cdfi);
      else
        di_info += 0.0;

      if( my_debug_out )
      {
        debug_out() << "Ti       = " << Ti << endl
                    << "T_org/Ti = " << T_ratio << endl
                    << "cdfi     = " << cdfi << endl
                    << "sum_cdfi = " << di_info.sum() << endl;
      }
    }
    else
      di_info += 0.0;

    double pi  = di_info.mean();
    double var = di_info.variance();
    double n   = (double)di_info.count();
    double m   = n*pi*pi / var;

    if( my_debug_out )
    {
      debug_out() << "n        = " << n << endl
                  << "pi       = " << pi << endl
                  << "m        = " << m << endl;
      debug_out() << "var      = " << var << endl;
    }

    if( i > min_replicates && finite(m) && m > precision )
      break;
  }

  rep_count = di_info.count();

  return di_info.mean();
}

double
two_level_score_test::get_random_normal()
{
  static int iset = 0;
  static double gset;

  if( iset == 0 )
  {
    double u1, u2, s;

    do
    {
      u1 = 2.0 * my_random_generator.uniform_real() - 1.0;
      u2 = 2.0 * my_random_generator.uniform_real() - 1.0;

      s = u1 * u1 + u2 * u2;
    }
    while ( s > 1.0 );

    s = sqrt( (-2.0 * log(s) ) / s );

    gset = u1 * s;

    iset = 1;

    return u2 * s;
  }
  else
  {
    iset = 0;
    return gset;
  }
}   

void
two_level_score_test::correction_tests()
{
  matrix sigma_identity = eye<double>(6);

  print_matrix(sigma_identity, cout, "sigma_identity");

  matrix U_test1;
  U_test1.resize(6,1, 1.0);

  print_matrix(U_test1, cout, "U1");

  double correction1 = compute_correction(U_test1, sigma_identity);

  cout << "correction1 = " << correction1 << endl;

  matrix U_test2;
  U_test2.resize(6,1, -1.0);

  print_matrix(U_test2, cout, "U2");

  double correction2 = compute_correction(U_test2, sigma_identity);

  cout << "correction2 = " << correction2 << endl;

  matrix U_test3;
  U_test3.resize(6,1, 1.0);
  U_test3(0,0) = -1.0;
  U_test3(3,0) = 10.0;

  print_matrix(U_test3, cout, "U3");

  double correction3 = compute_correction(U_test3, sigma_identity);

  cout << "correction3 = " << correction3 << endl;

  return;
}

void
two_level_score_test::get_B_matrix(size_t n, const trimatrix& b, const matrix& S, matrix& B)
{
  trimatrix Wi;
  my_score_2.get_tri_kron(b, Wi);

  matrix C;
  my_score_2.triXZ(Wi, S, C);

#if 0
  //print_trimatrix(b, cout, "b");
  //print_trimatrix(Wi, cout, "Wi");
  //print_matrix(S, cout, "S");
  print_matrix(C, cout, "C");
#endif

  // compute B(given C);
  size_t t_count = my_model.get_trait_count();
  size_t c = 0;
  for( size_t i = 0; i < n; ++i )
  {
    for( size_t j = 0; j < n; ++j, ++c )
    {
      size_t r = 0;
      for( size_t t1 = 0; t1 < t_count; ++t1 )
      {
        for( size_t t2 = t1; t2 < t_count; ++t2, ++r )
        {
          size_t col = t_count * i + t1;
          size_t row = t_count * j + t2;

          //cout << "r = " << r << ", (" << i << ", " << j << ", " << t1 << ", " << t2 << ")  ";
          //cout << "col = " << col << ", row = " << row << ", B(" << r << ", " << c << ") += C(" << n*t_count*row + col<< ", 0) ";

          B(r, c) += C(n*t_count*row + col, 0);

          if( t1 != t2 )
          {
            col = t_count * i + t2;
            row = t_count * j + t1;

            //cout << "col = " << col << ", row = " << row << ", B(" << r << ", " << c << ") += C(" << n*t_count*row + col<< ", 0) ";
            B(r, c) += C(n*t_count*row + col, 0);
          }
          //cout << endl;
        }
      }
    }
  }

  return;
}

void
two_level_score_test::compute_covariance_matrix()
{
  if( my_null_IBD_covariances.size() )
    return;

  for( size_t sp = 0; sp < my_data.size(); ++sp )
  {
    covariance_calculator cov_cal;

    cov_cal.compute(my_data[sp].members[0]->subpedigree());

    trimatrix cov_matrix;
    cov_cal.get_covariance_matrix(my_data[sp].members, cov_matrix);
    my_null_IBD_covariances.push_back(cov_matrix);
  }

  return;
}

void
two_level_score_test::compute_covariance_matrix_given_m()
{
  if( !my_null_IBD_covariances.size() )
    compute_covariance_matrix();

  my_IBD_covariances.resize(0);

  size_t mi = (size_t) -1;

  for( size_t pi = 0; pi < my_model.get_ped_parameter_count(); ++pi )
  {
    if(    my_model.get_ped_parameter(pi).test_variable
        && my_model.get_ped_parameter(pi).type == independent_variable::MARKER )
    {
      mi = my_model.get_ped_parameter(pi).markers[0].marker_index;
    }
  }

  for( size_t sp = 0; sp < my_data.size(); ++sp )
  {
    trimatrix new_cov_matrix = my_null_IBD_covariances[sp];

    ibd_state_info istate_info;

    if( my_pairs->get_ibd_state(my_data[sp].members[0]->subpedigree(), istate_info) )
    {
      covariance_calculator cov_cal;

      trimatrix cov_matrix;
      cov_cal.get_covariance_matrix(my_data[sp].members, istate_info, mi, cov_matrix);

      for( size_t r = 0; r < cov_matrix.size(); ++ r )
      {
        for( size_t c = 0; c < cov_matrix.size(); ++ c )
        {
          new_cov_matrix(r, c) = my_null_IBD_covariances[sp](r, c) - cov_matrix(r, c);
        }
      }
    }

    my_IBD_covariances.push_back(new_cov_matrix);
  }

  return;
}

} // end of namespace RELPAL
} // end of namespace SAGE
