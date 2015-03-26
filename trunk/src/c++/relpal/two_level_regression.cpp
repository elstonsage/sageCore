#include "relpal/two_level_test.h"

namespace SAGE   {
namespace RELPAL {

two_level_regression::two_level_regression(cerrorstream& err)
                    : two_level_base(err)
{}

bool
two_level_regression::do_two_level_regression()
{
  if( my_debug_out )
  {
    debug_out() << endl << "do_two_level_regression..." << endl;
    //my_model.dump_model(*my_pairs, debug_out());
  }

  vector< pair<trimatrix, trimatrix> > vs;
  vector<matrix> residuals;

  return get_null_weight(vs, residuals);
}

bool
two_level_regression::get_null_weight(vector<pair <trimatrix, trimatrix> >& null_ws,
                                      vector<matrix>& residuals)
{
  compute_covariate_means();

  if( my_debug_out )
  {
    debug_out() << endl << "# get_null_weight().." << endl;
    debug_out() << "my_data.size() = " << my_data.size() << endl;
    //my_model.dump_model(*my_pairs, debug_out());
  }

  // Initial variances
  size_t t_count = my_model.get_trait_count();
  size_t p_count = my_model.get_ped_parameter_count();

  vector<double> previous_variances;
  my_variances.resize(p_count, 0.0);

  for( size_t r = p_count - ((t_count+1)*t_count)/2; r < p_count ; ++r )
  {
    if( my_model.get_ped_parameter(r).t1 == my_model.get_ped_parameter(r).t2 )
      my_variances[r] = 1.0;
  }

  if( my_debug_out )
  {
    debug_out() << "initial variances:" << endl;
    print_variances(debug_out());
  }

  // Iterrations
  bool iterate = true;

  size_t max_iterations = 80;

  my_gls_1.reset(my_model.get_ind_parameter_count());
  my_gls_2.reset(my_model.get_ped_parameter_count());

  vector<double> previous_good_variances = my_variances;

  for( size_t i = 0; iterate && i < max_iterations; ++i )
  {
    if( my_debug_out )
      debug_out() << endl << "iteration No. " << i << endl;

    vector<matrix> ys, xs;
    null_ws.resize(0);

    // 1. First level
    //
    do_individual_level(ys, xs, null_ws);

    if( !my_gls_1.Variance || !my_gls_1.Variance.rows() || !my_gls_1.Variance.cols() )
    {
      errors << priority(error)
             << "Failed at the first level estimation.  Can't perform analysis!"
             << endl;

      return false;
    }

    // 2. Get residuals from 1st level to go to 2nd level
    //
    residuals.resize(0);

    if( !build_ind_residuals(ys, xs, residuals) )
    {
      errors << priority(error)
             << "Residuals from first level have only two distinct values.  Can't do more!"
             << endl;

      return false;
    }

    // 3. Second level
    //
    do_pedigree_level(null_ws, residuals);

    if( !my_gls_2.Variance || !my_gls_2.Variance.rows() || !my_gls_2.Variance.cols() )
    {
      errors << priority(error)
             << "Failed at the second level estimation.  Can't perform analysis!"
             << endl;

      return false;
    }

    for( size_t b = 0; b < my_variances.size(); ++b )
      my_variances[b] = my_gls_2.beta(b,0);

    if( i )
    {
      // 4. Check positive semidefinite
      //
      bool pos_sem = check_positive_semidefinite();

      if( pos_sem )
      {
        // 5.1 Check convergency
        //
        double max_diff = find_max_diff(previous_good_variances, my_variances);

        if( max_diff < 0.01 )
          iterate = false;
        else
        {
          for( size_t p = 0; p < my_variances.size(); ++p )
            previous_good_variances[p] = my_variances[p];
        }
      }
      else
      {
        for( size_t j = 1; j < 20 && !pos_sem; ++j )
        {
          if( my_debug_out )
            debug_out() << "Check No. " << j << endl;

          for( size_t p = 0; p < my_variances.size(); ++p )
            my_variances[p] = 0.5 * (my_variances[p] + previous_good_variances[p]);

          if( my_debug_out )
          {
            print_variances(debug_out());
          }

          pos_sem = check_positive_semidefinite();
        }

        if( !pos_sem )
        {
          errors << priority(warning)
                 << "Variance components out of bound!" << endl;

          iterate = false;
        }
        else
        {
          // 5.2 Check convergency
          //
          double max_diff = find_max_diff(previous_good_variances, my_variances);

          if( max_diff < 0.01 )
            iterate = false;
          else
          {
            for( size_t p = 0; p < my_variances.size(); ++p )
              previous_good_variances[p] = my_variances[p];
          }
        }
      }
    }
    else
    {
      for( size_t p = 0; p < my_variances.size(); ++p )
        previous_good_variances[p] = my_variances[p];
    }
  }

  if( iterate )
  {
    errors << priority(warning)
           << "May not be fully converged!"
           << endl;
    //return false;
  }

  return true;
}

void
two_level_regression::do_individual_level(vector<matrix>&    ys,
                                          vector<matrix>&    xs,
                                          vector< pair<trimatrix, trimatrix> >& vs)
{
  // first level
  //
  size_t t_count = my_model.get_trait_count();
  size_t p_count = my_model.get_ind_parameter_count();

  my_gls_1.reset(p_count);

  for( size_t sp = 0; sp < my_data.size(); ++sp )
  {
    const vector< mem_pointer >&    subped_members = my_data[sp].members;
    const TriangleMatrix< size_t >& subped_pairs   = my_data[sp].member_to_pair;

    size_t n  = subped_members.size();
    size_t nt = n * t_count;

    matrix y;
    y.resize_nofill(nt, 1);

    matrix x;
    x.resize_nofill(nt, p_count);

    get_ind_x_y(subped_members, x, y);

    trimatrix v;
    v.resize(nt);

    for( size_t i = 0; i < n; ++i )
    {
      mem_pointer mem_i = subped_members[i];

      for( size_t t1 = 0; t1 < t_count; ++t1 )
      {
        for( size_t j = i; j < n; ++j )
        {
          mem_pointer mem_j = subped_members[j];

          size_t pair_index = subped_pairs(i,j);

          size_t t2 = 0;
          if( i == j )
            t2 = t1;

          for( ; t2 < t_count; ++t2 )
          {
            double par_val = 0.0;

            if( my_debug_out && sp < 10 )
              debug_out() << i << "," << j << ":" << t1 << "," << t2;

            for( size_t p = 0; p < my_model.get_ped_parameter_count(); ++p )
            {
              double a_par_val = get_ped_parameter_value(t1, t2, p, pair_index, mem_i, mem_j);
              par_val += (a_par_val * my_variances[p]);

              if( my_debug_out && sp < 10 )
                debug_out() << "  " << a_par_val << "*" << my_variances[p];
            }

            v(i*t_count+t1, j*t_count+t2) = par_val;

            if( my_debug_out && sp < 10 )
              debug_out() << "  v(" << i*t_count+t1 << "," << j*t_count+t2 << ") = " << par_val << endl;
          }
        }
      }
    }

    trimatrix vi;
    my_gls_1.add_block(y, v, x, vi);

    ys.push_back(y);
    xs.push_back(x);
    vs.push_back(make_pair(v, vi));

    if( my_debug_out && sp < 10 )
    {
      debug_out() << "subped " << sp << ", name = "
                  << subped_members[0]->pedigree()->name() << " "
                  << subped_members[0]->subpedigree()->name() << endl;
      debug_out() << "Y size : " << nt << " by 1" << endl; 
      print_matrix_first10(y, debug_out(), "Y");
      print_matrix_first10(x, debug_out(), "X");
      print_trimatrix_first10(v, debug_out(), "V");
      print_trimatrix_first10(vi, debug_out(), "V^-1");
    }
  }

  my_gls_1.compute();

  if( my_debug_out )
  {
    debug_out() << "from ind level regression:" << endl;
    print_matrix(my_gls_1.beta, debug_out(), "beta");
    print_matrix_first10(my_gls_1.Variance, debug_out(), "variances");
  }

  return;
}

void
two_level_regression::do_pedigree_level(const vector< pair<trimatrix, trimatrix> >& vs,
                                        const vector<matrix>& resids)
{
  // second level
  //
  size_t t_count = my_model.get_trait_count();
  size_t p_count = my_model.get_ped_parameter_count();

  my_gls_2.reset(p_count);

  for( size_t sp = 0; sp < my_data.size(); ++sp )
  {
    const vector< mem_pointer >&    subped_members = my_data[sp].members;
    const TriangleMatrix< size_t >& subped_pairs   = my_data[sp].member_to_pair;

    // Construct second level matrix.
    //
    size_t n        = subped_members.size();
    size_t nt       = resids[sp].rows();
    size_t row_size = nt*nt;

    matrix r;
    r.resize_nofill(row_size, 1);

    matrix z;
    z.resize_nofill(row_size, p_count);

    for( size_t i = 0, k = 0; i < nt; ++i )
    {
      for( size_t j = 0; j < nt; ++j, ++k )
      {
        double resid_val1 = resids[sp](i, 0);
        double resid_val2 = resids[sp](j, 0);

        r(k, 0) = resid_val1 * resid_val2;
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
            for( size_t p = 0; p < p_count; ++p )
            {
              double a_par_val = get_ped_parameter_value(t1, t2, p, pair_index, mem_i, mem_j);
              z(k, p) = a_par_val;
            }
          }
        }
      }
    }

    if( my_debug_out && sp < 10 )
    {
      debug_out() << "subped " << sp << ", name = "
                  << subped_members[0]->pedigree()->name() << " "
                  << subped_members[0]->subpedigree()->name() << endl;
      debug_out() << "R size : " << row_size << " by 1" << endl; 
      print_matrix_first10(r, debug_out(), "R");
      print_matrix_first10(z, debug_out(), "Z");
      print_trimatrix_first10(vs[sp].first, debug_out(), "V");
      print_trimatrix_first10(vs[sp].second, debug_out(), "V^-1");
    }

    my_gls_2.add_block_kron(r, vs[sp].second, z);
  }

  my_gls_2.compute();

  if( my_debug_out )
  {
    debug_out() << "from ped level regression:" << endl;
    print_matrix(my_gls_2.beta, debug_out(), "beta");
    print_matrix(my_gls_2.Variance, debug_out(), "variances");

    for( size_t r = 0, i = 0; r < p_count ; ++i )
    {
      double tv_sum = 0.0;
      for( size_t t = r; t < (i+1)*((t_count+1)*t_count)/2; ++t )
        if( my_model.get_ped_parameter(t).t1 == my_model.get_ped_parameter(t).t2 )
          tv_sum += my_gls_2.beta(t,0);

      debug_out() << "trait " << i << " variance component sum = " << tv_sum << endl;
      r += ((t_count+1)*t_count)/2;
    }
  }

  return;
}

bool
two_level_regression::check_positive_semidefinite()
{
  if( my_debug_out )
    debug_out() << "check_positive_semidefinite()... ";

  bool pos_sem = true;

  size_t t_count = my_model.get_trait_count();

  for( size_t sp = 0; sp < my_data.size(); ++sp )
  {
    const vector< mem_pointer >&    subped_members = my_data[sp].members;
    const TriangleMatrix< size_t >& subped_pairs   = my_data[sp].member_to_pair;

    size_t n  = subped_members.size();
    size_t nt = n * t_count;

    matrix b;
    b.resize_nofill(nt, nt);

    for( size_t i = 0; i < n; ++i )
    {
      mem_pointer mem_i = subped_members[i];

      for( size_t t1 = 0; t1 < t_count; ++t1 )
      {
        for( size_t j = i; j < n; ++j )
        {
          mem_pointer mem_j = subped_members[j];

          size_t pair_index = subped_pairs(i,j);

          size_t t2 = 0;
          if( i == j )
            t2 = t1;

          for( ; t2 < t_count; ++t2 )
          {
            double par_val = 0.0;

            //cout << i << "," << j << ":" << t1 << "," << t2;
            for( size_t p = 0; p < my_model.get_ped_parameter_count(); ++p )
            {
              double a_par_val = get_ped_parameter_value(t1, t2, p, pair_index, mem_i, mem_j);
              par_val += (a_par_val * my_variances[p]);
              //cout << "  " << a_par_val << "*" << my_variances[p];
            }

            b(i*t_count+t1, j*t_count+t2) = b(j*t_count+t2, i*t_count+t1) = par_val;
          }
        }
      }
    }

    vector<double> w;
    int info = Eigen(b, w);

    if( w[0] < 0.0 || info != 0 )
    {
      pos_sem = false;

      if( my_debug_out && sp < 10 )
      {
        debug_out() << endl << "eigenvalue for sped " << sp << endl;
        debug_out() << "info = " << info << endl;
        debug_out() << "eigenvalue : ";
        for( size_t i = 0; i < b.rows(); ++i )
          debug_out() << w[i] << " ";
        debug_out() << endl;
        print_matrix(b, debug_out(), "eigenvector");
      }
    }
  }

  if( my_debug_out )
  {
    if( pos_sem )
      debug_out() << " good" << endl;
    else
      debug_out() << " bad" << endl;
  }

  return pos_sem;
}

double
two_level_regression::find_max_diff(const vector<double>& prev_var, const vector<double>& curr_var) const
{
  double max_diff = NE_INF;

  for( size_t p = 0; p < prev_var.size(); ++p )
  {
    double prev = prev_var[p];
    double curr = curr_var[p];

    max_diff = max(max_diff, fabs(prev-curr)/(0.5*(curr+prev)));

    //cout << "curr = " << curr << ", prev = " << prev << ", max_diff = " << max_diff << endl;
  }

  return max_diff;
}

} // end of namespace RELPAL
} // end of namespace SAGE
