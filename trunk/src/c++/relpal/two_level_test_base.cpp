#include "relpal/two_level_test.h"

namespace SAGE   {
namespace RELPAL {

two_level_base::two_level_base(cerrorstream& err)
              : my_data_out(NULL), my_debug_out(NULL), errors(err)
{}

void
two_level_base::compute_covariate_means()
{
  for( size_t p = 0; p < my_model.get_ind_parameter_count(); ++p )
  {
    for( size_t c = 0; c < my_model.get_ind_parameter(p).covariates.size(); ++c )
    {
      my_model.get_ind_parameter(p).covariates[c].info.clear();
    }
  }

  for( size_t p = 0; p < my_model.get_ped_parameter_count(); ++p )
  {
    for( size_t c = 0; c < my_model.get_ped_parameter(p).covariates.size(); ++c )
    {
      my_model.get_ped_parameter(p).covariates[c].info.clear();
    }
  }

  for( size_t sp = 0; sp < my_data.size(); ++sp )
  {
    const vector< mem_pointer >& subped_members = my_data[sp].members;

    for( size_t i = 0; i < subped_members.size(); ++i )
    {
      mem_pointer mem = subped_members[i];

      for( size_t p = 0; p < my_model.get_ind_parameter_count(); ++p )
      {
        for( size_t c = 0; c < my_model.get_ind_parameter(p).covariates.size(); ++c )
        {
          covariate_type& cov = my_model.get_ind_parameter(p).covariates[c];

          cov.info += my_pairs->trait(*mem, cov.covariate_index);
        }
      }

      for( size_t p = 0; p < my_model.get_ped_parameter_count(); ++p )
      {
        for( size_t c = 0; c < my_model.get_ped_parameter(p).covariates.size(); ++c )
        {
          covariate_type& cov = my_model.get_ped_parameter(p).covariates[c];

          cov.info += my_pairs->trait(*mem, cov.covariate_index);
        }
      }
    }
  }

#if 0

  for( size_t p = 0; p < my_model.get_ind_parameter_count(); ++p )
  {
    for( size_t c = 0; c < my_model.get_ind_parameter(p).covariates.size(); ++c )
    {
      out << "ind covariate info:" << endl
           << " size = " << my_model.get_ind_parameter(p).covariates[c].info.count() << endl
           << " mean = " << my_model.get_ind_parameter(p).covariates[c].info.mean() << endl
           << " min  = " << my_model.get_ind_parameter(p).covariates[c].info.min() << endl
           << " max  = " << my_model.get_ind_parameter(p).covariates[c].info.max() << endl
           << " var  = " << my_model.get_ind_parameter(p).covariates[c].info.variance() << endl;
    }
  }
  out << endl;

  for( size_t p = 0; p < my_model.get_ped_parameter_count(); ++p )
  {
    for( size_t c = 0; c < my_model.get_ped_parameter(p).covariates.size(); ++c )
    {
      out << "ped covariate info:" << endl
           << " size = " << my_model.get_ped_parameter(p).covariates[c].info.count() << endl
           << " mean = " << my_model.get_ped_parameter(p).covariates[c].info.mean() << endl
           << " min  = " << my_model.get_ped_parameter(p).covariates[c].info.min() << endl
           << " max  = " << my_model.get_ped_parameter(p).covariates[c].info.max() << endl
           << " var  = " << my_model.get_ped_parameter(p).covariates[c].info.variance() << endl;
    }
  }
  out << endl;
#endif

  return;
}

void
two_level_base::get_ind_x_y(const vector< mem_pointer >& subped_members,
                            matrix& x, matrix& y) const
{
  size_t t_count = my_model.get_trait_count();
  size_t p_count = my_model.get_ind_parameter_count();

  for( size_t i = 0; i < subped_members.size(); ++i )
  {
    mem_pointer mem = subped_members[i];

    for( size_t t = 0; t < t_count; ++t )
    {
      y(i*t_count+t, 0) = get_trait_value(t, mem);

      for( size_t p = 0; p < p_count; ++p )
      {
        x(i*t_count+t, p) = get_ind_parameter_value(t, p, mem);
      }
    }
  }

  return;
}

bool
two_level_base::build_ind_residuals(const vector<matrix>& ys,
                                    const vector<matrix>& xs,
                                          vector<matrix>& resids)
{
  if( !my_gls_1.beta || !my_gls_1.Variance )

  resids.resize(0);

  for( size_t sp = 0; sp < my_data.size(); ++sp )
  {
    // Get residuals from individual level.
    //
    matrix y = ys[sp];
    matrix x = xs[sp];

    matrix residuals;
    my_gls_1.build_residuals(y, x, residuals);

    resids.push_back(residuals);
#if 0
    print_matrix(residuals, cout, "original residuals");
#endif
  }

  if( my_model.get_analysis_options().transform_residuals )
  {
    size_t t_count = my_model.get_trait_count();

    // Check the variability of residual values. At least 2 different value needed to go on.
    //
    vector< set<double> > r_sets;
    r_sets.resize(t_count);

    for( size_t sp = 0; sp < my_data.size(); ++sp )
    {
      const vector< mem_pointer >& subped_members = my_data[sp].members;

      for( size_t i = 0, nt = 0; i < subped_members.size(); ++i )
        for( size_t t = 0; t < t_count; ++t, ++nt )
          r_sets[t].insert(resids[sp](nt, 0));
    }

    for( size_t t = 0; t < t_count; ++t )
    {
      if( r_sets[t].size() < 3 )
        return false;
    }

    // Split multi-trait to uni-trait to normalize marginal distribution.
    //
    vector< SampleInfo >         r_sample_infos;
    vector< vector<resid_info> > r_loc_infos;

    r_sample_infos.resize(t_count);
    r_loc_infos.resize(t_count);

    for( size_t sp = 0; sp < my_data.size(); ++sp )
    {
      const vector< mem_pointer >& subped_members = my_data[sp].members;

      for( size_t i = 0, nt = 0; i < subped_members.size(); ++i )
      {
        for( size_t t = 0; t < t_count; ++t, ++nt )
        {
          double r_val = resids[sp](nt, 0);

          r_sample_infos[t] += r_val;
          r_loc_infos[t].push_back(resid_info(sp, i, r_val));
        }
      }
    }

    // For each trait, obtain empirical quantile of rank,
    //  then inverse_normal transform.
    //  tr = Inver_normal_CDF( (rank_of_residual - 1/3) / (n + 1/3) )
    //
    for( size_t t = 0; t < t_count; ++t )
    {
#if 0
      cout << "trait " << t+1 << " residual_set size = " << r_sets[t].size() << endl;
      cout << "trait " << t+1 << " residual info:" << endl;
      cout << "  size = " << r_sample_infos[t].count() << endl
           << "  mean = " << r_sample_infos[t].mean() << endl
           << "  var  = " << r_sample_infos[t].variance() << endl
           << "  min  = " << r_sample_infos[t].min() << endl
           << "  max  = " << r_sample_infos[t].max() << endl
           << "  skewness = " << r_sample_infos[t].skewness() << endl
           << "  kurtosis = " << r_sample_infos[t].kurtosis() << endl << endl;
#endif

      // 1. Sort.
      //
      sort(r_loc_infos[t].begin(), r_loc_infos[t].end(), resid_info_less());

      // 2. Find rank & assign rank value.
      //
      size_t n        = r_loc_infos[t].size();
      size_t pi       = 0;
      double previous = r_loc_infos[t][0].res_val;

      for( size_t i = 0; i < n; ++i )
      {
        double current = r_loc_infos[t][i].res_val;

        if( current != previous )
        {
          SampleInfo rank_info;
          for( size_t j = pi; j < i; ++j )
            rank_info += (j+1);

          for( size_t j = pi; j < i; ++j )
            r_loc_infos[t][j].rank_val = rank_info.mean();

          pi = i;
        }

        previous = current;
      }

      SampleInfo rank_info;
      for( size_t j = pi; j < n; ++j )
        rank_info += (j+1);

      for( size_t j = pi; j < n; ++j )
        r_loc_infos[t][j].rank_val = rank_info.mean();

      // 3. Inverse normal transform of rank value.
      //
      SampleInfo tr_info;

      for( size_t i = 0; i < n; ++i )
      {
        double ra  = r_loc_infos[t][i].rank_val;
        double tr  = inv_normal_cdf((ra-(1./3.))/(n+(1./3.)));

        tr_info                    += tr;
      }
#if 0
      cout << "trait " << t+1 << " tr info:" << endl;
      cout << "  size = " << tr_info.count() << endl
           << "  mean = " << tr_info.mean() << endl
           << "  var  = " << tr_info.variance() << endl
           << "  min  = " << tr_info.min() << endl
           << "  max  = " << tr_info.max() << endl
           << "  skewness = " << tr_info.skewness() << endl
           << "  kurtosis = " << tr_info.kurtosis() << endl << endl;
#endif

      // 4. Multiply by variance ratio.
      //
      rank_info.clear();
      for( size_t i = 0; i < n; ++i )
      {
        size_t sp  = r_loc_infos[t][i].sp_id;
        size_t ii  = r_loc_infos[t][i].ind_id;
        double ra  = r_loc_infos[t][i].rank_val;
        double tr  = inv_normal_cdf((ra-(1./3.))/(n+(1./3.)));
        double inr = tr * (r_sample_infos[t].standard_deviation() / tr_info.standard_deviation());

        resids[sp](ii*t_count+t, 0) = inr;
        r_loc_infos[t][i].inr_val   = inr;
        rank_info                  += inr;
      }
#if 0
      for( size_t i = 0; i < n; ++i )
      {
        cout << i << ":" << r_loc_infos[t][i].sp_id
                  << " " << r_loc_infos[t][i].ind_id
                  << " " << r_loc_infos[t][i].res_val
                  << " " << r_loc_infos[t][i].rank_val
                  << " " << r_loc_infos[t][i].inr_val
                  << endl;
      }
      cout << endl;
#endif
#if 0
      cout << "trait " << t+1 << " inverse normal rank info:" << endl;
      cout << "  size = " << rank_info.count() << endl
           << "  mean = " << rank_info.mean() << endl
           << "  var  = " << rank_info.variance() << endl
           << "  min  = " << rank_info.min() << endl
           << "  max  = " << rank_info.max() << endl
           << "  skewness = " << rank_info.skewness() << endl
           << "  kurtosis = " << rank_info.kurtosis() << endl << endl;
#endif
    }

#if 0
    for( size_t sp = 0; sp < my_data.size(); ++sp )
    {
      matrix residuals = resids[sp];

      print_matrix(residuals, cout, "new residuals");
    }
#endif
  }

  return true;
}

double
two_level_base::get_trait_value(size_t t, mem_pointer mem) const
{
  return my_pairs->trait(*mem, my_model.get_trait(t).trait_index);
}

double
two_level_base::get_ind_parameter_value(size_t t, size_t p, mem_pointer mem) const
{
  double param_val = 0.0;

  if( my_model.get_ind_parameter(p).t1 == t )
  {
    if( my_model.get_ind_parameter(p).type == independent_variable::INTERCEPT )
    {
      param_val =  1.0;
    }
    else if( my_model.get_ind_parameter(p).type == independent_variable::COVARIATE )
    {
      size_t cc = my_model.get_ind_parameter(p).covariates[0].covariate_index;
      double cov_mean = my_model.get_ind_parameter(p).covariates[0].info.mean();

      param_val =  my_pairs->trait(*mem, cc) - cov_mean;
    }
  }

  return param_val;
}

double
two_level_base::get_ped_parameter_value(size_t t1, size_t t2, size_t p, size_t pair_index,
                                        mem_pointer mem1,
                                        mem_pointer mem2) const
{
  if(    (    my_model.get_ped_parameter(p).t1 == t1
           && my_model.get_ped_parameter(p).t2 == t2 )
      || (    my_model.get_ped_parameter(p).t1 == t2
           && my_model.get_ped_parameter(p).t2 == t1 ) )
  {
    switch( my_model.get_ped_parameter(p).type )
    {
      case independent_variable::RANDOM_ERR:
        if( mem1 == mem2 ) return 1.0;
        else return 0.0;
        break;

      case independent_variable::COMMON_ENV:
        return 1.0;
        break;

      case independent_variable::POLYGENIC:
        if( mem1 == mem2 ) return 1.0;
        else return my_pairs->prior_avg_share(mem1->pedigree()->index(), mem1->index(), mem2->index());
        break;

      case independent_variable::MARKER:
        return   get_marker_sharing_value(p, pair_index, mem1, mem2)
               - my_pairs->prior_avg_share(mem1->pedigree()->index(), mem1->index(), mem2->index());
        break;

      case independent_variable::COVARIATE:
        return get_covariate_value(p, pair_index, mem1, mem2);
        break;

      case independent_variable::MM_INTER:
      case independent_variable::MC_INTER:
      case independent_variable::CC_INTER:
        return get_interaction_value(p, pair_index, mem1, mem2);
        break;

      default:
        return 0.0;
        break;
    }
  }

  return 0.0;
}

double
two_level_base::get_marker_sharing_value(size_t p, size_t pair_index,
                                         mem_pointer mem1,
                                         mem_pointer mem2) const
{
  double marker_sharing = 1.0;

  //cout << "mem1 = " << mem1->name() << ", mem2 = " << mem2->name();

  if( mem1 != mem2 )
  {
    if( pair_index != (size_t)-1 )
    {
      const vector<marker_type>& markers = my_model.get_ped_parameter(p).markers;
      for( size_t m = 0; m < markers.size(); ++m )
      {
        if( markers[m].effect == marker_type::DOMINANCE )
          marker_sharing *= my_pairs->prob_share(pair_index, markers[m].marker_index, 2);
        else
          marker_sharing *= my_pairs->avg_share(pair_index, markers[m].marker_index);
      }
    }
    else
    {
      const vector<marker_type>& markers = my_model.get_ped_parameter(p).markers;
      for( size_t m = 0; m < markers.size(); ++m )
      {
        if( markers[m].effect == marker_type::DOMINANCE )
          marker_sharing *= my_pairs->prior_prob_share(mem1->pedigree()->index(), mem1->index(), mem2->index(), 2);
        else
          marker_sharing *= my_pairs->prior_avg_share(mem1->pedigree()->index(), mem1->index(), mem2->index());
      }
    }
  }

  //cout << " returning mar_val = " << marker_sharing
  //     << " for param " << my_model.get_ped_parameter(p).name(*my_pairs) << endl;

  return marker_sharing;
}

double
two_level_base::get_covariate_value(size_t p, size_t pair_index,
                                    mem_pointer mem1,
                                    mem_pointer mem2) const
{
  double covariate_value = 1.0;

  //cout << "mem1 = " << mem1->name() << ", mem2 = " << mem2->name();

  const vector<covariate_type>& covariates = my_model.get_ped_parameter(p).covariates;
  for( size_t c = 0; c < covariates.size(); ++c )
  {
    double cov_mean = covariates[c].info.mean();

    double t1 = my_pairs->trait(*mem1, covariates[c].covariate_index);
    double t2 = my_pairs->trait(*mem2, covariates[c].covariate_index);

    //cout << " c = " << c << "mean = " << cov_mean << ", t1 = " << t1 << ", t2 = " << t2;

    covariate_value *= ((t1 - cov_mean)*(t2 - cov_mean)); break;
  }

  //cout << " returning cov_val = " << covariate_value
  //     << " for param " << my_model.get_ped_parameter(p).name(*my_pairs) << endl;

  return covariate_value;
}

double
two_level_base::get_interaction_value(size_t p, size_t pair_index,
                                      mem_pointer mem1,
                                      mem_pointer mem2) const
{
  return   get_marker_sharing_value(p, pair_index, mem1, mem2)
         * get_covariate_value(p, pair_index, mem1, mem2);
}

void
two_level_base::get_V(const matrix& w, matrix& v) const
{
  for( size_t r1 = 0, r = 0; r1 < w.rows(); ++r1)
  {
    for( size_t r2 = r1; r2 < w.rows(); ++r2, ++r )
    {
      for( size_t c1 = 0, c = 0; c1 < w.cols(); ++c1 )
      {
        for( size_t c2 = c1; c2 < w.cols(); ++c2, ++c )
        {
          //cout << "  (" << r << "," << c
          //     << ") = (" << r1 << "," << c1 << ")(" << r2 << "," << c2 << "),";
          v(r,c) = 2.0 * w(r1,c1) * w(r2,c2);
        }
      }
      //cout << endl;
    }    
  }

  return;
}

void
two_level_base::print_variances(ostream &out) const
{
  double sum = 0.0;
  out << "variance values:" << endl;
  for( size_t p = 0; p < my_variances.size(); ++p )
  {
    out << p << " : " << my_variances[p] << endl;
    sum += my_variances[p];
  }
  out << "var sum = " << sum << endl;
}

} // end of namespace RELPAL
} // end of namespace SAGE
