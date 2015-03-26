#include "relpal/output.h"

namespace SAGE   {
namespace RELPAL {

relpal_analysis::relpal_analysis(cerrorstream& err)
               : errors(err)
{}

bool
relpal_analysis::run_analysis(relative_pairs& r_pairs,
                              relpal_parser&  r_parser,
                              bool            multipoint_ibd,
                              ostream&        sum,
                              ostream&        det,
                              ostream&        exp)
{
  if( !sum ) return false;

  my_pairs    = &r_pairs;
  my_reg_type = r_parser.get_regression_type();

  build_models(r_parser);

#if 0
  cout << "number of models = " << my_analysis_models.size() << endl;
  cout << "number of null models = " << my_null_models.size() << endl;
  cout << "filter model: " << endl;
  my_filter_model.dump_model(r_pairs, cout);
#endif

  // Output
  //
  typedef boost::shared_ptr<relpal_outfile> ro_ptr;
  ro_ptr rel_output;

  rel_output = ro_ptr( new relpal_outfile(sum, det, exp, *this) );

  ostream* dat_out = NULL;
  ostream* deb_out = NULL;
  ostream* res_out = NULL;

  ofstream res_file;

  if( get_output_options().residual_out && my_null_models.size() < 2 )
  {
    string res_filename = "relpal.resid";
    res_file.open( res_filename.c_str() );
    res_out = &res_file;
  }

  // Start analysis
  //
  two_level_score_test current_score(errors);

  for( size_t i = 0; i < my_analysis_models.size(); ++i )
  {
    my_current_result     = analysis_result();
    my_current_model      = my_analysis_models[i];
    my_current_null_model = my_null_models[0];

    if( my_null_models.size() > 1 )
      my_current_null_model = my_null_models[i];

    ofstream dat_file;
    ofstream deb_file;

    string   dat_filename;
    string   deb_filename;

    if( get_output_options().data_out )
    {
      dat_filename = long2str(i+1, 1) + "_" + my_current_model.get_name() + ".dat";
      dat_file.open( dat_filename.c_str() );
      dat_out = &dat_file;

      current_score.set_data_out(dat_out);

      dat_file << "*** current test model ***" << endl;
      my_current_model.dump_model(r_pairs, *dat_out);
      dat_file << "*** current null model ***" << endl;
      my_current_null_model.dump_model(r_pairs, *dat_out);
    }

    if( get_output_options().debug_out )
    {
      deb_filename = long2str(i+1, 1) + "_" + my_current_model.get_name() + ".debug";
      deb_file.open( deb_filename.c_str() );
      deb_out = &deb_file;

      current_score.set_debug_out(deb_out);

      deb_file << "*** current test model ***" << endl;
      my_current_model.dump_model(r_pairs, *deb_out);
      deb_file << "*** current null model ***" << endl;
      my_current_null_model.dump_model(r_pairs, *deb_out);
    }

    if( get_output_options().residual_out && my_null_models.size() > 1 )
    {
      string res_filename = long2str(i+1, 1) + "_" + my_current_model.get_name() + ".resid";
      res_file.open( res_filename.c_str() );
      res_out = &res_file;
    }

    bool first_level_test = (my_current_model == my_current_null_model);

    bool print_result = true;

    if( !multipoint_ibd || first_level_test )
    {
      //cout << "Singlepoint IBD, update the ped structure..." << endl;

      build_by_pedigree(r_pairs, my_current_model, *dat_out);

      if( my_current_data.size() <= my_current_model.get_ped_parameter_count() )
      {
        cout << "Not enough pedigree to do analysis!" << endl;
        return false;
      }

      current_score.set_pairs(my_pairs);
      current_score.set_data(my_current_data);

      current_score.compute_null_weight(my_current_null_model,
                                        my_current_result.H0_result);

      if( current_score.has_valid_null() && get_output_options().residual_out )
        dump_residuals(current_score.get_residuals(), *res_out);
    }
    else
    {
      if( i == 0 )
      {
        //cout << "Multipoint IBD, construct the default ped structure & compute the default null weight..." << endl;

        build_by_pedigree(r_pairs, my_filter_model, *dat_out);

        if( my_current_data.size() <= my_current_model.get_ped_parameter_count() )
        {
          cout << "Not enough pedigree to do analysis!" << endl;
          return false;
        }

        current_score.set_pairs(my_pairs);
        current_score.set_data(my_current_data);

        current_score.compute_null_weight(my_current_null_model,
                                          my_current_result.H0_result);

        if( current_score.has_valid_null() && get_output_options().residual_out )
          dump_residuals(current_score.get_residuals(), *res_out);
      }
      else if( my_null_models.size() > 1 )
      {
        current_score.compute_null_weight(my_current_null_model,
                                          my_current_result.H0_result);

        if( current_score.has_valid_null() && get_output_options().residual_out )
          dump_residuals(current_score.get_residuals(), *res_out);
      }
    }

    if( first_level_test )
    {
      cout << "  Performing first level test............." << flush;

      if( current_score.has_valid_null() )
        cout << "done." << endl << flush;
      else
      {
        cout << "failed." << endl << flush;
        print_result = false;
      }
    }
    else
    {
      cout << "  Performing second level score test......" << flush;

      if( current_score.has_valid_null() )
      {
        current_score.set_model(my_current_model);

        if( current_score.do_two_level_score_test(my_current_result.score_result) )
        {
          cout << "done." << endl << flush;

          if( my_current_model.get_trait_count() == 1 )
          {
            cout << "  Performing second level regression test." << flush;

            print_result = do_two_level_regression(deb_out);
            cout << "done." << endl << flush;
          }
        }
        else
        {
          cout << "failed." << endl << flush;
          print_result = false;
        }

      }
      else
      {
        cout << "failed." << endl << flush;
        print_result = false;
      }
    }

    // Print results.
    //
    if( i == 0 )
      rel_output->print_header(first_level_test);

    rel_output->print_results(i + 1, first_level_test, print_result,
                              current_score.has_reliable_score());

    if( i == my_analysis_models.size()-1 )
      rel_output->print_footer(first_level_test);
  }

  return true;
}

void
relpal_analysis::build_by_pedigree(relative_pairs& r_pairs, const regression_model& f_model, ostream& out)
{
  pair_filter current_filter = make_filter(f_model);

  pair_map_type current_pair_map;

  for( size_t rp = 0; rp < r_pairs.pair_count(); ++rp )
  {
    rel_pair a_pair(rp, &r_pairs);
#if 0
    cout << "pair No " << rp << " : "
         << a_pair.rels().pair.first->pedigree()->name() << " ("
         << a_pair.rels().pair.first->name()             << ","
         << a_pair.rels().pair.second->name()            << ") ";
    cout << a_pair.rels().type << " ";
    cout << a_pair.is_fsib_pair() << " " << a_pair.is_hsib_pair() << " ";
#endif
    if( get_data_options().use_pairs == FSIB && !a_pair.is_fsib_pair() )
    {
      //cout << "not a full sib pair, skip..." << endl;
      continue;
    }

    if( get_data_options().use_pairs == HSIB && !a_pair.is_hsib_pair() )
    {
      //cout << "not a half sib pair, skip..." << endl;
      continue;
    }

    if(     get_data_options().use_pairs == SIB
        && !(a_pair.is_fsib_pair() || a_pair.is_hsib_pair()) )
    {
      //cout << "not a sib pair, skip..." << endl;
      continue;
    }

    if( current_filter.valid(a_pair) )
    {
      FPED::SubpedigreeConstPointer sp = a_pair.rels().pair.first->subpedigree();

      current_pair_map[sp].members.insert(a_pair.rels().pair.first);
      current_pair_map[sp].members.insert(a_pair.rels().pair.second);

      id_pair pair1 = make_pair(a_pair.rels().pair.first, a_pair.rels().pair.second);
      id_pair pair2 = make_pair(a_pair.rels().pair.second, a_pair.rels().pair.first);

      current_pair_map[sp].member_to_pair[pair1] = rp;
      current_pair_map[sp].member_to_pair[pair2] = rp;

      //cout << "valid pair, added..." << endl;
    }
    //else
    //  cout << "invalid pair, skip..." << endl;
  }

  my_current_data.resize(0);
  my_member_count    = 0;
  my_pair_count      = 0;
  my_fsib_pair_count = 0;
  my_hsib_pair_count = 0;

  pair_map_const_iterator pi = current_pair_map.begin();
  for( ; pi != current_pair_map.end(); ++pi )
  {
    const member_set&             subped_members = pi->second.members;
    const map< id_pair, size_t >& subped_map     = pi->second.member_to_pair;

    subpedigree_data a_subped_data;

    a_subped_data.member_to_pair.resize(subped_members.size(), (size_t)-1);

    size_t valid_pair_cnt      = 0;
    size_t valid_fsib_pair_cnt = 0;
    size_t valid_hsib_pair_cnt = 0;

    member_set_const_iterator ind_i = subped_members.begin();
    for( size_t i = 0; ind_i != subped_members.end(); ++ind_i, ++i )
    {
      a_subped_data.members.push_back(*ind_i);

      member_set_const_iterator ind_j = ind_i;
      for( size_t j = i; ind_j != subped_members.end(); ++ind_j, ++j )
      {
        size_t pair_index = (size_t)-1;

        id_pair pair_ij = make_pair(*ind_i, *ind_j);

        map< id_pair, size_t >::const_iterator mp = subped_map.find(pair_ij);
        if( mp != subped_map.end() )
        {
          pair_index = mp->second;
          ++valid_pair_cnt;

          if( is_fsib(*ind_i, *ind_j) )
            ++valid_fsib_pair_cnt;
          else if( is_hsib(*ind_i, *ind_j) )
            ++valid_hsib_pair_cnt;
        }
        
        a_subped_data.member_to_pair(i,j) = pair_index;
      }
    }

    my_current_data.push_back(a_subped_data);
    my_member_count    += a_subped_data.members.size();
    my_pair_count      += valid_pair_cnt;
    my_fsib_pair_count += valid_fsib_pair_cnt;
    my_hsib_pair_count += valid_hsib_pair_cnt;
  }

  if( get_output_options().data_out )
    dump_by_pedigree(out);

  return;
}

pair_filter
relpal_analysis::make_filter(const regression_model& r_model)
{
  pair_filter current_filter;

  for( size_t t = 0; t < r_model.get_trait_count(); ++t )
  {
    current_filter.add_trait(r_model.get_trait(t).trait_index);
  }

  for( size_t p = 0; p < r_model.get_ind_parameter_count(); ++p )
  {
    const independent_variable& param = r_model.get_ind_parameter(p);

    for( size_t c = 0; c < param.covariates.size(); ++c )
      current_filter.add_covariate(param.covariates[c].covariate_index);
  }

  for( size_t p = 0; p < r_model.get_ped_parameter_count(); ++p )
  {
    const independent_variable& param = r_model.get_ped_parameter(p);

    for( size_t m = 0; m < param.markers.size(); ++m )
      current_filter.add_marker(param.markers[m].marker_index);

    for( size_t c = 0; c < param.covariates.size(); ++c )
      current_filter.add_covariate(param.covariates[c].covariate_index);
  }

  for( size_t s = 0; s < my_data_opt.subsets.size(); ++s )
  {
    current_filter.add_subset(my_data_opt.subsets[s].subset_index, (size_t)2);
  }
  
  return current_filter;
}

bool
relpal_analysis::do_two_level_regression(ostream* out)
{
  two_level_regression current_reg(errors);

  current_reg.set_pairs(my_pairs);
  current_reg.set_model(my_current_model);
  current_reg.set_data(my_current_data);

  current_reg.set_debug_out(out);

  if( !current_reg.do_two_level_regression() )
    return false;

  my_current_result.H1_result.ind_beta     = current_reg.get_gls_1().beta;
  my_current_result.H1_result.ind_variance = current_reg.get_gls_1().Variance;

  my_current_result.H1_result.ped_beta     = current_reg.get_gls_2().beta;
  my_current_result.H1_result.ped_variance = current_reg.get_gls_2().Variance;

  return true;
}

void
relpal_analysis::dump_by_pedigree(ostream &out) const
{
  out << "subped data:" << endl;
  out << "  total member count = " << my_member_count
      << "  total pair count   = " << my_pair_count
      << endl << endl;

  for( size_t s = 0; s < my_current_data.size(); ++s )
  {
    const vector< mem_pointer >&     members = my_current_data[s].members;
    const TriangleMatrix< size_t >&  pairs   = my_current_data[s].member_to_pair;

    out << "subped " << s << ", name = " << members[0]->subpedigree()->name() << endl;

    out << "pairs: " << endl;
    for( size_t r = 0; r < pairs.size(); ++r )
    {
      for( size_t c = r; c < pairs.size(); ++c )
      {
        size_t pair_index = pairs(r,c);
      
        if( pair_index != (size_t)-1 )
        {
          out << " " << pair_index
              << " " << my_pairs->rels(pair_index).pair.first->pedigree()->name()
              << "(" << my_pairs->rels(pair_index).pair.first->name()
              << "," << my_pairs->rels(pair_index).pair.second->name()
              << ")";
        }
      }
      out << endl;
    }

    out << "members: ";
    for( size_t i = 0; i < members.size(); ++i )
      out << " " << members[i]->name();
    out << endl;

    out << "member to pair: " << endl;
    for( size_t r = 0; r < pairs.size(); ++r )
    {
      for( size_t c = 0; c < pairs.size(); ++c )
      {
        if( pairs(r,c) == (size_t)-1 )
          out << " " << "*";
        else
          out << " " << pairs(r,c);
      }
      out << endl;
    }
    out << endl;
  }
  out << endl;

  return;
}

void
relpal_analysis::dump_residuals(const vector<matrix>& res, ostream &out) const
{
  size_t t_count = my_current_model.get_trait_count();

  out << "ped_id\tind_id";
  for( size_t t = 0; t < t_count; ++t )
    out << "\t" << my_current_model.get_trait(t).name(*get_pairs());
  out << endl;

  for( size_t s = 0; s < my_current_data.size(); ++s )
  {
    const vector< mem_pointer >& members = my_current_data[s].members;
    const matrix&                a_res   = res[s];

    for( size_t i = 0; i < members.size(); ++i )
    {
      out << members[0]->pedigree()->name() << "\t" << members[i]->name();

      for( size_t t = 0; t < t_count; ++t )
        out << "\t" << fp(a_res(i*t_count+t, 0), 10, 8);

      out << endl;
    }
  }

  return;
}

void
relpal_analysis::build_models(const relpal_parser& r_parser)
{
  my_analysis_models.resize(0);
  my_null_models.resize(0);

  add_traits(my_filter_model, r_parser);
  add_ind_level_effects(my_filter_model, r_parser); 
  add_ped_level_null_effects(my_filter_model, r_parser, true); 
  add_other_options(my_filter_model, r_parser); 

  if( r_parser.get_ped_test_interactions().size() )
  {
    //cout << "Building a interaction test..." << endl;

    for( size_t i = 0; i < r_parser.get_ped_test_interactions().size(); ++i )
    {
      regression_model current_model;
      regression_model null_model;

      add_traits(current_model, r_parser);
      add_traits(null_model, r_parser);

      add_a_ped_interaction(current_model, r_parser.get_ped_test_interactions()[i], true);
      add_a_ped_interaction(my_filter_model, r_parser.get_ped_test_interactions()[i], true);

      for( size_t m = 0; m < r_parser.get_ped_test_interactions()[i].markers.size(); ++m )
        add_a_marker(null_model, r_parser.get_ped_test_interactions()[i].markers[m], false);

      add_ind_level_effects(current_model, r_parser); 
      add_ind_level_effects(null_model,    r_parser); 

      add_ped_level_null_effects(current_model, r_parser, false); 
      add_ped_level_null_effects(null_model,    r_parser, false); 

      add_other_options(current_model, r_parser); 
      add_other_options(null_model,    r_parser); 

      my_analysis_models.push_back(current_model);
      my_null_models.push_back(null_model);
    }
  }
  else
  {
    if( my_reg_type == STZM || my_reg_type == MTZM )
    {
      //cout << "Building a zero_marker test..." << endl;

      if( r_parser.get_ped_test_covariates().size() )
      {
        regression_model current_model;
        regression_model null_model;

        add_traits(current_model, r_parser);
        add_traits(null_model, r_parser);

        add_a_ped_covariate(current_model, r_parser.get_ped_test_covariates()[0], true);
        add_a_ped_covariate(my_filter_model, r_parser.get_ped_test_covariates()[0], true);

        add_ind_level_effects(current_model, r_parser); 
        add_ind_level_effects(null_model,    r_parser); 

        add_ped_level_null_effects(current_model, r_parser, false); 
        add_ped_level_null_effects(null_model,    r_parser, false); 

        add_other_options(current_model, r_parser); 
        add_other_options(null_model,    r_parser); 

        my_analysis_models.push_back(current_model);
        my_null_models.push_back(null_model);
      }
      else
      {
        if( r_parser.get_ind_batch_covariates().size() )
        {
          for( size_t i = 0; i < r_parser.get_ind_batch_covariates().size(); ++i )
          {
            regression_model current_model;
            regression_model null_model;

            add_traits(current_model, r_parser);
            add_traits(null_model, r_parser);

            add_ind_level_effects(current_model, r_parser); 
            add_ind_level_effects(null_model,    r_parser); 

            add_a_ind_covariate(current_model, r_parser.get_ind_batch_covariates()[i]);
            add_a_ind_covariate(null_model, r_parser.get_ind_batch_covariates()[i]);
            add_a_ind_covariate(my_filter_model, r_parser.get_ind_batch_covariates()[i]);

            add_ped_level_null_effects(current_model, r_parser, false); 
            add_ped_level_null_effects(null_model,    r_parser, false); 

            add_other_options(current_model, r_parser); 
            add_other_options(null_model,    r_parser); 

            current_model.set_name(my_pairs->trait_name(r_parser.get_ind_batch_covariates()[i].covariate_index));

            my_analysis_models.push_back(current_model);
            my_null_models.push_back(null_model);
          }
        }
        else
        {
          regression_model current_model;
          regression_model null_model;

          add_traits(current_model, r_parser);
          add_traits(null_model, r_parser);

          add_ind_level_effects(current_model, r_parser); 
          add_ind_level_effects(null_model,    r_parser); 

          add_ped_level_null_effects(current_model, r_parser, false); 
          add_ped_level_null_effects(null_model,    r_parser, false); 

          add_other_options(current_model, r_parser); 
          add_other_options(null_model,    r_parser); 

          current_model.set_name("first_level");

          my_analysis_models.push_back(current_model);
          my_null_models.push_back(null_model);
        }
      }
    }
    else if( my_reg_type == STMM || my_reg_type == MTMM )
    {
      //cout << "Building a multiple_marker test..." << endl;

      if(    r_parser.get_test_markers().size() == 1
          && r_parser.get_null_markers().size() == my_pairs->marker_count() - 1 )
      {
        add_a_marker(my_filter_model, r_parser.get_test_markers()[0], true);

        for( size_t i = 0; i < r_parser.get_null_markers().size(); ++i )
        {
          regression_model current_model;
          regression_model null_model;

          add_traits(current_model, r_parser);
          add_traits(null_model, r_parser);

          add_a_marker(current_model, r_parser.get_test_markers()[0], true);

          add_a_marker(current_model, r_parser.get_null_markers()[i], false);
          add_a_marker(null_model, r_parser.get_null_markers()[i], false);
          add_a_marker(my_filter_model, r_parser.get_null_markers()[i], false);

          add_ind_level_effects(current_model, r_parser); 
          add_ind_level_effects(null_model,    r_parser); 

          add_ped_level_null_effects(current_model, r_parser, false); 
          add_ped_level_null_effects(null_model,    r_parser, false); 

          add_other_options(current_model, r_parser); 
          add_other_options(null_model,    r_parser); 

          my_analysis_models.push_back(current_model);
          my_null_models.push_back(null_model);
        }
      }
      else
      {
        for( size_t i = 0; i < r_parser.get_test_markers().size(); ++i )
        {
          regression_model current_model;
          regression_model null_model;

          add_traits(current_model, r_parser);
          add_traits(null_model, r_parser);

          add_a_marker(current_model, r_parser.get_test_markers()[i], true);
          add_a_marker(my_filter_model, r_parser.get_test_markers()[i], true);

          add_ind_level_effects(current_model, r_parser); 
          add_ind_level_effects(null_model,    r_parser); 

          add_ped_level_null_effects(current_model, r_parser, true); 
          add_ped_level_null_effects(null_model,    r_parser, true); 

          add_other_options(current_model, r_parser); 
          add_other_options(null_model,    r_parser); 

          my_analysis_models.push_back(current_model);
          my_null_models.push_back(null_model);
        }
      }
    }
    else
    {
      //cout << "Building a single_marker test..." << endl;

      regression_model null_model;

      add_traits(null_model, r_parser);
      add_ind_level_effects(null_model, r_parser); 
      add_ped_level_null_effects(null_model, r_parser, true); 
      add_other_options(null_model, r_parser); 
      my_null_models.push_back(null_model);

      for( size_t i = 0; i < r_parser.get_test_markers().size(); ++i )
      {
        regression_model current_model;

        add_traits(current_model, r_parser);

        add_a_marker(current_model, r_parser.get_test_markers()[i], true);
        add_a_marker(my_filter_model, r_parser.get_test_markers()[i], true);

        add_ind_level_effects(current_model, r_parser); 
        add_ped_level_null_effects(current_model, r_parser, true); 
        add_other_options(current_model, r_parser); 

        my_analysis_models.push_back(current_model);
      }
    }
  }

  set_data_options(r_parser);
  set_output_options(r_parser);

  return;
}

void
relpal_analysis::add_traits(regression_model&    a_model,
                            const relpal_parser& r_parser)
{
  //======================
  // Dependent variable
  //======================
  // trait
  for( size_t t = 0; t < r_parser.get_traits().size(); ++t )
    a_model.add_trait( r_parser.get_traits()[t] );

  return;
}

void
relpal_analysis::add_ind_level_effects(regression_model&    a_model,
                                       const relpal_parser& r_parser)
{
  //======================
  // Independent variable
  //======================
  //-------------------
  // Individual level
  //-------------------

  // intercept
  for( size_t t = 0; t < a_model.get_trait_count(); ++t )
    a_model.add_intercept(t);

  // covariates
  add_ind_covariates(a_model, r_parser);

  // interactions
  add_ind_interactions(a_model, r_parser);

  return;
}

void
relpal_analysis::add_ped_level_null_effects(regression_model&    a_model,
                                            const relpal_parser& r_parser,
                                            bool                 add_null_marker)
{
  //======================
  // Independent variable
  //======================
  //-------------------
  // Pedigree level
  //-------------------

  // null marker
  if( add_null_marker )
    for( size_t m = 0; m < r_parser.get_null_markers().size(); ++m )
      add_a_marker(a_model, r_parser.get_null_markers()[m], false);

  // null covariates
  for( size_t c = 0; c < r_parser.get_ped_null_covariates().size(); ++c )
    add_a_ped_covariate(a_model, r_parser.get_ped_null_covariates()[c], false);

  // null interactions
  for( size_t i = 0; i < r_parser.get_ped_null_interactions().size(); ++i )
    add_a_ped_interaction(a_model, r_parser.get_ped_null_interactions()[i], false);

  // polygenic
  for( size_t t1 = 0; t1 < a_model.get_trait_count(); ++t1 )
  {
    for( size_t t2 = t1; t2 < a_model.get_trait_count(); ++t2 )
    {
      a_model.add_polygenic_variance(t1, t2);
    }
  }

  // random_error
  for( size_t t1 = 0; t1 < a_model.get_trait_count(); ++t1 )
  {
    for( size_t t2 = t1; t2 < a_model.get_trait_count(); ++t2 )
    {
      a_model.add_random_err_variance(t1, t2);
    }
  }

  return;
}

void
relpal_analysis::add_other_options(regression_model&    a_model,
                                   const relpal_parser& r_parser)
{
  //======================
  // Other oprions
  //======================
  a_model.set_analysis_options(r_parser.get_analysis_options());
  a_model.set_pvalue_options(r_parser.get_pvalue_options());

  a_model.validate();

  return;
}

void
relpal_analysis::add_ind_covariates(regression_model& a_model, const relpal_parser& r_parser)
{
  for( size_t c = 0; c < r_parser.get_ind_covariates().size(); ++c )
  {
    add_a_ind_covariate(a_model, r_parser.get_ind_covariates()[c]);
  }

  return;
}

void
relpal_analysis::add_a_ind_covariate(regression_model& a_model, const covariate_type& cov)
{
  if( cov.adj_trait_index == (size_t)-1 )
  {
    for( size_t t = 0; t < a_model.get_trait_count(); ++t )
    {
      independent_variable param;

      param.type = independent_variable::COVARIATE;
      param.t1   = param.t2 = t;
      param.covariates.push_back(cov);

      a_model.add_ind_parameter(param);
    }
  }
  else
  {
    bool   adj_t_exist = false;
    size_t adj_t       = (size_t)-1;

    for( size_t t = 0; t < a_model.get_trait_count(); ++t )
    {
      if( a_model.get_trait(t).trait_index == cov.adj_trait_index )
      {
        adj_t_exist = true;
        adj_t       = t;
        break;
      }
    }
    
    if( adj_t_exist )
    {
      independent_variable param;

      param.type = independent_variable::COVARIATE;
      param.t1   = param.t2 = adj_t;
      param.covariates.push_back(cov);

      a_model.add_ind_parameter(param);
    }
    else
      errors << priority(error)
             << "Trait to be adjusted by covariate not in the model.  Skipping..." << endl;
  }

  return;
}

void
relpal_analysis::add_ind_interactions(regression_model& a_model, const relpal_parser& r_parser)
{
  for( size_t i = 0; i < r_parser.get_ind_interactions().size(); ++i )
  {
    for( size_t t = 0; t < r_parser.get_traits().size(); ++t )
    {
      independent_variable param;

      param.type = independent_variable::CC_INTER;
      param.t1   = param.t2 = t;
      param.covariates.push_back(r_parser.get_ind_interactions()[i].covariates[0]);
      param.covariates.push_back(r_parser.get_ind_interactions()[i].covariates[1]);

      a_model.add_ind_parameter(param);
    }
  }

  return;
}

void
relpal_analysis::add_a_marker(regression_model& a_model, const marker_type& mar, bool test_variable)
{
  for( size_t t1 = 0; t1 < a_model.get_trait_count(); ++t1 )
  {
    for( size_t t2 = t1; t2 < a_model.get_trait_count(); ++t2 )
    {
      independent_variable param;

      param.type = independent_variable::MARKER;
      param.t1   = t1;
      param.t2   = t2;
      param.markers.push_back(mar);
      param.test_variable = test_variable;
      param.valid         = true;

      a_model.add_ped_parameter(param);
    }
  }

  if( test_variable )
    a_model.set_name(my_pairs->marker_name(mar.marker_index));

  return;
}

void
relpal_analysis::add_a_ped_covariate(regression_model& a_model, const covariate_type& cov, bool test_variable)
{
  for( size_t t1 = 0; t1 < a_model.get_trait_count(); ++t1 )
  {
    for( size_t t2 = t1; t2 < a_model.get_trait_count(); ++t2 )
    {
      independent_variable param;

      param.type = independent_variable::COVARIATE;
      param.t1   = t1;
      param.t2   = t2;
      param.covariates.push_back(cov);
      param.test_variable = test_variable;
      param.valid         = true;

      a_model.add_ped_parameter(param);
    }
  }

  if( test_variable )
    a_model.set_name(my_pairs->trait_name(cov.covariate_index));

  return;
}

void
relpal_analysis::add_a_ped_interaction(regression_model& a_model, const interaction_type& inter, bool test_variable)
{
  for( size_t t1 = 0; t1 < a_model.get_trait_count(); ++t1 )
  {
    for( size_t t2 = t1; t2 < a_model.get_trait_count(); ++t2 )
    {
      independent_variable param;

      if( inter.type == interaction_type::MM )
        param.type = independent_variable::MM_INTER;
      else if( inter.type == interaction_type::MC )
        param.type = independent_variable::MC_INTER;
      else
        param.type = independent_variable::CC_INTER;

      param.t1   = t1;
      param.t2   = t2;

      for( size_t c = 0; c < inter.covariates.size(); ++c )
        param.covariates.push_back(inter.covariates[c]);

      for( size_t m = 0; m < inter.markers.size(); ++m )
        param.markers.push_back(inter.markers[m]);

      param.test_variable = test_variable;
      param.valid         = true;

      a_model.add_ped_parameter(param);
    }
  }

  if( test_variable )
    a_model.set_name("interaction");

  if( inter.markers.size() )
  {
    for( size_t m = 0; m < inter.markers.size(); ++m )
    {
      add_a_marker(a_model, inter.markers[m], false);
    }
  }

  return;
}


void
relpal_analysis::set_output_options(const relpal_parser& r_parser)
{
  size_t param_name_max = DEFAULT_PARAM_NAME_MAX;

  for( size_t i = 0; i < my_analysis_models.size(); ++i )
  {
    regression_model& a_model = my_analysis_models[i];

    // Get the maximum independent variables name length.
    //
    for( size_t i = 0; i < a_model.get_ind_parameter_count(); ++i )
    {
      string param_name = a_model.get_ind_parameter(i).name(*my_pairs);

      param_name_max = std::max(param_name.size(), param_name_max);
    }

    for( size_t i = 0; i < a_model.get_ped_parameter_count(); ++i )
    {
      string param_name = a_model.get_ped_parameter(i).name(*my_pairs);

      param_name_max = std::max(param_name.size(), param_name_max);
    }
  }

  my_output_opt = r_parser.get_output_options();
  my_output_opt.param_name_max = param_name_max;

  return;
}


} // end of namespace RELPAL
} // end of namespace SAGE
