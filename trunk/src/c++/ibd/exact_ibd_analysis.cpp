//==========================================================================
//  File:       exact_ibd_analysis.cpp
//
//  Author:     Qing Sun
//
//  History:    initial coding of the file                          05/02/99
//              Updated to new libraries                         yjs Jul. 03
//
//  Notes:      Exact multipoint IBD generator
//
//  Copyright (c) 1999 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "ibd/exact_ibd_analysis.h"
#include "ibd/basic_storage_ibd.h"

namespace SAGE
{

exact_ibd_analysis::exact_ibd_analysis(cerrorstream& e, ostream& output, bool verbose)
                  : my_ldata(e, verbose), my_verbose(verbose),
                    my_output(output), errors(e)
{
  my_ibds  = NULL;
  my_dots  = NULL;

  my_built = false;
  my_valid = false;
}

exact_ibd_analysis::~exact_ibd_analysis()
{
  if( my_ibds )
    delete my_ibds;

  if( my_dots )
    delete my_dots;
}

bool
exact_ibd_analysis::build(long max_markers, long max_bits,
                          bool mp, bool sp, bool intervals)
{
  if( sp )
    my_ldata.request_resource(mpoint_likelihood_data::SINGLE_POINT);

  if( mp )
  {
    my_ldata.request_resource(mpoint_likelihood_data::MP_COMBINED);

    if( intervals )
      my_ldata.request_resource(mpoint_likelihood_data::MP_SEPARATE);
  }

  if( !my_ldata.build(max_markers, max_bits) )
  {
    errors << SAGE::priority(SAGE::error)
           << "Unable to allocate memory for analysis.  "
           << "Please reduce the maximum cutoff." << endl;

    return false;
  }

  if( !my_ldata.built() )
  {
    errors << SAGE::priority(SAGE::error)
           << "Unable to perform analysis." << endl;

    return false;
  }

  temp1 = Likelihood_Vector(max_bits);
  temp2 = temp1;

  if( !temp1.is_valid() || !temp2.is_valid() )
  {
    errors << SAGE::priority(SAGE::error)
           << "Unable to allocate memory for analysis.  "
           << "Please reduce the maximum cutoff." << endl;

    return false;
  }

  return my_built = true;
}

bool
exact_ibd_analysis::set_pedigree(const meiosis_map& pm, const pedigree_region& pr)
{
  my_meiosis_map = pm;

  my_ped_region = pr;

  my_region = my_ped_region.get_region();

#if 0
  cout << endl
       << "Meiosis_Map test on Pedigree: " << my_meiosis_map.get_pedigree()->name() << endl;

  cout << my_meiosis_map.get_subpedigree()->member_count() << endl;
  cout << my_meiosis_map.founder_count() << endl;
  cout << my_meiosis_map.nonfounder_count() << endl;
  cout << my_meiosis_map.founder_mask() << endl;
  cout << my_meiosis_map.nonfounder_mask() << endl;

  for( int count = 0; count < my_meiosis_map.get_subpedigree()->member_count(); count++ )
  {
    cout << setw(3) << my_meiosis_map.member(count)->name() << ' ';

    cout << setw(3) << (signed) my_meiosis_map.mother_meiosis(count) << ' '
         << setw(3) << (signed) my_meiosis_map.father_meiosis(count) << ' ';

    cout << setw(3) << (signed) my_meiosis_map.mother_index(count) << ' '
         << setw(3) << (signed) my_meiosis_map.father_index(count) << ' ';

    cout << setw(8) << (signed) my_meiosis_map.mother_mask(count) << ' '
         << setw(8) << (signed) my_meiosis_map.father_mask(count);

    cout << endl;
  }

  for( int count = 0; count < my_meiosis_map.founder_count(); ++count )
    cout << (signed) my_meiosis_map.mask(count) << endl;

  for( size_t m = 0; m < my_ped_region.inheritance_model_count(); ++m )
    cout << my_ped_region[m].name() << ":" << my_ped_region[m].gmodel().name() << endl;

#endif

  return true;
}

bool
exact_ibd_analysis::build_ibds(bool use_intervals)
{
  if( my_ibds )
    delete my_ibds;

  my_ibds = new basic_storage_ibd(my_meiosis_map, my_region, use_intervals);

  if( use_intervals )
  {
    my_params.ibd_count = my_region.point_count();
  }
  else
  {
    my_params.ibd_count = my_region.locus_count();
  }

  return true;
}

size_t
exact_ibd_analysis::add_pair(mem_pointer i1, mem_pointer i2, pair_type pt)
{
  if( !my_ibds )
    return false;

  return my_ibds->add_pair(i1, i2, pt);
}

bool
exact_ibd_analysis::compute(const string& title, bool sp, bool intervals, bool i_state)
{
  if( my_dots )
    delete my_dots;

  my_dots = new text_dot_formatter(cout);

  if( sp )
    intervals = false;

  my_params.generate_ibd_intervals = intervals;

  if( !my_ibds->built() )
  {
    my_ibds->build();

    assert(my_ibds->built());
  }

  my_ldata.set_markers(my_ped_region);
  my_ldata.set_meiosis_map(my_meiosis_map);

  //  my_ldata.set_prefix(title);

  if( !built() || !my_ldata.valid() )
    return my_valid = false;

  size_t ibd = 0;

  size_t mc = my_region.locus_count();
  size_t pc = my_region.point_count();

#if 0
  cout << "exact_ibd_analysis::compute().." << endl;
  cout << "mc = " << mc << endl;
  cout << "pc = " << pc << endl;
#endif

  // Generate all our likelihood vectors

  if( sp )
    my_ldata.single_point_vector(0);

  else
    my_ldata.multi_point_vector(0);

  if( my_verbose )
  {
    my_dots->set_prefix_width(43);
  }

  if( my_params.generate_ibd_intervals )
  {
    my_ldata. left_sided_vector(0);
    my_ldata.right_sided_vector(0);

    if( my_verbose )
      my_dots->set_trigger_count(pc);
  }
  else
  {
    if( my_verbose )
      my_dots->set_trigger_count(mc);
  }

  if( sp )
  {
    if( my_dots )
      my_dots->set_prefix("        Generating Singlepoint IBDs");

    compute_single_point(my_ldata.single_point_vector(0), ibd, i_state);
  }
  else
  {
    if( my_dots )
      my_dots->set_prefix("        Generating Multipoint IBDs");

    compute_single_point(my_ldata.multi_point_vector(0), ibd, i_state);
  }

  // Generate our stored variables
  ++ibd;

  if( my_dots ) my_dots->trigger();

  for(size_t m = 0; m < mc - 1; ++m )
  {
    if( my_params.generate_ibd_intervals )
    {
      for(size_t pt = 1; pt < my_region.locus(m).interval_point_count(1); ++pt )
      {
        compute_two_point(my_region, ibd, pt, m, i_state);
        ++ibd;

        if( my_dots ) my_dots->trigger();
      }
    }

    if( sp )
      compute_single_point(my_ldata.single_point_vector(m+1), ibd, i_state);
    else
      compute_single_point(my_ldata.multi_point_vector(m+1), ibd, i_state);

    ++ibd;

    if( my_dots ) my_dots->trigger();
  }

#if 0
  vector<double>  f0, f1, f2;

  my_ibds->get_ibd(0, f0, f1, f2);

  for( size_t i = 0; i < f0.size(); ++i )
    cout << f0[i] << ", " << f1[i] << ", " << f2[i] << endl;
#endif

  return my_valid = true;
}

//
//-------------------------------------------------------------------
//

void
exact_ibd_analysis::compute_single_point(const lvector& mkr, long m, bool ibd_state_out)
{
#if 0
  cout << "exact_ibd_analysis::compute_single_point(" << m << ")" << endl;
#endif

#if 0
  // - Print inheritance vector bits likelihoods.
  //
  inheritance_vector  local_ivec(my_meiosis_map);
  int  max = local_ivec.num_equivalence_classes();
  for( inheritance_vector::storage_type s = 0; s < max; s++ )
  {
    local_ivec.set_equivalence_class(s);
    int  nf = local_ivec.get_nonfounders();

    cout << "\nnf = " << nf << ", Non-founder bits:" << endl;
    for(int p = 0; p < 64; p++)
    {
      cout << nf % 2;
      nf >>= 1;
    }
    cout << "\nLikelihood  " << mkr[s] << endl;
  }
#endif

  my_ibd_analysis.compute(mkr, my_meiosis_map, ibd_state_out);

  a_marker_ibd_state current_ibd_state(0);

  if( ibd_state_out )
  {
    vector<double> lv_prob;
    my_ibd_analysis.get_lvec_probability(lv_prob);

    vector<size_t> pair_is(my_ibds->pair_count(), (size_t)-1);

    for( size_t j = 0; j < lv_prob.size(); ++j )
      current_ibd_state.push_back(make_pair(lv_prob[j], pair_is));
  }

  for( size_t i = 0; i < my_ibds->pair_count(); ++i )
  {
    id_pair p = my_ibds->get_pair(i);

    // Why are these pairs, anyway?
    if(    p.first->subindex()  >= my_meiosis_map.get_pedigree()->member_count()
        || p.second->subindex() >= my_meiosis_map.get_pedigree()->member_count() )
    { 
      DEBUG(cout << "Skipping invalid pair " << p.first->name() << ", " << p.second->name() << endl;)
      continue;
    }

    double s0 = my_ibd_analysis.prob_share(p.first, p.second, 0);
    double s2 = my_ibd_analysis.prob_share(p.first, p.second, 2);

    // Added for maternal & paternal bit split - yjs
    //
    double s1m = my_ibd_analysis.sib_prob_share(p.first, p.second, 0);
    double s1p = my_ibd_analysis.sib_prob_share(p.first, p.second, 2);
    double s1m_s1p = numeric_limits<double>::quiet_NaN();

    if( finite(s1m) && finite(s1p) )
    {
      s1m_s1p = s1m-s1p;

      if( s1m_s1p < 0. && s1m_s1p > -std::numeric_limits<double>::epsilon() ) 
        s1m_s1p = 0.;
    }

    if( finite(s0) && finite(s2) )
    {
      my_ibds->set_ibd(i, m, s0, s1m_s1p, s2);
    }

    if( ibd_state_out )
    {
      vector<size_t> i_state;
      my_ibd_analysis.get_pair_ibd_state(p.first, p.second, i_state);

      for( size_t j = 0; j < current_ibd_state.size(); ++j )
        current_ibd_state[j].second[i] = i_state[j];
    }
  }

  if( ibd_state_out )
  {
    my_ibds->set_ibd_state(m, current_ibd_state);
  }
}

void
exact_ibd_analysis::compute_two_point(region_type& r, long index, long pt, long m, bool ibd_state_out)
{
#if 0
  cout << "exact_ibd_analysis::compute_two_point(" << index << "," << pt << "," << m << ")" << endl;
#endif

  temp1 = my_ldata. left_sided_vector(m    );
  temp2 = my_ldata.right_sided_vector(m + 1);

#if 0
  cout << "temp1 lvector storage : " << endl;
  lvector::const_iterator li;
  for( li = temp1.begin(); li != temp1.end(); ++li )
    cout << *li << "    ";
  cout << endl;
  cout << "     log_scale = " << temp1.log_scale() << endl;

  cout << "temp2 lvector storage : " << endl;
  for( li = temp2.begin(); li != temp2.end(); ++li )
    cout << *li << "    ";
  cout << endl;
  cout << "     log_scale = " << temp2.log_scale() << endl;
#endif

  temp1(&my_meiosis_map, my_region.locus(m).point_theta(pt));
  temp2(&my_meiosis_map, my_region.locus(m+1).point_theta(pt - my_region.locus(m).interval_point_count(1)));

#if 0
  cout << "temp1 lvector storage : " << endl;
//  lvector::const_iterator li;
  for( li = temp1.begin(); li != temp1.end(); ++li )
    cout << *li << "    ";
  cout << endl;
  cout << "     log_scale = " << temp1.log_scale() << endl;

  cout << "temp2 lvector storage : " << endl;
  for( li = temp2.begin(); li != temp2.end(); ++li )
    cout << *li << "    ";
  cout << endl;
  cout << "     log_scale = " << temp2.log_scale() << endl;
#endif

  temp1 *= temp2;

#if 0
  cout << "temp1 lvector storage : " << endl;
//  lvector::const_iterator li;
  for( li = temp1.begin(); li != temp1.end(); ++li )
    cout << *li << "    ";
  cout << endl;
  cout << "     log_scale = " << temp1.log_scale() << endl;
#endif

  my_ibd_analysis.compute(temp1, my_meiosis_map, ibd_state_out);

  a_marker_ibd_state current_ibd_state(0);

  if( ibd_state_out )
  {
    vector<double> lv_prob;
    my_ibd_analysis.get_lvec_probability(lv_prob);

    vector<size_t> pair_is(my_ibds->pair_count(), (size_t)-1);

    for( size_t j = 0; j < lv_prob.size(); ++j )
      current_ibd_state.push_back(make_pair(lv_prob[j], pair_is));
  }

  for( size_t i = 0; i < my_ibds->pair_count(); ++i )
  {
    id_pair p = my_ibds->get_pair(i);

    // Why are these pairs, anyway?
    if(    p.first->subindex()  >= my_meiosis_map.get_pedigree()->member_count()
        || p.second->subindex() >= my_meiosis_map.get_pedigree()->member_count() )
    { 
      DEBUG(cout << "Skipping invalid pair " << p.first->name() << ", " << p.second->name() << endl;)
      continue;
    }

    double s0 = my_ibd_analysis.prob_share(p.first, p.second, 0);
    double s2 = my_ibd_analysis.prob_share(p.first, p.second, 2);

    // Added for maternal & paternal bit split - yjs
    //
    double s1m = my_ibd_analysis.sib_prob_share(p.first, p.second, 0);
    double s1p = my_ibd_analysis.sib_prob_share(p.first, p.second, 2);
    double s1m_s1p = numeric_limits<double>::quiet_NaN();

    if( finite(s1m) && finite(s1p) )
    {
      s1m_s1p = s1m-s1p;

      if( s1m_s1p < 0. && s1m_s1p > -std::numeric_limits<double>::epsilon() ) 
        s1m_s1p = 0.;
    }

    if( finite(s0) && finite(s2) )
    {
      //cout << "set ibd with " << i << "," << index << "," << s0 << "," << s1m-s1p << "," << s2 << endl;
      my_ibds->set_ibd(i, index, s0, s1m_s1p, s2);
    }

    if( ibd_state_out )
    {
      vector<size_t> i_state;
      my_ibd_analysis.get_pair_ibd_state(p.first, p.second, i_state);

      for( size_t j = 0; j < current_ibd_state.size(); ++j )
        current_ibd_state[j].second[i] = i_state[j];
    }
  }

  if( ibd_state_out )
  {
    my_ibds->set_ibd_state(index, current_ibd_state);
  }
}

} // end of namespace SAGE
