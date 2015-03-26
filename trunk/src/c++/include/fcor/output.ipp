// ---------------------------------------------------------------------------
// Inline Implementation of FcorViewBase
// ---------------------------------------------------------------------------

inline size_t
FcorView::trait_count() const
{
  return my_pairsetdata->get_parser()->get_trait_count();
}

inline const analysis_option_type&
FcorView::get_analysis_options() const
{
  return my_pairsetdata->get_parser()->get_analysis_options();
}

inline string
FcorView::get_trait_name(size_t t) const
{
  const vector<name_index_type>& trait = my_pairsetdata->get_parser()->get_trait_list();

  return trait[t].first;
}

inline double
FcorView::eff_count(double cor, double se) const
{
  if( SAGE::isnan(cor) || SAGE::isnan(se) )
    return std::numeric_limits<double>::quiet_NaN();

  if( fabs(se) < std::numeric_limits<double>::epsilon() )
    return std::numeric_limits<double>::quiet_NaN();

  double N = 1.0 + ( 1.0 - cor * cor ) * ( 1.0 - cor * cor ) / ( se * se );

  return ( N + sqrt( N * N + 22.0 * cor * cor * ( 1.0 - cor * cor ) / ( se * se ) ) ) / 2.0;
}
