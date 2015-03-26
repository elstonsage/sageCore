//---------------------------------------------------------------------------
// Inline Implementation of genibd_analysis
//---------------------------------------------------------------------------

inline const genibd_parameters*
genibd_analysis::get_parameters() const
{
  return my_parameters;
}

inline const RefMultiPedigree*
genibd_analysis::get_multipedigree() const
{
  return my_multipedigree;
}

inline bool
genibd_analysis::split_pedigree(const meiosis_map& mm)
{
  if( !mm.is_x_linked() && my_parameters->allow_family_splitting() == NO )
    return false;

  if( my_parameters->allow_family_splitting() == ALWAYS )
    return true;

  if( !mm.is_x_linked() && my_parameters->allow_simulation() != NO )
    return false;

  size_t max_bits = mm.bit_count();

  if( my_parameters->is_multipoint() || !my_parameters->allow_loops() )
    if( max_bits > my_parameters->max_exact_size() )
      return true;

  return false;
}

inline bool
genibd_analysis::use_single_point(const meiosis_map& mm)
{
  if( my_parameters->is_multipoint() )
    return false;

  if( !mm.is_x_linked() && my_parameters->allow_simulation() == ALWAYS )
    return false;

  if( mm.family_count() != mm.founder_count() - 1 ) //loop exist
    return false;

  return true;
}

inline bool
genibd_analysis::use_exact(const meiosis_map& mm)
{
  if( !mm.is_x_linked() && my_parameters->allow_simulation() == ALWAYS )
    return false;

  if(   !my_parameters->is_multipoint()                // singlepoint,
      && my_parameters->pair_category() != ALL         //  not all pair type,
      && mm.family_count() == mm.founder_count() - 1 ) //  no loop, don't need to use exact.
    return false;

  if( mm.bit_count() > my_parameters->max_exact_size() )
    return false;

  return true;
}

inline bool
genibd_analysis::use_simulation(const meiosis_map& mm)
{
  if( !mm.is_x_linked() && my_parameters->allow_simulation() == ALWAYS )
    return true;

  if( my_parameters->allow_simulation() == NO )
    return false;

  if(   !my_parameters->is_multipoint()                // no loop single,
      && mm.family_count() == mm.founder_count() - 1 ) //  no need for simulation
    return false;

  size_t ped_size = mm.bit_count();

  if( ped_size > my_parameters->max_exact_size() )
    return true;

  return false;
}

inline bool
genibd_analysis::check_pedigree(const meiosis_map& mm,
                                const region_type& r,
                                analysis_data&     data)
{
  if( use_exact(mm) )
  {
    if( my_parameters->is_multipoint() )
    {
      data.allow_exact_multi = true;

      if( my_parameters->scan_interval() )
        data.intervals = true;
    }
    else
      data.allow_exact_single = true;

    size_t ped_size = mm.bit_count();

    if( data.max_ped_size < ped_size )
      data.max_ped_size = ped_size;

    if( data.max_loci < r.locus_count() )
      data.max_loci = r.locus_count();
  }
  else if( use_simulation(mm) )
    data.allow_sim = true;

  else if( use_single_point(mm) )
    data.allow_single = true;

  else
    return false;

  return true;
}
