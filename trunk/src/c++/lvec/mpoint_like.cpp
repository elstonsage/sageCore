#include "lvec/mpoint_like.h"

#if 0
#define DEBUG(x) x
#endif

namespace SAGE
{

mpoint_likelihood_data::mpoint_likelihood_data(SAGE::cerrorstream& e, bool verbose)
    : err(e), my_size(0), dots(NULL),
      my_build(false), my_valid(false), my_verbose(verbose)
{
  for(int i = 0; i < 3; ++i)
  {
    ref_count[i] = 0;
    compute[i]   = false;
  }

  bad_lvector.set_valid(false);
}

mpoint_likelihood_data::~mpoint_likelihood_data()
{
  if(dots) delete dots;
}

bool mpoint_likelihood_data::build(size_t max_markers, size_t max_size)
{
  if(built() && max_size <= my_size && max_markers <= my_mcount) return true;

  if(max_size < 1 || max_markers < 1) return true;

  my_size   = max_size;
  my_mcount = max_markers;

  // If temp too small
  if(temp.bit_capacity() <= my_size)
  {
    temp = lvector(my_size);

    if(!temp.is_valid()) return false;
  }

  for(int i = 0; i < 3; ++i)
  {
    if(ref_count[i])
    {
      bool b = false;
      switch(i)
      {
        case 0 : b = build_single_point(); break;
        case 1 : b = build_separate();     break;
        case 2 : b = build_combined();     break;
      }
      if(!b)
      {
        err << SAGE::priority(SAGE::error) << "Unable to allocate memory for ";

        switch(i)
        {
          case 0 : err << "generating single point likelihood vectors.";   break;
          case 1 : err << "generating left and right likelihood vectors."; break;
          case 2 : err << "generating multipoint likelihood vectors.";     break;
        }

        err << "  Please reduce the size of your pedigrees." << endl;

        exit(7);
      }
    }
  }

  return (my_build = true);
}

bool mpoint_likelihood_data::build_single_point()
{
  single_point.resize(my_mcount);
  for(size_t i = 0; i < my_mcount; ++i)
  {
    single_point[i] = temp;

#if 0
  cout << "single_point[" << i << "] lvector storage : " << endl;
  lvector::const_iterator li;
  for( li = single_point[i].l.begin(); li != single_point[i].l.end(); ++li )
    cout << *li << "	";
  cout << endl;
  cout << "	log_scale = " << single_point[i].l.log_scale() << endl;
#endif
      
    // Error allocating the memory
    if(!single_point[i].l.is_valid()) return false;
  }      

  return true;
}

bool mpoint_likelihood_data::build_separate()
{
  if(right.size()) return true;

  right.resize(my_mcount);
  left .resize(my_mcount);

  for(size_t i = 0; i < my_mcount; ++i)
  {
    right[i].l = left[i].l = temp;

#if 0
  cout << "right[" << i << "] lvector storage : " << endl;
  lvector::const_iterator li;
  for( li = right[i].l.begin(); li != right[i].l.end(); ++li )
    cout << *li << "	";
  cout << endl;
  cout << "	log_scale = " << right[i].l.log_scale() << endl;

  cout << "left[" << i << "] lvector storage : " << endl;
  for( li = left[i].l.begin(); li != left[i].l.end(); ++li )
    cout << *li << "	";
  cout << endl;
  cout << "	log_scale = " << left[i].l.log_scale() << endl;
#endif
      
    // Error allocating the memory
    if(!right[i].l.is_valid() || !left[i].l.is_valid()) return false;
  }      

  return true;
}

bool mpoint_likelihood_data::build_combined()
{
  if(multi_point.size()) return true;

  multi_point.resize(my_mcount);

  for(size_t i = 0; i < my_mcount; ++i)
  {
    multi_point[i].l = temp;
      
#if 0
  cout << "multi_point[" << i << "] lvector storage : " << endl;
  lvector::const_iterator li;
  for( li = multi_point[i].l.begin(); li != multi_point[i].l.end(); ++li )
    cout << *li << "	";
  cout << endl;
  cout << "	log_scale = " << multi_point[i].l.log_scale() << endl;
#endif

    // Error allocating the memory
    if(!multi_point[i].l.is_valid()) return false;
  }      

  return true;
}

void mpoint_likelihood_data::compute_single_point(vector<lvector_profile>& v)
{
  if(compute[0]) return;

  if(!dots) instantiate_dot_formatter();

  if(my_verbose) 
  {
    dots->set_prefix("        Generating Marker Likelihoods");
  }

  for(int i = 0; i < lvector_count(); ++i)
  {
    SAGE::MLOCUS::inheritance_model& pm = my_markers[i];
    
#if 0
  cout << i << ":" << pm.gmodel().name() << ", lvec size = " << v[i].l.size() << endl;
#endif

    clock_type start_time = clock();

    if(my_markers.model_informative(i))
    {
      v[i].l(&my_mm, pm);
#if 0
  cout << "  total = " << v[i].l.total()
       << "  scale = " << v[i].l.log_scale()
       << "  storage : ";
  for( lvector::const_iterator li = v[i].l.begin(); li != v[i].l.end(); ++li )
    if( *li ) cout << *li << "	";
  cout << endl;
#endif

      if(!v[i].l.total() && v[i].l.bit_count())
      {
        err << SAGE::priority(SAGE::warning) << "Pedigree '"
            << my_markers.get_subpedigree().pedigree()->name()
            << "' has zero likelihood at marker '"
            << my_region.locus(i).locus()->name()
            << "'.  No information will be used at this marker." << endl;

        v[i].l.flatten(my_mm.nonfounder_meiosis_count());
      }
      else
        v[i].l.normalize();
    }
    else
    {
      // If the marker is inconsistent or otherwise bad, we just turn it
      // off.
      v[i].l.flatten(my_mm.nonfounder_meiosis_count());
    }

    v[i].self_time = clock() - start_time;
    v[i].inherited_time = 0;
    v[i].operation_time = 0;
    v[i].copy_time = 0;

#if 0
  cout << "  total = " << v[i].l.total()
       << "  scale = " << v[i].l.log_scale()
       << "  storage : ";
  for( lvector::const_iterator li = v[i].l.begin(); li != v[i].l.end(); ++li )
    if( *li ) cout << *li << "	";
  cout << endl;
#endif

    if(my_verbose) dots->trigger();   
  }

  compute[0] = true;
}

void mpoint_likelihood_data::compute_separate()
{
  if(compute[1]) return;

  if(!dots) instantiate_dot_formatter();

  vector<lvector_profile>* sp;

  // If single point available
  if(compute[0] || ref_count[0])
  {
    compute_single_point(single_point);
    sp = &single_point;
  }
  else
  {
    compute_single_point(right);
    sp = &right;
  }
  
  if(my_verbose)
  {
    dots->set_prefix("        Generating Multipoint Information");
  }

  clock_type start_time = clock();

  left[0] = (*sp)[0];

#if 0
  cout << "  total = " << left[0].l.total()
       << "  scale = " << left[0].l.log_scale()
       << "  left[0] storage : ";
  for( lvector::const_iterator li = left[0].l.begin(); li != left[0].l.end(); ++li )
    if( *li ) cout << *li << "	";
  cout << endl;
#endif

  left[0].copy_time = clock() - start_time;
  left[0].inherited_time = (*sp)[0].self_time;
  left[0].operation_time = 0;

  left[0].self_time = left[0].copy_time +
                      left[0].inherited_time +
                      left[0].operation_time;

  for(long i = 1; i < lvector_count(); ++i)
  {
    if(my_verbose && (i % 2)) dots->trigger();

    double d = my_region.locus(i).locus_distance(-1);

#if 0
  cout << endl << "distance = " << d << endl;
#endif

    start_time = clock();

    // Set it to the previous vector
    left[i].l = left[i-1].l;

#if 0
  cout << "  total = " << left[i].l.total()
       << "  scale = " << left[i].l.log_scale()
       << "  **left[" << i << "] storage : ";
  for( lvector::const_iterator li = left[i].l.begin(); li != left[i].l.end(); ++li )
    if( *li ) cout << *li << "	";
  cout << endl;
#endif

    left[i].copy_time = clock() - start_time;

    left[i].inherited_time = left[i-1].self_time + (*sp)[i].self_time;

    start_time = clock();
  
    // multiply by the transition matrix
    left[i].l(&my_mm, my_region.map()->rec_frac(d) );

#if 0
  cout << "  total = " << left[i].l.total()
       << "  scale = " << left[i].l.log_scale()
       << "  **left[" << i << "] storage : ";
  for( lvector::const_iterator li = left[i].l.begin(); li != left[i].l.end(); ++li )
    if( *li ) cout << *li << "	";
  cout << endl;
#endif
    
    // multiply by the single point vector at the target
    left[i].l *= (*sp)[i].l;

#if 0
  cout << "  total = " << left[i].l.total()
       << "  scale = " << left[i].l.log_scale()
       << "  **left[" << i << "] storage : ";
  for( lvector::const_iterator li = left[i].l.begin(); li != left[i].l.end(); ++li )
    if( *li ) cout << *li << "	";
  cout << endl;
#endif

    left[i].operation_time = clock() - start_time;

    left[i].self_time = left[i].inherited_time + left[i].copy_time +
                        left[i].operation_time;
  }
    
  start_time = clock();

  right[lvector_count() - 1] = (*sp)[lvector_count() - 1];

  right[lvector_count() - 1].copy_time = clock() - start_time;
  right[lvector_count() - 1].inherited_time = (*sp)[lvector_count()-1].self_time;

  if(my_verbose && (lvector_count() % 2)) dots->trigger();
  
  for(long i = lvector_count() - 2; i >= 0; --i)
  {
    if(my_verbose && !(i % 2)) dots->trigger();

    start_time = clock();

    right[i].l = (*sp)[i].l;

    right[i].copy_time = clock() - start_time;

    right[i].inherited_time = right[i+1].self_time + (*sp)[i].self_time;

    double d = my_region.locus(i).locus_distance(1);

    start_time = clock();

    temp = right[i+1].l;

    temp(&my_mm, my_region.map()->rec_frac(d) );

    right[i].l *= temp;

#if 0
  cout << "  total = " << right[i].l.total()
       << "  scale = " << right[i].l.log_scale()
       << "  ++right[" << i << "] storage : ";
  for( lvector::const_iterator li = right[i].l.begin(); li != right[i].l.end(); ++li )
    if( *li ) cout << *li << "	";
  cout << endl;
#endif

    right[i].operation_time = clock() - start_time;

    right[i].self_time = right[i].inherited_time + right[i].copy_time +
                         right[i].operation_time;
  }

  for(long i = 0; i < lvector_count(); ++i)
  {
    left[i].l.normalize();
    right[i].l.normalize();
  }

  compute[1] = true;
}

void mpoint_likelihood_data::compute_combined()
{
  if(compute[2]) return;

  if(!dots) instantiate_dot_formatter();

  if(compute[1] || ref_count[1])
  {
    compute_separate();

    if( my_verbose )
    {
      dots->set_prefix("        Generating Multipoint Combined Info");
      dots->trigger();
    }

    clock_type start_time = clock();

    multi_point[0].l = right[0].l;

    multi_point[0].copy_time = clock() - start_time;
    multi_point[0].inherited_time = right[0].self_time;
    multi_point[0].operation_time = 0;
    multi_point[0].self_time = multi_point[0].copy_time +
                               multi_point[0].inherited_time +
                               multi_point[0].operation_time;

    for(int i = 1; i < lvector_count() - 1; ++i)
    {
      if(my_verbose) dots->trigger();

      start_time = clock();

      multi_point[i].l = left[i].l;

      multi_point[i].copy_time = clock() - start_time;
      
      multi_point[i].inherited_time = left[i].self_time +
                                      right[i+1].self_time;

      double d = my_region.locus(i).locus_distance(1);

      start_time = clock();

      temp = right[i+1].l;

      temp (&my_mm, my_region.map()->rec_frac(d) );
      
      multi_point[i].l *= temp;

      multi_point[i].operation_time = clock() - start_time;
      multi_point[i].self_time = multi_point[i].copy_time +
                                 multi_point[i].inherited_time +
                                 multi_point[i].operation_time;
    }

    size_t t = lvector_count() - 1;

    start_time = clock();

    multi_point[t].l = left[t].l;
    
    multi_point[t].copy_time = clock() - start_time;
    multi_point[t].inherited_time = left[t].self_time;
    multi_point[t].operation_time = 0;
    multi_point[t].self_time = multi_point[t].copy_time +
                               multi_point[t].inherited_time +
                               multi_point[t].operation_time;

    if(my_verbose) dots->trigger();

    compute[2] = true;

    return;
  }

  compute_single_point(single_point);

  if(my_verbose)
  {
    dots->set_prefix("        Generating Multipoint Combined Info");
  }

  clock_type start_time = clock();

#if 0
  cout << "0:" << my_markers[0].gmodel().name()
       << ", m lvec size = " << multi_point[0].l.size() << endl;
#endif

  multi_point[0] = single_point[0];

  multi_point[0].copy_time = clock() - start_time;
  multi_point[0].inherited_time = single_point[0].self_time;
  multi_point[0].operation_time = 0;
  multi_point[0].self_time = multi_point[0].copy_time +
                             multi_point[0].inherited_time +
                             multi_point[0].operation_time;

#if 0
  cout << "  total = " << multi_point[0].l.total()
       << "  log_scale = " << multi_point[0].l.log_scale()
       << "  m storage : ";
  for( lvector::const_iterator li = multi_point[0].l.begin(); li != multi_point[0].l.end(); ++li )
    if( *li ) cout << *li << "	";
  cout << endl;
#endif

  for(long i = 1; i < lvector_count(); ++i)
  {
    if(my_verbose && (i % 2)) dots->trigger();

#if 0
  cout << i << ":" << my_markers[i].gmodel().name()
       << ", m lvec size = " << multi_point[i].l.size() << endl;
#endif

    double d = my_region.locus(i).locus_distance(-1);
    
    start_time = clock();

    multi_point[i] = multi_point[i-1];

#if 0
  cout << "  total = " << multi_point[i].l.total()
       << "  log_scale = " << multi_point[i].l.log_scale()
       << "  m storage 1 : ";
  for( lvector::const_iterator li = multi_point[i].l.begin(); li != multi_point[i].l.end(); ++li )
    if( *li ) cout << *li << "	";
  cout << endl;
#endif

    multi_point[i].copy_time = clock() - start_time;

    multi_point[i].inherited_time = multi_point[i-1].self_time;

    start_time = clock();

    multi_point[i].l (&my_mm, my_region.map()->rec_frac(d) );

#if 0
  cout << "  total = " << multi_point[i].l.total()
       << "  log_scale = " << multi_point[i].l.log_scale()
       << "  m storage 2 : ";
  for( lvector::const_iterator li = multi_point[i].l.begin(); li != multi_point[i].l.end(); ++li )
    if( *li ) cout << *li << "	";
  cout << endl;
#endif
    
    multi_point[i].l *= single_point[i].l;

    multi_point[i].operation_time = clock() - start_time;

#if 0
  cout << "  total = " << multi_point[i].l.total()
       << "  log_scale = " << multi_point[i].l.log_scale()
       << "  m storage 3 : ";
  for( lvector::const_iterator li = multi_point[i].l.begin(); li != multi_point[i].l.end(); ++li )
    if( *li ) cout << *li << "	";
  cout << endl;
#endif
  }

  for(long i = 0; i < lvector_count(); ++i)
    multi_point[i].l.normalize();

  temp = single_point[lvector_count() - 1].l;

  if(my_verbose && lvector_count() % 2) dots->trigger();

  for( long i = lvector_count() - 2; i >= 0; --i )
  {
#if 0
  cout << i << ":" << my_markers[i].gmodel().name()
       << ", m lvec size = " << multi_point[i].l.size() << endl;
#endif
    if(my_verbose && !(i % 2)) dots->trigger();

    double d = my_region.locus(i).locus_distance(1);

    start_time = clock();

    temp(&my_mm, my_region.map()->rec_frac(d));
 
    temp.normalize();
#if 0
  cout << "  total = " << temp.total()
       << "  log_scale = " << temp.log_scale()
       << "  temp storage : ";
  for( lvector::const_iterator li = temp.begin(); li != temp.end(); ++li )
    if( *li ) cout << *li << "	";
  cout << endl;
#endif

    multi_point[i].l *= temp;

    if(i != 0) temp *= single_point[i].l;

#if 0
  cout << "  total = " << multi_point[i].l.total()
       << "  log_scale = " << multi_point[i].l.log_scale()
       << "  m storage 4 : ";
  for( lvector::const_iterator li = multi_point[i].l.begin(); li != multi_point[i].l.end(); ++li )
    if( *li ) cout << *li << "	";
  cout << endl;
#endif
  }

  compute[2] = true;
}

bool mpoint_likelihood_data::valid()
{
  if(my_valid) return true;

  DEBUG( cout << "VALID: " << my_region.valid() << ' ' << my_region.locus_count()
              << my_markers.marker_count() << ' ' << my_mm.build() << endl; )

  if(!my_region.valid() || my_region.locus_count() != my_markers.inheritance_model_count() 
                        || my_markers.inheritance_model_count() == 0)
    return false;

  DEBUG( cout << my_mm.founder_count() << ' ' << my_mm.nonfounder_count()
              << endl; )

  DEBUG( cout << my_mm.nonfounder_meiosis_count() << ' ' << my_size << ' '
              << my_mm.get_subpedigree()->name() << ' ' << my_markers.get_subpedigree().name() << endl; )

  if(my_mm.nonfounder_meiosis_count() > my_size)
    return false;

  if(my_region.locus_count() > my_mcount) return false;

  return my_valid = true;
}

void mpoint_likelihood_data::dump_times(ostream& o) const
{
  o << "Timing for Vectors:" << endl;

  o << "Single Point: " << endl;

  o << "Vector\tSelf\tCopy\tOp\tIn" << endl;

  for(size_t i = 0; i < single_point.size(); ++i)
  {
    cout << i << '\t' << ((double) single_point[i].self_time)/CLOCKS_PER_SEC << '\t'
                      << ((double) single_point[i].copy_time)/CLOCKS_PER_SEC << '\t'
                      << ((double) single_point[i].operation_time)/CLOCKS_PER_SEC << '\t'
                      << ((double) single_point[i].inherited_time)/CLOCKS_PER_SEC << endl;
  }

  o << "Left: " << endl;

  o << "Vector\tSelf\tCopy\tOp\tIn" << endl;

  for(size_t i = 0; i < left.size(); ++i)
  {
    cout << i << '\t' << ((double) left[i].self_time)/CLOCKS_PER_SEC << '\t'
                      << ((double) left[i].copy_time)/CLOCKS_PER_SEC << '\t'
                      << ((double) left[i].operation_time)/CLOCKS_PER_SEC << '\t'
                      << ((double) left[i].inherited_time)/CLOCKS_PER_SEC << endl;
  }
  o << "Right: " << endl;

  o << "Vector\tSelf\tCopy\tOp\tIn" << endl;

  for(size_t i = 0; i < right.size(); ++i)
  {
    cout << i << '\t' << ((double) right[i].self_time)/CLOCKS_PER_SEC << '\t'
                      << ((double) right[i].copy_time)/CLOCKS_PER_SEC << '\t'
                      << ((double) right[i].operation_time)/CLOCKS_PER_SEC << '\t'
                      << ((double) right[i].inherited_time)/CLOCKS_PER_SEC << endl;
  }
  o << "Multi: " << endl;

  o << "Vector\tSelf\tCopy\tOp\tIn" << endl;

  for(size_t i = 0; i < multi_point.size(); ++i)
  {
    cout << i << '\t' << ((double) multi_point[i].self_time)/CLOCKS_PER_SEC << '\t'
                      << ((double) multi_point[i].copy_time)/CLOCKS_PER_SEC << '\t'
                      << ((double) multi_point[i].operation_time)/CLOCKS_PER_SEC << '\t'
                      << ((double) multi_point[i].inherited_time)/CLOCKS_PER_SEC << endl;
  }
}

} // End of SAGE namespace
