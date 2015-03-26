//
//  Genomic Description - Useful class defining the regions of interest
//     with markers attacked.
//
//  Author: Geoff Wedig (wedig@darwin.cwru.edu)
//
//  History  0.01 gcw Initial Interface design  Jun 15, 1998
//           0.02 gcw Moved out of mlod         Aug  5, 1998
//           0.03 gcw Redefined interface and   Aug 20, 1998
//                    made more robust
//                djb Updated to use rped       Mar 23, 2001
//
//  Copyright (c) R. C. Elston

#include <iomanip>
#include "rped/genome_description.h"

namespace SAGE {
namespace RPED {


//============================================================================
// IMPLEMENTATION:  genome_description
//============================================================================

//============================================================================
//
//  Constructor
//
//============================================================================
genome_description::genome_description(const RefMPedInfo& rmp, cerrorstream& errors) : 
	my_rmp           (rmp), 
	my_map_func      (),
	distance         (2.0), 
	interval         (2), 
	absolute         (true), 
	method           (false),
	my_built         (false),
	my_frozen        (false),
	my_errors        (errors) 
{}
  
//============================================================================
//
//  Destructor
//
//============================================================================
genome_description::~genome_description()
{ }
  
//============================================================================
//
//  add_locus(...)
//
//============================================================================
bool 
genome_description::add_locus(long m, double pos)
{
  if(my_frozen) return false;

  my_built = false;

  long n = loci.size();

  if(n > 0 && (int) regions.size() - 1 == loci[n-1].region
           && pos < loci[n-1].location)
    return false;

  loci.push_back(locus_struct());

  loci[n].locus      = m;
  loci[n].location   = pos;
  loci[n].region     = regions.size() - 1;
  loci[n].linfo      = &(my_rmp.marker_info(m));

  return true;
}

//============================================================================
//
//  add_region(...)
//
//============================================================================
size_t 
genome_description::add_region(std::string name, bool x)
{
  if(my_frozen) return (size_t) -1;

  long n = regions.size();

  if(!name.size())
  {
    name = string("region ");

    name += doub2str(n+1, 0, 0);

    my_errors << priority(information) << "Region missing name.  Setting"
              << " to '" + name + "'." << std::endl;
  }

  regions.push_back(region_struct());

  regions[n].name = name;
  regions[n].x_linked = x;

  // - Does region with this name already exist?
  //
  for(int i = 0; i < n; ++i)
    if(name == region_name(i))
    {
      my_errors << priority(warning) << "Region name '" + name + "' repeats.  "
                << "Skipping second region with this name." << std::endl;

      regions.pop_back();

      return (size_t) -1;
    }

  regions[n].locus_loc = loci.size();
 
  my_built = false;

  return n;
}
  
//============================================================================
//
//  build()
//
//============================================================================
bool 
genome_description::build()
{
  if(my_built) return true;

  if(my_frozen || my_map_func == NULL)
    return false;

  for(long i = locus_count() - 1; i >= 0; --i)
    if(!loci[i].linfo)
    {
      loci[i].linfo = NULL; 
    }

  // Get the point positions
  if(locus_count() > 1)
    if(method)     interval_build();
    else
      if(absolute) abs_distance_build();
      else         segment_distance_build();

  return my_built = true;
}

//============================================================================
//
//  interval_build()
//
//============================================================================
void 
genome_description::interval_build()
{
  long p = interval * (locus_count() - region_count()) + region_count();

  points.resize(p);

  size_t pt     = 0;
  size_t region = 0;

  regions[0].point_loc = 0;

  for(long i = 0; i < locus_count() - 1; ++i)
  {
    if(loci[i + 1].region != loci[i].region)
    {
      points[pt].location = loci[i].location;
      points[pt].locus    = i;
      points[pt].region   = region;

      loci[i].point_loc = pt;

      ++pt;
      ++region;

      regions[region].point_loc = pt;

      continue;
    }

    double delta = (loci[i + 1].location - loci[i].location) / interval;

    points[pt].location = loci[i].location;
    points[pt].locus    = i;
    points[pt].region   = region;

    loci[i].point_loc   = pt;

    ++pt;

    for(int j = 1; j < interval; ++j, ++pt)
    {
      points[pt].location = loci[i].location + j * delta;
      points[pt].locus    = -1;
      points[pt].region   = region;
    }
  }

  points[pt].location = loci[locus_count() - 1].location;
  points[pt].locus    = locus_count() - 1;
  points[pt].region   = region;

  loci[locus_count() - 1].point_loc = point_count() - 1;
}

//============================================================================
//
//  abs_distance_build()
//
//============================================================================
void 
genome_description::abs_distance_build()
{
  long max_points = locus_count();

  for(long i = 0; i < locus_count() - 1; ++i)
  {
    if(loci[i].region != loci[i+1].region)
      max_points += long((loci[i].location + 0.5) / distance);
  }

  max_points += long((loci[locus_count() - 1].location + 0.5) / distance);

  points.resize(max_points);

  double point  = distance + loci[0].location;
  size_t j      = 1;
  size_t region = 0;

  points[0].location = loci[0].location;
  points[0].locus    = 0;
  points[0].region   = 0;

  loci[0].point_loc  = 0;

  regions[0].point_loc = 0;

  for(long i = 1; i < locus_count(); ++i, ++j)
  {
    if(loci[i].region != loci[i-1].region)
    {
      points[j].location  = loci[i].location;
      points[j].locus     = i;
      points[j].region    = loci[i].region;
      loci[i].point_loc   = j;

      point = points[j].location + distance;

      ++region;

      regions[region].point_loc = j;

      continue;
    }

    while(point < loci[i].location - distance * 0.01)
    {
      points[j].location = point;
      points[j].locus    = -1;
      points[j++].region = loci[i].region;

      point += distance;
    }

    points[j].location = loci[i].location;
    points[j].locus    = i;
    points[j].region   = loci[i].region;

    loci[i].point_loc         = j;

    if(point < loci[i].location + distance * 0.01) point += distance;
  }

  points.resize(j);
}

//============================================================================
//
//  segment_distance_build()
//
//============================================================================
void 
genome_description::segment_distance_build()
{
  long max_points = locus_count();

  for(long i = 0; i < locus_count() - 1; ++i)
  {
    if(loci[i].region != loci[i+1].region)
      max_points += long((loci[i].location + 0.5) / distance);
  }

  max_points += long((loci[locus_count() - 1].location + 0.5) / distance);

  points.resize(max_points);

  double point  = distance;
  size_t j      = 1;
  size_t region = 0;

  points[0].location = loci[0].location;
  points[0].locus    = 0;
  points[0].region   = 0;
  loci[0].point_loc  = 0;

  regions[0].point_loc = 0;

  for(long i = 1; i < locus_count(); ++i, ++j)
  {
    if(loci[i].region != loci[i-1].region)
    {
      point = distance;

      points[j].location = 0.0;
      points[j].locus    = i;
      points[j].region   = loci[i].region;
      loci[i].point_loc  = j;

      ++region;

      regions[region].point_loc = j;

      continue;
    }

    while(point < loci[i].location - distance * 0.01)
    {
      points[j].location = point;
      points[j].locus  = -1;

      ++j;

      point += distance;
    }

    points[j].location = loci[i].location;
    points[j].locus    = i;
    loci[i].point_loc  = j;

    point = loci[i].location + distance;
  }

  points.resize(j);
}

//============================================================================
// IMPLEMENTATION:  LSFgenome_description
//============================================================================
//
LSFgenome_description::LSFgenome_description(
	const RefMPedInfo       & rmp, 
	      cerrormultistream   errors, 
	      LSFBase           * params, 
	      bool                multipoint)
	: 
	genome_description(rmp, errors)
{
  if(!params || !params->List())
  {
    // GCW - Add default behavior to make single all marker single point?
    my_errors << priority(error) << "Unable to find marker map.  Unable to "
              << "continue." << std::endl;
    exit(1);
  }

//cout << "genome=\"simulate\", map=\"KOSAMBI\"" << endl
//     << "{" << endl;

  // - Map type.
  //
  LSFList::iterator i;
  MapTypeShPtr m;

  if(multipoint)
  {
    AttrVal v = attr_value(params, "MAP");

    if( v.has_value() )
    {
      string s = toUpper(v.String());

      if(s == "KOSAMBI")
      {
        m = MapTypeShPtr(new Kosambi());
      }
      else
      {
        if(s != "HALDANE")
          my_errors << priority(warning) << "Mapping function incorrectly "
                    << "specified.  Using Haldane Mapping Function." << std::endl;
      }
    }
    else
    {
      if(multipoint)
        my_errors << priority(information) << "Mapping function not "
                  << "specified.  Using Haldane Mapping Function." << std::endl;
    }
  }

  if(!m) m = MapTypeShPtr(new Haldane());

  set_mapping_function(m);

  // - Regions.
  //
  for(i=params->List()->begin(); i != params->List()->end(); ++i)
  {
    if(!*i || toUpper((*i)->name()) != "REGION") 
    {
      continue;
    }
    else
    {
      init_region(*i, multipoint);
    }
  }

//cout << "}" << endl;

  if(!region_count())
  {
    cerr << "No valid regions!  Unable to continue." << endl;

    exit(1);
  }
}

bool LSFgenome_description::init_region(LSFBase* region, bool multipoint)
{
  Haldane haldane;
  Kosambi kosambi;

  string s = attr_value(region,0).String();

  bool x = false;
  
  // Check for x linkage attribute
  if(region->attrs())
  {
    x = region->attrs()->has_attr("x_linked");

    if( !x )
      x = region->attrs()->has_attr("x");
  }

  // if the name itself shows DX...
  if( !x && s.size() )
  {
    if(    toUpper(s) == "X"
        || toUpper(s) == "CHROMOSOME X" || toUpper(s) == "CHROMOSOME_X"
        || toUpper(s) == "CHR X" || toUpper(s) == "CHR_X" )
      x = true;
  }

  if(add_region(s, x) == (size_t) -1) return false;

  if(!region->List())
  {
    my_errors << priority(error) << "Region '" + region_name(region_count()-1)
              << "' is empty.  Skipping." << std::endl;

    regions.pop_back();

    return false;
  }

  s = region_name(region_count() - 1);

//cout << "  region=\"" << s << "\"" << endl
//     << "  {" << endl;

  double pos = 0.0;

  bool marker   = false; // has seen a marker
  bool mgood    = false; // has seen a *valid* marker.
  bool distance = true;  // Have seen a distance after this marker
  bool skip     = false;

  string last_marker;

  LSFList::const_iterator i;
  AttrVal v;

  // - Markers.
  //
  for(i = region->List()->begin(); i != region->List()->end(); ++i)
  {
    if(!*i) continue;

    v=attr_value(*i,"MARKER", 0);
    if( v.has_value() )
    {
      // No distance between markers.
      if(!distance && multipoint)
      {
        my_errors << priority(warning) << "No distance given after marker "
                  << last_marker + " in region " + s + ".  Markers assumed to be at a "
                  << "distance of 0." << std::endl;
      }

      last_marker = v.String();

      size_t m = my_rmp.marker_find(v.String());

      if(m == (size_t)(-1))
      {
        skip = true;

        my_errors << priority(warning) <<  "Marker '" + v.String() 
                  << "' not found.  Marker will not be used for analysis." << std::endl;
      }
      else
      {
        // If this is the first good marker in the region.
        if(!mgood)
        {
          mgood = true;
          pos = 0;

          // If we've already seen invalid markers, we want to
          // print a warning about where the markers will be placed.
          if(marker)
            my_errors << priority(warning) << "Marker " + last_marker
                      << " is the first valid marker in region " + s + ".  It will be "
                      << "set to a position of 0." << std::endl;
        }
        
//cout << "    marker=\"" << last_marker << "\", position=" << pos << endl;
        add_locus(m, pos);
      }
      marker   = true;
      distance = false;

      continue;
    }

    v=attr_value(*i,"MISSING_MARKER", 0);
    if( v.has_value() && marker)
    {
      distance = false;
      continue;
    }

    if(!marker)
    {
      if(!skip)
      {
        skip = true;

        my_errors << priority(information) << "Parameters before the first "
                  << "marker in region " + s + " in the Genome Map File will be "
                  << "ignored." << std::endl;
      }

      continue;
    }

    if(toUpper((*i)->name()) == "THETA")
    {
      if(distance)
      {
         my_errors << priority(warning) << "Multiple distances appear "
                   << "after marker " + last_marker + " in region " + s + " in "
                   << "Genome Map File." << std::endl;
      }

      v=attr_value(*i, 0);

      if(v.Real() < 0.0 || v.Real() >= 0.5)
      {
        my_errors << priority(warning) << "Recombination fraction '"
                  << v.String() + "' outside of bounds.\n  Skipping." << std::endl;
      }
      else if( v.has_value() )
      {
        pos += map()->distance(v.Real());
      }

      distance = true;
    }
    else if(toUpper((*i)->name()) == "DISTANCE")
    {
      if(distance)
      {
         my_errors << priority(warning) <<"Multiple distances appear "
                   << "after marker " + last_marker + " in region " + s + " in "
                   << "Genome Map File." << std::endl;
      }

      v=attr_value(*i, 0);
      if( v.has_value() )
      {
        if(v.Real() < 0.0)
        {
          my_errors << priority(warning) << "Distance less than 0.0 ignored." << std::endl;

          continue;
        }

        bool h = has_attr(region, "HALDANE");
        bool k = has_attr(region, "KOSAMBI");

        if(!(h ^ k))
          pos += v.Real();
        else if(h && !k)
          pos += map()->distance(haldane.rec_frac(v.Real()));
        else 
          pos += map()->distance(kosambi.rec_frac(v.Real()));
      }

      distance = true;
    }
    else
    {
      my_errors << priority(warning) << "Parameter "
                << (*i)->name() + " in region " + s + " in Genome Map File is "
                << "unknown.  Parameter will be ignored." << std::endl;

      continue;
    }
  }

//cout << "  }" << endl;

  if(!marker)
  {
    my_errors << priority(error) << "Region '" + region_name(region_count()-1)
              << "' is empty.  Skipping." << std::endl;

    regions.pop_back();

    return false;
  }

  return true;
}

} // End namespace RPED
} // End namespace SAGE
