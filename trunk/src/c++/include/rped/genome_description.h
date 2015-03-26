#ifndef _GENOME_DESCRIPTION_H_
#define _GENOME_DESCRIPTION_H_

#include <vector>
#include "boost/shared_ptr.hpp"
#include "LSF/LSF.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "rped/rped.h"
#include "rped/mapping.h"

namespace SAGE {
namespace RPED {

/** \brief Provide access to information read from genome description file
  *
  * \par Review of genetic theory
  *
  * To understand the genome_description class, it is necessary to review the basic genetic
  * theory of how a genome is structured.
  *
  * First, recall that a genome consists of a linear sequence of nucleotide base pairs. This
  * sequence can be broken down into constituent regions (generally synonymous with chromosomes).
  * Within the regions, any single base pair can be uniquely identified as the physical distance
  * (measured in centimorgans) from the start of the region. Thus, the first pair (or point, if you will)
  * is located at position 0.0.
  *
  * Furthermore, in addition to the unique point-location system, there are also \b markers; a marker
  * is a sequence of base pairs beginning at a known location on a genomic region. Consequently, a marker
  * can be uniquely identified by its location on the region. Otherwise put, a marker's location is given by
  * its point.
  *
  * Although the distance between two points on a region can be expressed in a strict physical form
  * (centimorgans), it is also possible to interpret that distance as a recombination probability. That is,
  * given some genetic distance x and a mapping function f(x) (where f(x)'s range is 0.0 to 1.0), we can
  * interpret f(x) as the recombination fraction between two given points on a genome.
  *
  * \par genome_description and assorted classes
  *
  * The genome_description class tracks a sequence of regions.
  *
  * Each region is described by a single region_type instance (accessible through the genome_description).
  *
  * Within each region, a base pair location is defined by the point_struct, and a single locus is defined
  * by the locus_struct.
  *
  *
  *
  *
  *
  *
  *
  * Provide access to information read from genome description file.
  */
class genome_description : public LSFBase
{
public:

  /// Object describing a locus (MLOCUS::inheritance_model -> RPED::RefMarkerInfo)
  typedef RefMarkerInfo locus_info;

  /// Mapping function object
  typedef Mapping_Function map_type;

  /// Shared ptr to the mapping function
  typedef boost::shared_ptr<Mapping_Function> MapTypeShPtr;

  ///
  /// Double limit
  typedef numeric_limits<double> dlimit;

  class region_type;
  class locus_type;

  /** \brief Provide access to locus information via the locus_type class
    *
    * Provide access to locus information via the locus_type class.
    */
  class locus_iterator
  {
  public:
    friend class genome_description;
    friend class region_type;
    
    /** \name Typedefs
      *
      * Typedef-ed local to this class
      */
    //@{

      typedef genome_description   genome;
      typedef locus_iterator       iterator;
      typedef locus_type           reference;
      typedef size_t               difference_type;
      typedef size_t               size_type;

    //@}

    /** \name Standard iterator interface
      *
      * Standard iterator interface
      */
    //@{

      locus_iterator();
      locus_iterator(const iterator& b);

      reference       operator*() const;

      iterator& operator++();
      iterator  operator++(int);
      iterator& operator--();
      iterator  operator--(int);

      bool operator==(const iterator& x) const;
      bool operator!=(const iterator& x) const;

    //@}
    
  protected:
    // To be defined
    locus_iterator(long, genome_description&);

    long m;
    
    genome_description* gd;
  };

  ///
  /// Note that currently these iterators for the genome and the region are
  /// the same, but this need not be the case if the data representation
  /// should change.  Also note there is not currently iterators for
  /// regions.
  typedef locus_iterator       iterator;

  /** \brief A region of markers where each interval < 0.5 r.f.
    *
    * A region of markers where each interval < 0.5 r.f.  For example, a region may be a chromosome in a 
    * genome scan.
    */
  class region_type
  {
  public:

    /// Locus iterator is the default iterator
    typedef locus_iterator iterator;
  
    friend class genome_description;
    friend class locus_type;
    
    /// @name Constructors/operators
    //@{

      ///
      /// Constructor.
      region_type();
    
      ///
      /// Assignment operator.
      /// \param other The region_type to copy
      region_type& operator=(const region_type & other);

      ///
      /// Equality operator.
      /// \param other The region_type to compare
      bool operator==(const region_type & other) const;

      ///
      /// Inequality operator.
      /// \param other The region_type to compare
      bool operator!=(const region_type & other) const;

    //@}

    /// @name Basic information
    //@{

      ///
      /// Returns \c true if this region is valid, \c false if it is not.
      bool valid() const;

      ///
      /// Returns \c true if this region is x-linked, \c false if it is not.
      bool is_x_linked() const;
    
      ///
      /// Returns the name of this region.
      ///
      /// Please note that the name of a region is generally assigned when it is
      /// added to a genome_description (via genome_description::add_region()).
      /// In that case, add_region() takes as one of its arguments the intended
      /// name for that region.
      ///
      /// If the name argument is not specified, add_region will automatically generate
      /// a name, as well as a warning message.
      const std::string& name() const;

      ///
      /// Sets the name of this region.
      /// \param name The name of this region
      /// \returns The name of this region (newly updated)
      const std::string& name(const std::string & name);

      ///
      /// Returns the index number of this region (since the genome description stores a vector
      /// of region_type's, that vector can be indexed).
      size_t index() const;
    
      ///
      /// Returns the length of this region (in number of points)
      double length() const { return locus_location(locus_count()-1) - locus_location(0); }

    //@}


    /// @name Locus functions
    //@{

      ///
      /// Returns the indicated locus.
      /// \param l The index number of the requested locus
      locus_type locus(size_t l) const;
      
      ///
      /// Returns the number of locii in this region
      size_t locus_count() const;
      
      ///
      /// Returns the exact genomic location of the indicated locus (in number of points)
      /// \param l The index number of the requested locus
      double locus_location(size_t l) const;

      ///
      /// Returns a non-const begin iterator for the list of locii in this region.
      iterator   begin() const;
      
      ///
      /// Returns a non-const end iterator for the list of locii in this region.
      iterator   end() const;

    //@}

    /// @name Point functions
    //@{

      ///
      /// Returns the number of distinct points in this region.
      size_t point_count() const;
      
      ///
      /// Returns the exact location of the indicated point.
      /// \param p The index number of the requested point
      double point_location(size_t p) const;

    //@}

    /// @name Mapping function
    //@{

      ///
      /// Returns a pointer to the map function associated with this region.
      MapTypeShPtr map() const;

    //@}


  protected:

    region_type(long m, genome_description& v);
    
    long m;

    genome_description* gd;
  };

  /** \brief Provide access to locus information
    *
    * Provide access to locus information
    */
  class locus_type
  {
  public:

    friend class genome_description;
    friend class region_type;
    friend class locus_iterator;
    
    bool operator==(const locus_type&);

    const std::string& name() const;

    const locus_info* locus()      const;

    size_t      marker_index() const { return mpp.loci[m].locus; }

  // Position of the locus.

    double location() const;
    
    size_t       region_index() const;
    region_type  region()       const;

    // returns the location of the point/locus relative to this locus or beginning of region.  If
    // n == 0, returns the current locus location.  Location is the location
    // on the region, distance is the distance from this locus
    // If the point/locus is on a different region, returns infinity
    // for distances and 0.5 for recombinations
    double  point_location(long p) const;  // centimorgans
    double  point_distance(long p) const;  // centimorgans
    double  point_theta(long p)    const;  // recombination fraction

    double  locus_location(long m) const;  // centimorgans
    double  locus_distance(long m) const;  // centimorgans
    double  locus_theta(long m)    const;  // recombination fraction

    // returns the number of points between this and the m-th locus, not
    // including this locus.  If the locus is on a different region,
    // returns infinity.
    size_t interval_point_count(long m) const;

    // returns the number of items previous to or after the current
    // locus on this region
    size_t  point_prev_count() const;
    size_t  point_next_count() const;

    size_t  locus_prev_count() const;
    size_t  locus_next_count() const;
    
    // return the next/previous locus to this one in the genome.
    locus_type  next_locus() const;
    locus_type  prev_locus() const;

    // index of this point/locus on the region.
    size_t  point_index()  const;
    size_t  locus_index() const;

  protected:

    locus_type(long m, genome_description& v);
    
    // Data members.
    long                   m;
    genome_description&    mpp;
  };

  // - GENOME_DESCRIPTION CLASS CONTINUED ...
  //

  friend class region_type;
  friend class locus_type;

  /// @name Constructor/destructor
  //@{

    ///
    /// Constructor.
    /// \param rmp ???
    /// \param errors ???
    genome_description(const RefMPedInfo& rmp, cerrorstream& errors = sage_cerr);
  
    ///
    /// Destructor.
    ~genome_description();
  
  //@}

  /// @name Input specification
  //@{

    ///
    /// The input of loci into the genome_description is kludged and should be
    /// fixed when the port is done! - GCW
    /// \param n The id number of the locus
    /// \param position The absolute genomic location of the locus
    /// \retval true Locus was added successfully
    /// \retval true Locus was \b not added successfully
    bool add_locus(long n, double position);
  
    ///
    /// Adds a region to this genome.
    /// \param name The name of the new region
    /// \param x_linked Boolean value indicating whether this region is x-linked
    /// \returns The index number of the newly added region
    size_t add_region(std::string name, bool x_linked = false);
  
    ///
    /// Order and build the data structures.  
    /// Note that build currently blows away any
    /// information about regions (name, etc) that might have existed before.
    /// This is due to the restructuring necessary for the build.  If
    /// possible, this should be eventually fixed.
    /// \retval true Data structures were built successfully
    /// \retval false Data structures were \b not built successfully
    bool build();

    ///
    /// Precondition: built() = true;  if so, returns true, otherwise returns
    /// false
    bool freeze();

    ///
    /// returns boolean from previous build(), or false if markers have been
    /// added since previous build.
    bool built()  const;

    ///
    /// returns true if freeze() has returned successfully.
    bool frozen() const;

  //@}

  /// @name Positioning
  //@{

    ///
    /// Set the conversion function from centimorgans to recombination fraction
    /// and vice-versa/retrieve mapping function
    bool set_mapping_function(MapTypeShPtr);

    ///
    /// Returns a pointer to the mapping function associated with this genome.
    MapTypeShPtr map() const;

    ///
    /// Set the number of points between markers
    bool set_interval_count(long);

    ///
    /// Set distance between points.  If absolute distancing is specified,
    /// distances are in increments from the first locus on the region.
    /// If absolute is not specified, distances are in increments from the
    /// previous locus.
    bool set_scan_distance(double, bool _absolute = true);

    ///
    /// NOT YET IMPLEMENTED!
    /// Set the number of centimorgans to either side of each region in
    /// which we're interested.  If interval counting is specified, there will
    /// be that many points in this interval.  If scan distance is specified,
    /// points will be allocated based on increments in this area.  Note that
    /// points previous to the first marker will have negative locations
    /// GCW (Should that be changed?)
    bool set_edge_distance(double);

  //@}

  /// @name Output specification
  //@{

    ///
    /// Returns the number of points/loci/regions.
    /// point_count() and region_count() return -1 if built() == false;
    long point_count()   const;

    ///
    /// Returns the total number of locii.
    long locus_count() const;

    ///
    /// Returns the total number of regions.
    long region_count() const;

    ///
    /// Returns the requested locus.
    /// \param i The index number of the requested locus
    locus_type locus(size_t i);

    ///
    /// Returns the requested region.
    /// \param i The index number of the requested region
    region_type region(size_t i);

    ///
    /// Returns the requested region.
    /// \param name The name of the requested region
    region_type  region(const std::string & name);

    ///
    /// Returns the name of the requested region.
    /// \param i The id number of the requested region
    const std::string & region_name(size_t i) const;

    ///
    /// Assigns the name of the indicated region.
    /// \param i The index number of the requested region
    /// \param name The new name for the requested region
    /// \returns The new name of the requested region ???
    const std::string & region_name(size_t i, const std::string & name);

  /// @name Iterators
  //@{

    ///
    /// Returns a non-const begin iterator for this object's list of locii.
    iterator  begin();

    ///
    /// Returns a non-const end iterator for this object's list of locii.
    iterator  end();

  //@}

protected:
  void  interval_build();
  void  abs_distance_build();
  void  segment_distance_build();

// Member structures
  
  struct region_struct           // region extends from 1st locus to 1st locus of next region.
  {
    size_t       locus_loc;      // loci vector index of first locus in region.
    size_t       point_loc;      // points vector index of point corresponding to 1st locus in region.
    std::string  name;
    bool         x_linked;
  };

  struct locus_struct
  {
    long               locus;         // loci vector index of this locus.
    double             location;      // distance in cm from first locus in region.
    long               point_loc;     // points vector index of point corresponding to this locus.
    long               region;        // region to which locus belongs.
    const locus_info*  linfo;
  };

  struct point_struct
  {
    point_struct();
    
    double  location;     // distance in cm from first locus in region.
    long    locus;        // set to -1 if not a locus, loci vector index if a locus.
    long    region;
  };

// Data members.

  vector<region_struct>    regions;
  vector<locus_struct>     loci;
  vector<point_struct>     points;

  const RefMPedInfo&       my_rmp;
  MapTypeShPtr             my_map_func;
  
  double   distance;
  int      interval;
  bool     absolute;
  bool     method;

  bool   my_built;
  bool   my_frozen;
  
  cerrorstream   my_errors;
};


/** \brief Genome description object with parameter file parsing abilities
  *
  * the 'real' genome_description.  Calling constructor initiates
  * reading of genome desription file and populating of class
  * member containers.
  */
class LSFgenome_description : public genome_description
{
  public:

  /// @name Constructor
  //@{

    ///
    /// Constructor.
    /// \param rmp The RefMultiPedigree whose data will be used
    /// \param errors The errorstream to which error messages will be sent
    /// \param params An LSFBase pointer to the "GENOME" block in the parameter file
    /// \param multipoint Boolean value indicating whether or not use multipoint option ???
    LSFgenome_description(const RefMPedInfo& rmp, cerrormultistream errors,
                          LSFBase* params, bool multipoint = true);

  //@}

  private:
    bool init_region(LSFBase*, bool multipoint);
};

} // End namespace RPED
} // End namespace SAGE

#include "rped/genome_description.ipp"

#endif
