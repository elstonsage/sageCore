#ifndef __MAPPING_H
#define __MAPPING_H

//
// Mapping functions between recombination fraction and centimorgan distances
// Eventually should be moved in with the marker locus classes.

#include <limits>
#include <math.h>
#include "globals/config.h"

namespace SAGE {
namespace RPED {

/** \brief Base class for developping a recombination mapping function.
  * \par A brief review of recombination math
  *
  * Given that the distance between two locii can be expressed in pure physical
  * form (centimorgans), and that distance can be used to calculate the probability
  * of a recombination between those two locii, a mapping function can carry out
  * that distance-to-probability conversion.
  *
  * \par Virtual interface
  *
  * The (pure) virtual interface consists of two functions:
  *
  * distance() - Converts a recombination fraction into a distance (in centimorgans)
  *
  * rec_frac() - Converts a distance (in centimorgans) into a recombination fraction
  */
class Mapping_Function
{
public:

  // Required for making the compiler happy.
  virtual ~Mapping_Function() { }

  /// @name Required virtual interface
  //@{
  
    ///
    /// Calculates the physical distance between two locii given their recombination fraction.
    /// \param rf The recombination fraction
    virtual double distance(double rf) const = 0;
    
    ///
    /// Calculates the recombination fraction between two locii given their physical distance.
    /// \param d The distance
    virtual double rec_frac(double d)  const = 0;
    
  //@}
};

/** \brief Haldane mapping function
  *
  * \par Recombination calculation
  *
  * Given a distance \c d, the recombination fraction = \f$ 0.5 * ( 1.0 - exp ( -0.02 * d ) ) \f$
  *
  * \par Distance calculation
  *
  * Given a recombination fraction \c rf, the distance = \f$ -50.0 * log ( 1.0 - 2.0 * rf ) \f$
  */
class Haldane : public Mapping_Function
{
public:

  ///
  /// Calculates the physical distance between two locii given their recombination fraction.
  /// \param rf The recombination fraction
  inline virtual double distance(double rf) const;

  ///
  /// Calculates the recombination fraction between two locii given their physical distance.
  /// \param d The distance
  inline virtual double rec_frac(double d)  const;
};

/** \brief Kosambi mapping function
  *
  * \par Recombination calculation
  *
  * Given a distance \c d, the recombination fraction = \f$ tanh ( 0.02 * d ) \over 2.0 \f$
  *
  * \par Distance calculation
  *
  * Given a recombination fraction \c rf, the distance = \f$ 25.0 * log ( { 1.0 + 2.0 * rf \over 1.0 - 2.0 * rf } ) \f$
  */
class Kosambi : public Mapping_Function
{
public:

  ///  
  /// Calculates the physical distance between two locii given their recombination fraction.      
  /// \param rf The recombination fraction     
  inline virtual double distance(double rf) const;

  ///  
  /// Calculates the recombination fraction between two locii given their physical distance.      
  /// \param d The distance          
  inline virtual double rec_frac(double d)  const;
};

// ================
// Inline functions
// ================

inline double Haldane::distance(double rf) const
{
  if(rf >= 0.5) return numeric_limits<double>::infinity();

  return -50.0 * log(1.0 - 2.0 * rf);
}

inline double Haldane::rec_frac(double d) const
{
  if(!finite(d)) return 0.5;

  return 0.5 * (1.0 - exp(-0.02 * d));
}

inline double Kosambi::distance(double rf) const
{
  if(rf >= 0.5) return numeric_limits<double>::infinity();

  return 25.0 * log( (1.0 + 2.0 * rf) / (1.0 - 2.0 * rf) ); 
}

inline double Kosambi::rec_frac(double d) const
{
  if(!finite(d)) return 0.5;

  return tanh(0.02 * d) / 2.0;
}

} // End namespace RPED
} // End namespace SAGE

#endif
