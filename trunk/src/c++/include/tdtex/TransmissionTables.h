#ifndef TDTEX_TABLES
#define TDTEX_TABLES

#include <string>
#include <vector>
#include "boost/shared_ptr.hpp"
#include "mlocus/genotype.h"
#include "numerics/fmatrix.h"
#include "tdtex/Transmission.h"

namespace SAGE  {
namespace TDTEX {

/** \class TransmissionTable
  * \brief Base class for storing information in a  table.
  *
  * A TransmissionTable is basically a fancy interface for a matrix. The matrix is intended
  * to stores \b counts of transmission occurences. The particular type of transmission
  * is given by the derived types of this class (AlleleTable and GenotypeTable).
  *
  */
class TransmissionTable
{
public:

  /// @name Typedefs
  //@{
  
    typedef FortranMatrix<size_t> Matrix;

    /// A vector where the i'th element should be remapped (for alphabetic ordering) as vector[i]
    typedef std::vector<size_t> ReorderingVector;

    struct ReorderedMatrix
    {
      typedef std::vector<std::string> HeadingVector;
      
      Matrix        matrix;
      HeadingVector headings;
    };

  //@}

  /// @name Constructors
  //@{

    ///
    /// Constructor
    /// \param n The size of the table (an n-by-n matrix)
    explicit TransmissionTable(size_t n);
    
    ///
    /// Destructor.
    virtual ~TransmissionTable() {}
    
  //@}

  /// @name Basic information
  //@{
  
    ///
    /// ?
    size_t  max_width(size_t min_width = 5) const;

    ///
    /// Returns a const reference to the matrix of counts.
    const Matrix& get_counts() const { return my_counts; }

    ///
    /// Returns a non-const reference to the matrix of counts.
    Matrix& get_counts() { return my_counts; }

    ///
    /// Reorders the counts matrix so that the indices are sorted according to their
    /// corresponding names (headings).
    ReorderedMatrix generate_reordered_matrix() const;

    void populate_reordering_vector(ReorderingVector & v) const;

  //@}

  /// @name Virtual interface (optional in derived class)
  //@{

    virtual void resize(size_t n);

  //@}
  
  /// @name Pure virtual interface (\b required in derived class)
  //@{

    ///
    /// Clone function.
    virtual TransmissionTable * clone() const = 0;

    ///
    /// Count up the transmissions for the given TransmissionList.
    virtual void count_transmissions(const TransmissionList& xmit) = 0;

    ///
    /// Indicates whether the given TransmissionList is informative.
    virtual bool informative(const TransmissionList& xmit) = 0;

    ///
    /// The name of the column heading for the n'th column.
    virtual string heading(size_t n) const = 0;

    ///
    /// The name of the units represented by the cells.
    virtual string units() const = 0;

  //@}

  // Debugging
  void dump() const;

protected:

  /// @name Protected copy constructor & operator=
  //@{
  
    ///
    /// Copy constructor.
    TransmissionTable(const TransmissionTable & other);
    
    ///
    /// Operator=
    TransmissionTable& operator=(const TransmissionTable & other);

  //@}

  Matrix my_counts;

};

/// A pointer to a table (since objects will be storing TransmissionTable-derivations).
typedef boost::shared_ptr<TransmissionTable> TransmissionTablePtr;

/** \class AlleleTable
  * \brief Describes a TransmissionTable whose index is based on a specific allele
  */
class AlleleTransmissionTable : public TransmissionTable
{
public:

  /// @name Constructors
  //@{
  
    ///
    /// \param gmodel The genotype model to use as the basis for this table's construction
    explicit AlleleTransmissionTable(const MLOCUS::genotype_model& gmodel);
    
    ///
    /// Copy constructor.
    AlleleTransmissionTable(const AlleleTransmissionTable & other);
    
    ///
    /// Assignment operator.
    AlleleTransmissionTable& operator=(const AlleleTransmissionTable & other);

    virtual ~AlleleTransmissionTable() {}
    
  //@}

  /// @name Required virtual interface
  //@{

    virtual AlleleTransmissionTable * clone() const;
    virtual void   count_transmissions(const TransmissionList& xmit);
    virtual bool   informative(const TransmissionList& xmit);
    virtual string heading(size_t n) const;
    virtual string units() const;

  //@}

  private:

    MLOCUS::genotype_model my_gmodel;
};

/** \class GenotypeTransmissionTable
  * \brief Describes a TransmissionTable whose index is based on a specific genotype
  */
class GenotypeTransmissionTable : public TransmissionTable
{
public:

  /// @name Constructors
  //@{
  
    explicit GenotypeTransmissionTable(const MLOCUS::genotype_model& gmodel);
  
    GenotypeTransmissionTable(const GenotypeTransmissionTable & other);

    GenotypeTransmissionTable& operator=(const GenotypeTransmissionTable & other);

    virtual ~GenotypeTransmissionTable() {}

  //@}

  virtual GenotypeTransmissionTable * clone() const;
  virtual void   count_transmissions(const TransmissionList& xmit);
  virtual bool   informative(const TransmissionList& xmit);
  virtual string heading(size_t n) const;
  virtual string units() const;

private:
  MLOCUS::genotype_model my_gmodel;
};

} // End namespace TDTEX
} // End namespace SAGE

#endif
