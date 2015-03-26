#ifndef ASSOC_MATRIXDEFS_H
#define ASSOC_MATRIXDEFS_H

#include "numerics/fmatrix.h"
#include "maxfunapi/maxfunapi.h"
#include "output/Output.h"
#include <set>
#include <list>
#include <vector>

namespace SAGE  {
namespace ASSOC {

// A vector of names for terms-of-integration.
typedef vector<string> TermNameVector;

// Represents a pair of values: a coefficient of an integration term, 
// and the index number of that term itself.
struct CoeffPair
{
  CoeffPair() : coefficient(0), term_idx(0) {}
  CoeffPair(double _coefficient, size_t _term_idx) : coefficient(_coefficient), term_idx(_term_idx) {}
  CoeffPair(const CoeffPair& other) : coefficient(other.coefficient), term_idx(other.term_idx) {}
  CoeffPair& operator=(const CoeffPair& other) { if(this != &other) { coefficient = other.coefficient; term_idx = other.term_idx; } return *this; }

  // The coefficient of the integration term.
  double coefficient;

  // The index of the corresponding integration term (note that it is not the responsibility of the 
  // CoeffPair to keep track of those indices).
  size_t term_idx;
};

inline
ostream& operator <<(ostream& out, const CoeffPair& pair)
{
  out << "coefficient " << pair.coefficient << ", " 
      << "term index " << pair.term_idx;
      
  return  out;
}

typedef vector<CoeffPair>  CoeffPairVector;

// A vector of indices of terms.
typedef vector<size_t>  IdxVector;

// A set of indices.
typedef set<size_t>  IdxSet;

// Represents the contents of a single phi factor.
//
// or any invocation of phi (the normal density function) of the form
//
// \f[ \phi ( \beta_1 v_1 + \beta_2 v_2 + ... + \beta_k v_k, \sigma^2 ) \f]
//
// you can represent its contents as a sum of products (of two terms), along with a variance.
// That's exactly what this class represents: one invocation of the normal density function.
//
struct PointDensity
{
  PointDensity() : name(""), coeff_pairs(0), var_idxs(0) { }
  PointDensity(const string& _name, CoeffPairVector _coeff_pairs, IdxVector _var_idxs) : name(_name), coeff_pairs(_coeff_pairs), var_idxs(_var_idxs) { }
  PointDensity(const PointDensity& other) : name(other.name), coeff_pairs(other.coeff_pairs), var_idxs(other.var_idxs) { }
  PointDensity& operator=(const PointDensity& other) { if(this != &other) { name = other.name; coeff_pairs = other.coeff_pairs; var_idxs = other.var_idxs; } return *this; }
    
  void dump(const TermNameVector& term_names) const;

  // Data members
  string  name;
  CoeffPairVector  coeff_pairs;

  // The vector of integration terms (which are added together to create the variance).
  IdxVector var_idxs;
};

typedef vector<PointDensity>  PointDensities;
typedef vector<vector<double> >  Matrix;

// Represents the contents of a PointDensity, but in matrix form.
//
// Basically, this class corresponds to the A matrix (explained in both the Wedig and Gross documentation).
class PointDensityMatrix
{
  public:
    PointDensityMatrix() { }
    explicit PointDensityMatrix(const std::string& name) : my_name(name) { }
    PointDensityMatrix(const PointDensityMatrix& other) : my_name(other.my_name), my_idxs(other.my_idxs), my_matrix(other.my_matrix) { }
    PointDensityMatrix& operator=(const PointDensityMatrix& other) 
    { 
      if(this != &other) 
      { 
        my_name   = other.my_name;
        my_idxs   = other.my_idxs; 
        my_matrix = other.my_matrix; 
      } 
      
      return *this; 
    }

    // Ordering operator.
    bool operator< (const PointDensityMatrix& other)
    {
      return my_name < other.my_name;
    }

    // Returns the vector of indices; for every index i in the IdxVector, i represents
    // the column/row in the Matrix, and the vector's value the index number
    // of the corresponding term.
    IdxVector& get_idxs() { return my_idxs; }
    
    // Returns the vector of indices; for every index i in the IdxVector, i represents
    // the column/row in the Matrix, and the vector's value the index number
    // of the corresponding term.
    const IdxVector& get_idxs() const { return my_idxs; }
    
    // Returns the matrix itself. To get information about which terms correspond
    // to its rows/columns, use get_idxs().
    FortranMatrix<double>& get_matrix() { return my_matrix; }

    // Returns the matrix itself. To get information about which terms correspond
    // to its rows/columns, use get_idxs().
    const FortranMatrix<double>& get_matrix() const { return my_matrix; }

    // Sets the user-friendly name.
    void set_name(const string& name) { my_name = name; }
    
    // Gets the user-friendly name.
    const string& get_name() const { return my_name; }

    // Debugging:
    void dump() const;

    // Debugging:
    void dump(const TermNameVector& term_names) const;

  private:
    std::string  my_name;
    IdxVector    my_idxs;
    FortranMatrix<double>  my_matrix;
};


// Creates a PointDensityMatrix, using the ruleset in a PointDensity and the given standard deviation.
PointDensityMatrix generate_pdm(const PointDensity& pd);

// Adds up the variance components indicated in the PointDensity (taken from the mgr) and returns the square root of the total.
double calculate_variance(const PointDensity& pd, const MAXFUN::ParameterMgr& mgr);

// A vector of AvailablePDM's.
typedef vector<PointDensityMatrix> PDMVector;

// A vector of pairs, where:
//   the first element is an H matrix
//   the second element is a reverse lookup for idxs (index is the term-of-integration, value is the index within the matrix)
typedef vector<pair<PointDensityMatrix, IdxVector> >  PDMWithLookupVector;

// A vector of IdxSet's, where for every index i:
//   i is the index number of a term-of-integration
//   the corresponding value is a list of all the PDMs that have non-zero coefficients for that term.
typedef vector<IdxSet> TermLookupTable;

} 
} 

#endif
