#ifndef FPMM_SUB_MODEL_H
#define FPMM_SUB_MODEL_H

#include "numerics/binomial_dist.h"
#include "maxfunapi/Submodel.h"
#include "error/internal_error.h"
#include "app/aparser.h"

namespace SAGE   {
namespace SEGREG {

/// @name fpmm sub-model constants
//@{
extern const std::string  FPMM_NAME;
extern const double       FPMM_DEFAULT_FREQ;
extern const size_t       FPMM_DEFAULT_LOCI;
extern const size_t       FPMM_MAX_LOCI;
extern const double       FPMM_DEFAULT_VALUE;              // 1;
extern const double       FPMM_EPSILON;
extern const double       FPMM_LB;
extern const double       FPMM_UB;
extern const bool         FPMM_DEFAULT_FIXED;
//@}

//----------------------------------------------------------------------------
//  Class:    FPMMSubmodel
//                                                                          
//----------------------------------------------------------------------------
//
class FPMMSubmodel : public MAXFUN::Submodel
{
  public:
    
  /// @name Constructor / destructor
  //@{

    FPMMSubmodel(cerrorstream& errors = sage_cerr);
    FPMMSubmodel(const FPMMSubmodel& other);
    FPMMSubmodel&  operator=(const FPMMSubmodel& other);
    virtual ~FPMMSubmodel();
    
  //@}

  /// @name Setup
  //@{

    bool  set(model_input var, double freq, size_t loci); 

    bool is_complete() const;
    
  //@}
  
  /// @name Accessors
  //@{
  
    virtual string  name() const;
    string          option_description() const;
    double          variance() const;
    double          frequency() const;
    size_t          loci() const;
    size_t          max_pgt() const;

    double          mean(size_t polygeno) const;
    double          pop_freq(size_t polygeno) const;
    double          pop_freq(size_t polygeno, size_t spouse_pg, double corr) const;
    
  //@}
    
  /// \internal
  /// @name Debugging
  //@{
  
    void  dump(std::ostream& out) const;
    
  //@}

  
  protected:
    virtual int update();
    
    // Ancillary functions.
    bool  input_meets_constraints(model_input& input);
    
    // Data members. 
    double   my_variance;     // Variance of polygenic variable.
    double   my_frequency;    // Allele frequency of polygenic loci.
    size_t   my_max_pgt;      // Max polygenotypes  2 * loci + 1

    vector<double> my_means; // Polygenic means calculated using D11
};

//=================================================================
//      INLINE FUNCTIONS
//=================================================================

//============================================================================
// IMPLEMENTATION:  FPMMSubmodel
//============================================================================
//
inline  
FPMMSubmodel::FPMMSubmodel
      (cerrorstream& errors)
    : MAXFUN::Submodel(errors)
{
  //lint -e{534}
  set(model_input(FPMM_DEFAULT_VALUE, FPMM_DEFAULT_FIXED), FPMM_DEFAULT_FREQ,
                  FPMM_DEFAULT_LOCI);
}

inline
FPMMSubmodel::FPMMSubmodel
      (const FPMMSubmodel& other)
   : MAXFUN::Submodel(other)
{
  my_variance  = other.my_variance;
  my_frequency = other.my_frequency;
  my_max_pgt   = other.my_max_pgt;
  my_means     = other.my_means;
}

inline FPMMSubmodel&
FPMMSubmodel::operator=
        (const FPMMSubmodel& other)
{
  if(this != &other)
  {
    MAXFUN::Submodel::operator=(other);

    my_variance  = other.my_variance;
    my_frequency = other.my_frequency;
    my_max_pgt   = other.my_max_pgt;
    my_means     = other.my_means;
  }
  
  return *this;
}

inline
FPMMSubmodel::~FPMMSubmodel()
{}


inline string
FPMMSubmodel::option_description() const
{
  return //lint -e(713)
         long2str(loci()) + " loci, " +
         "frequency = " + doub2str(my_frequency);
}

inline string  
FPMMSubmodel::name() const
{
  return FPMM_NAME; 
}

inline double
FPMMSubmodel::variance() const
{
  return my_variance;
}

inline double
FPMMSubmodel::mean(size_t p) const
{
  return my_means[p];
}

inline double
FPMMSubmodel::frequency() const
{
  return my_frequency;
}

inline size_t  
FPMMSubmodel::max_pgt() const
{
  return my_max_pgt;
}

inline size_t  
FPMMSubmodel::loci() const
{
  return (my_max_pgt - 1) / 2;
}

inline bool
FPMMSubmodel::is_complete() const
{
  return finite(my_variance) && my_variance > FPMM_LB + FPMM_EPSILON;
}

// - Write sub-model values in LSF readable format.
//
inline void
FPMMSubmodel::dump(std::ostream& out) const
{
  int  old_precision = out.precision();
  out.precision(DUMP_PRECISION);

  bool  var_fixed = my_parameters[0].initial_type == MAXFUN::Parameter::FIXED;

  assert(! SAGE::isnan(my_variance));

  out << "# " << name() << "\n"
      << "fpmm\n" 
      << "{\n"
      << "  loci=" << loci() << "\n"
      << "  freq=" << my_frequency << "\n"
      << "  var, val=" << my_variance << ", fixed=" << std::boolalpha << var_fixed << "\n"
      << "}" << std::noboolalpha << std::endl;
      
  out.precision(old_precision);
}

inline double          
FPMMSubmodel::pop_freq(size_t s) const
{
  double p     = frequency();
  size_t l     = loci() * 2;

  return NUMERICS::bin_prob(l,p,s);
}

inline double          
FPMMSubmodel::pop_freq(size_t i2, size_t s2, double corr) const
{
  double p          = frequency();
  size_t l          = 2 * loci();

  double sfrac      = (double) s2 / (double) l;

  double p_bar = p + corr * (sfrac - p);

  return NUMERICS::bin_prob(l, p_bar, i2);
}

} // End namespace SEGREG
} // End namespace SAGE

#endif


