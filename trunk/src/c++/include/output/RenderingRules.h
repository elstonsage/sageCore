#ifndef OUTPUT_RENDERING_RULES_H
#define OUTPUT_RENDERING_RULES_H

#include "output/Element.h"
#include "util/StringUtils.h"
#include "util/RegexUtils.h"
#include <string>
#include <set>

namespace SAGE {
namespace OUTPUT {

/// \brief Describes rules for how to render a floating point value
///
/// ALLOWABLE CHILD ELEMENTS: [None]
///
/// \par Introduction
///
/// The RenderingRules object allows you to specify how you want
/// a floating point (double) value rendered. In particular, you can
/// specify numeric format (scientific vs. fixed notation), precision,
/// You can also specify 'threshold' values, below or above which 
/// the rendered version will simply indicated greater-than or less-than.
///
/// \par Numeric formatting & precision
///
/// The number format can be set to either DEFAULT, FIXED, or SCIENTIFIC. If it is
/// FIXED, numbers will be rendered with fixed-width formatting. If it is
/// SCIENTIFIC, numbers will be rendered in scientific notation. If it is
/// DEFAULT, numbers will be rendered according to the governing logic of the
/// system. Generally speaking, this means that numbers lower than 1e-5 will be
/// rendered as scientific notation.
///
/// Numbers will always be rendered with the given precision (default is 6).
/// 
/// \par Threshold
///
/// You can set lower and upper threshold values for rendering a number. Consider
/// a p-value, for instance. Below a certain value, you may want it simply reported as 
/// "< 0.0001". To accomplish this, you can add as many lower thresholds as you want,
/// and the render() function will figure out whether a given number lies below one of
/// the thresholds.
///
/// For example:
///
/// \code
/// SAGE::OUTPUT::RenderingRules r;
/// r.addLowerThreshold(0.001);
/// r.addLowerThreshold(0.0001);
/// r.addLowerThreshold(0.00001);
///
/// std::cout << r.render(0.002)    << std::endl;
/// std::cout << r.render(0.0002)   << std::endl;
/// std::cout << r.render(0.00002)  << std::endl;
/// std::cout << r.render(0.000002) << std::endl;
/// \endcode
///
/// The above example will print:
/// \verbatim
/// 0.002
/// < 0.001
/// < 0.0001
/// < 0.00001
/// \endverbatim
///
/// \par XML Attributes
///
/// "NUMBER_FORMAT" (String - "SCIENTIFIC" / "FIXED") - The type of numeric formatting
///
/// "PRECISION" (Integer >= 0) - The number of digits to the right of the decimal place to report
///
/// "LOWER_THRESHOLD" ( ... ) -
///
/// "UPPER_THRESHOLD" ( ...) -
///
class RenderingRules : public Element<RenderingRules, AttrMgr<HasNumberFormat, HasPrecision, HasUpperThreshold, HasLowerThreshold, HasStripZeroes> >
{
public:

  typedef Element<RenderingRules, AttrMgr<HasNumberFormat, HasPrecision, HasUpperThreshold, HasLowerThreshold, HasStripZeroes> > Base;

  /// @name Constructors
  //@{

    ///
    /// Constructor.
    /// \param e The number format
    /// \param p The precision
    RenderingRules(HasNumberFormat::NumberFormat e = DEFAULT, int p = 6, bool strip_zeroes = false) :
      Base("RenderingRules")
    {
      setNumberFormat (e);
      setPrecision    (p);
      setStripZeroes  (false);
    }

    ///
    /// Copy constructor.
    RenderingRules(const RenderingRules & other) : Base(other) {}
    
    ///
    /// Assignment operator.
    RenderingRules& operator= (const RenderingRules & other)
    {
      Base::operator=(other);
      return *this;
    }
    
  //@}
  
  /// @name Rendering a value
  //@{
  
    ///
    /// Returns \c true if the given number falls below a lower
    /// threshold or above an upper threshold.
    bool exceedsThreshold(double v) const;

    template<class VALUE_TYPE>
    std::string render(const HasValue<VALUE_TYPE> & val) const
    {
      /* You have tried to render something other than HasValue<double / int / string> */ BOOST_STATIC_ASSERT((sizeof(VALUE_TYPE) == (size_t)-1));
      
      return "";
    }
    
  //@}
  
};

extern const RenderingRules default_rules;

//==================
// INLINE FUNCTIONS
//==================

//===================================
//
//  render(...) STRING
//
//===================================
template<>
inline std::string 
RenderingRules::render(const HasValue<std::string> & val) const
{
  return val.getValue();
}        

//===================================
//
//  render(...) INT
//
//===================================
template<>
inline std::string 
RenderingRules::render(const HasValue<int> & val) const
{
  // Check for rational number:
  
  if(val.getVType() != HasValue<int>::RATIONAL)
    return val.vtypeToStr();

  // Create object to store output:

  std::ostringstream s;

  double v = (double)val.getValue();

  // Check for lower threshold:

  for(HasLowerThreshold::ThresholdSet::const_iterator i = getLowerThresholds().begin(); i != getLowerThresholds().end(); ++i)
  {
    if(v < *i)
    {
      std::ostringstream t; 
      t << std::scientific << std::setprecision(2) << *i;

      t.str(UTIL::StringUtils::stripLeadingExpZero(t.str()));
      
      return "< " + t.str();
    }
  }

  // Check for upper threshold:

  if(getUpperThresholds().size())
  {
    HasUpperThreshold::ThresholdSet::const_iterator i = getUpperThresholds().end();

    while(1)
    {
      i--;

      if(v > *i)
      {
        std::ostringstream t; 
        t << std::scientific << std::setprecision(2) << *i;

        t.str(UTIL::StringUtils::stripLeadingExpZero(t.str()));
       
        return "> " + t.str();
      }
      
      if(i == getUpperThresholds().begin())
        break;
    }
  }

  // No threshold to report, do the normal thing:

  s << (v >= 0 ? " " : "") << v;

  return s.str();
}

//===================================
//
//  render(...) DOUBLE
//
//===================================
template<>
inline std::string 
RenderingRules::render(const HasValue<double> & val) const
{
  // Check for rational number:
  
  if(val.getVType() != HasValue<double>::RATIONAL)
    return val.vtypeToStr();

  // Create object to store output:

  std::ostringstream s;

  double v = val.getValue();

  // Check for lower threshold:

  for(HasLowerThreshold::ThresholdSet::const_iterator i = getLowerThresholds().begin(); i != getLowerThresholds().end(); ++i)
  {
    if(v < *i)
    {
      std::ostringstream t; 
      t << std::scientific << std::setprecision(2) << *i;

      t.str(UTIL::StringUtils::stripLeadingExpZero(t.str()));
      
      return "< " + t.str();
    }
  }

  // Check for upper threshold:

  if(getUpperThresholds().size())
  {
    HasUpperThreshold::ThresholdSet::const_iterator i = getUpperThresholds().end();

    while(1)
    {
      i--;

      if(v > *i)
      {
        std::ostringstream t; 
        t << std::scientific << std::setprecision(2) << *i;

        t.str(UTIL::StringUtils::stripLeadingExpZero(t.str()));
      
        return "> " + t.str();
      }
      
      if(i == getUpperThresholds().begin())
        break;
    }
  }

  // No threshold to report, do the normal thing:

  HasNumberFormat::NumberFormat number_format = getNumberFormat();
  int                           precision     = getPrecision();
  
  // If set to "default", let's find out what default looks like:
  
  if(getNumberFormat() == DEFAULT)
  {
    std::ostringstream temp_str;

    temp_str << v;
  
    number_format = temp_str.str().find("e") == std::string::npos ? FIXED : SCIENTIFIC;
    precision     = number_format            == SCIENTIFIC        ? 2     : 6;
  }

  // Now go ahead and generate correct output:

  s << std::setprecision(precision)
    << (v >= 0 ? " " : "")
    << (number_format == SCIENTIFIC ? std::scientific : std::fixed)
    << (v == 0 ? 0 : v); // This line is designed to correctly process a negative 0 (-0.0); it strips the unnecessary negative sign.

  // Replace e+0XX with e+XX:
  
  if(number_format == SCIENTIFIC)
  {
    s.str(UTIL::StringUtils::stripLeadingExpZero(s.str()));
  }
  
  // Strip trailing zeroes if necessary:

  std::ostringstream s2;

  if(number_format == SCIENTIFIC || ((number_format == FIXED) && (getStripZeroes() == false)))
  {
    s2 << s.str();
  }
  else
  {
    size_t decimal_idx = 0;

    for( ; decimal_idx < s.str().length(); ++decimal_idx)
      if(s.str()[decimal_idx] == '.')
        break;

    if(decimal_idx == s.str().length())
    {
      s2 << s.str();
    }
    else
    {
      size_t last_idx = s.str().length() - 1;

      for( ; last_idx > decimal_idx; --last_idx)
        if(s.str()[last_idx] != '0')
          break;

      if(s.str()[last_idx] == '.')
        last_idx--;
        
      s2 << s.str().substr(0, last_idx + 1);
    }
  }

  return s2.str();
}

} // End namespace OUTPUT
} // End namespace SAGE

#endif
