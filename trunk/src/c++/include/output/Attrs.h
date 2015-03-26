#ifndef OUTPUT_ATTRS_H
#define OUTPUT_ATTRS_H

#include "numerics/isnan.h"
#include "output/AttrMgr.h"
#include "util/RegexUtils.h"
#include <cmath>
#include <set>

namespace SAGE   {
namespace OUTPUT {

CREATE_STRING_ATTRIBUTE(Title, "TITLE")

CREATE_STRING_ATTRIBUTE(Filename, "")

CREATE_BOOL_ATTRIBUTE(FileExists, false)

CREATE_STRING_ATTRIBUTE(Tooltip, "")

CREATE_STRING_ATTRIBUTE(HelpText, "")

CREATE_ENUM_ATTRIBUTE3(GraphType, LINE_GRAPH, LINE_GRAPH, SCATTER_PLOT, BAR_CHART)

CREATE_ENUM_ATTRIBUTE2(XColumnType, NUMERICAL, NUMERICAL, CATEGORICAL)

CREATE_ENUM_ATTRIBUTE3(BulletType, UNSPECIFIED, UNSPECIFIED, NUMBERED, BULLETED)

CREATE_STRING_ATTRIBUTE(XColumn, "")

CREATE_ENUM_ATTRIBUTE3(NumberFormat, DEFAULT, DEFAULT, SCIENTIFIC, FIXED)

CREATE_INT_ATTRIBUTE(Precision, 6)

CREATE_INT_ATTRIBUTE(ByteCount, 0)

CREATE_BOOL_ATTRIBUTE(StripZeroes, false)

CREATE_BOOL_ATTRIBUTE(Visible, true)

CREATE_STRING_ATTRIBUTE(GroupName, "")

CREATE_ENUM_ATTRIBUTE3(Justification, LEFT, LEFT, CENTER, RIGHT)

CREATE_STRING_ATTRIBUTE(UnavailableCode, "Unavailable")

/// \brief YCOLUMNS attribute
///
class HasYColumns : public AttrBase
{
public:

  /// A list of column names
  typedef std::vector<std::string> ColumnList;

  HasYColumns           ()                          : AttrBase("YColumns")        { my_y_columns.clear();                            }
  HasYColumns           (const HasYColumns & other) : AttrBase(other)             { my_y_columns = other.my_y_columns;               }
  HasYColumns& operator=(const HasYColumns & other) { AttrBase::operator=(other);   my_y_columns = other.my_y_columns; return *this; }

  /// @name Querying configuration
  //@{

    ///
    /// Returns the number of y-columns to draw.
    size_t getYColumnCnt() const { return my_y_columns.size(); }

    ///
    /// Returns a const begin iterator to the list of y-columns.
    ColumnList::const_iterator beginYColumn() const { return my_y_columns.begin(); }

    ///
    /// Returns a const end iterator to the list of y-columns.
    ColumnList::const_iterator endYColumn() const { return my_y_columns.end(); }

  //@}
  
  /// @name Setting configuration
  //@{
  
    ///
    /// Add another y-column to the graph.
    /// \param y_column The name of the y-column to draw
    /// \returns \c true if successful, \c false otherwise (if the column has already been added)
    bool addYColumn(const std::string & y_column)
    {
      for(size_t i = 0; i < my_y_columns.size(); ++i)
        if(my_y_columns[i] == y_column)
          return false;
        
      my_y_columns.push_back(y_column);

      std::ostringstream s;

      for(size_t i = 0; i < my_y_columns.size(); ++i)
      {
        s << (i ? "," : "") << my_y_columns[i];
      }

      setAttrValue(s.str());

      return true;
    }
  
  //@}

private:

  ColumnList my_y_columns;  
};

/// HasValue
///
template<typename VALUE_TYPE>
class HasValue : public AttrBase
{
public:

  enum VType { STRING, RATIONAL, NON_NUMBER, INF, NEG_INF };

  typedef VALUE_TYPE val_type;

  HasValue() : AttrBase("Value") { }

  HasValue(const HasValue & other) : AttrBase(other) { my_v = other.my_v; my_t = other.my_t; }

  HasValue& operator=(const HasValue & other) { AttrBase::operator=(other); my_v = other.my_v; my_t = other.my_t; return *this; }

  void setValue(const VALUE_TYPE & v);

  const VALUE_TYPE & getValue() const { return my_v; }

  VType getVType() const { return my_t; };

  std::string vtypeToStr() const
  {
    switch(getVType())
    {
      case STRING     : return "String";
      case RATIONAL   : return "Rational";
      case NON_NUMBER : return "Non-number";
      case INF        : return "Infinity";
      case NEG_INF    : return "-Infinity";

      default         : return "";
    }
  }

private:

  std::string valToStr() const;

  void calculateVType();

  val_type my_v;
  VType    my_t;
};

/// Specialize valToStr:

template<> inline std::string HasValue<double>::valToStr() const
{
  std::ostringstream s; s << std::scientific << std::setprecision(10) << my_v; 

  // Replace e0XX with eXX: (this is for windows formatting purposes):
  
  std::string pattern = "(.+e[\\+-])0(\\d\\d)";
    
  std::vector<std::string> results = UTIL::RegexUtils::matchPattern(pattern, s.str());
    
  if(results.size() == 3)
  {
    s.str(results[1] + results[2]);
  }
  
  return s.str();
}

template<> inline std::string HasValue<int>::valToStr() const
{
  std::ostringstream s; s << my_v; return s.str();
}

template<> inline std::string HasValue<std::string>::valToStr() const
{
  return my_v;
}

/// Specialize calculateType:

template<> inline void HasValue<double>::calculateVType()
{
  my_t = finite(my_v) ? RATIONAL : (SAGE::isnan(my_v) ? NON_NUMBER : (my_v < 0 ? NEG_INF : INF));
}

template<> inline void HasValue<int>::calculateVType()
{
  my_t = finite(my_v) ? RATIONAL : (SAGE::isnan(my_v) ? NON_NUMBER : (my_v < 0 ? NEG_INF : INF));
}

template<> inline void HasValue<std::string>::calculateVType()
{
  my_t = STRING;
}

/// Specialize setValue:

template<> inline void HasValue<double>::setValue(const double & v)
{
  my_v = v;
  calculateVType();
  setAttrValue(getVType() == RATIONAL ? valToStr() : vtypeToStr());
}

template<> inline void HasValue<int>::setValue(const int & v)
{
  my_v = v;
  calculateVType();
  setAttrValue(getVType() == RATIONAL ? valToStr() : vtypeToStr());
}

template<> inline void HasValue<std::string>::setValue(const std::string & v)
{
  my_v = v;
  calculateVType();
  setAttrValue(valToStr());
}

// HasUpperThreshold
class HasUpperThreshold : public AttrBase
{
public:

  typedef std::set<double> ThresholdSet;

  HasUpperThreshold() : AttrBase("UpperThresholds") 
  {
    my_upper_thresholds.clear(); 
  }

  HasUpperThreshold(const HasUpperThreshold & other) : AttrBase(other) 
  {
    my_upper_thresholds = other.my_upper_thresholds;
  }

  HasUpperThreshold& operator=(const HasUpperThreshold & other) 
  { 
    AttrBase::operator=(other); 
    my_upper_thresholds = other.my_upper_thresholds;
    return *this; 
  }
      
  ///
  /// Adds another upper threshold.
  void addUpperThreshold(double t)
  {
    if(my_upper_thresholds.count(t))
      return;
  
    my_upper_thresholds.insert(t);
  
    std::ostringstream uppers_str;

    for(ThresholdSet::const_iterator i = my_upper_thresholds.begin(); i != my_upper_thresholds.end(); ++i)
      uppers_str << (i == my_upper_thresholds.begin() ? "" : ",") << *i;

    setAttrValue(uppers_str.str());
  }
    
  ///
  /// Returns the set of upper thresholds.
  const ThresholdSet & getUpperThresholds() const { return my_upper_thresholds; }

private:

  ThresholdSet my_upper_thresholds;

};

// HasLowerThreshold
class HasLowerThreshold : public AttrBase
{
public:

  typedef std::set<double> ThresholdSet;

  HasLowerThreshold() : AttrBase("LowerThresholds") 
  {
    my_lower_thresholds.clear(); 
  }

  HasLowerThreshold(const HasLowerThreshold & other) : AttrBase(other) 
  {
    my_lower_thresholds = other.my_lower_thresholds;
  }

  HasLowerThreshold& operator=(const HasLowerThreshold & other) 
  { 
    AttrBase::operator=(other); 
    my_lower_thresholds = other.my_lower_thresholds;
    return *this; 
  }
      
  ///
  /// Adds another lower threshold.
  void addLowerThreshold(double t)
  {
    if(my_lower_thresholds.count(t))
      return;
  
    my_lower_thresholds.insert(t);
  
    std::ostringstream lowers_str;

    for(ThresholdSet::const_iterator i = my_lower_thresholds.begin(); i != my_lower_thresholds.end(); ++i)
      lowers_str << (i == my_lower_thresholds.begin() ? "" : ",") << *i;

    setAttrValue(lowers_str.str());
  }
    
  ///
  /// Returns the set of lower thresholds.
  const ThresholdSet & getLowerThresholds() const { return my_lower_thresholds; }

private:

  ThresholdSet my_lower_thresholds;

};

extern const int screen_width;

} // End namespace OUTPUT
} // End namespace SAGE

#endif
