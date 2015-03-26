//=====================================================
//
//  File:	Field.ipp
//
//  Author:	Stephen Gross
//
// Copyright 2004 R.C. Elston
//=====================================================

namespace SAGE {
namespace SAMPLING {

//=====================================================
//  INLINE FUNCTIONS for FieldID, GroupID
//=====================================================

inline          FieldID::FieldID   (const FieldID & other) { id = other.id;               }
inline FieldID& FieldID::operator= (const FieldID & other) { id = other.id; return *this; }
inline          GroupID::GroupID   (const GroupID & other) { id = other.id;               }
inline GroupID& GroupID::operator= (const GroupID & other) { id = other.id; return *this; }

//=====================================================
//  INLINE FUNCTIONS for Field
//=====================================================

// Accessors
inline       double                Field::getOrigValue         (size_t i) const { return my_orig_values[i];        }
inline       double                Field::getReplacedValue     (size_t i) const { return my_replaced_values[i];        }
inline       double                Field::getAdjValue          (size_t i) const { return my_adj_values[i];         }
inline       bool                  Field::isAdjValuePresent    (size_t i) const { return my_adj_value_presents[i]; }
inline const std::string         & Field::getOrigName          ()         const { return my_orig_name;             }
inline const std::string         & Field::getFieldName         ()         const { return my_field_name;            }
inline const std::string         & Field::getGroupName         ()         const { return my_group_name;            }
inline const FieldID             & Field::getFieldID           ()         const { return my_field_id;              }
inline const GroupID             & Field::getGroupID           ()         const { return my_group_id;              }
inline       double                Field::getMean              ()         const { return my_sample_info.mean();              }
inline       double                Field::getVariance          ()         const { return my_sample_info.variance();          }
inline       double                Field::getStdev             ()         const { return my_sample_info.standard_deviation();             }
inline const std::vector<double> & Field::getAllOrigValues     ()         const { return my_orig_values;            }
inline const std::vector<double> & Field::getAllReplacedValues ()         const { return my_replaced_values;            }
inline const std::vector<double> & Field::getAllAdjValues      ()         const { return my_adj_values;            }
inline       unsigned long         Field::getInputFlags        ()         const { return my_input_flags;       }
inline       bool                  Field::getUserDefined       ()         const { return my_user_defined;      }
inline const SampleInfo          & Field::getSampleInfo        ()         const { return my_sample_info; }

// Mutators
inline int Field::setOrigName    (const std::string  & orig_name)  { my_orig_name    = orig_name;    return 0; }
inline int Field::setFieldName   (const std::string  & field_name) { my_field_name   = field_name;   return 0; }
inline int Field::setGroupName   (const std::string  & group_name) { my_group_name   = group_name;   return 0; }
inline int Field::setFieldID     (const FieldID & id)         { my_field_id     = id;           return 0; }
inline int Field::setGroupID     (const GroupID & id)         { my_group_id     = id;           return 0; }
inline int Field::setInputFlags  (unsigned long flags)        { my_input_flags  = flags;        return 0; }
inline int Field::setUserDefined (bool user_defined)          { my_user_defined = user_defined; return 0; }

inline int Field::setSampleInfo (const SampleInfo & sinfo) { my_sample_info = sinfo; return 0; }

inline       std::vector<double> & Field::getAllOrigValues  ()               { return my_orig_values;            }
inline       std::vector<double> & Field::getAllAdjValues   ()               { return my_adj_values;            }


} // End namespace SAMPLING
} // End namespace SAGE
